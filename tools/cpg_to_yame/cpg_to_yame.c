#define _GNU_SOURCE

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

#define INDEX_MAGIC "CPGYIDX1"
#define INDEX_VERSION 1U
#define CHROM_NAME_BYTES 64U
#define YAME_CDSIG UINT64_C(266563789635)
#define YAME_FMT3_MAX_ZERO_RUN 16382U

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t chrom_count;
    uint64_t row_count;
    uint64_t source_size;
    int64_t source_mtime;
    uint64_t table_offset;
    uint64_t data_offset;
    uint64_t reserved[2];
} index_header_t;

typedef struct {
    char name[CHROM_NAME_BYTES];
    uint64_t row_start;
    uint64_t count;
    uint64_t positions_offset;
} index_chrom_t;

typedef struct {
    char name[CHROM_NAME_BYTES];
    uint32_t *positions;
    size_t count;
    size_t capacity;
    uint64_t row_start;
} build_chrom_t;

typedef struct {
    uint32_t row;
    uint32_t methylated;
    uint32_t unmethylated;
} call_t;

typedef struct {
    const index_chrom_t *disk;
    const uint32_t *positions;
    uint64_t search_hint;
    uint32_t last_call_row;
    int have_last_call;
    call_t *calls;
    size_t call_count;
    size_t call_capacity;
} chrom_state_t;

typedef struct {
    uint8_t *data;
    size_t size;
    size_t capacity;
} byte_buffer_t;

typedef struct {
    char *name;
    double sum;
    uint64_t count;
} conversion_stat_t;

typedef struct {
    int present;
    int chrom_index;
    int64_t beg;
    int64_t end;
    char ref;
    double beta;
    int64_t depth;
} meth_record_t;

typedef struct {
    int fd;
    size_t size;
    void *mapping;
    const index_header_t *header;
    const index_chrom_t *chroms;
    chrom_state_t *states;
} mapped_index_t;

typedef struct {
    mapped_index_t *index;
    chrom_state_t *chrom_states;
    int owns_chrom_states;
    meth_record_t previous;
    conversion_stat_t *conversion;
    size_t conversion_count;
    size_t conversion_capacity;
    uint64_t ignored;
} pack_state_t;

typedef struct {
    char *name;
    pack_state_t pack;
    char *cg_path;
    char *conversion_path;
    char *cg_temporary;
    char *conversion_temporary;
} batch_sample_t;

_Static_assert(sizeof(index_header_t) == 72, "unexpected index header layout");
_Static_assert(sizeof(index_chrom_t) == 88, "unexpected chromosome layout");
_Static_assert(sizeof(call_t) == 12, "unexpected call layout");

static void usage(FILE *stream)
{
    fprintf(stream,
            "Usage:\n"
            "  cpg_to_yame index REF_CPG_BED[.gz] CPG_INDEX\n"
            "  cpg_to_yame rows CPG_INDEX\n"
            "  cpg_to_yame pack -i CPG_INDEX -o OUTPUT.cg "
            "[-c CONVERSION.tsv] [CONTEXT_BED|-]\n"
            "  cpg_to_yame pack-vcf -i CPG_INDEX -o OUTPUT.cg "
            "[-c CONVERSION.tsv] [-@ THREADS] INPUT.vcf[.gz]\n"
            "  cpg_to_yame pack-vcf-batch -i CPG_INDEX -o OUTPUT_DIR "
            "[-@ THREADS] INPUT.vcf[.gz]\n");
}

static int checked_write(FILE *stream, const void *data, size_t size,
                         const char *label)
{
    if (size == 0 || fwrite(data, 1, size, stream) == size) {
        return 0;
    }
    fprintf(stderr, "Cannot write %s: %s\n", label, strerror(errno));
    return -1;
}

static int reserve_positions(build_chrom_t *chrom, size_t needed)
{
    if (needed <= chrom->capacity) {
        return 0;
    }
    size_t capacity = chrom->capacity ? chrom->capacity : 1U << 20;
    while (capacity < needed) {
        if (capacity > SIZE_MAX / 2) {
            errno = ENOMEM;
            return -1;
        }
        capacity *= 2;
    }
    uint32_t *positions = realloc(chrom->positions,
                                  capacity * sizeof(*positions));
    if (positions == NULL) {
        return -1;
    }
    chrom->positions = positions;
    chrom->capacity = capacity;
    return 0;
}

static int split_bed_prefix(char *line, char **chrom, char **beg, char **end)
{
    *chrom = line;
    char *tab1 = strchr(line, '\t');
    if (tab1 == NULL) {
        return -1;
    }
    *tab1 = '\0';
    *beg = tab1 + 1;
    char *tab2 = strchr(*beg, '\t');
    if (tab2 == NULL) {
        return -1;
    }
    *tab2 = '\0';
    *end = tab2 + 1;
    char *tail = strpbrk(*end, "\t\r\n");
    if (tail != NULL) {
        *tail = '\0';
    }
    return 0;
}

static int parse_u32(const char *text, uint32_t *value)
{
    char *end = NULL;
    errno = 0;
    unsigned long long parsed = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed > UINT32_MAX) {
        return -1;
    }
    *value = (uint32_t)parsed;
    return 0;
}

static void free_build_chroms(build_chrom_t *chroms, size_t count)
{
    if (chroms == NULL) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        free(chroms[i].positions);
    }
    free(chroms);
}

static int build_index(const char *reference_path, const char *output_path)
{
    gzFile input = gzopen(reference_path, "rb");
    if (input == NULL) {
        fprintf(stderr, "Cannot open %s\n", reference_path);
        return 1;
    }

    build_chrom_t *chroms = NULL;
    size_t chrom_count = 0;
    size_t chrom_capacity = 0;
    uint64_t row_count = 0;
    char line[65536];
    uint64_t line_number = 0;
    int status = 0;

    while (gzgets(input, line, (int)sizeof(line)) != NULL) {
        line_number++;
        size_t length = strlen(line);
        if (length == sizeof(line) - 1 && line[length - 1] != '\n'
                && !gzeof(input)) {
            fprintf(stderr, "Reference line too long at line %" PRIu64 "\n",
                    line_number);
            status = 1;
            break;
        }
        if (length == 0 || line[0] == '\n' || line[0] == '\r') {
            continue;
        }

        char *chrom_name = NULL;
        char *beg_text = NULL;
        char *end_text = NULL;
        uint32_t beg = 0;
        uint32_t end = 0;
        if (split_bed_prefix(line, &chrom_name, &beg_text, &end_text) != 0
                || parse_u32(beg_text, &beg) != 0
                || parse_u32(end_text, &end) != 0 || end <= beg) {
            fprintf(stderr, "Malformed CpG reference at line %" PRIu64 "\n",
                    line_number);
            status = 1;
            break;
        }
        if (strlen(chrom_name) >= CHROM_NAME_BYTES) {
            fprintf(stderr, "Chromosome name too long at line %" PRIu64 "\n",
                    line_number);
            status = 1;
            break;
        }

        if (chrom_count == 0
                || strcmp(chroms[chrom_count - 1].name, chrom_name) != 0) {
            for (size_t i = 0; i < chrom_count; ++i) {
                if (strcmp(chroms[i].name, chrom_name) == 0) {
                    fprintf(stderr,
                            "Chromosome %s occurs in multiple reference blocks\n",
                            chrom_name);
                    status = 1;
                    break;
                }
            }
            if (status != 0) {
                break;
            }
            if (chrom_count == chrom_capacity) {
                size_t next = chrom_capacity ? chrom_capacity * 2 : 32;
                build_chrom_t *grown = realloc(chroms,
                                              next * sizeof(*grown));
                if (grown == NULL) {
                    perror("realloc");
                    status = 1;
                    break;
                }
                chroms = grown;
                memset(chroms + chrom_capacity, 0,
                       (next - chrom_capacity) * sizeof(*chroms));
                chrom_capacity = next;
            }
            build_chrom_t *chrom = &chroms[chrom_count++];
            snprintf(chrom->name, sizeof(chrom->name), "%s", chrom_name);
            chrom->row_start = row_count;
        }

        build_chrom_t *chrom = &chroms[chrom_count - 1];
        if (chrom->count > 0 && beg <= chrom->positions[chrom->count - 1]) {
            fprintf(stderr,
                    "CpG positions are not strictly increasing at line %" PRIu64
                    "\n", line_number);
            status = 1;
            break;
        }
        if (reserve_positions(chrom, chrom->count + 1) != 0) {
            perror("realloc");
            status = 1;
            break;
        }
        chrom->positions[chrom->count++] = beg;
        row_count++;
    }

    if (status == 0 && !gzeof(input)) {
        int zerr = Z_OK;
        const char *message = gzerror(input, &zerr);
        if (zerr != Z_OK && zerr != Z_STREAM_END) {
            fprintf(stderr, "Cannot read %s: %s\n", reference_path, message);
            status = 1;
        }
    }
    if (gzclose(input) != Z_OK) {
        fprintf(stderr, "Cannot close %s\n", reference_path);
        status = 1;
    }
    if (status != 0 || row_count == 0 || chrom_count == 0) {
        free_build_chroms(chroms, chrom_count);
        return 1;
    }

    struct stat source_stat;
    if (stat(reference_path, &source_stat) != 0) {
        perror(reference_path);
        free_build_chroms(chroms, chrom_count);
        return 1;
    }
    if (chrom_count > UINT32_MAX) {
        fprintf(stderr, "Too many chromosomes\n");
        free_build_chroms(chroms, chrom_count);
        return 1;
    }

    index_header_t header = {0};
    memcpy(header.magic, INDEX_MAGIC, sizeof(header.magic));
    header.version = INDEX_VERSION;
    header.chrom_count = (uint32_t)chrom_count;
    header.row_count = row_count;
    header.source_size = (uint64_t)source_stat.st_size;
    header.source_mtime = (int64_t)source_stat.st_mtime;
    header.table_offset = sizeof(header);
    header.data_offset = sizeof(header) + chrom_count * sizeof(index_chrom_t);

    index_chrom_t *disk_chroms = calloc(chrom_count, sizeof(*disk_chroms));
    if (disk_chroms == NULL) {
        perror("calloc");
        free_build_chroms(chroms, chrom_count);
        return 1;
    }
    uint64_t offset = header.data_offset;
    for (size_t i = 0; i < chrom_count; ++i) {
        snprintf(disk_chroms[i].name, sizeof(disk_chroms[i].name), "%s",
                 chroms[i].name);
        disk_chroms[i].row_start = chroms[i].row_start;
        disk_chroms[i].count = chroms[i].count;
        disk_chroms[i].positions_offset = offset;
        offset += chroms[i].count * sizeof(uint32_t);
    }

    char *temporary_path = NULL;
    if (asprintf(&temporary_path, "%s.tmp.%ld", output_path, (long)getpid()) < 0) {
        perror("asprintf");
        free(disk_chroms);
        free_build_chroms(chroms, chrom_count);
        return 1;
    }
    FILE *output = fopen(temporary_path, "wb");
    if (output == NULL) {
        perror(temporary_path);
        free(temporary_path);
        free(disk_chroms);
        free_build_chroms(chroms, chrom_count);
        return 1;
    }

    status = checked_write(output, &header, sizeof(header), "index header")
        || checked_write(output, disk_chroms,
                         chrom_count * sizeof(*disk_chroms), "index table");
    for (size_t i = 0; status == 0 && i < chrom_count; ++i) {
        status = checked_write(output, chroms[i].positions,
                               chroms[i].count * sizeof(uint32_t),
                               "index positions");
    }
    if (status == 0 && fflush(output) != 0) {
        perror("fflush");
        status = 1;
    }
    if (status == 0 && fsync(fileno(output)) != 0) {
        perror("fsync");
        status = 1;
    }
    if (fclose(output) != 0) {
        perror("fclose");
        status = 1;
    }
    if (status == 0 && rename(temporary_path, output_path) != 0) {
        perror("rename");
        status = 1;
    }
    if (status != 0) {
        unlink(temporary_path);
    } else {
        fprintf(stderr,
                "Indexed %" PRIu64 " CpGs across %zu chromosomes (%" PRIu64
                " bytes)\n", row_count, chrom_count, offset);
    }

    free(temporary_path);
    free(disk_chroms);
    free_build_chroms(chroms, chrom_count);
    return status != 0;
}

static void free_chrom_states(const mapped_index_t *index,
                              chrom_state_t *states)
{
    if (states != NULL && index->header != NULL) {
        for (uint32_t i = 0; i < index->header->chrom_count; ++i) {
            free(states[i].calls);
        }
    }
    free(states);
}

static void close_index(mapped_index_t *index)
{
    free_chrom_states(index, index->states);
    if (index->mapping != NULL && index->mapping != MAP_FAILED) {
        munmap(index->mapping, index->size);
    }
    if (index->fd >= 0) {
        close(index->fd);
    }
    *index = (mapped_index_t){.fd = -1};
}

static chrom_state_t *clone_chrom_states(const mapped_index_t *index)
{
    chrom_state_t *states = calloc(index->header->chrom_count, sizeof(*states));
    if (states == NULL) {
        return NULL;
    }
    for (uint32_t i = 0; i < index->header->chrom_count; ++i) {
        states[i].disk = &index->chroms[i];
        states[i].positions =
            (const uint32_t *)((const uint8_t *)index->mapping
                               + index->chroms[i].positions_offset);
    }
    return states;
}

static int open_index(const char *path, mapped_index_t *index)
{
    *index = (mapped_index_t){.fd = -1};
    index->fd = open(path, O_RDONLY);
    if (index->fd < 0) {
        perror(path);
        return -1;
    }
    struct stat st;
    if (fstat(index->fd, &st) != 0 || st.st_size < (off_t)sizeof(index_header_t)) {
        fprintf(stderr, "Invalid CpG index: %s\n", path);
        close_index(index);
        return -1;
    }
    index->size = (size_t)st.st_size;
    index->mapping = mmap(NULL, index->size, PROT_READ, MAP_SHARED, index->fd, 0);
    if (index->mapping == MAP_FAILED) {
        perror("mmap");
        close_index(index);
        return -1;
    }
    index->header = index->mapping;
    if (memcmp(index->header->magic, INDEX_MAGIC, 8) != 0
            || index->header->version != INDEX_VERSION
            || index->header->table_offset != sizeof(index_header_t)
            || index->header->chrom_count == 0
            || index->header->row_count == 0) {
        fprintf(stderr, "Unsupported or corrupt CpG index: %s\n", path);
        close_index(index);
        return -1;
    }
    uint64_t table_end = index->header->table_offset
        + (uint64_t)index->header->chrom_count * sizeof(index_chrom_t);
    if (table_end > index->size || index->header->data_offset < table_end
            || index->header->data_offset > index->size) {
        fprintf(stderr, "Corrupt CpG index table: %s\n", path);
        close_index(index);
        return -1;
    }
    index->chroms = (const index_chrom_t *)((const uint8_t *)index->mapping
                                            + index->header->table_offset);
    index->states = calloc(index->header->chrom_count, sizeof(*index->states));
    if (index->states == NULL) {
        perror("calloc");
        close_index(index);
        return -1;
    }
    for (uint32_t i = 0; i < index->header->chrom_count; ++i) {
        const index_chrom_t *chrom = &index->chroms[i];
        uint64_t end = chrom->positions_offset + chrom->count * sizeof(uint32_t);
        if (chrom->name[CHROM_NAME_BYTES - 1] != '\0' || chrom->count == 0
                || chrom->positions_offset < index->header->data_offset
                || end > index->size
                || chrom->row_start + chrom->count > index->header->row_count) {
            fprintf(stderr, "Corrupt chromosome entry in %s\n", path);
            close_index(index);
            return -1;
        }
        index->states[i].disk = chrom;
        index->states[i].positions =
            (const uint32_t *)((const uint8_t *)index->mapping
                               + chrom->positions_offset);
    }
    return 0;
}

static int find_chromosome(const mapped_index_t *index, const char *name)
{
    for (uint32_t i = 0; i < index->header->chrom_count; ++i) {
        if (strcmp(index->chroms[i].name, name) == 0) {
            return (int)i;
        }
    }
    return -1;
}

static int find_position(chrom_state_t *state, uint32_t target,
                         uint32_t *row)
{
    uint64_t lo = state->search_hint;
    uint64_t count = state->disk->count;
    if (lo > count) {
        lo = count;
    }
    if (lo < count && state->positions[lo] == target) {
        *row = (uint32_t)lo;
        state->search_hint = lo + 1;
        return 1;
    }
    if (lo < count && state->positions[lo] > target) {
        return 0;
    }

    uint64_t hi = lo < count ? lo + 1 : count;
    uint64_t step = 1;
    while (hi < count && state->positions[hi] < target) {
        lo = hi + 1;
        if (step > count / 2) {
            hi = count;
        } else {
            step *= 2;
            hi = lo + step;
            if (hi > count) {
                hi = count;
            }
        }
    }
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo) / 2;
        if (state->positions[mid] < target) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    state->search_hint = lo;
    if (lo < count && state->positions[lo] == target) {
        if (lo > UINT32_MAX) {
            return -1;
        }
        *row = (uint32_t)lo;
        state->search_hint = lo + 1;
        return 1;
    }
    return 0;
}

static int append_call(chrom_state_t *state, uint32_t row,
                       uint32_t methylated, uint32_t unmethylated)
{
    if (state->have_last_call && row <= state->last_call_row) {
        fprintf(stderr, "CpG calls are not strictly increasing within %s\n",
                state->disk->name);
        return -1;
    }
    if (state->call_count == state->call_capacity) {
        size_t next = state->call_capacity ? state->call_capacity * 2 : 4096;
        call_t *grown = realloc(state->calls, next * sizeof(*grown));
        if (grown == NULL) {
            perror("realloc");
            return -1;
        }
        state->calls = grown;
        state->call_capacity = next;
    }
    state->calls[state->call_count++] = (call_t){
        .row = row,
        .methylated = methylated,
        .unmethylated = unmethylated,
    };
    state->last_call_row = row;
    state->have_last_call = 1;
    return 0;
}

static int reserve_bytes(byte_buffer_t *buffer, size_t additional)
{
    if (additional > SIZE_MAX - buffer->size) {
        errno = ENOMEM;
        return -1;
    }
    size_t needed = buffer->size + additional;
    if (needed <= buffer->capacity) {
        return 0;
    }
    size_t capacity = buffer->capacity ? buffer->capacity : 1U << 20;
    while (capacity < needed) {
        if (capacity > SIZE_MAX / 2) {
            capacity = needed;
            break;
        }
        capacity *= 2;
    }
    uint8_t *grown = realloc(buffer->data, capacity);
    if (grown == NULL) {
        return -1;
    }
    buffer->data = grown;
    buffer->capacity = capacity;
    return 0;
}

static int append_little_endian(byte_buffer_t *buffer, uint64_t value,
                                size_t width)
{
    if (reserve_bytes(buffer, width) != 0) {
        return -1;
    }
    for (size_t i = 0; i < width; ++i) {
        buffer->data[buffer->size++] = (uint8_t)(value & 0xffU);
        value >>= 8;
    }
    return 0;
}

static int append_zero_run(byte_buffer_t *buffer, uint64_t count)
{
    while (count > 0) {
        uint16_t run = count > YAME_FMT3_MAX_ZERO_RUN
            ? YAME_FMT3_MAX_ZERO_RUN : (uint16_t)count;
        if (append_little_endian(buffer, (uint64_t)run << 2, 2) != 0) {
            return -1;
        }
        count -= run;
    }
    return 0;
}

static int append_value(byte_buffer_t *buffer, uint64_t methylated,
                        uint64_t unmethylated)
{
    if (methylated < 7 && unmethylated < 7) {
        return append_little_endian(buffer,
                                    (methylated << 5) | (unmethylated << 2) | 1,
                                    1);
    }
    if (methylated < 127 && unmethylated < 127) {
        return append_little_endian(buffer,
                                    (methylated << 9) | (unmethylated << 2) | 2,
                                    2);
    }
    const uint64_t maximum = (UINT64_C(1) << 31) - 1;
    while (methylated >= maximum || unmethylated >= maximum) {
        methylated >>= 1;
        unmethylated >>= 1;
    }
    return append_little_endian(buffer,
                                (methylated << 33) | (unmethylated << 2) | 3,
                                8);
}

static int split_fields(char *line, char **fields, int capacity)
{
    int count = 0;
    char *cursor = line;
    while (count < capacity) {
        fields[count++] = cursor;
        char *tab = strchr(cursor, '\t');
        if (tab == NULL) {
            break;
        }
        *tab = '\0';
        cursor = tab + 1;
    }
    return count;
}

static int parse_i64(const char *text, int64_t *value)
{
    char *end = NULL;
    errno = 0;
    long long parsed = strtoll(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        return -1;
    }
    *value = (int64_t)parsed;
    return 0;
}

static int parse_double_value(const char *text, double *value)
{
    char *end = NULL;
    errno = 0;
    double parsed = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(parsed)) {
        return -1;
    }
    *value = parsed;
    return 0;
}

static int add_conversion(conversion_stat_t **stats, size_t *count,
                          size_t *capacity, const char *name, double beta)
{
    for (size_t i = 0; i < *count; ++i) {
        if (strcmp((*stats)[i].name, name) == 0) {
            (*stats)[i].sum += beta;
            (*stats)[i].count++;
            return 0;
        }
    }
    if (*count == *capacity) {
        size_t next = *capacity ? *capacity * 2 : 8;
        conversion_stat_t *grown = realloc(*stats, next * sizeof(*grown));
        if (grown == NULL) {
            return -1;
        }
        *stats = grown;
        *capacity = next;
    }
    char *copy = strdup(name);
    if (copy == NULL) {
        return -1;
    }
    (*stats)[*count] = (conversion_stat_t){
        .name = copy,
        .sum = beta,
        .count = 1,
    };
    (*count)++;
    return 0;
}

static int compare_conversion(const void *left, const void *right)
{
    const conversion_stat_t *a = left;
    const conversion_stat_t *b = right;
    return strcmp(a->name, b->name);
}

static int report_record(chrom_state_t *chrom_states,
                         const meth_record_t *record, double g_beta,
                         int64_t g_depth, uint64_t *ignored)
{
    double c_beta = 0.0;
    int64_t c_depth = 0;
    int64_t beg = record->beg;
    if (record->ref == 'C') {
        c_beta = record->beta;
        c_depth = record->depth;
    } else {
        g_beta = record->beta;
        g_depth = record->depth;
        beg--;
    }
    int64_t coverage = c_depth + g_depth;
    if (beg < 0 || coverage < 1 || coverage > UINT32_MAX) {
        fprintf(stderr, "Invalid merged CpG coverage or coordinate\n");
        return -1;
    }

    int64_t retained = (int64_t)rintf((float)(c_beta * c_depth))
        + (int64_t)rintf((float)(g_beta * g_depth));
    double merged_beta = retained / (double)coverage;
    char rounded_beta[64];
    if (snprintf(rounded_beta, sizeof(rounded_beta), "%.3f", merged_beta)
            >= (int)sizeof(rounded_beta)) {
        return -1;
    }
    double beta_after_text = 0.0;
    if (parse_double_value(rounded_beta, &beta_after_text) != 0) {
        return -1;
    }
    int64_t methylated = (int64_t)floor(beta_after_text * coverage + 0.5);
    if (methylated < 0 || methylated > coverage) {
        fprintf(stderr, "Invalid methylated count after strand merge\n");
        return -1;
    }

    chrom_state_t *state = &chrom_states[record->chrom_index];
    uint32_t local_row = 0;
    int found = find_position(state, (uint32_t)beg, &local_row);
    if (found < 0) {
        return -1;
    }
    if (found == 0) {
        (*ignored)++;
        return 0;
    }
    return append_call(state, local_row, (uint32_t)methylated,
                       (uint32_t)(coverage - methylated));
}

static int normalize_beta(double input, double *output)
{
    char text[64];
    if (snprintf(text, sizeof(text), "%.3f", input) >= (int)sizeof(text)) {
        return -1;
    }
    return parse_double_value(text, output);
}

static int consume_context_record(pack_state_t *state, const char *chrom,
                                  int64_t beg, int64_t end, char ref,
                                  const char *cx, const char *context2,
                                  double beta, int64_t depth)
{
    if (add_conversion(&state->conversion, &state->conversion_count,
                       &state->conversion_capacity, context2, beta) != 0) {
        perror("conversion statistics");
        return -1;
    }
    if (strcmp(cx, "CG") != 0) {
        return 0;
    }

    int chrom_index = find_chromosome(state->index, chrom);
    if (chrom_index < 0) {
        if (state->previous.present
                && report_record(state->chrom_states, &state->previous, 0.0, 0,
                                 &state->ignored) != 0) {
            return -1;
        }
        state->previous.present = 0;
        state->ignored++;
        return 0;
    }
    meth_record_t current = {
        .present = 1,
        .chrom_index = chrom_index,
        .beg = beg,
        .end = end,
        .ref = ref,
        .beta = beta,
        .depth = depth,
    };
    int is_pair = state->previous.present && state->previous.ref == 'C'
        && current.ref == 'G'
        && state->previous.chrom_index == current.chrom_index
        && current.beg == state->previous.beg + 1
        && current.end == state->previous.end + 1;
    if (is_pair) {
        if (report_record(state->chrom_states, &state->previous, current.beta,
                          current.depth, &state->ignored) != 0) {
            return -1;
        }
        state->previous.present = 0;
    } else {
        if (state->previous.present
                && report_record(state->chrom_states, &state->previous, 0.0, 0,
                                 &state->ignored) != 0) {
            return -1;
        }
        state->previous = current;
    }
    return 0;
}

static int finish_context_records(pack_state_t *state)
{
    if (state->previous.present
            && report_record(state->chrom_states, &state->previous, 0.0, 0,
                             &state->ignored) != 0) {
        return -1;
    }
    state->previous.present = 0;
    return 0;
}

static void free_pack_state(pack_state_t *state)
{
    for (size_t i = 0; i < state->conversion_count; ++i) {
        free(state->conversion[i].name);
    }
    free(state->conversion);
    state->conversion = NULL;
    state->conversion_count = 0;
    state->conversion_capacity = 0;
    if (state->owns_chrom_states) {
        free_chrom_states(state->index, state->chrom_states);
        state->chrom_states = NULL;
        state->owns_chrom_states = 0;
    }
}

static int write_conversion(const char *path, conversion_stat_t *stats,
                            size_t count)
{
    if (path == NULL) {
        return 0;
    }
    if (count > 1) {
        qsort(stats, count, sizeof(*stats), compare_conversion);
    }
    FILE *output = fopen(path, "w");
    if (output == NULL) {
        perror(path);
        return -1;
    }
    int status = 0;
    for (size_t i = 0; i < count; ++i) {
        if (fprintf(output, "%s\t%.6f\t%" PRIu64 "\n", stats[i].name,
                    stats[i].sum, stats[i].count) < 0) {
            status = -1;
            break;
        }
    }
    if (fclose(output) != 0) {
        status = -1;
    }
    if (status != 0) {
        fprintf(stderr, "Cannot write conversion summary %s\n", path);
    }
    return status;
}

static int write_yame(const char *path, const mapped_index_t *index,
                      const chrom_state_t *chrom_states)
{
    byte_buffer_t payload = {0};
    uint64_t next_row = 0;
    int status = 0;
    for (uint32_t i = 0; status == 0 && i < index->header->chrom_count; ++i) {
        const chrom_state_t *state = &chrom_states[i];
        for (size_t j = 0; j < state->call_count; ++j) {
            const call_t *call = &state->calls[j];
            uint64_t row = state->disk->row_start + call->row;
            if (row < next_row
                    || append_zero_run(&payload, row - next_row) != 0
                    || append_value(&payload, call->methylated,
                                    call->unmethylated) != 0) {
                status = -1;
                break;
            }
            next_row = row + 1;
        }
    }
    if (status == 0 && (next_row > index->header->row_count
            || append_zero_run(&payload,
                               index->header->row_count - next_row) != 0)) {
        status = -1;
    }
    if (status != 0) {
        perror("building YAME payload");
        free(payload.data);
        return -1;
    }

    BGZF *output = bgzf_open(path, "w");
    if (output == NULL) {
        fprintf(stderr, "Cannot open YAME output %s\n", path);
        free(payload.data);
        return -1;
    }
    uint64_t signature = YAME_CDSIG;
    uint8_t format = '3';
    uint64_t payload_size = payload.size;
    if (bgzf_write(output, &signature, sizeof(signature)) != sizeof(signature)
            || bgzf_write(output, &format, sizeof(format)) != sizeof(format)
            || bgzf_write(output, &payload_size, sizeof(payload_size))
                != sizeof(payload_size)
            || bgzf_write(output, payload.data, payload.size)
                != (ssize_t)payload.size) {
        fprintf(stderr, "Cannot write YAME output %s\n", path);
        status = -1;
    }
    if (bgzf_close(output) != 0) {
        fprintf(stderr, "Cannot close YAME output %s\n", path);
        status = -1;
    }
    free(payload.data);
    return status;
}

static int info_value(const char *info, const char *key,
                      const char **value, size_t *length)
{
    size_t key_length = strlen(key);
    const char *field = info;
    while (*field != '\0') {
        const char *end = strchr(field, ';');
        if (end == NULL) {
            end = field + strlen(field);
        }
        const char *equals = memchr(field, '=', (size_t)(end - field));
        if (equals != NULL && (size_t)(equals - field) == key_length
                && memcmp(field, key, key_length) == 0) {
            *value = equals + 1;
            *length = (size_t)(end - equals - 1);
            return 1;
        }
        field = *end == ';' ? end + 1 : end;
    }
    return 0;
}

static int format_index(const char *format, const char *key)
{
    size_t key_length = strlen(key);
    const char *field = format;
    int index = 0;
    while (*field != '\0') {
        const char *end = strchr(field, ':');
        if (end == NULL) {
            end = field + strlen(field);
        }
        if ((size_t)(end - field) == key_length
                && memcmp(field, key, key_length) == 0) {
            return index;
        }
        if (*end == '\0') {
            break;
        }
        field = end + 1;
        index++;
    }
    return -1;
}

static int sample_value(const char *sample, int wanted,
                        const char **value, size_t *length)
{
    const char *field = sample;
    for (int index = 0; ; ++index) {
        const char *end = field;
        while (*end != '\0' && *end != ':' && *end != '\t') {
            end++;
        }
        if (index == wanted) {
            *value = field;
            *length = (size_t)(end - field);
            return 1;
        }
        if (*end != ':') {
            return 0;
        }
        field = end + 1;
    }
}

static int parse_i64_span(const char *text, size_t length, int64_t *value)
{
    if (length == 0 || length >= 64) {
        return -1;
    }
    char buffer[64];
    memcpy(buffer, text, length);
    buffer[length] = '\0';
    return parse_i64(buffer, value);
}

static int parse_double_span(const char *text, size_t length, double *value)
{
    if (length == 0 || length >= 64) {
        return -1;
    }
    char buffer[64];
    memcpy(buffer, text, length);
    buffer[length] = '\0';
    return parse_double_value(buffer, value);
}

static int pack_vcf(const char *index_path, const char *output_path,
                    const char *conversion_path, const char *input_path,
                    int threads)
{
    mapped_index_t index = {.fd = -1};
    if (open_index(index_path, &index) != 0) {
        return 1;
    }
    htsFile *input = hts_open(input_path, "r");
    if (input == NULL) {
        fprintf(stderr, "Cannot open VCF input %s\n", input_path);
        close_index(&index);
        return 1;
    }
    if (threads > 1 && hts_set_threads(input, threads) != 0) {
        fprintf(stderr, "Cannot enable %d VCF reader threads\n", threads);
        hts_close(input);
        close_index(&index);
        return 1;
    }

    kstring_t line = {0, 0, NULL};
    uint64_t line_number = 0;
    pack_state_t pack = {
        .index = &index,
        .chrom_states = index.states,
    };
    int status = 0;
    int read_status = 0;
    char cached_format[256] = "";
    int cached_cv_index = -1;
    int cached_bt_index = -1;

    while ((read_status = hts_getline(input, '\n', &line)) >= 0) {
        line_number++;
        if (line.l == 0 || line.s[0] == '#') {
            continue;
        }
        char *fields[10] = {0};
        int field_count = split_fields(line.s, fields, 10);
        if (field_count < 10 || fields[3][0] == '\0') {
            fprintf(stderr, "Malformed VCF record at line %" PRIu64 "\n",
                    line_number);
            status = 1;
            break;
        }

        const char *cx_value = NULL;
        const char *n5_value = NULL;
        size_t cx_length = 0;
        size_t n5_length = 0;
        if (!info_value(fields[7], "CX", &cx_value, &cx_length)) {
            continue;
        }
        char ref = fields[3][0];
        if (ref != 'C' && ref != 'G') {
            continue;
        }
        info_value(fields[7], "N5", &n5_value, &n5_length);

        if (strcmp(cached_format, fields[8]) != 0) {
            size_t format_length = strlen(fields[8]);
            if (format_length >= sizeof(cached_format)) {
                fprintf(stderr, "VCF FORMAT is too long at line %" PRIu64 "\n",
                        line_number);
                status = 1;
                break;
            }
            memcpy(cached_format, fields[8], format_length + 1);
            cached_cv_index = format_index(fields[8], "CV");
            cached_bt_index = format_index(fields[8], "BT");
        }
        if (cached_cv_index < 0 || cached_bt_index < 0) {
            continue;
        }
        const char *coverage_text = NULL;
        const char *beta_text = NULL;
        size_t coverage_length = 0;
        size_t beta_length = 0;
        int64_t coverage = 0;
        double beta = 0.0;
        if (!sample_value(fields[9], cached_cv_index, &coverage_text,
                          &coverage_length)
                || parse_i64_span(coverage_text, coverage_length, &coverage) != 0
                || coverage < 1) {
            continue;
        }
        if (!sample_value(fields[9], cached_bt_index, &beta_text, &beta_length)
                || parse_double_span(beta_text, beta_length, &beta) != 0
                || normalize_beta(beta, &beta) != 0) {
            fprintf(stderr, "Missing BT value at covered VCF line %" PRIu64
                    "\n", line_number);
            status = 1;
            break;
        }

        int64_t position = 0;
        if (parse_i64(fields[1], &position) != 0 || position < 1
                || cx_length == 0 || cx_length >= 16) {
            fprintf(stderr, "Invalid cytosine VCF record at line %" PRIu64
                    "\n", line_number);
            status = 1;
            break;
        }
        char cx[16];
        memcpy(cx, cx_value, cx_length);
        cx[cx_length] = '\0';
        char context2[3] = "NN";
        if (n5_value != NULL && n5_length == 5) {
            context2[0] = n5_value[2];
            context2[1] = n5_value[3];
        }
        if (consume_context_record(&pack, fields[0], position - 1, position,
                                   ref, cx, context2, beta, coverage) != 0) {
            status = 1;
            break;
        }
    }
    if (status == 0 && read_status < -1) {
        fprintf(stderr, "Cannot read VCF input %s\n", input_path);
        status = 1;
    }
    if (status == 0 && finish_context_records(&pack) != 0) {
        status = 1;
    }
    free(line.s);
    if (hts_close(input) != 0) {
        fprintf(stderr, "Cannot close VCF input %s\n", input_path);
        status = 1;
    }

    if (status == 0 && write_conversion(conversion_path, pack.conversion,
                                        pack.conversion_count) != 0) {
        status = 1;
    }
    if (status == 0 && write_yame(output_path, &index, pack.chrom_states) != 0) {
        status = 1;
    }
    if (status != 0) {
        unlink(output_path);
        if (conversion_path != NULL) {
            unlink(conversion_path);
        }
    }
    free_pack_state(&pack);
    close_index(&index);
    return status;
}

static int split_all_fields(char *line, char ***fields, size_t *capacity)
{
    size_t count = 1;
    for (const char *cursor = line; *cursor != '\0'; ++cursor) {
        if (*cursor == '\t') {
            count++;
        }
    }
    if (count > INT_MAX) {
        return -1;
    }
    if (count > *capacity) {
        char **grown = realloc(*fields, count * sizeof(*grown));
        if (grown == NULL) {
            return -1;
        }
        *fields = grown;
        *capacity = count;
    }
    char *cursor = line;
    for (size_t i = 0; i < count; ++i) {
        (*fields)[i] = cursor;
        char *tab = strchr(cursor, '\t');
        if (tab == NULL) {
            break;
        }
        *tab = '\0';
        cursor = tab + 1;
    }
    return (int)count;
}

static int valid_sample_name(const char *name)
{
    size_t length = strlen(name);
    if (length == 0 || length > 200
            || strcmp(name, ".") == 0 || strcmp(name, "..") == 0) {
        return 0;
    }
    for (const unsigned char *cursor = (const unsigned char *)name;
         *cursor != '\0'; ++cursor) {
        if (!(('a' <= *cursor && *cursor <= 'z')
                || ('A' <= *cursor && *cursor <= 'Z')
                || ('0' <= *cursor && *cursor <= '9')
                || *cursor == '_' || *cursor == '-' || *cursor == '.')) {
            return 0;
        }
    }
    return 1;
}

static void free_batch_samples(batch_sample_t *samples, size_t count)
{
    if (samples == NULL) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        free_pack_state(&samples[i].pack);
        free(samples[i].name);
        free(samples[i].cg_path);
        free(samples[i].conversion_path);
        if (samples[i].cg_temporary != NULL) {
            unlink(samples[i].cg_temporary);
        }
        if (samples[i].conversion_temporary != NULL) {
            unlink(samples[i].conversion_temporary);
        }
        free(samples[i].cg_temporary);
        free(samples[i].conversion_temporary);
    }
    free(samples);
}

static batch_sample_t *initialize_batch_samples(mapped_index_t *index,
                                                char **names, size_t count,
                                                const char *output_dir)
{
    if (count == 0 || count > 4096) {
        fprintf(stderr, "Invalid VCF sample count: %zu\n", count);
        return NULL;
    }
    batch_sample_t *samples = calloc(count, sizeof(*samples));
    if (samples == NULL) {
        return NULL;
    }
    for (size_t i = 0; i < count; ++i) {
        if (!valid_sample_name(names[i])) {
            fprintf(stderr, "Unsafe VCF sample name: %s\n", names[i]);
            free_batch_samples(samples, count);
            return NULL;
        }
        for (size_t j = 0; j < i; ++j) {
            if (strcmp(names[i], names[j]) == 0) {
                fprintf(stderr, "Duplicate VCF sample name: %s\n", names[i]);
                free_batch_samples(samples, count);
                return NULL;
            }
        }
        samples[i].name = strdup(names[i]);
        samples[i].pack = (pack_state_t){
            .index = index,
            .chrom_states = clone_chrom_states(index),
            .owns_chrom_states = 1,
        };
        if (samples[i].name == NULL || samples[i].pack.chrom_states == NULL
                || asprintf(&samples[i].cg_path, "%s/%s.cg", output_dir,
                            names[i]) < 0
                || asprintf(&samples[i].conversion_path, "%s/%s.conv.tsv",
                            output_dir, names[i]) < 0
                || asprintf(&samples[i].cg_temporary, "%s.tmp.%ld",
                            samples[i].cg_path, (long)getpid()) < 0
                || asprintf(&samples[i].conversion_temporary, "%s.tmp.%ld",
                            samples[i].conversion_path, (long)getpid()) < 0) {
            perror("batch sample allocation");
            free_batch_samples(samples, count);
            return NULL;
        }
    }
    return samples;
}

static int pack_vcf_batch(const char *index_path, const char *output_dir,
                          const char *input_path, int threads)
{
    struct stat output_stat;
    if (stat(output_dir, &output_stat) != 0 || !S_ISDIR(output_stat.st_mode)) {
        fprintf(stderr, "Batch output directory does not exist: %s\n",
                output_dir);
        return 1;
    }
    mapped_index_t index = {.fd = -1};
    if (open_index(index_path, &index) != 0) {
        return 1;
    }
    htsFile *input = hts_open(input_path, "r");
    if (input == NULL) {
        fprintf(stderr, "Cannot open VCF input %s\n", input_path);
        close_index(&index);
        return 1;
    }
    if (threads > 1 && hts_set_threads(input, threads) != 0) {
        fprintf(stderr, "Cannot enable %d VCF reader threads\n", threads);
        hts_close(input);
        close_index(&index);
        return 1;
    }

    kstring_t line = {0, 0, NULL};
    char **fields = NULL;
    size_t field_capacity = 0;
    uint64_t line_number = 0;
    batch_sample_t *samples = NULL;
    size_t sample_count = 0;
    int status = 0;
    int read_status = 0;
    char cached_format[256] = "";
    int cached_cv_index = -1;
    int cached_bt_index = -1;

    while ((read_status = hts_getline(input, '\n', &line)) >= 0) {
        line_number++;
        if (line.l == 0 || (line.s[0] == '#' && line.s[1] == '#')) {
            continue;
        }
        int field_count = split_all_fields(line.s, &fields, &field_capacity);
        if (field_count < 0) {
            perror("VCF fields");
            status = 1;
            break;
        }
        if (line.s[0] == '#') {
            if (strcmp(fields[0], "#CHROM") != 0 || field_count < 10
                    || samples != NULL) {
                fprintf(stderr, "Malformed VCF sample header at line %" PRIu64
                        "\n", line_number);
                status = 1;
                break;
            }
            sample_count = (size_t)field_count - 9;
            samples = initialize_batch_samples(&index, fields + 9,
                                               sample_count, output_dir);
            if (samples == NULL) {
                status = 1;
                break;
            }
            continue;
        }
        if (samples == NULL || field_count < 9 + (int)sample_count
                || fields[3][0] == '\0') {
            fprintf(stderr, "Malformed VCF record at line %" PRIu64 "\n",
                    line_number);
            status = 1;
            break;
        }

        const char *cx_value = NULL;
        const char *n5_value = NULL;
        size_t cx_length = 0;
        size_t n5_length = 0;
        if (!info_value(fields[7], "CX", &cx_value, &cx_length)) {
            continue;
        }
        char ref = fields[3][0];
        if (ref != 'C' && ref != 'G') {
            continue;
        }
        info_value(fields[7], "N5", &n5_value, &n5_length);

        if (strcmp(cached_format, fields[8]) != 0) {
            size_t format_length = strlen(fields[8]);
            if (format_length >= sizeof(cached_format)) {
                fprintf(stderr, "VCF FORMAT is too long at line %" PRIu64 "\n",
                        line_number);
                status = 1;
                break;
            }
            memcpy(cached_format, fields[8], format_length + 1);
            cached_cv_index = format_index(fields[8], "CV");
            cached_bt_index = format_index(fields[8], "BT");
        }
        if (cached_cv_index < 0 || cached_bt_index < 0) {
            continue;
        }
        int64_t position = 0;
        if (parse_i64(fields[1], &position) != 0 || position < 1
                || cx_length == 0 || cx_length >= 16) {
            fprintf(stderr, "Invalid cytosine VCF record at line %" PRIu64
                    "\n", line_number);
            status = 1;
            break;
        }
        char cx[16];
        memcpy(cx, cx_value, cx_length);
        cx[cx_length] = '\0';
        char context2[3] = "NN";
        if (n5_value != NULL && n5_length == 5) {
            context2[0] = n5_value[2];
            context2[1] = n5_value[3];
        }

        for (size_t i = 0; i < sample_count; ++i) {
            const char *coverage_text = NULL;
            const char *beta_text = NULL;
            size_t coverage_length = 0;
            size_t beta_length = 0;
            int64_t coverage = 0;
            double beta = 0.0;
            if (!sample_value(fields[9 + i], cached_cv_index, &coverage_text,
                              &coverage_length)
                    || parse_i64_span(coverage_text, coverage_length,
                                      &coverage) != 0
                    || coverage < 1) {
                continue;
            }
            if (!sample_value(fields[9 + i], cached_bt_index, &beta_text,
                              &beta_length)
                    || parse_double_span(beta_text, beta_length, &beta) != 0
                    || normalize_beta(beta, &beta) != 0) {
                fprintf(stderr,
                        "Missing BT value for %s at covered VCF line %" PRIu64
                        "\n", samples[i].name, line_number);
                status = 1;
                break;
            }
            if (consume_context_record(&samples[i].pack, fields[0], position - 1,
                                       position, ref, cx, context2, beta,
                                       coverage) != 0) {
                status = 1;
                break;
            }
        }
        if (status != 0) {
            break;
        }
    }
    if (status == 0 && read_status < -1) {
        fprintf(stderr, "Cannot read VCF input %s\n", input_path);
        status = 1;
    }
    if (status == 0 && samples == NULL) {
        fprintf(stderr, "VCF sample header not found in %s\n", input_path);
        status = 1;
    }
    for (size_t i = 0; status == 0 && i < sample_count; ++i) {
        if (finish_context_records(&samples[i].pack) != 0) {
            status = 1;
        }
    }
    free(fields);
    free(line.s);
    if (hts_close(input) != 0) {
        fprintf(stderr, "Cannot close VCF input %s\n", input_path);
        status = 1;
    }

    for (size_t i = 0; status == 0 && i < sample_count; ++i) {
        if (write_conversion(samples[i].conversion_temporary,
                             samples[i].pack.conversion,
                             samples[i].pack.conversion_count) != 0
                || write_yame(samples[i].cg_temporary, &index,
                              samples[i].pack.chrom_states) != 0) {
            status = 1;
        }
    }
    for (size_t i = 0; status == 0 && i < sample_count; ++i) {
        if (rename(samples[i].cg_temporary, samples[i].cg_path) != 0
                || rename(samples[i].conversion_temporary,
                          samples[i].conversion_path) != 0) {
            perror("renaming batch output");
            status = 1;
        }
    }

    free_batch_samples(samples, sample_count);
    close_index(&index);
    return status;
}

static int pack_context(const char *index_path, const char *output_path,
                        const char *conversion_path, const char *input_path)
{
    mapped_index_t index = {.fd = -1};
    if (open_index(index_path, &index) != 0) {
        return 1;
    }
    FILE *input = stdin;
    if (input_path != NULL && strcmp(input_path, "-") != 0) {
        input = fopen(input_path, "r");
        if (input == NULL) {
            perror(input_path);
            close_index(&index);
            return 1;
        }
    }

    char *line = NULL;
    size_t line_capacity = 0;
    ssize_t line_length = 0;
    uint64_t line_number = 0;
    pack_state_t pack = {
        .index = &index,
        .chrom_states = index.states,
    };
    int status = 0;

    while ((line_length = getline(&line, &line_capacity, input)) >= 0) {
        line_number++;
        while (line_length > 0
                && (line[line_length - 1] == '\n'
                    || line[line_length - 1] == '\r')) {
            line[--line_length] = '\0';
        }
        if (line_length == 0) {
            continue;
        }
        char *fields[10] = {0};
        int field_count = split_fields(line, fields, 10);
        double beta = 0.0;
        int64_t beg = 0;
        int64_t end = 0;
        int64_t depth = 0;
        if (field_count < 9 || parse_i64(fields[1], &beg) != 0
                || parse_i64(fields[2], &end) != 0
                || parse_double_value(fields[7], &beta) != 0
                || parse_i64(fields[8], &depth) != 0
                || beg < 0 || end <= beg || depth < 1
                || fields[3][1] != '\0'
                || (fields[3][0] != 'C' && fields[3][0] != 'G')) {
            fprintf(stderr, "Malformed context BED record at line %" PRIu64
                    "\n", line_number);
            status = 1;
            break;
        }
        if (consume_context_record(&pack, fields[0], beg, end, fields[3][0],
                                   fields[4], fields[5], beta, depth) != 0) {
            status = 1;
            break;
        }
    }
    if (status == 0 && ferror(input)) {
        perror("context BED input");
        status = 1;
    }
    if (status == 0 && finish_context_records(&pack) != 0) {
        status = 1;
    }
    if (input != stdin && fclose(input) != 0) {
        status = 1;
    }
    free(line);

    if (status == 0 && write_conversion(conversion_path, pack.conversion,
                                        pack.conversion_count) != 0) {
        status = 1;
    }
    if (status == 0 && write_yame(output_path, &index, pack.chrom_states) != 0) {
        status = 1;
    }
    if (status != 0) {
        unlink(output_path);
        if (conversion_path != NULL) {
            unlink(conversion_path);
        }
    }
    free_pack_state(&pack);
    close_index(&index);
    return status;
}

int main(int argc, char **argv)
{
    setlocale(LC_NUMERIC, "C");
    if (argc >= 2 && strcmp(argv[1], "index") == 0) {
        if (argc != 4) {
            usage(stderr);
            return 2;
        }
        return build_index(argv[2], argv[3]);
    }
    if (argc >= 2 && strcmp(argv[1], "rows") == 0) {
        if (argc != 3) {
            usage(stderr);
            return 2;
        }
        mapped_index_t index = {.fd = -1};
        if (open_index(argv[2], &index) != 0) {
            return 1;
        }
        printf("%" PRIu64 "\n", index.header->row_count);
        close_index(&index);
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "pack") == 0) {
        const char *index_path = NULL;
        const char *output_path = NULL;
        const char *conversion_path = NULL;
        int option = 0;
        optind = 2;
        while ((option = getopt(argc, argv, "i:o:c:h")) != -1) {
            switch (option) {
            case 'i': index_path = optarg; break;
            case 'o': output_path = optarg; break;
            case 'c': conversion_path = optarg; break;
            case 'h': usage(stdout); return 0;
            default: usage(stderr); return 2;
            }
        }
        if (index_path == NULL || output_path == NULL || argc - optind > 1) {
            usage(stderr);
            return 2;
        }
        const char *input_path = optind < argc ? argv[optind] : "-";
        return pack_context(index_path, output_path, conversion_path, input_path);
    }
    if (argc >= 2 && strcmp(argv[1], "pack-vcf") == 0) {
        const char *index_path = NULL;
        const char *output_path = NULL;
        const char *conversion_path = NULL;
        int threads = 1;
        int option = 0;
        optind = 2;
        while ((option = getopt(argc, argv, "i:o:c:@:h")) != -1) {
            switch (option) {
            case 'i': index_path = optarg; break;
            case 'o': output_path = optarg; break;
            case 'c': conversion_path = optarg; break;
            case '@':
                threads = atoi(optarg);
                if (threads < 1) {
                    fprintf(stderr, "VCF reader threads must be positive\n");
                    return 2;
                }
                break;
            case 'h': usage(stdout); return 0;
            default: usage(stderr); return 2;
            }
        }
        if (index_path == NULL || output_path == NULL || argc - optind != 1) {
            usage(stderr);
            return 2;
        }
        return pack_vcf(index_path, output_path, conversion_path, argv[optind],
                        threads);
    }
    if (argc >= 2 && strcmp(argv[1], "pack-vcf-batch") == 0) {
        const char *index_path = NULL;
        const char *output_dir = NULL;
        int threads = 1;
        int option = 0;
        optind = 2;
        while ((option = getopt(argc, argv, "i:o:@:h")) != -1) {
            switch (option) {
            case 'i': index_path = optarg; break;
            case 'o': output_dir = optarg; break;
            case '@':
                threads = atoi(optarg);
                if (threads < 1) {
                    fprintf(stderr, "VCF reader threads must be positive\n");
                    return 2;
                }
                break;
            case 'h': usage(stdout); return 0;
            default: usage(stderr); return 2;
            }
        }
        if (index_path == NULL || output_dir == NULL || argc - optind != 1) {
            usage(stderr);
            return 2;
        }
        return pack_vcf_batch(index_path, output_dir, argv[optind], threads);
    }
    usage(argc == 2 && strcmp(argv[1], "-h") == 0 ? stdout : stderr);
    return argc == 2 && strcmp(argv[1], "-h") == 0 ? 0 : 2;
}
