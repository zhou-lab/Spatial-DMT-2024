/* Stream a coordinate-sorted BAM into per-cell BAMs using one shared BGZF pool. */
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#define COORD_VALUES 256
#define CELL_SLOTS (COORD_VALUES * COORD_VALUES)
#define MALFORMED_SLOT CELL_SLOTS
#define UNTAGGED_SLOT (CELL_SLOTS + 1)
#define OUTPUT_SLOTS (CELL_SLOTS + 2)

typedef struct {
    samFile *file;
    char *path;
} output_t;

static void usage(FILE *stream) {
    fprintf(stream,
            "Usage: split_bam_by_cb_hts [options] INPUT.bam OUTPUT_DIR\n"
            "Options:\n"
            "  -@, --threads INT  shared HTSlib compression/decompression threads [8]\n"
            "  -t, --tag TAG      two-character auxiliary tag [CB]\n"
            "  -b, --buckets INT  route cells into this many BAM buckets instead\n"
            "  -h, --help         show this help\n");
}

static int raise_fd_limit(rlim_t target) {
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) != 0) return -1;
    if (limit.rlim_cur >= target) return 0;
    limit.rlim_cur = target < limit.rlim_max ? target : limit.rlim_max;
    return setrlimit(RLIMIT_NOFILE, &limit);
}

static int base4(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

static int decode_cb(const char *cb, char label[16]) {
    if (!cb || strlen(cb) != 8) {
        strcpy(label, "MALFORMED");
        return MALFORMED_SLOT;
    }

    unsigned x = 0, y = 0;
    for (int i = 0; i < 8; i++) {
        int digit = base4(cb[i]);
        if (digit < 0) {
            strcpy(label, "MALFORMED");
            return MALFORMED_SLOT;
        }
        if (i < 4) x = x * 4 + (unsigned)digit;
        else y = y * 4 + (unsigned)digit;
    }
    x++;
    y++;
    snprintf(label, 16, "%02u%02u", x, y);
    return (int)((x - 1) * COORD_VALUES + (y - 1));
}

static samFile *open_output(output_t *output, const char *out_dir,
                            const char *label, const sam_hdr_t *header,
                            htsThreadPool *pool) {
    size_t path_len = strlen(out_dir) + strlen(label) + 6;
    output->path = malloc(path_len);
    if (!output->path) return NULL;
    snprintf(output->path, path_len, "%s/%s.bam", out_dir, label);

    output->file = sam_open(output->path, "wb");
    if (!output->file) return NULL;
    if (pool->pool && hts_set_opt(output->file, HTS_OPT_THREAD_POOL, pool) < 0) {
        return NULL;
    }
    if (sam_hdr_write(output->file, header) < 0) return NULL;
    return output->file;
}

int main(int argc, char **argv) {
    int threads = 8;
    int buckets = 0;
    const char *tag = "CB";
    static const struct option long_options[] = {
        {"threads", required_argument, NULL, '@'},
        {"tag", required_argument, NULL, 't'},
        {"buckets", required_argument, NULL, 'b'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0},
    };

    int option;
    while ((option = getopt_long(argc, argv, "@:t:b:h", long_options, NULL)) != -1) {
        switch (option) {
            case '@': threads = atoi(optarg); break;
            case 't': tag = optarg; break;
            case 'b': buckets = atoi(optarg); break;
            case 'h': usage(stdout); return 0;
            default: usage(stderr); return 2;
        }
    }
    if (argc - optind != 2 || threads < 0 || buckets < 0
            || buckets > CELL_SLOTS || strlen(tag) != 2) {
        usage(stderr);
        return 2;
    }

    const char *input_path = argv[optind];
    const char *out_dir = argv[optind + 1];
    if (mkdir(out_dir, 0775) != 0 && errno != EEXIST) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot create %s: %s\n",
                out_dir, strerror(errno));
        return 1;
    }
    if (raise_fd_limit(65535) != 0) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot raise file limit: %s\n",
                strerror(errno));
        return 1;
    }

    int status = 1;
    samFile *input = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *record = NULL;
    output_t *outputs = calloc(OUTPUT_SLOTS, sizeof(*outputs));
    htsThreadPool pool = {NULL, 0};
    uint64_t n_reads = 0, n_tagged = 0, n_untagged = 0, n_malformed = 0;
    uint64_t n_cells = 0;

    if (!outputs) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot allocate output table\n");
        goto cleanup;
    }
    if (threads > 0) {
        pool.pool = hts_tpool_init(threads);
        if (!pool.pool) {
            fprintf(stderr, "split_bam_by_cb_hts: cannot create %d threads\n", threads);
            goto cleanup;
        }
        pool.qsize = threads * 2;
    }

    input = sam_open(input_path, "rb");
    if (!input) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot open %s\n", input_path);
        goto cleanup;
    }
    if (pool.pool && hts_set_opt(input, HTS_OPT_THREAD_POOL, &pool) < 0) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot attach input thread pool\n");
        goto cleanup;
    }
    header = sam_hdr_read(input);
    if (!header) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot read input header\n");
        goto cleanup;
    }
    record = bam_init1();
    if (!record) {
        fprintf(stderr, "split_bam_by_cb_hts: cannot allocate BAM record\n");
        goto cleanup;
    }

    int read_status;
    while ((read_status = sam_read1(input, header, record)) >= 0) {
        int slot;
        char label[16];
        uint8_t *aux = bam_aux_get(record, tag);
        n_reads++;
        if (!aux || !bam_aux2Z(aux)) {
            slot = UNTAGGED_SLOT;
            strcpy(label, "UNTAGGED");
            n_untagged++;
        } else {
            n_tagged++;
            slot = decode_cb(bam_aux2Z(aux), label);
            if (slot == MALFORMED_SLOT) n_malformed++;
        }

        if (buckets > 0) {
            slot %= buckets;
            snprintf(label, sizeof(label), "bucket_%05d", slot);
        }

        if (!outputs[slot].file) {
            if (!open_output(&outputs[slot], out_dir, label, header, &pool)) {
                fprintf(stderr, "split_bam_by_cb_hts: cannot open output %s\n",
                        outputs[slot].path ? outputs[slot].path : label);
                goto cleanup;
            }
            n_cells++;
        }
        if (sam_write1(outputs[slot].file, header, record) < 0) {
            fprintf(stderr, "split_bam_by_cb_hts: write failed for %s\n",
                    outputs[slot].path);
            goto cleanup;
        }
    }
    if (read_status < -1) {
        fprintf(stderr, "split_bam_by_cb_hts: failed while reading %s\n", input_path);
        goto cleanup;
    }
    status = 0;

cleanup:
    if (record) bam_destroy1(record);
    if (input && sam_close(input) < 0) status = 1;
    if (outputs) {
        for (int i = 0; i < OUTPUT_SLOTS; i++) {
            if (outputs[i].file && sam_close(outputs[i].file) < 0) {
                fprintf(stderr, "split_bam_by_cb_hts: close failed for %s\n",
                        outputs[i].path);
                status = 1;
            }
            free(outputs[i].path);
        }
    }
    if (pool.pool) hts_tpool_destroy(pool.pool);
    if (header) sam_hdr_destroy(header);
    free(outputs);

    if (status == 0) {
        fprintf(stderr,
                "reads_total=%" PRIu64 "  tagged=%" PRIu64
                "  untagged=%" PRIu64 "  malformed=%" PRIu64
                "  cells=%" PRIu64 "\n",
                n_reads, n_tagged, n_untagged, n_malformed, n_cells);
    }
    return status;
}
