/* Summarize selected YAME format-3 records into fixed CpG-ordinal bins.
 *
 * Usage: summarize_cg_bins CG RECORD_LIST BIN_WIDTH > OUTPUT.bin
 *
 * Binary output (little endian):
 *   uint64 magic = 0x314E494247434D53 ("SMCGBIN1")
 *   uint64 selected_record_count
 *   uint64 bin_width
 *   repeated selected_record_count times:
 *     uint64 physical_record
 *     uint64 number_of_bins
 *     uint32 methylated[number_of_bins]
 *     uint32 unmethylated[number_of_bins]
 *
 * Bin coordinates are ordinals in the YAME reference CpG vector.  They are
 * deliberately described as ordinal tiles, not genomic VMRs.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAGIC 0x314E494247434D53ULL
#define MAX_RECORDS 20000
#define MAX_BINS 100000

static uint64_t le(const unsigned char *p, int n) {
    uint64_t v = 0;
    for (int j = 0; j < n; ++j) v |= (uint64_t)p[j] << (8*j);
    return v;
}

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 4) die("usage: summarize_cg_bins CG RECORD_LIST BIN_WIDTH");
    uint64_t bin_width = strtoull(argv[3], NULL, 10);
    if (bin_width == 0) die("BIN_WIDTH must be positive");

    unsigned char keep[MAX_RECORDS] = {0};
    uint64_t selected = 0, requested;
    FILE *ids = fopen(argv[2], "r");
    if (!ids) die("record list missing");
    while (fscanf(ids, "%llu", (unsigned long long *)&requested) == 1) {
        if (requested >= MAX_RECORDS) die("record outside supported range");
        if (!keep[requested]) {
            keep[requested] = 1;
            ++selected;
        }
    }
    fclose(ids);

    fwrite(&(uint64_t){MAGIC}, sizeof(uint64_t), 1, stdout);
    fwrite(&selected, sizeof(uint64_t), 1, stdout);
    fwrite(&bin_width, sizeof(uint64_t), 1, stdout);

    uint32_t *methylated = calloc(MAX_BINS, sizeof(uint32_t));
    uint32_t *unmethylated = calloc(MAX_BINS, sizeof(uint32_t));
    if (!methylated || !unmethylated) die("bin allocation failed");

    gzFile input = gzopen(argv[1], "rb");
    if (!input) die("CG file missing");
    unsigned char header[17];
    uint64_t record = 0, emitted = 0;
    int got;
    while ((got = gzread(input, header, 17)) != 0) {
        if (got != 17 || le(header, 8) != 266563789635ULL || header[8] != '3')
            die("invalid YAME format-3 record");
        uint64_t bytes = le(header + 9, 8);
        if (bytes > 1073741824ULL) die("record too large");
        unsigned char *payload = malloc(bytes ? bytes : 1);
        if (!payload || gzread(input, payload, (unsigned)bytes) != (int)bytes)
            die("truncated YAME record");

        uint64_t position = 0, i = 0;
        if (record < MAX_RECORDS && keep[record]) {
            memset(methylated, 0, MAX_BINS * sizeof(uint32_t));
            memset(unmethylated, 0, MAX_BINS * sizeof(uint32_t));
        }
        while (i < bytes) {
            int type = payload[i] & 3;
            int width = type == 1 ? 1 : type == 3 ? 8 : 2;
            if (i + (uint64_t)width > bytes) die("truncated encoded value");
            uint64_t value = le(payload + i, width), m = 0, u = 0;
            i += width;
            if (type == 0) {
                position += value >> 2;
                continue;
            }
            if (type == 1) {
                m = value >> 5;
                u = (value >> 2) & 7;
            } else if (type == 2) {
                m = value >> 9;
                u = (value >> 2) & 127;
            } else {
                m = value >> 33;
                u = (value >> 2) & 2147483647ULL;
            }
            if (record < MAX_RECORDS && keep[record] && m + u) {
                uint64_t bin = position / bin_width;
                if (bin >= MAX_BINS) die("increase MAX_BINS");
                if (methylated[bin] > UINT32_MAX - m ||
                    unmethylated[bin] > UINT32_MAX - u)
                    die("bin count exceeds uint32");
                methylated[bin] += (uint32_t)m;
                unmethylated[bin] += (uint32_t)u;
            }
            ++position;
        }
        if (record < MAX_RECORDS && keep[record]) {
            uint64_t bins = (position + bin_width - 1) / bin_width;
            if (bins > MAX_BINS) die("increase MAX_BINS");
            fwrite(&record, sizeof(uint64_t), 1, stdout);
            fwrite(&bins, sizeof(uint64_t), 1, stdout);
            fwrite(methylated, sizeof(uint32_t), bins, stdout);
            fwrite(unmethylated, sizeof(uint32_t), bins, stdout);
            ++emitted;
        }
        free(payload);
        ++record;
    }
    if (!gzeof(input)) die("compressed stream read failed");
    gzclose(input);
    free(methylated);
    free(unmethylated);
    if (emitted != selected) die("not all requested records were emitted");
    return 0;
}
