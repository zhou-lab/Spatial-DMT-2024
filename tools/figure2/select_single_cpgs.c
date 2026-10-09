/* Select broadly covered and typical single CpGs from a multi-record YAME
 * format-3 file, then emit their M/U counts for every record.
 *
 * Usage: select_single_cpgs CG VALUES.tsv SITES.tsv HISTOGRAM.tsv [RECORD_LIST]
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAX_SITES 30000000U
#define MAX_RECORDS 20000U

static uint64_t le(const unsigned char *p, int n) {
    uint64_t value = 0;
    for (int i = 0; i < n; ++i) value |= (uint64_t)p[i] << (8 * i);
    return value;
}

static void die(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static int compare_u32(const void *left, const void *right) {
    uint32_t a = *(const uint32_t *)left;
    uint32_t b = *(const uint32_t *)right;
    return (a > b) - (a < b);
}

static uint32_t median_depth_at_breadth(
    const uint16_t *covered,
    const uint32_t *depth,
    uint64_t number_of_sites,
    uint16_t breadth,
    uint64_t class_size
) {
    uint32_t *depths = malloc(class_size * sizeof(uint32_t));
    if (!depths) die("class-depth allocation failed");
    uint64_t j = 0;
    for (uint64_t site = 0; site < number_of_sites; ++site)
        if (covered[site] == breadth) depths[j++] = depth[site];
    if (j != class_size) die("breadth-class size mismatch");
    qsort(depths, class_size, sizeof(uint32_t), compare_u32);
    uint32_t median = depths[(class_size - 1) / 2];
    free(depths);
    return median;
}

static void decode_value(uint64_t value, int type, uint64_t *m, uint64_t *u) {
    if (type == 1) {
        *m = value >> 5;
        *u = (value >> 2) & 7U;
    } else if (type == 2) {
        *m = value >> 9;
        *u = (value >> 2) & 127U;
    } else {
        *m = value >> 33;
        *u = (value >> 2) & 2147483647ULL;
    }
}

int main(int argc, char **argv) {
    if (argc != 5 && argc != 6)
        die("usage: select_single_cpgs CG VALUES.tsv SITES.tsv HISTOGRAM.tsv [RECORD_LIST]");

    unsigned char *keep = calloc(MAX_RECORDS, 1);
    if (!keep) die("record-mask allocation failed");
    int use_filter = argc == 6;
    if (use_filter) {
        FILE *record_list = fopen(argv[5], "r");
        if (!record_list) die("record list missing");
        uint64_t record_id;
        while (fscanf(record_list, "%" SCNu64, &record_id) == 1) {
            if (record_id >= MAX_RECORDS) die("record outside supported range");
            keep[record_id] = 1;
        }
        fclose(record_list);
    }

    uint16_t *covered = calloc(MAX_SITES, sizeof(uint16_t));
    uint32_t *depth = calloc(MAX_SITES, sizeof(uint32_t));
    uint32_t *methylated = calloc(MAX_SITES, sizeof(uint32_t));
    if (!covered || !depth || !methylated) die("site-stat allocation failed");

    gzFile input = gzopen(argv[1], "rb");
    if (!input) die("CG file missing");
    unsigned char header[17];
    uint64_t records = 0, selected_records = 0, number_of_sites = 0;
    int got;
    while ((got = gzread(input, header, 17)) != 0) {
        if (got != 17 || le(header, 8) != 266563789635ULL || header[8] != '3')
            die("invalid YAME format-3 record");
        uint64_t bytes = le(header + 9, 8);
        if (bytes > 1073741824ULL) die("record too large");
        unsigned char *payload = malloc(bytes ? bytes : 1);
        if (!payload || gzread(input, payload, (unsigned)bytes) != (int)bytes)
            die("truncated YAME record");

        uint64_t position = 0;
        for (uint64_t i = 0; i < bytes;) {
            int type = payload[i] & 3;
            int width = type == 1 ? 1 : type == 3 ? 8 : 2;
            if (i + (uint64_t)width > bytes) die("truncated encoded value");
            uint64_t value = le(payload + i, width), m = 0, u = 0;
            i += (uint64_t)width;
            if (type == 0) {
                position += value >> 2;
                continue;
            }
            if (position >= MAX_SITES) die("increase MAX_SITES");
            decode_value(value, type, &m, &u);
            uint64_t cov = m + u;
            if ((!use_filter || keep[records]) && cov) {
                if (covered[position] == UINT16_MAX) die("pixel coverage count overflow");
                if (depth[position] > UINT32_MAX - cov ||
                    methylated[position] > UINT32_MAX - m)
                    die("site aggregate count overflow");
                covered[position]++;
                depth[position] += (uint32_t)cov;
                methylated[position] += (uint32_t)m;
            }
            position++;
        }
        if (records == 0) number_of_sites = position;
        else if (position != number_of_sites) die("inconsistent logical site count");
        free(payload);
        if (!use_filter || keep[records]) selected_records++;
        records++;
        if (records > MAX_RECORDS) die("increase MAX_RECORDS");
    }
    if (!gzeof(input)) die("compressed stream read failed");
    gzclose(input);
    if (!records || !selected_records || !number_of_sites) die("empty input or selection");

    uint64_t *histogram = calloc(selected_records + 1, sizeof(uint64_t));
    if (!histogram) die("histogram allocation failed");
    uint64_t observed_sites = 0;
    uint16_t maximum_breadth = 0;
    for (uint64_t site = 0; site < number_of_sites; ++site) {
        histogram[covered[site]]++;
        if (covered[site]) observed_sites++;
        if (covered[site] > maximum_breadth) maximum_breadth = covered[site];
    }

    uint64_t median_rank = (observed_sites - 1) / 2;
    uint64_t top_decile_midpoint_rank = (uint64_t)(0.95 * (observed_sites - 1));
    uint64_t cumulative = 0;
    uint16_t median_breadth = 0, top_decile_breadth = 0;
    for (uint64_t breadth = 1; breadth <= selected_records; ++breadth) {
        cumulative += histogram[breadth];
        if (!median_breadth && cumulative > median_rank) {
            median_breadth = (uint16_t)breadth;
        }
        if (!top_decile_breadth && cumulative > top_decile_midpoint_rank)
            top_decile_breadth = (uint16_t)breadth;
        if (median_breadth && top_decile_breadth) break;
    }

    uint64_t typical_class_size = histogram[median_breadth];
    uint64_t top_decile_class_size = histogram[top_decile_breadth];
    uint32_t median_typical_depth = median_depth_at_breadth(
        covered, depth, number_of_sites, median_breadth, typical_class_size
    );
    uint32_t median_top_decile_depth = median_depth_at_breadth(
        covered, depth, number_of_sites, top_decile_breadth, top_decile_class_size
    );

    uint64_t high_site = UINT64_MAX, top_decile_site = UINT64_MAX;
    uint64_t typical_site = UINT64_MAX;
    uint32_t high_depth = 0;
    uint64_t closest_typical_depth_difference = UINT64_MAX;
    uint64_t closest_top_decile_depth_difference = UINT64_MAX;
    for (uint64_t site = 0; site < number_of_sites; ++site) {
        if (covered[site] == maximum_breadth && depth[site] > high_depth) {
            high_site = site;
            high_depth = depth[site];
        }
        if (covered[site] == median_breadth) {
            uint64_t difference = depth[site] > median_typical_depth
                ? (uint64_t)depth[site] - median_typical_depth
                : (uint64_t)median_typical_depth - depth[site];
            if (difference < closest_typical_depth_difference) {
                closest_typical_depth_difference = difference;
                typical_site = site;
            }
        }
        if (covered[site] == top_decile_breadth) {
            uint64_t difference = depth[site] > median_top_decile_depth
                ? (uint64_t)depth[site] - median_top_decile_depth
                : (uint64_t)median_top_decile_depth - depth[site];
            if (difference < closest_top_decile_depth_difference) {
                closest_top_decile_depth_difference = difference;
                top_decile_site = site;
            }
        }
    }
    if (high_site == UINT64_MAX || top_decile_site == UINT64_MAX ||
        typical_site == UINT64_MAX)
        die("site selection failed");

    FILE *site_output = fopen(argv[3], "w");
    FILE *histogram_output = fopen(argv[4], "w");
    if (!site_output || !histogram_output) die("cannot open summary outputs");
    fprintf(site_output,
            "site_class\tordinal0\tcovered_pixels\ttotal_pixels\ttotal_depth\t"
            "aggregate_methylated\taggregate_beta\n");
    uint64_t selected[3] = {high_site, top_decile_site, typical_site};
    const char *labels[3] = {
        "broad_high_depth", "top_decile_representative", "typical_observed"
    };
    for (int k = 0; k < 3; ++k) {
        uint64_t site = selected[k];
        fprintf(site_output, "%s\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%u\t%.10f\n",
                labels[k], site, covered[site], selected_records, depth[site], methylated[site],
                depth[site] ? (double)methylated[site] / depth[site] : 0.0);
    }
    fclose(site_output);
    fprintf(histogram_output, "covered_pixels\tsite_count\n");
    for (uint64_t breadth = 0; breadth <= selected_records; ++breadth)
        if (histogram[breadth])
            fprintf(histogram_output, "%" PRIu64 "\t%" PRIu64 "\n",
                    breadth, histogram[breadth]);
    fclose(histogram_output);

    FILE *value_output = fopen(argv[2], "w");
    if (!value_output) die("cannot open pixel-value output");
    fprintf(value_output,
            "record\thigh_ordinal0\thigh_m\thigh_u\t"
            "top10_ordinal0\ttop10_m\ttop10_u\t"
            "typical_ordinal0\ttypical_m\ttypical_u\n");
    input = gzopen(argv[1], "rb");
    if (!input) die("cannot reopen CG file");
    uint64_t record = 0;
    while ((got = gzread(input, header, 17)) != 0) {
        if (got != 17 || le(header, 8) != 266563789635ULL || header[8] != '3')
            die("invalid YAME record on second pass");
        uint64_t bytes = le(header + 9, 8);
        unsigned char *payload = malloc(bytes ? bytes : 1);
        if (!payload || gzread(input, payload, (unsigned)bytes) != (int)bytes)
            die("truncated record on second pass");
        uint64_t position = 0;
        uint64_t selected_m[3] = {0, 0, 0}, selected_u[3] = {0, 0, 0};
        for (uint64_t i = 0; i < bytes;) {
            int type = payload[i] & 3;
            int width = type == 1 ? 1 : type == 3 ? 8 : 2;
            uint64_t value = le(payload + i, width), m = 0, u = 0;
            i += (uint64_t)width;
            if (type == 0) {
                position += value >> 2;
                continue;
            }
            decode_value(value, type, &m, &u);
            if (position == high_site) {
                selected_m[0] = m;
                selected_u[0] = u;
            }
            if (position == top_decile_site) {
                selected_m[1] = m;
                selected_u[1] = u;
            }
            if (position == typical_site) {
                selected_m[2] = m;
                selected_u[2] = u;
            }
            position++;
        }
        if (!use_filter || keep[record]) {
            fprintf(value_output,
                    "%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
                    "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
                    "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                record, high_site, selected_m[0], selected_u[0],
                    top_decile_site, selected_m[1], selected_u[1],
                    typical_site, selected_m[2], selected_u[2]);
        }
        free(payload);
        record++;
    }
    if (!gzeof(input)) die("second compressed-stream read failed");
    gzclose(input);
    fclose(value_output);

    free(histogram);
    free(covered);
    free(depth);
    free(methylated);
    free(keep);
    return 0;
}
