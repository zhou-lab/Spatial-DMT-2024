#define _GNU_SOURCE

/*
 * Single-sample context-aware implementation of the merge behavior in
 * BISCUIT's MIT-licensed src/mergecg.c.
 *
 * Copyright (c) 2016-2020 Wanding Zhou
 * Copyright (c) 2021-2026 Jacob Morrison
 */

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char *chrom;
    int64_t beg;
    int64_t end;
    char ref;
    double beta;
    int64_t depth;
    int present;
} meth_record_t;

static void usage(FILE *stream)
{
    fprintf(stream, "Usage: mergecg_context [CONTEXT_BED|-]\n");
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

static int parse_double(const char *text, double *value)
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

static int assign_record(meth_record_t *record, char **fields, int line_number)
{
    int64_t beg = 0;
    int64_t end = 0;
    int64_t depth = 0;
    double beta = 0.0;

    if ((fields[3][0] != 'C' && fields[3][0] != 'G') || fields[3][1] != '\0'
            || parse_i64(fields[1], &beg) != 0
            || parse_i64(fields[2], &end) != 0
            || parse_double(fields[7], &beta) != 0
            || parse_i64(fields[8], &depth) != 0
            || beg < 0 || end <= beg || depth < 1) {
        fprintf(stderr, "Malformed context BED record at line %d.\n", line_number);
        return -1;
    }

    char *chrom = strdup(fields[0]);
    if (chrom == NULL) {
        perror("strdup");
        return -1;
    }
    free(record->chrom);
    record->chrom = chrom;
    record->beg = beg;
    record->end = end;
    record->ref = fields[3][0];
    record->beta = beta;
    record->depth = depth;
    record->present = 1;
    return 0;
}

static int report_record(const meth_record_t *record,
                         double g_beta, int64_t g_depth)
{
    double c_beta = 0.0;
    int64_t c_depth = 0;
    int64_t beg = record->beg;
    int64_t end = record->end;

    if (record->ref == 'C') {
        c_beta = record->beta;
        c_depth = record->depth;
        end++;
    } else {
        g_beta = record->beta;
        g_depth = record->depth;
        if (beg == 0) {
            fprintf(stderr, "Cannot extend a G-strand call before position zero.\n");
            return -1;
        }
        beg--;
    }

    int64_t coverage = c_depth + g_depth;
    if (coverage < 1) {
        return 0;
    }

    float c_retained = rintf(c_beta * c_depth);
    float g_retained = rintf(g_beta * g_depth);
    int64_t retained = (int64_t)c_retained + (int64_t)g_retained;
    double merged_beta = retained / (double)coverage;

    if (printf("%s\t%" PRId64 "\t%" PRId64 "\t%1.3f\t%" PRId64,
               record->chrom, beg, end, merged_beta, coverage) < 0) {
        return -1;
    }
    if (c_depth == 0) {
        fputs("\tC:.:0", stdout);
    } else {
        printf("\tC:%1.3f:%" PRId64, c_beta, c_depth);
    }
    if (g_depth == 0) {
        fputs(",G:.:0\n", stdout);
    } else {
        printf(",G:%1.3f:%" PRId64 "\n", g_beta, g_depth);
    }
    return ferror(stdout) ? -1 : 0;
}

int main(int argc, char **argv)
{
    if (argc > 2 || (argc == 2 && strcmp(argv[1], "-h") == 0)) {
        usage(argc > 2 ? stderr : stdout);
        return argc > 2 ? 2 : 0;
    }

    FILE *input = stdin;
    if (argc == 2 && strcmp(argv[1], "-") != 0) {
        input = fopen(argv[1], "r");
        if (input == NULL) {
            perror(argv[1]);
            return 1;
        }
    }

    setvbuf(stdout, NULL, _IOFBF, 1U << 20);
    char *line = NULL;
    size_t line_capacity = 0;
    ssize_t line_length;
    int line_number = 0;
    int status = 0;
    meth_record_t previous = {0};

    while ((line_length = getline(&line, &line_capacity, input)) >= 0) {
        line_number++;
        while (line_length > 0
               && (line[line_length - 1] == '\n' || line[line_length - 1] == '\r')) {
            line[--line_length] = '\0';
        }
        if (line_length == 0) {
            continue;
        }

        char *fields[10] = {0};
        int field_count = split_fields(line, fields, 10);
        if (field_count < 9) {
            fprintf(stderr, "Expected nine context BED fields at line %d.\n", line_number);
            status = 1;
            break;
        }
        if (strcmp(fields[4], "CG") != 0) {
            continue;
        }

        meth_record_t current = {0};
        if (assign_record(&current, fields, line_number) != 0) {
            status = 1;
            break;
        }

        int is_pair = previous.present
            && previous.ref == 'C'
            && current.ref == 'G'
            && strcmp(previous.chrom, current.chrom) == 0
            && current.beg == previous.beg + 1
            && current.end == previous.end + 1;
        if (is_pair) {
            if (report_record(&previous, current.beta, current.depth) != 0) {
                status = 1;
            }
            free(previous.chrom);
            previous = (meth_record_t){0};
            free(current.chrom);
            if (status != 0) {
                break;
            }
        } else {
            if (previous.present && report_record(&previous, 0.0, 0) != 0) {
                free(current.chrom);
                status = 1;
                break;
            }
            free(previous.chrom);
            previous = current;
        }
    }

    if (status == 0 && ferror(input)) {
        perror("context BED input");
        status = 1;
    }
    if (status == 0 && previous.present
            && report_record(&previous, 0.0, 0) != 0) {
        status = 1;
    }
    if (fflush(stdout) != 0) {
        perror("stdout");
        status = 1;
    }

    free(previous.chrom);
    free(line);
    if (input != stdin && fclose(input) != 0) {
        perror("context BED input");
        status = 1;
    }
    return status;
}
