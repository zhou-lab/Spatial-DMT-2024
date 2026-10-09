/* Summarize YAME format-3 records without expanding uncovered CpGs.
 * Build: cc -O3 summarize_cg_stream.c -lz -o summarize_cg_stream
 * Output is in physical record order; callers must validate record names.
 * Format specification: https://github.com/zhou-lab/YAME (format3.c).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>

static uint64_t le(const unsigned char *p, int n) {
    uint64_t v = 0;
    for (int j = 0; j < n; j++) v |= (uint64_t)p[j] << (8*j);
    return v;
}
static void fail(const char *s) { fprintf(stderr, "%s\n", s); exit(1); }
int main(int argc, char **argv) {
    if (argc != 2) fail("Supply one YAME .cg path");
    gzFile f = gzopen(argv[1], "rb");
    if (!f) fail("Cannot open .cg");
    unsigned char h[17]; uint64_t record = 0;
    puts("record\tcovered_cpgs\tsite_mean_beta\tcall_weighted_beta\tmethylated_calls\ttotal_calls");
    int got;
    while ((got = gzread(f, h, 17)) != 0) {
        if (got != 17 || le(h,8) != 266563789635ULL || h[8] != '3')
            fail("Invalid or unsupported YAME record header");
        uint64_t n = le(h+9,8), i = 0, covered = 0, tm = 0, tc = 0;
        if (n > 1073741824ULL) fail("Record too large");
        unsigned char *p = malloc(n ? n : 1);
        if (!p || gzread(f,p,(unsigned)n) != (int)n) fail("Truncated record");
        double sum = 0;
        while (i < n) {
            int t = p[i]&3, width = t==1 ? 1 : t==3 ? 8 : 2;
            if (i + width > n) fail("Truncated format-3 value");
            uint64_t v = le(p+i,width), m=0, u=0; i += width;
            if (!t) continue;
            if (t==1) { m=v>>5; u=(v>>2)&7; }
            else if (t==2) { m=v>>9; u=(v>>2)&127; }
            else { m=v>>33; u=(v>>2)&2147483647ULL; }
            if (m+u) { covered++; sum += (double)m/(m+u); tm+=m; tc+=m+u; }
        }
        printf("%llu\t%llu\t%.10f\t%.10f\t%llu\t%llu\n",
            (unsigned long long)record++, (unsigned long long)covered,
            covered ? sum/covered : NAN, tc ? (double)tm/tc : NAN,
            (unsigned long long)tm, (unsigned long long)tc);
        free(p);
    }
    if (!gzeof(f)) fail("Compressed stream read failed");
    if (gzclose(f) != Z_OK) fail("Compressed stream close failed");
    return 0;
}
