#!/usr/bin/env python3
"""Generate a combined HTML diagnostic report from Spatial-DMT DNA and RNA QC."""

import argparse
import html
import os


DIAGNOSTIC_CSS = """
  .diagnostic-section {
    --ink: #172026;
    --muted: #5c6870;
    --line: #d8dee3;
    --paper: #ffffff;
    --panel-bg: #f5f7f8;
    --bc: #d9f0ff;
    --linker: #ffe6c7;
    --tn5: #eadcff;
    --genomic: #dff3df;
    --miss: #c83737;
    --wm: #2e6f9e;
    color: var(--ink);
  }
  .diagnostic-section h3 {
    margin: 0 0 12px;
    font-size: 18px;
    line-height: 1.15;
  }
  .diagnostic-panel {
    background: var(--paper);
    border: 1px solid var(--line);
    border-radius: 8px;
    padding: 18px;
    margin: 16px 0;
  }
  .diagnostic-section .muted {
    color: var(--muted);
  }
  .diagnostic-section .structure-grid {
    display: grid;
    grid-template-columns: 0.8fr 3fr 0.8fr 3fr 1.9fr 2fr;
    gap: 6px;
    align-items: stretch;
  }
  .diagnostic-section .dna-callout-row {
    display: grid;
    grid-template-columns: 0.8fr 3fr 0.8fr 3fr 1.9fr 2fr;
    gap: 6px;
    margin-top: 7px;
    align-items: center;
  }
  .diagnostic-section .seg {
    border: 1px solid rgba(0,0,0,0.18);
    border-radius: 6px;
    padding: 10px 8px;
    min-height: 48px;
    text-align: center;
    font-weight: 700;
  }
  .diagnostic-section .seg small {
    display: block;
    margin-top: 4px;
    font-weight: 500;
    color: var(--muted);
  }
  .diagnostic-section .bc { background: var(--bc); }
  .diagnostic-section .linker { background: var(--linker); }
  .diagnostic-section .tn5 { background: var(--tn5); }
  .diagnostic-section .genomic { background: var(--genomic); }
  .diagnostic-section .rna-structure-grid,
  .diagnostic-section .rna-callout-row {
    display: grid;
    grid-template-columns: 2.2fr 0.8fr 3fr 0.8fr 3fr 1fr 2fr;
    gap: 6px;
  }
  .diagnostic-section .rna-callout-row {
    margin-top: 7px;
    align-items: center;
  }
  .diagnostic-section .bracket {
    border: 2px solid var(--miss);
    border-top: 0;
    border-radius: 0 0 7px 7px;
    min-height: 18px;
    position: relative;
  }
  .diagnostic-section .bracket::after {
    content: attr(data-label);
    position: absolute;
    left: 50%;
    top: 20px;
    transform: translateX(-50%);
    white-space: nowrap;
    background: var(--paper);
    color: var(--miss);
    font-weight: 700;
    padding: 0 6px;
  }
  .diagnostic-section .wm-bracket {
    border-color: var(--wm);
  }
  .diagnostic-section .wm-bracket::after {
    color: var(--wm);
  }
  .diagnostic-section table {
    width: 100%;
    border-collapse: collapse;
    background: var(--paper);
  }
  .diagnostic-section th,
  .diagnostic-section td {
    border-bottom: 1px solid var(--line);
    padding: 9px 10px;
    text-align: left;
    vertical-align: top;
  }
  .diagnostic-section thead th {
    color: var(--muted);
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0;
  }
  .diagnostic-section tbody th {
    width: 220px;
  }
  .diagnostic-section tbody tr.sf-stage th {
    padding-left: 28px;
  }
  .diagnostic-section tbody tr.sf-stage th::before {
    content: "↳ ";
    color: var(--muted);
  }
  .diagnostic-section tbody tr.stage-2 th { padding-left: 44px; }
  .diagnostic-section tbody tr.stage-3 th { padding-left: 60px; }
  .diagnostic-section tbody tr.summary-row {
    background: var(--panel-bg);
  }
  .diagnostic-section tbody tr.summary-row th {
    font-weight: 700;
  }
  .diagnostic-section tbody tr.overlap-row {
    color: var(--muted);
  }
  .diagnostic-section tbody tr {
    cursor: help;
  }
  .diagnostic-section tbody tr:hover {
    background: #f7fbff;
  }
  @media (max-width: 800px) {
    .diagnostic-section table {
      table-layout: fixed;
      font-size: 11px;
    }
    .diagnostic-section th,
    .diagnostic-section td {
      padding: 6px 4px;
      overflow-wrap: anywhere;
    }
    .diagnostic-section tbody th {
      width: auto;
    }
    .diagnostic-section .structure-grid,
    .diagnostic-section .dna-callout-row,
    .diagnostic-section .rna-structure-grid,
    .diagnostic-section .rna-callout-row {
      grid-template-columns: 1fr;
    }
    .diagnostic-section .dna-callout-row .bracket,
    .diagnostic-section .rna-callout-row .bracket {
      grid-column: 1 !important;
    }
    .diagnostic-section .bracket::after {
      position: static;
      transform: none;
      display: block;
      text-align: center;
      margin-top: 18px;
    }
  }
"""


def parse_kv(path):
    data = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            key, value = parts
            try:
                data[key] = int(value)
            except ValueError:
                try:
                    data[key] = float(value)
                except ValueError:
                    data[key] = value
    return data


def parse_unmatched(path):
    data = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["section", "key", "value"]:
            raise ValueError(f"Unexpected unmatched diagnostics header in {path}")
        for line in handle:
            section, key, value = line.rstrip("\n").split("\t")
            try:
                parsed = int(value)
            except ValueError:
                try:
                    parsed = float(value)
                except ValueError:
                    parsed = value
            data[(section, key)] = parsed
    return data


def parse_bbduk(path):
    data = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#Total\t"):
                data["total_reads"] = int(line.rstrip("\n").split("\t")[1])
            elif line.startswith("#Matched\t"):
                data["matched_pairs"] = int(line.rstrip("\n").split("\t")[1])
    return data


def parse_barcode_stats(path):
    data = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 2:
                data[fields[0]] = int(fields[1])
    return data


def pct(num, denom, digits=1):
    if denom in (None, 0) or num is None:
        return "NA"
    return f"{100.0 * num / denom:.{digits}f}%"


def fmt_int(value):
    if value is None:
        return "NA"
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return html.escape(str(value))


def miss_description(name, count, total, sf_total=None):
    frac_total = (count / total) if total else 0.0
    frac_sf = (count / sf_total) if sf_total else 0.0
    if name == "linker1":
        base = "Expected linker1 after 8 bp BC1 was not uniquely found within <=2 edits."
        if frac_sf > 0.30:
            return base + " Large first-anchor loss; check BC1/linker1 chemistry and R2 cycle quality."
        return base + " Linker1 loss is present but not dominant among structure failures."
    if name == "linker2":
        base = "Expected linker2 after 8 bp BC2 was not uniquely found within <=2 edits."
        if frac_sf > 0.10:
            return base + " BC2/linker2 loss is noticeable; this can indicate indels around BC2."
        return base + " Linker2-specific failure is relatively small."
    if name == "tn5":
        base = "Tn5 mosaic-end plus flank after linker2 was not uniquely found within <=1 edit."
        if frac_sf > 0.20:
            return base + " Tn5-anchor loss is substantial; short inserts or readthrough are plausible."
        return base + " Tn5-anchor failure is present but not the main miss class."
    if name == "wm":
        base = "BC1+BC2 sequence exists structurally but was not in the whitelist after Hamming-1 rescue."
        if frac_total > 0.05:
            return base + " Whitelist misses are non-trivial; inspect barcode quality or possible Hamming-2 rescue."
        return base + " Whitelist misses are modest relative to total reads."
    return ""


def row(label, requirement, value, frac_failure, frac_total, description, class_name=""):
    class_attr = f' class="{html.escape(class_name, quote=True)}"' if class_name else ""
    tooltip = f"{requirement}. {description}" if requirement else description
    return (
        f'<tr{class_attr} title="{html.escape(tooltip, quote=True)}">'
        f"<th>{html.escape(label)}</th>"
        f"<td>{fmt_int(value)}</td>"
        f"<td>{html.escape(frac_failure)}</td>"
        f"<td>{html.escape(frac_total)}</td>"
        "</tr>"
    )


def genomic_total(trim):
    """DNA read pairs excluding the barcode85/lambda bucket."""
    if trim.get("Reads_detected") is not None:
        return trim["Reads_detected"]
    total = trim.get("Reads_total", 0)
    for key in ("Reads_barcode85", "Reads_written_barcode85"):
        if trim.get(key) is not None:
            return total - trim[key]
    return total


def build_sections(sample, trim, unmatched, rna=None, protocol="spatialdmt"):
    total = genomic_total(trim)
    wl_unmatched = trim.get("WL_unmatched", 0)
    sf_reads = trim.get(
        "Unmatched_SF",
        unmatched.get(
            ("fastq_stats", "structure_fail_n_reads"),
            unmatched.get(("align_stats", "structure_fail_n_reads"), 0),
        ),
    )
    wm_reads = trim.get(
        "Unmatched_WM",
        unmatched.get(
            ("fastq_stats", "wl_miss_n_reads"),
            unmatched.get(("align_stats", "wl_miss_n_reads"), 0),
        ),
    )

    sf_linker1 = unmatched.get(("sf_anchor_fail", "linker1_fail_n"), 0)
    sf_linker2 = unmatched.get(("sf_anchor_fail", "linker2_fail_n"), 0)
    sf_tn5 = unmatched.get(("sf_anchor_fail", "tn5_me_fail_n"), 0)
    sf_all_three_absent = unmatched.get(("sf_anchor_overlap", "all_three_absent_n"), 0)
    wm_count = wl_unmatched or wm_reads
    failure_reads = sf_reads + wm_count

    miss_rows = "\n".join([
        row("Total read pairs", "", total, "—", "100.0%",
            "Genomic DNA read pairs entering trim/tag, excluding the barcode85/lambda bucket.", "summary-row"),
        row("Total missed pairs", "Structure failures plus barcode whitelist misses", failure_reads,
            "100.0%", pct(failure_reads, total),
            "All DNA read pairs that did not enter a designed spatial barcode.", "summary-row"),
        row("SF pairs", "Any structural stage failed", sf_reads, pct(sf_reads, failure_reads), pct(sf_reads, total),
            "All structure-fail reads. Each read appears under exactly one first-failure stage below."),
        row("Stage 1: linker1 miss", "Start; linker1 expected after BC1", sf_linker1,
            pct(sf_linker1, failure_reads), pct(sf_linker1, total),
            miss_description("linker1", sf_linker1, total, sf_reads), "sf-stage stage-1"),
        row("Stage 2: linker2 miss", "Reached after linker1 is found; BC1 is extractable", sf_linker2,
            pct(sf_linker2, failure_reads), pct(sf_linker2, total),
            miss_description("linker2", sf_linker2, total, sf_reads), "sf-stage stage-2"),
        row("Stage 3: Tn5 ME miss", "Reached after linker1/linker2 are found; BC1/BC2 are extractable", sf_tn5,
            pct(sf_tn5, failure_reads), pct(sf_tn5, total),
            miss_description("tn5", sf_tn5, total, sf_reads), "sf-stage stage-3"),
        row("All three anchors absent", "Independent anchor search; overlapping subset of Stage 1", sf_all_three_absent,
            pct(sf_all_three_absent, failure_reads), pct(sf_all_three_absent, total),
            "Linker1, linker2, and Tn5 ME each had zero fuzzy matches. This row overlaps Stage 1 and must not be added to the first-failure rows.",
            "sf-stage stage-1 overlap-row"),
        row("WM: barcode whitelist miss", "All structure stages passed", wm_count,
            pct(wm_count, failure_reads), pct(wm_count, total),
            miss_description("wm", wm_count, total)),
    ])

    rna_html = ""
    if rna:
        if protocol == "smcseq":
            rna_stage1 = "Published motif1 at R2 position 2, allowing <=1 substitution."
            rna_stage2 = "Published 30 bp motif2 after BC1, allowing <=1 edit."
            rna_stage3 = (
                "Published 49 bp motif3 after BC2, allowing <=2 substitutions "
                "and <=1 insertion/deletion."
            )
            rna_structure = """
      <div class="seg linker">1 bp + motif1<small>discarded + 14 bp</small></div>
      <div class="seg bc">BC1<small>8 bp</small></div>
      <div class="seg linker">motif2<small>30 bp, fuzzy</small></div>
      <div class="seg bc">BC2<small>8 bp</small></div>
      <div class="seg linker">motif3<small>49 bp, fuzzy</small></div>
      <div class="seg tn5">UMI<small>10 bp</small></div>
      <div class="seg genomic">remaining R2<small>discarded</small></div>
"""
        else:
            rna_stage1 = "Expected 22 bp primer in R2 within Hamming distance 2."
            rna_stage2 = "Expected 30 bp linker1 in R2 within Hamming distance 3."
            rna_stage3 = "Expected 30 bp linker2 in R2 within Hamming distance 3."
            rna_structure = """
      <div class="seg linker">primer<small>22 bp, Hamming &lt;=2</small></div>
      <div class="seg bc">BC2<small>8 bp</small></div>
      <div class="seg linker">linker1<small>30 bp, Hamming &lt;=3</small></div>
      <div class="seg bc">BC1<small>8 bp</small></div>
      <div class="seg linker">linker2<small>30 bp, Hamming &lt;=3</small></div>
      <div class="seg tn5">UMI<small>10 bp</small></div>
      <div class="seg genomic">remaining R2<small>not used by STARsolo</small></div>
"""
        primer, linker1, linker2, barcodes = rna
        rna_total = primer["total_reads"] // 2
        primer_pass = primer["matched_pairs"]
        linker1_pass = linker1["matched_pairs"]
        linker2_pass = linker2["matched_pairs"]
        primer_miss = max(rna_total - primer_pass, 0)
        linker1_miss = max(primer_pass - linker1_pass, 0)
        linker2_miss = max(linker1_pass - linker2_pass, 0)
        structure_fail = primer_miss + linker1_miss + linker2_miss
        wl_miss = barcodes.get("noNoWLmatch", 0)
        rna_failure = structure_fail + wl_miss

        rna_miss_rows = "\n".join([
            row("Total read pairs", "", rna_total, "—", "100.0%",
                "All raw RNA read pairs entering sequential structure filtering.", "summary-row"),
            row("Total missed pairs", "Structure failures plus barcode whitelist misses", rna_failure,
                "100.0%", pct(rna_failure, rna_total),
                "All RNA read pairs removed during structure filtering or lacking a STARsolo whitelist match.",
                "summary-row"),
            row("RNA: structure fail", "Any sequential RNA structure stage failed", structure_fail,
                pct(structure_fail, rna_failure), pct(structure_fail, rna_total),
                "All RNA read pairs removed by sequential primer, linker1, or linker2 filtering.",
                "summary-row"),
            row("Stage 1: primer miss", "Start", primer_miss, pct(primer_miss, rna_failure), pct(primer_miss, rna_total),
                rna_stage1, "sf-stage stage-1"),
            row("Stage 2: linker1 miss", "Reached only after primer passes", linker1_miss,
                pct(linker1_miss, rna_failure), pct(linker1_miss, rna_total),
                rna_stage2, "sf-stage stage-2"),
            row("Stage 3: linker2 miss", "Reached only after primer and linker1 pass", linker2_miss,
                pct(linker2_miss, rna_failure), pct(linker2_miss, rna_total),
                rna_stage3, "sf-stage stage-3"),
            row("RNA: barcode whitelist miss", "All RNA structure stages passed", wl_miss,
                pct(wl_miss, rna_failure), pct(wl_miss, rna_total),
                "Structure-passing RNA read pairs whose cell barcode had no STARsolo whitelist match."),
        ])

        rna_html = f"""
  <section class="diagnostic-panel">
    <h3>RNA R2 Structure</h3>
    <div class="rna-structure-grid" aria-label="RNA R2 structure">
{rna_structure}
    </div>
    <div class="rna-callout-row"><div class="bracket" style="grid-column: 1 / 2;" data-label="Primer miss"></div></div>
    <div class="rna-callout-row"><div class="bracket" style="grid-column: 3 / 4;" data-label="Linker1 miss"></div></div>
    <div class="rna-callout-row"><div class="bracket" style="grid-column: 5 / 6;" data-label="Linker2 miss"></div></div>
    <div class="rna-callout-row"><div class="bracket wm-bracket" style="grid-column: 2 / 5;" data-label="WL: BC1 + BC2 not in whitelist"></div></div>
  </section>

  <section class="diagnostic-panel">
    <h3>RNA Failure Breakdown</h3>
    <table>
      <thead><tr>
        <th>Failure class</th>
        <th title="RNA read-pair count for this miss class.">Count</th>
        <th title="Fraction of all RNA failures: structure failures plus barcode whitelist misses.">Fraction of total failures</th>
        <th title="Fraction of all raw RNA read pairs entering primer filtering for this sample.">Fraction of total</th>
      </tr></thead>
      <tbody>{rna_miss_rows}</tbody>
    </table>
  </section>

"""

    return f"""
<section class="diagnostic-section">
  <h2>DNA/RNA Structure Diagnostic</h2>

  <section class="diagnostic-panel">
    <h3>DNA R2 Structure</h3>
    <div class="structure-grid" aria-label="R2 structure">
      <div class="seg bc">BC1<small>8 bp</small></div>
      <div class="seg linker">linker1<small>30 bp, fuzzy &lt;=2</small></div>
      <div class="seg bc">BC2<small>8 bp</small></div>
      <div class="seg linker">linker2<small>30 bp, fuzzy &lt;=2</small></div>
      <div class="seg tn5">Tn5 ME + flank<small>19 bp, fuzzy &lt;=1</small></div>
      <div class="seg genomic">genomic<small>alignment insert</small></div>
    </div>
    <div class="dna-callout-row"><div class="bracket" style="grid-column: 2 / 3;" data-label="Linker1 miss"></div></div>
    <div class="dna-callout-row"><div class="bracket" style="grid-column: 4 / 5;" data-label="Linker2 miss"></div></div>
    <div class="dna-callout-row"><div class="bracket" style="grid-column: 5 / 6;" data-label="Tn5 ME miss"></div></div>
    <div class="dna-callout-row"><div class="bracket wm-bracket" style="grid-column: 1 / 4;" data-label="WL: BC1 + BC2 not in whitelist"></div></div>
  </section>

  <section class="diagnostic-panel">
    <h3>DNA Failure Breakdown</h3>
    <table>
      <thead><tr>
        <th>Failure class</th>
        <th title="Read-pair count for this miss class.">Count</th>
        <th title="Fraction of all DNA failures: structure failures plus barcode whitelist misses.">Fraction of total failures</th>
        <th title="Fraction of genomic DNA read pairs entering trim/tag, excluding the barcode85/lambda bucket.">Fraction of total</th>
      </tr></thead>
      <tbody>
        {miss_rows}
      </tbody>
    </table>
  </section>
  {rna_html}
</section>
"""


def build_html(sample, trim, unmatched, rna=None, protocol="spatialdmt"):
    sections = build_sections(sample, trim, unmatched, rna, protocol)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(sample)} DNA and RNA Diagnostic</title>
<style>
  body {{
    margin: 0;
    background: #f5f7f8;
    color: #172026;
    font: 14px/1.45 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  }}
  main {{
    max-width: 1180px;
    margin: 0 auto;
    padding: 28px;
  }}
  h1 {{
    margin: 0 0 18px;
    font-size: 26px;
    line-height: 1.15;
  }}
{DIAGNOSTIC_CSS}
</style>
</head>
<body>
<main>
  <h1>{html.escape(sample)} DNA and RNA Diagnostic</h1>
  {sections}
</main>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--protocol", choices=["spatialdmt", "smcseq"], default="spatialdmt")
    parser.add_argument("--trim-stats", required=True)
    parser.add_argument("--unmatched-diagnostics", required=True)
    parser.add_argument("--rna-primer-stats")
    parser.add_argument("--rna-linker1-stats")
    parser.add_argument("--rna-linker2-stats")
    parser.add_argument("--rna-barcode-stats")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    trim = parse_kv(args.trim_stats)
    unmatched = parse_unmatched(args.unmatched_diagnostics)
    rna_paths = [
        args.rna_primer_stats,
        args.rna_linker1_stats,
        args.rna_linker2_stats,
        args.rna_barcode_stats,
    ]
    if any(rna_paths) and not all(rna_paths):
        parser.error("All RNA input arguments must be supplied together")
    rna = None
    if all(rna_paths):
        rna = (
            parse_bbduk(args.rna_primer_stats),
            parse_bbduk(args.rna_linker1_stats),
            parse_bbduk(args.rna_linker2_stats),
            parse_barcode_stats(args.rna_barcode_stats),
        )
    html_text = build_html(args.sample, trim, unmatched, rna, args.protocol)

    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write(html_text)
    print(args.output)


if __name__ == "__main__":
    main()
