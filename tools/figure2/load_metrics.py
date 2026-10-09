"""Discovery, harmonization, and provenance for Figure 2."""

from __future__ import annotations

import csv
import logging
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class SampleRecord:
    sample: str
    species: str
    tissue: str
    modality: str
    revised_root: Path | None
    legacy_root: Path | None

    @property
    def context(self) -> str:
        return f"{self.species} {self.tissue}"


DEFAULT_SAMPLES = (
    SampleRecord(
        "SpMETSLE14D", "mouse", "embryo", "DNA",
        Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
        Path("/scr1/users/suc3/tmp"),
    ),
    SampleRecord(
        "SpMETSLE14DM", "mouse", "embryo", "DNA",
        Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
        Path("/scr1/users/suc3/tmp"),
    ),
    SampleRecord(
        "SpMETSLE16D", "mouse", "embryo", "DNA+RNA",
        Path("/scr1/users/suc3/20260528_runNewSpatialPipe"),
        Path("/scr1/users/suc3/tmp"),
    ),
    SampleRecord(
        "SpMETSLB16D", "mouse", "brain", "DNA+RNA",
        Path("/scr1/users/suc3/20260704_runSpMET_latest_mm10"),
        Path("/scr1/users/suc3/tmp"),
    ),
    SampleRecord(
        "SpMETSLHSK03D", "human", "skin", "DNA",
        Path("/scr1/users/suc3/20260624_runNewSpatialPipe_hg38"),
        Path("/scr1/users/suc3/tmp"),
    ),
    SampleRecord(
        "SpMETSLHSK05D", "human", "skin", "DNA",
        Path("/scr1/users/suc3/20260706_runSpMETSLHSK_more_hg38"),
        Path("/scr1/users/suc3/tmp"),
    ),
)


METRIC_META = {
    "raw_read_pairs": ("Raw pairs", "read pairs", "higher"),
    "pass_fraction": ("Passed", "% input", "higher"),
    "mapping_rate": ("Mapped", "% spatial fragments", "higher"),
    "unique_fragment_fraction": ("Unique", "% mapped", "higher"),
    "duplication_rate": ("Duplicates", "% mapped", "lower"),
    "conversion_efficiency": ("Conversion", "%", "higher"),
    "unique_cpgs": ("Unique CpGs", "sites", "higher"),
    "median_cpgs_per_pixel": ("CpGs/pixel", "sites", "higher"),
    "covered_pixel_fraction": ("Covered pixels", "% designed", "higher"),
    "median_unique_fragments": ("Median fragments", "fragments/pixel", "higher"),
    "rna_unique_mapping_rate": ("RNA mapped", "% reads", "higher"),
    "median_rna_umis": ("RNA UMIs", "UMIs/pixel", "higher"),
    "median_rna_genes": ("RNA genes", "genes/pixel", "higher"),
}


def load_sample_config(path: Path | None) -> list[SampleRecord]:
    if path is None:
        return list(DEFAULT_SAMPLES)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(path, sep=sep)
    required = {"sample", "species", "tissue", "modality", "revised_root", "legacy_root"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Sample config lacks columns: {sorted(missing)}")
    records = []
    for row in frame.to_dict("records"):
        records.append(
            SampleRecord(
                row["sample"], row["species"], row["tissue"], row["modality"],
                Path(row["revised_root"]) if pd.notna(row["revised_root"]) else None,
                Path(row["legacy_root"]) if pd.notna(row["legacy_root"]) else None,
            )
        )
    return records


def _first(paths: Iterable[Path]) -> Path | None:
    return next((p for p in paths if p.is_file()), None)


def revised_sample_dir(record: SampleRecord) -> Path | None:
    if record.revised_root is None:
        return None
    path = record.revised_root / "pipeline_output" / record.sample
    return path if path.is_dir() else None


def find_curated_multiqc_metrics(record: SampleRecord) -> Path | None:
    if record.revised_root is None:
        return None
    return _first(
        (
            record.revised_root / "final_output" / record.sample
            / "multiqc_report_data" / "multiqc_curated_metrics.txt",
            record.revised_root / "pipeline_output" / record.sample / "qc"
            / "biscuit_qc_merged" / f"{record.sample}_curated_metrics_mqc.tsv",
        )
    )


def find_trim_stats(record: SampleRecord) -> tuple[Path | None, str, bool]:
    sample_dir = revised_sample_dir(record)
    if sample_dir:
        found = _first(sorted((sample_dir / "trim").glob("*trim_stats*")))
        if found:
            return found, "revised", False
    if record.legacy_root:
        found = _first(sorted((record.legacy_root / "trim" / record.sample).glob("*_stats.txt")))
        if found:
            LOG.warning("Using legacy trim stats for %s: %s", record.sample, found)
            return found, "legacy_fallback", True
    return None, "missing", False


def find_spatial_barcode_counts(record: SampleRecord) -> tuple[Path | None, str, bool]:
    sample_dir = revised_sample_dir(record)
    if sample_dir:
        found = _first(sorted((sample_dir / "qc" / "table").glob("*_spatial_barcode_counts.tsv")))
        if found:
            return found, "revised", False
    if record.legacy_root:
        found = _first(sorted((record.legacy_root / "trim" / record.sample).glob("*_spatial_barcode_counts.tsv")))
        if found:
            LOG.warning("Using legacy spatial counts for %s: %s", record.sample, found)
            return found, "legacy_fallback", True
    return None, "missing", False


def find_unmatched_diagnostics(record: SampleRecord) -> Path | None:
    sample_dir = revised_sample_dir(record)
    path = sample_dir / "qc" / "unmatched_diagnostics.tsv" if sample_dir else None
    return path if path and path.is_file() else None


def find_cg_file(record: SampleRecord) -> tuple[Path | None, str, bool]:
    if record.revised_root:
        candidates = (
            record.revised_root / "final_output" / record.sample / f"{record.sample}.cg",
            record.revised_root / "pipeline_output" / record.sample / "pileup" / f"{record.sample}.cg",
        )
        found = _first(candidates)
        if found:
            return found, "recomputed_revised", False
    return None, "missing", False


def find_rna_matrix(record: SampleRecord) -> Path | None:
    sample_dir = revised_sample_dir(record)
    if not sample_dir:
        return None
    return _first(
        (
            sample_dir / "rna_STAR" / "Solo.out" / "GeneFull" / "raw" / "matrix.mtx",
            record.revised_root / "final_output" / record.sample / "Solo.out"
            / "GeneFull" / "raw" / "matrix.mtx",
        )
    )


def find_legacy_metric(record: SampleRecord, pattern: str) -> Path | None:
    if record.legacy_root is None:
        return None
    return _first(sorted(record.legacy_root.glob(pattern)))


def read_kv(path: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                try:
                    result[fields[0]] = float(fields[1])
                except ValueError:
                    continue
    return result


def read_custom_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def read_star_summary(path: Path) -> dict[str, float]:
    result = {}
    with path.open(newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2:
                try:
                    result[row[0]] = float(row[1])
                except ValueError:
                    continue
    return result


def read_biscuit_conversion(path: Path) -> float:
    lines = [line.rstrip("\n") for line in path.open() if line.strip()]
    columns = lines[-2].split("\t")
    values = lines[-1].split("\t")
    parsed = dict(zip(columns, values))
    residual = []
    for base in ("CA", "CC", "CT"):
        try:
            value = float(parsed[base])
        except (KeyError, ValueError):
            continue
        if value >= 0:
            residual.append(value)
    return 100 * (1 - float(np.mean(residual))) if residual else np.nan


def _metric_row(
    record: SampleRecord,
    metric: str,
    value: float,
    source: Path,
    source_type: str,
    definition_version: str,
    is_fallback: bool = False,
    panel: str = "A-C",
) -> dict:
    unit = METRIC_META.get(metric, (metric, "", "higher"))[1]
    return {
        "sample": record.sample,
        "species": record.species,
        "tissue": record.tissue,
        "modality": record.modality,
        "pipeline": "legacy" if "legacy" in source_type else "revised",
        "metric": metric,
        "value": float(value) if value is not None else np.nan,
        "unit": unit,
        "source_file": str(source),
        "source_type": source_type,
        "definition_version": definition_version,
        "is_fallback": bool(is_fallback),
        "panel": panel,
    }


def run_yame_summary(
    record: SampleRecord,
    cache_dir: Path,
    force: bool = False,
    skip_expensive: bool = False,
) -> tuple[pd.DataFrame | None, Path | None, str, bool]:
    cg, source_type, fallback = find_cg_file(record)
    if cg is None:
        return None, None, "missing", False
    cache = cache_dir / f"{record.sample}_yame_summary.tsv"
    if cache.is_file() and not force:
        return pd.read_csv(cache, sep="\t"), cg, source_type, fallback
    if skip_expensive:
        LOG.warning("Skipping uncached yame summary for %s", record.sample)
        return None, cg, source_type, fallback
    LOG.info("Computing per-pixel CpG summary for %s", record.sample)
    try:
        with cache.open("w") as out:
            subprocess.run(["yame", "summary", str(cg)], stdout=out, check=True)
    except BaseException:
        cache.unlink(missing_ok=True)
        raise
    return pd.read_csv(cache, sep="\t"), cg, source_type, fallback


def run_yame_aggregate(
    record: SampleRecord,
    cache_dir: Path,
    force: bool = False,
    skip_expensive: bool = False,
) -> tuple[pd.DataFrame | None, Path | None, str, bool]:
    """Sum indexed cells and summarize the union of covered CpGs."""
    cg, source_type, fallback = find_cg_file(record)
    if cg is None:
        return None, None, "missing", False
    cache = cache_dir / f"{record.sample}_yame_aggregate.tsv"
    if cache.is_file() and not force:
        return pd.read_csv(cache, sep="\t"), cg, source_type, fallback
    if skip_expensive:
        LOG.warning("Skipping uncached aggregate yame summary for %s", record.sample)
        return None, cg, source_type, fallback
    LOG.info("Computing aggregate CpG coverage for %s", record.sample)
    rowop = subprocess.Popen(
        ["yame", "rowop", "-o", "musum", str(cg)],
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        with cache.open("w") as out:
            subprocess.run(
                ["yame", "summary", "-"],
                stdin=rowop.stdout,
                stdout=out,
                check=True,
                text=True,
            )
        if rowop.stdout:
            rowop.stdout.close()
        if rowop.wait() != 0:
            raise subprocess.CalledProcessError(rowop.returncode, rowop.args)
    except BaseException:
        rowop.kill()
        cache.unlink(missing_ok=True)
        raise
    return pd.read_csv(cache, sep="\t"), cg, source_type, fallback


def collect_sample_metrics(
    record: SampleRecord,
    cache_dir: Path,
    force: bool = False,
    skip_expensive: bool = False,
    include_yame: bool = False,
) -> tuple[list[dict], pd.DataFrame | None]:
    rows: list[dict] = []
    trim_path, trim_type, trim_fallback = find_trim_stats(record)
    trim = read_kv(trim_path) if trim_path else {}
    raw = trim.get("Reads_total", trim.get("Reads_detected", np.nan))
    passed = trim.get("Reads_passed", np.nan)
    if trim_path:
        rows.append(_metric_row(record, "raw_read_pairs", raw, trim_path, trim_type, "trim_stats_v2", trim_fallback))
        rows.append(_metric_row(record, "pass_fraction", 100 * passed / raw if raw else np.nan, trim_path, trim_type, "Reads_passed/Reads_total_or_detected", trim_fallback))

    spatial_path, spatial_type, spatial_fallback = find_spatial_barcode_counts(record)
    spatial = read_custom_tsv(spatial_path) if spatial_path else None
    if spatial is not None and {"dna_reads", "mapped_reads", "dup_reads"} <= set(spatial):
        dna = pd.to_numeric(spatial["dna_reads"], errors="coerce")
        mapped = pd.to_numeric(spatial["mapped_reads"], errors="coerce")
        dup = pd.to_numeric(spatial["dup_reads"], errors="coerce")
        total_dna, total_mapped, total_dup = dna.sum(), mapped.sum(), dup.sum()
        covered = dna > 0
        values = {
            "mapping_rate": 100 * total_mapped / total_dna if total_dna else np.nan,
            "duplication_rate": 100 * total_dup / total_mapped if total_mapped else np.nan,
            "unique_fragment_fraction": 100 * (total_mapped - total_dup) / total_mapped if total_mapped else np.nan,
            "covered_pixel_fraction": 100 * covered.mean(),
            "median_unique_fragments": float((mapped - dup)[covered].median()) if covered.any() else np.nan,
        }
        for metric, value in values.items():
            rows.append(_metric_row(record, metric, value, spatial_path, spatial_type, "revised_spatial_bam_fragments_v1", spatial_fallback))

    curated_path = find_curated_multiqc_metrics(record)
    if curated_path:
        curated = read_custom_tsv(curated_path)
        if not curated.empty:
            one = curated.iloc[0]
            if "cpg_coverage_yame" in one:
                rows.append(_metric_row(record, "unique_cpgs", float(one["cpg_coverage_yame"]) * 1e6, curated_path, "revised", "yame_aggregate_N_query_v1"))
            contexts = [one.get(f"{b}_retention_biscuitqc", np.nan) for b in ("CA", "CC", "CT")]
            if np.isfinite(contexts).any():
                rows.append(_metric_row(record, "conversion_efficiency", 100 * (1 - np.nanmean(contexts)), curated_path, "revised", "1-mean_CA_CC_CT_retention_v1"))

    if not any(row["metric"] == "conversion_efficiency" for row in rows):
        sample_dir = revised_sample_dir(record)
        conv = _first(sorted((sample_dir / "qc" / "biscuit_qc").glob("*totalBaseConversionRate.txt"))) if sample_dir else None
        if conv:
            rows.append(
                _metric_row(
                    record,
                    "conversion_efficiency",
                    read_biscuit_conversion(conv),
                    conv,
                    "recomputed_revised",
                    "1-mean_CA_CC_CT_retention_v1",
                )
            )

    sample_dir = revised_sample_dir(record)
    if sample_dir:
        rna_table = _first(sorted((sample_dir / "qc" / "table").glob("*_rna_spatial_counts.tsv")))
        summary = sample_dir / "rna_STAR" / "Solo.out" / "GeneFull" / "Summary.csv"
        if rna_table:
            rna = pd.read_csv(rna_table, sep="\t")
            covered = pd.to_numeric(rna["rna_umis"], errors="coerce") > 0
            for metric, col in (("median_rna_umis", "rna_umis"), ("median_rna_genes", "rna_genes")):
                value = float(pd.to_numeric(rna.loc[covered, col], errors="coerce").median()) if covered.any() else np.nan
                rows.append(_metric_row(record, metric, value, rna_table, "revised", "STARsolo_raw_covered_pixels_v1"))
        if summary.is_file():
            star = read_star_summary(summary)
            value = star.get("Reads Mapped to Genome: Unique", np.nan)
            rows.append(_metric_row(record, "rna_unique_mapping_rate", 100 * value, summary, "revised", "STARsolo_GeneFull_Summary_v1"))

    yame = None
    if include_yame:
        yame, cg_path, cg_type, cg_fallback = run_yame_summary(record, cache_dir, force, skip_expensive)
        if yame is not None and not yame.empty:
            n_query = pd.to_numeric(yame["N_query"], errors="coerce")
            rows.append(_metric_row(record, "median_cpgs_per_pixel", float(n_query[n_query > 0].median()), cg_path, cg_type, "yame_per_index_N_query_v1", cg_fallback))
            if not any(row["metric"] == "unique_cpgs" for row in rows):
                aggregate, agg_path, agg_type, agg_fallback = run_yame_aggregate(
                    record, cache_dir, force, skip_expensive
                )
                if aggregate is not None and not aggregate.empty:
                    rows.append(
                        _metric_row(
                            record,
                            "unique_cpgs",
                            float(aggregate.iloc[0]["N_query"]),
                            agg_path,
                            agg_type,
                            "yame_musum_union_N_query_v1",
                            agg_fallback,
                        )
                    )
                else:
                    LOG.warning("Aggregate unique CpGs unavailable for %s", record.sample)
    return rows, yame


def read_legacy_whitelist(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None)
    return pd.DataFrame({"barcode": frame.iloc[:, 3], "x": frame.iloc[:, 1], "y": frame.iloc[:, 2]})


def _count_legacy_bam(args: tuple[str, int, int]) -> tuple[int, int, int, int]:
    path, x, y = args
    # One primary mapped read1 represents one sequenced fragment. Flag mask 2308
    # excludes unmapped, secondary, and supplementary records while retaining
    # duplicate records for the mapped denominator.
    mapped = int(
        subprocess.check_output(
            ["samtools", "view", "-c", "-f", "64", "-F", "2308", path], text=True
        ).strip()
    )
    dup = int(
        subprocess.check_output(
            ["samtools", "view", "-c", "-f", "1088", "-F", "2308", path], text=True
        ).strip()
    )
    return x, y, mapped, dup


def legacy_pixel_counts(
    record: SampleRecord,
    whitelist_path: Path,
    cache_dir: Path,
    force: bool = False,
    skip_expensive: bool = False,
    workers: int = 8,
) -> tuple[pd.DataFrame | None, Path | None]:
    cache = cache_dir / f"{record.sample}_legacy_pixel_counts.tsv"
    if cache.is_file() and not force:
        return pd.read_csv(cache, sep="\t"), cache
    if skip_expensive:
        LOG.warning("Skipping uncached legacy BAM recomputation for %s", record.sample)
        return None, None
    if record.legacy_root is None:
        return None, None
    bam_dir = record.legacy_root / "bam" / record.sample
    if not bam_dir.is_dir():
        return None, None
    whitelist = read_legacy_whitelist(whitelist_path).set_index("barcode")
    tasks = []
    for path in sorted(bam_dir.glob("*.bam")):
        if path.stem in whitelist.index:
            row = whitelist.loc[path.stem]
            tasks.append((str(path), int(row["x"]), int(row["y"])))
    LOG.info("Recomputing legacy per-pixel mapped fragments from %d BAMs", len(tasks))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        result = list(pool.map(_count_legacy_bam, tasks, chunksize=20))
    frame = pd.DataFrame(result, columns=["x", "y", "mapped_reads", "dup_reads"])
    frame["unique_fragments"] = frame["mapped_reads"] - frame["dup_reads"]
    frame.to_csv(cache, sep="\t", index=False)
    return frame, cache


def records_frame(records: list[SampleRecord]) -> pd.DataFrame:
    rows = []
    for record in records:
        row = asdict(record)
        row["revised_root"] = str(row["revised_root"] or "")
        row["legacy_root"] = str(row["legacy_root"] or "")
        rows.append(row)
    return pd.DataFrame(rows)
