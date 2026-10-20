import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class QcReportingTest(unittest.TestCase):
    def test_yame_sample_summaries_runs_rowop_once(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fake = root / "yame"
            count = root / "rowop.count"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "if [[ $1 == rowop ]]; then\n"
                f"  echo rowop >> {count}\n"
                "  printf 'aggregate-stream\\n'\n"
                "elif [[ $1 == summary ]]; then\n"
                "  cat >/dev/null\n"
                "  printf 'QFile\\tQuery\\tMFile\\tMask\\tN_univ\\tN_query\\tN_mask\\tN_overlap\\tLog2OddsRatio\\tBeta\\tDepth\\n'\n"
                "  printf 'sample\\tall\\tNA\\tall\\t100\\t42\\t100\\t42\\t0\\t0.61\\t3.5\\n'\n"
                "else exit 2; fi\n"
            )
            fake.chmod(0o755)
            cg = root / "sample.cg"
            cm = root / "mask.cm"
            aggregate = root / "aggregate.tsv"
            chromhmm = root / "chromhmm.tsv"
            cg.touch()
            cm.touch()
            subprocess.run(
                [
                    str(REPO / "tools" / "yame_sample_summaries.sh"),
                    "--yame", str(fake),
                    "--cg", str(cg),
                    "--aggregate-output", str(aggregate),
                    "--chromhmm-cm", str(cm),
                    "--chromhmm-output", str(chromhmm),
                ],
                check=True,
            )
            self.assertEqual(count.read_text().splitlines(), ["rowop"])
            self.assertEqual(aggregate.read_text(), chromhmm.read_text())

    def test_parse_precomputed_yame_summary(self):
        curated = load_module(
            "write_curated_metrics_mqc_summary",
            REPO / "tools" / "write_curated_metrics_mqc.py",
        )
        text = (
            "QFile\tQuery\tMFile\tMask\tN_univ\tN_query\tN_mask\t"
            "N_overlap\tLog2OddsRatio\tBeta\tDepth\n"
            "sample\tall\tNA\tall\t100\t42\t100\t42\t0\t0.61\t3.5\n"
        )
        self.assertEqual(curated.parse_yame_summary(text), (42.0, 3.5, 0.61))

    def test_linker_yield_uses_raw_structure_pass_denominator(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            trim = root / "trim.tsv"
            spatial = root / "spatial.tsv"
            chromhmm = root / "chromhmm.tsv"
            trim.write_text(
                "Reads_total\t100\n"
                "Reads_passed\t60\n"
                "Unmatched_SF\t20\n"
                "WL_unmatched\t10\n"
            )
            spatial.write_text("barcode\tx\ty\tdna_reads\n0101\t1\t1\t55\n")
            chromhmm.write_text("")
            subprocess.run(
                [
                    sys.executable,
                    str(REPO / "tools" / "write_stats_mqc.py"),
                    "--trim-stats", str(trim),
                    "--spatial-counts", str(spatial),
                    "--chromhmm-summary", str(chromhmm),
                    "--sample", "fixture",
                    "--output-dir", str(root),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            output = root / "fixture_stats_mqc.tsv"
            lines = [line for line in output.read_text().splitlines()
                     if line and not line.startswith("#")]
            row = next(csv.DictReader(lines, delimiter="\t"))
            self.assertEqual(row["total read pairs"], "100")
            self.assertEqual(row["both linkers"], "80")
            self.assertEqual(row["with designed barcodes"], "60")
            self.assertEqual(row["both linkers / total"], "80.00")
            self.assertEqual(row["designed / linkers"], "75.00")

    def test_failure_html_partitions_post_parse_drop(self):
        diagnostic = load_module(
            "generate_dna_r2_diagnostic_html",
            REPO / "tools" / "generate_dna_r2_diagnostic_html.py",
        )
        trim = {
            "Reads_total": 100,
            "Reads_passed": 60,
            "Unmatched_SF": 20,
            "Unmatched_WM": 10,
        }
        rendered = diagnostic.build_sections("fixture", trim, {})
        self.assertIn("Post-parse/read-length drop", rendered)
        self.assertIn("All pairs not emitted on the matched track", rendered)
        self.assertIn("<td>10</td>", rendered)

    def test_dupsifter_counts_both_duplicate_mates(self):
        curated = load_module(
            "write_curated_metrics_mqc",
            REPO / "tools" / "write_curated_metrics_mqc.py",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            stats = Path(tmpdir) / "dupsifter.stat"
            stats.write_text(
                "individual reads processed: 100\n"
                "both reads marked as duplicates: 10\n"
                "forward strand marked as duplicates: 5\n"
                "reverse strand marked as duplicates: 7\n"
            )
            self.assertAlmostEqual(curated.parse_dupsifter(stats), 0.32)


if __name__ == "__main__":
    unittest.main()
