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
