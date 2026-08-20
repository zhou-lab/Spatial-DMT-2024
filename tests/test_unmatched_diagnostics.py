import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pysam


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "tools" / "unmatched_diagnostics.py"


class UnmatchedDiagnosticsCliTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.fastq = self.root / "unmatched_R2.fq"
        self.out = self.root / "diagnostics.tsv"
        self.mqc = self.root / "summary_mqc.tsv"
        self.records = [
            ("read1_SF", "A" * 120),
            (
                "read2_SF",
                "GGTGTAGTGGGTTTGGAGG" + "A" * 111,
            ),
            ("read3_WM", "C" * 120),
            ("read4_WM", "T" * 120),
        ]
        with self.fastq.open("w") as handle:
            for name, sequence in self.records:
                handle.write(
                    f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                )

    def tearDown(self):
        self.tempdir.cleanup()

    def run_cli(self, *extra_args):
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--r2-fastq",
                str(self.fastq),
                "--protocol",
                "smcseq",
                "--sample",
                "fixture",
                "--out",
                str(self.out),
                "--mqc-summary",
                str(self.mqc),
                "--threads",
                "2",
                *extra_args,
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    def parse_long_tsv(self):
        rows = {}
        with self.out.open() as handle:
            self.assertEqual(handle.readline(), "section\tkey\tvalue\n")
            for line in handle:
                section, key, value = line.rstrip("\n").split("\t")
                rows[(section, key)] = value
        return rows

    def test_fastq_only_mode_omits_alignment_metrics(self):
        self.run_cli()
        rows = self.parse_long_tsv()
        self.assertEqual(
            rows[("fastq_stats", "structure_fail_n_reads")], "2")
        self.assertEqual(rows[("fastq_stats", "wl_miss_n_reads")], "2")
        self.assertEqual(rows[("sf_anchor_fail", "primer_fail_n")], "1")
        self.assertEqual(rows[("sf_anchor_fail", "linker1_fail_n")], "1")
        self.assertFalse(any(section == "align_stats" for section, _ in rows))

        mqc_text = self.mqc.read_text()
        table_header = next(
            line for line in mqc_text.splitlines()
            if line.startswith("Sample\t")
        )
        self.assertIn("sf_n_reads", table_header)
        self.assertIn("wm_n_reads", table_header)
        self.assertNotIn("mapping_rate", mqc_text)

    def test_optional_bam_adds_alignment_metrics(self):
        bam = self.root / "unmatched.bam"
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
        }
        with pysam.AlignmentFile(bam, "wb", header=header) as handle:
            for index, (name, sequence) in enumerate(self.records):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = sequence
                read.query_qualities = pysam.qualitystring_to_array(
                    "I" * len(sequence))
                if name.endswith("_SF"):
                    read.flag = 77
                else:
                    read.flag = 73 | (1024 if index == 3 else 0)
                    read.reference_id = 0
                    read.reference_start = index * 10
                    read.mapping_quality = 60
                    read.cigarstring = f"{len(sequence)}M"
                handle.write(read)

        self.run_cli("--bam", str(bam))
        rows = self.parse_long_tsv()
        self.assertEqual(
            rows[("align_stats", "structure_fail_n_reads")], "2")
        self.assertEqual(
            rows[("align_stats", "structure_fail_n_mapped")], "0")
        self.assertEqual(rows[("align_stats", "wl_miss_n_reads")], "2")
        self.assertEqual(rows[("align_stats", "wl_miss_n_mapped")], "2")
        self.assertEqual(rows[("align_stats", "wl_miss_n_dup")], "1")
        self.assertIn("mapping_rate", self.mqc.read_text())


if __name__ == "__main__":
    unittest.main()
