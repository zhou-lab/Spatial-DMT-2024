import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "tools" / "write_pipeline_version.py"


class PipelineVersionTest(unittest.TestCase):
    def test_writes_release_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            version = root / "VERSION"
            output = root / "pipeline_version.tsv"
            version.write_text("2.0.1\n")
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--sample", "sample1",
                    "--version-file", str(version),
                    "--protocol", "spatialdmt",
                    "--reference", "mm10",
                    "--output", str(output),
                ],
                check=True,
            )
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["sample"], "sample1")
            self.assertEqual(row["pipeline_version"], "2.0.1")
            self.assertEqual(row["assignment"], "workflow_generated")

    def test_rejects_empty_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            version = root / "VERSION"
            output = root / "pipeline_version.tsv"
            version.write_text("\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--sample", "sample1",
                    "--version-file", str(version),
                    "--protocol", "smcseq",
                    "--reference", "mm10",
                    "--output", str(output),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
