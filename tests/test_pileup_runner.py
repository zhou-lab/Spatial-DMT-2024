import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "tools" / "run_biscuit_pileup_batches.sh"


class PileupRunnerTest(unittest.TestCase):
    def test_resolves_bare_tool_names_from_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bindir = root / "bin"
            bindir.mkdir()
            for name in ("fake-converter", "fake-yame", "fake-planner", "fake-python", "fake-biscuit"):
                tool = bindir / name
                tool.write_text("#!/usr/bin/env bash\nexit 0\n")
                tool.chmod(0o755)

            result = subprocess.run(
                [
                    str(RUNNER),
                    "--cell-list", str(root / "missing-cells.txt"),
                    "--bam-dir", str(root / "missing-bams"),
                    "--tmp-dir", str(root / "tmp"),
                    "--output-cg", str(root / "out.cg"),
                    "--output-index", str(root / "out.cg.idx"),
                    "--output-conversion", str(root / "conversion.tsv"),
                    "--reference", str(root / "missing.fa"),
                    "--cpg-converter", "fake-converter",
                    "--cpg-index", str(root / "missing.cpgidx"),
                    "--yame", "fake-yame",
                    "--planner", "fake-planner",
                    "--python", "fake-python",
                    "--biscuit", "fake-biscuit",
                    "--cpus", "2",
                    "--batch-size", "1",
                    "--max-batch-bytes", "0",
                    "--threads-per-batch", "1",
                    "--cpg-threads", "1",
                    "--window-step", "1000",
                ],
                env={**os.environ, "PATH": f"{bindir}:{os.environ['PATH']}"},
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("missing or empty input", result.stderr)
            self.assertNotIn("missing tool", result.stderr)
            self.assertNotIn("tool not found in PATH", result.stderr)


if __name__ == "__main__":
    unittest.main()
