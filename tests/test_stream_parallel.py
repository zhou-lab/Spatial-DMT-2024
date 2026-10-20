import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
LAUNCHER = REPO / "tools" / "stream_parallel.py"


class StreamParallelTest(unittest.TestCase):
    def test_preserves_complete_lines_and_block_order(self):
        payload = b"one\ntwo-two\nthree\nfour-four-four\nfive\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            worker = root / "worker.py"
            worker.write_text(
                "import pathlib, sys\n"
                "pathlib.Path(sys.argv[1]).write_bytes(sys.stdin.buffer.read())\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LAUNCHER),
                    "--jobs", "3",
                    "--block", "8",
                    "--",
                    sys.executable,
                    str(worker),
                    str(root / "block_{#}.txt"),
                ],
                input=payload,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr.decode())
            outputs = sorted(root.glob("block_*.txt"), key=lambda path: int(path.stem[6:]))
            self.assertGreater(len(outputs), 1)
            self.assertEqual(b"".join(path.read_bytes() for path in outputs), payload)
            self.assertTrue(all(path.read_bytes().endswith(b"\n") for path in outputs))

    def test_returns_worker_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            worker = Path(tmpdir) / "worker.py"
            worker.write_text(
                "import sys\n"
                "sys.stdin.buffer.read()\n"
                "raise SystemExit(7 if sys.argv[1] == '2' else 0)\n"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(LAUNCHER),
                    "--jobs", "1",
                    "--block", "2",
                    "--",
                    sys.executable,
                    str(worker),
                    "{#}",
                ],
                input=b"a\nb\nc\n",
                capture_output=True,
            )
            self.assertEqual(result.returncode, 7)
            self.assertIn(b"block 2 exited with status 7", result.stderr)


if __name__ == "__main__":
    unittest.main()
