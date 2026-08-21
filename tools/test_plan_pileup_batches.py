#!/usr/bin/env python3

from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plan_pileup_batches import Bam, balance_bams, batch_digest, write_if_changed


class BatchPlannerTests(unittest.TestCase):
    def test_balances_size_and_cell_limits(self) -> None:
        bams = [
            Bam(f"{index:04d}", Path(f"/{index:04d}.bam"), size, index)
            for index, size in enumerate(range(10, 0, -1), start=1)
        ]
        batches = balance_bams(bams, batch_size=4, max_batch_bytes=20)

        self.assertEqual(len(batches), 3)
        self.assertEqual(sum(map(len, batches)), 10)
        self.assertLessEqual(max(map(len, batches)), 4)
        totals = [sum(bam.size for bam in batch) for batch in batches]
        self.assertLessEqual(max(totals) - min(totals), 1)

    def test_single_large_bam_is_allowed(self) -> None:
        bams = [
            Bam("large", Path("/large.bam"), 100, 1),
            Bam("small", Path("/small.bam"), 1, 1),
        ]
        batches = balance_bams(bams, batch_size=16, max_batch_bytes=10)

        self.assertEqual(len(batches), 2)
        self.assertIn(["large"], [[bam.cell for bam in batch] for batch in batches])

    def test_digest_is_order_independent_and_sensitive(self) -> None:
        first = Bam("A", Path("/A.bam"), 10, 1)
        second = Bam("B", Path("/B.bam"), 20, 2)

        self.assertEqual(
            batch_digest([first, second], "settings"),
            batch_digest([second, first], "settings"),
        )
        self.assertNotEqual(
            batch_digest([first, second], "settings"),
            batch_digest([first, second], "changed"),
        )

    def test_write_if_changed_preserves_existing_inode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest"
            write_if_changed(path, "A\n")
            inode = path.stat().st_ino
            write_if_changed(path, "A\n")

            self.assertEqual(path.stat().st_ino, inode)
            self.assertEqual(path.read_text(encoding="utf-8"), "A\n")


if __name__ == "__main__":
    unittest.main()
