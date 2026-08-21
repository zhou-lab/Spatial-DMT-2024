#!/usr/bin/env python3
"""Build deterministic, size-balanced BAM manifests for BISCUIT pileup."""

from __future__ import annotations

import argparse
import hashlib
import heapq
import os
from dataclasses import dataclass
from pathlib import Path
import re
import tempfile


SAFE_CELL = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class Bam:
    cell: str
    path: Path
    size: int
    mtime_ns: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Partition cell BAMs into deterministic, size-balanced batches."
        )
    )
    parser.add_argument("--cell-list", required=True, type=Path)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest-list", required=True, type=Path)
    parser.add_argument("--batch-size", required=True, type=int)
    parser.add_argument(
        "--max-batch-bytes",
        type=int,
        default=0,
        help=(
            "Target maximum aggregate BAM bytes per batch; zero disables "
            "the byte limit. A single larger BAM remains valid."
        ),
    )
    parser.add_argument(
        "--signature",
        default="",
        help="Pileup settings included in each batch content key.",
    )
    return parser.parse_args()


def read_bams(cell_list: Path, bam_dir: Path) -> list[Bam]:
    cells: list[str] = []
    seen: set[str] = set()
    with cell_list.open(encoding="ascii") as stream:
        for line_number, raw in enumerate(stream, start=1):
            cell = raw.rstrip("\r\n")
            if not cell:
                continue
            if not SAFE_CELL.fullmatch(cell) or cell in {".", ".."}:
                raise ValueError(
                    f"unsafe cell name at {cell_list}:{line_number}: {cell!r}"
                )
            if cell in seen:
                raise ValueError(f"duplicate cell in {cell_list}: {cell}")
            seen.add(cell)
            cells.append(cell)

    if not cells:
        raise ValueError(f"no cells found in {cell_list}")

    bams: list[Bam] = []
    for cell in cells:
        path = (bam_dir / f"{cell}.bam").resolve()
        stat = path.stat()
        if not path.is_file() or stat.st_size == 0:
            raise ValueError(f"missing or empty cell BAM: {path}")
        bams.append(Bam(cell, path, stat.st_size, stat.st_mtime_ns))
    return bams


def balance_bams(
    bams: list[Bam], batch_size: int, max_batch_bytes: int
) -> list[list[Bam]]:
    if batch_size < 1:
        raise ValueError("batch size must be positive")
    if batch_size > 4096:
        raise ValueError("batch size cannot exceed 4096 VCF samples")
    if max_batch_bytes < 0:
        raise ValueError("maximum batch bytes cannot be negative")
    batch_count = (len(bams) + batch_size - 1) // batch_size
    if max_batch_bytes:
        total_size = sum(bam.size for bam in bams)
        byte_batches = (total_size + max_batch_bytes - 1) // max_batch_bytes
        batch_count = max(batch_count, byte_batches)
    batch_count = min(batch_count, len(bams))
    batches: list[list[Bam]] = [[] for _ in range(batch_count)]
    heap: list[tuple[int, int, int]] = [
        (0, 0, index) for index in range(batch_count)
    ]
    heapq.heapify(heap)

    for bam in sorted(bams, key=lambda item: (-item.size, item.cell)):
        total_size, count, index = heapq.heappop(heap)
        batches[index].append(bam)
        count += 1
        total_size += bam.size
        if count < batch_size:
            heapq.heappush(heap, (total_size, count, index))
    return batches


def batch_digest(batch: list[Bam], signature: str) -> str:
    digest = hashlib.sha256()
    digest.update(signature.encode("utf-8"))
    digest.update(b"\0")
    for bam in sorted(batch, key=lambda item: item.cell):
        digest.update(bam.cell.encode("ascii"))
        digest.update(b"\0")
        digest.update(str(bam.size).encode("ascii"))
        digest.update(b"\0")
        digest.update(str(bam.mtime_ns).encode("ascii"))
        digest.update(b"\0")
        digest.update(os.fsencode(bam.path))
        digest.update(b"\0")
    return digest.hexdigest()[:20]


def write_if_changed(path: Path, content: str) -> None:
    try:
        if path.read_text(encoding="utf-8") == content:
            return
    except FileNotFoundError:
        pass

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    bams = read_bams(args.cell_list, args.bam_dir)
    batches = balance_bams(bams, args.batch_size, args.max_batch_bytes)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifests: list[Path] = []
    totals: list[int] = []
    for batch in batches:
        digest = batch_digest(batch, args.signature)
        manifest = (args.output_dir / f"batch-{digest}.bams").resolve()
        content = "".join(f"{bam.path}\n" for bam in batch)
        write_if_changed(manifest, content)
        manifests.append(manifest)
        totals.append(sum(bam.size for bam in batch))

    manifest_content = "".join(f"{path}\n" for path in manifests)
    write_if_changed(args.manifest_list, manifest_content)
    print(
        f"planned {len(bams)} BAMs in {len(batches)} batches; "
        f"batch bytes min={min(totals)} max={max(totals)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
