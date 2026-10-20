#!/usr/bin/env python3
"""Run a command concurrently on newline-aligned blocks from standard input."""

from __future__ import annotations

import argparse
import os
import re
import signal
import subprocess
import sys
import threading
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait


_SIZE_RE = re.compile(r"^(\d+)([KMGT]?)(?:I?B)?$", re.IGNORECASE)


def parse_size(value: str) -> int:
    match = _SIZE_RE.fullmatch(value.strip())
    if not match:
        raise argparse.ArgumentTypeError(
            f"invalid block size {value!r}; use bytes or a K/M/G/T suffix"
        )
    number = int(match.group(1))
    suffix = match.group(2).upper()
    if number < 1:
        raise argparse.ArgumentTypeError("block size must be positive")
    return number * (1024 ** " KMGT".index(suffix))


def parse_env(values: list[str]) -> dict[str, str]:
    updates = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"invalid --env value {value!r}; expected KEY=VALUE")
        key, setting = value.split("=", 1)
        if not key:
            raise ValueError("environment variable name cannot be empty")
        updates[key] = setting
    return updates


def read_newline_block(stream, target_size: int) -> bytes:
    """Read at least target_size bytes without splitting the final line."""
    chunks = []
    remaining = target_size
    while remaining:
        chunk = stream.read(remaining)
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    if not chunks:
        return b""
    if chunks[-1].endswith(b"\n"):
        return b"".join(chunks)
    chunks.append(stream.readline())
    return b"".join(chunks)


class CommandRunner:
    def __init__(self, command: list[str], env: dict[str, str]) -> None:
        self.command = command
        self.env = env
        self._lock = threading.Lock()
        self._children: set[subprocess.Popen[bytes]] = set()
        self._stopping = threading.Event()

    def run(self, block_id: int, payload: bytes) -> int:
        command = [part.replace("{#}", str(block_id)) for part in self.command]
        with self._lock:
            if self._stopping.is_set():
                return 125
            process = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                env=self.env,
                start_new_session=True,
            )
            self._children.add(process)
        try:
            process.communicate(payload)
            return process.returncode
        finally:
            with self._lock:
                self._children.discard(process)

    def terminate_all(self) -> None:
        self._stopping.set()
        with self._lock:
            children = list(self._children)
        for process in children:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        for process in children:
            if process.poll() is None:
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass


def parse_args(argv: list[str]) -> tuple[argparse.Namespace, list[str]]:
    try:
        separator = argv.index("--")
    except ValueError:
        separator = -1
    if separator < 0 or separator == len(argv) - 1:
        raise ValueError("a worker command must follow --")

    parser = argparse.ArgumentParser(
        description=(
            "Feed newline-aligned stdin blocks to bounded concurrent workers. "
            "The token {#} in a command argument is replaced by the block number."
        )
    )
    parser.add_argument("--jobs", type=int, required=True)
    parser.add_argument("--block", type=parse_size, required=True)
    parser.add_argument("--env", action="append", default=[])
    args = parser.parse_args(argv[:separator])
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    return args, argv[separator + 1 :]


def completed_status(future: Future[int], block_id: int) -> int:
    try:
        status = future.result()
    except Exception as error:  # pragma: no cover - exercised by OS launch failures
        print(f"stream_parallel: block {block_id} failed to launch: {error}", file=sys.stderr)
        return 125
    if status:
        print(
            f"stream_parallel: block {block_id} exited with status {status}",
            file=sys.stderr,
        )
    return status


def run_stream(args: argparse.Namespace, command: list[str]) -> int:
    env = os.environ.copy()
    env.update(parse_env(args.env))
    runner = CommandRunner(command, env)
    active: dict[Future[int], int] = {}
    block_count = 0
    failed_status = 0
    exhausted = False

    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            while active or not exhausted:
                while not exhausted and len(active) < args.jobs:
                    payload = read_newline_block(sys.stdin.buffer, args.block)
                    if not payload:
                        exhausted = True
                        break
                    block_count += 1
                    future = executor.submit(runner.run, block_count, payload)
                    active[future] = block_count

                    for completed in [item for item in active if item.done()]:
                        block_id = active.pop(completed)
                        status = completed_status(completed, block_id)
                        if status and not failed_status:
                            failed_status = status
                    if failed_status:
                        break

                if failed_status:
                    runner.terminate_all()
                    break
                if active:
                    done, _ = wait(active, return_when=FIRST_COMPLETED)
                    for completed in done:
                        block_id = active.pop(completed)
                        status = completed_status(completed, block_id)
                        if status and not failed_status:
                            failed_status = status
                    if failed_status:
                        runner.terminate_all()
                        break
    except KeyboardInterrupt:
        runner.terminate_all()
        print("stream_parallel: interrupted", file=sys.stderr)
        return 130

    if failed_status:
        return failed_status if 0 < failed_status < 256 else 1
    print(f"stream_parallel: completed {block_count} blocks", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    try:
        args, command = parse_args(argv)
        return run_stream(args, command)
    except ValueError as error:
        print(f"stream_parallel: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
