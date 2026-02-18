#!/usr/bin/env python3
"""Unified verification pipeline for CLI + TUI."""

from __future__ import annotations

import argparse
import os
import pty
import select
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def run_cmd(cmd: list[str], cwd: Path) -> None:
    """Run a command and stream output."""
    pretty_cmd = " ".join(cmd)
    print(f"[check] Running: {pretty_cmd}", flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def run_cli_smoke(repo_root: Path) -> None:
    """Run required CLI smoke checks."""
    run_cmd(["uv", "run", "hedgehog", "--help"], repo_root)
    run_cmd(["uv", "run", "hedge", "--help"], repo_root)
    run_cmd(["uv", "run", "hedgehog", "run", "--help"], repo_root)
    run_cmd(["uv", "run", "hedgehog", "setup", "--help"], repo_root)
    run_cmd(["uv", "run", "hedgehog", "setup", "aizynthfinder", "--help"], repo_root)
    run_cmd(["uv", "run", "hedgehog", "version"], repo_root)


def build_tui(repo_root: Path) -> None:
    """Install and build TUI."""
    tui_dir = repo_root / "tui"
    run_cmd(["npm", "ci"], tui_dir)
    run_cmd(["npm", "run", "build"], tui_dir)


def run_tui_smoke(repo_root: Path, startup_timeout_sec: int = 25) -> None:
    """Run TUI in a pseudo-TTY and verify backend handshake."""
    if os.name != "posix":
        raise RuntimeError("TUI smoke requires a POSIX environment with PTY support.")

    tui_dir = repo_root / "tui"
    entry = tui_dir / "dist" / "index.js"
    if not entry.exists():
        raise FileNotFoundError(f"TUI entry point not found: {entry}")

    with tempfile.NamedTemporaryFile(
        prefix="hedgehog-tui-smoke-", suffix=".log", delete=False
    ) as tmp:
        log_path = Path(tmp.name)

    master_fd, slave_fd = pty.openpty()
    proc: subprocess.Popen[bytes] | None = None
    sent_quit = False
    started_at = time.monotonic()
    output = bytearray()

    env = os.environ.copy()
    env["HEDGEHOG_TUI_LOG"] = str(log_path)
    env["HEDGEHOG_TUI_DEBUG"] = "1"

    try:
        print("[check] Running: TUI smoke (PTY startup + graceful quit)", flush=True)
        proc = subprocess.Popen(
            ["node", str(entry)],
            cwd=tui_dir,
            stdin=slave_fd,
            stdout=slave_fd,
            stderr=slave_fd,
            env=env,
            start_new_session=True,
        )
        os.close(slave_fd)

        while True:
            if proc.poll() is not None:
                break

            now = time.monotonic()

            if not sent_quit and now - started_at >= 3:
                os.write(master_fd, b"q")
                sent_quit = True

            if now - started_at > startup_timeout_sec:
                raise TimeoutError(
                    f"TUI smoke timed out after {startup_timeout_sec} seconds."
                )

            ready, _, _ = select.select([master_fd], [], [], 0.25)
            if ready:
                try:
                    chunk = os.read(master_fd, 4096)
                except OSError:
                    # PTY can raise EIO when the child exits and the master side is drained.
                    if proc.poll() is not None:
                        break
                    continue
                if chunk:
                    output.extend(chunk)

        exit_code = proc.returncode
        if exit_code != 0:
            decoded = output.decode("utf-8", errors="ignore")
            raise RuntimeError(
                f"TUI exited with non-zero code {exit_code}. Captured output:\n{decoded}"
            )

        decoded = output.decode("utf-8", errors="ignore")
        if "HEDGEHOG PIPELINE ENGINE" not in decoded:
            raise RuntimeError("TUI banner was not observed during smoke run.")

        log_content = log_path.read_text(encoding="utf-8", errors="ignore")
        required_markers = [
            "Python backend ready",
            "Backend initialized",
            "Hedgehog TUI exited",
        ]
        missing = [marker for marker in required_markers if marker not in log_content]
        if missing:
            raise RuntimeError(
                "TUI smoke log is missing markers: " + ", ".join(missing)
            )

        print("[check] TUI smoke passed.", flush=True)
    finally:
        try:
            os.close(master_fd)
        except OSError:
            pass

        if proc and proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

        if log_path.exists():
            log_path.unlink()


def run_full_pipeline(repo_root: Path) -> None:
    """Run the production pipeline with the default config (all stages)."""
    run_cmd(["uv", "run", "hedgehog", "run"], repo_root)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Unified checks for CLI and TUI.",
    )
    parser.add_argument(
        "--mode",
        choices=["quick", "ci", "full"],
        default="quick",
        help="quick/ci: smoke checks, full: smoke checks + full production run.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]

    run_cli_smoke(repo_root)
    build_tui(repo_root)
    run_tui_smoke(repo_root)

    if args.mode == "full":
        print("[check] Running full production pipeline.", flush=True)
        run_full_pipeline(repo_root)
    else:
        print(f"[check] Full production run skipped in '{args.mode}' mode.", flush=True)

    print("[check] All requested checks passed.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
