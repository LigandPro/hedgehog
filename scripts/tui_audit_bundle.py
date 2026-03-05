#!/usr/bin/env python3
"""Build a TUI audit bundle: file hierarchy, screen/module map, and text captures.

This script creates a structured output folder for auditing Hedgehog TUI UI/UX.
It is intended for repeatable screen reviews after UI changes.
"""

from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import os
import pathlib
import pty
import re
import select
import subprocess
import time
from dataclasses import dataclass

ROOT = pathlib.Path(__file__).resolve().parents[1]
TUI_DIR = ROOT / "tui"
APP_TSX = TUI_DIR / "src" / "App.tsx"
TYPES_TS = TUI_DIR / "src" / "types" / "index.ts"

ANSI_RE = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
CTRL_RE = re.compile(r"[\x00-\x08\x0b-\x1f\x7f]")
TREE_SKIP_DIRS = {
    ".git",
    ".venv",
    ".mypy_cache",
    ".pytest_cache",
    "__pycache__",
    "node_modules",
    "dist",
    ".DS_Store",
}


@dataclass(frozen=True)
class CaptureStep:
    name: str
    keys: str
    wait_seconds: float = 0.9
    enter_delay_seconds: float = 0.0
    inter_key_delay_seconds: float = 0.0
    actions: tuple[tuple[str, float], ...] = ()


CAPTURE_STEPS: list[CaptureStep] = [
    CaptureStep("00_welcome", "", 1.2),
    CaptureStep("01_open_config_main", "c"),
    CaptureStep("02_back_to_home", "\x1b"),
    CaptureStep("03_open_history", "h"),
    CaptureStep("04_back_to_home", "\x1b"),
    CaptureStep("05_open_wizard_input", "n"),
    CaptureStep("06_wizard_stage_selection", "\r"),
    CaptureStep("07_wizard_review", "r"),
    CaptureStep("08_open_command_menu", "/"),
    CaptureStep("09_close_command_menu", "\x1b"),
    CaptureStep("10_open_command_menu_ctrl_p", "\x10"),
    # Select `/theme` from the already open command palette.
    CaptureStep(
        "11_open_theme_menu_from_commands",
        "",
        wait_seconds=1.0,
        inter_key_delay_seconds=0.25,
        actions=(
            ("\x1b", 0.6),
            ("/theme\r", 2.0),
        ),
    ),
    CaptureStep("12_close_theme_menu", "\x1b"),
    CaptureStep(
        "13_return_home",
        "",
        wait_seconds=1.3,
        inter_key_delay_seconds=0.25,
        actions=(
            ("/back\r", 2.0),
            ("/back\r", 2.0),
            ("/back\r", 2.2),
        ),
    ),
    CaptureStep("14_exit_app", "\x03"),
]


def parse_args() -> argparse.Namespace:
    default_out = (
        ROOT / ".thoughts" / "tui-screen-unification" / "live-audit" / "latest"
    )
    parser = argparse.ArgumentParser(description="Generate TUI audit bundle.")
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=default_out,
        help="Output directory for the audit bundle.",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Replace output directory if it exists.",
    )
    return parser.parse_args()


def ensure_dir(path: pathlib.Path, replace: bool = False) -> None:
    if replace and path.exists():
        for child in sorted(path.rglob("*"), reverse=True):
            if child.is_file() or child.is_symlink():
                child.unlink(missing_ok=True)
            elif child.is_dir():
                child.rmdir()
        path.rmdir()
    path.mkdir(parents=True, exist_ok=True)


def write_text(path: pathlib.Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def list_relative_files(base: pathlib.Path, glob_pattern: str) -> list[str]:
    return sorted(
        str(p.relative_to(ROOT)) for p in base.glob(glob_pattern) if p.is_file()
    )


def render_hierarchy() -> str:
    sections: list[tuple[str, list[str]]] = [
        (
            "TUI Components",
            list_relative_files(TUI_DIR / "src" / "components", "*.tsx"),
        ),
        ("TUI Screens", list_relative_files(TUI_DIR / "src" / "screens", "*.tsx")),
        (
            "TUI Wizard Screens",
            list_relative_files(TUI_DIR / "src" / "screens" / "wizard", "*.tsx"),
        ),
        ("TUI Theme", list_relative_files(TUI_DIR / "src" / "theme", "*.ts*")),
        ("TUI Hooks", list_relative_files(TUI_DIR / "src" / "hooks", "*.ts*")),
        ("TUI Store", list_relative_files(TUI_DIR / "src" / "store", "*.ts*")),
        (
            "Backend Handlers",
            list_relative_files(
                ROOT / "src" / "hedgehog" / "tui_backend" / "handlers", "*.py"
            ),
        ),
        (
            "Backend Core",
            list_relative_files(ROOT / "src" / "hedgehog" / "tui_backend", "*.py"),
        ),
    ]

    lines = [
        "# Hedgehog TUI Hierarchy",
        "",
        f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}",
        "",
    ]
    for title, entries in sections:
        lines.append(f"## {title}")
        if not entries:
            lines.append("- (none)")
        else:
            for item in entries:
                lines.append(f"- {item}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_repo_tree(base: pathlib.Path, max_depth: int = 5) -> str:
    lines: list[str] = []

    def walk(path: pathlib.Path, prefix: str, depth: int) -> None:
        if depth > max_depth:
            return
        try:
            entries = sorted(
                [p for p in path.iterdir() if p.name not in TREE_SKIP_DIRS],
                key=lambda p: (not p.is_dir(), p.name.lower()),
            )
        except PermissionError:
            lines.append(f"{prefix}[permission denied]")
            return

        for idx, entry in enumerate(entries):
            is_last = idx == len(entries) - 1
            branch = "└── " if is_last else "├── "
            lines.append(f"{prefix}{branch}{entry.name}{'/' if entry.is_dir() else ''}")
            if entry.is_dir() and depth < max_depth:
                extension = "    " if is_last else "│   "
                walk(entry, prefix + extension, depth + 1)

    lines.append(f"{base.name}/")
    walk(base, "", 1)
    return "\n".join(lines) + "\n"


def parse_screen_union() -> list[str]:
    text = TYPES_TS.read_text(encoding="utf-8")
    match = re.search(r"export type Screen\s*=\s*(.+?);", text, flags=re.S)
    if not match:
        return []
    body = match.group(1)
    screens = re.findall(r"'([^']+)'", body)
    return screens


def parse_screen_router_mapping() -> list[tuple[str, str]]:
    text = APP_TSX.read_text(encoding="utf-8")
    imports: dict[str, str] = {}
    for m in re.finditer(r"import\s+\{\s*([^}]+)\s*\}\s+from\s+'([^']+)';", text):
        names = [x.strip() for x in m.group(1).split(",")]
        path = m.group(2)
        for n in names:
            imports[n] = path

    mapping: list[tuple[str, str]] = []
    case_re = re.compile(r"case '([^']+)':\s*return\s+<([A-Za-z0-9_]+)\s*/>;", re.S)
    for m in case_re.finditer(text):
        screen, component = m.group(1), m.group(2)
        src = imports.get(component, "(unknown)")
        mapping.append((screen, f"{component} -> {src}"))
    return mapping


def render_screen_map() -> str:
    declared_screens = parse_screen_union()
    router_mapping = parse_screen_router_mapping()
    mapped = {screen for screen, _ in router_mapping}

    lines = [
        "# Screen and Module Map",
        "",
        f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Declared Screen IDs",
    ]
    for s in declared_screens:
        lines.append(f"- {s}")
    lines.extend(["", "## App Router Mapping"])
    for screen, target in router_mapping:
        lines.append(f"- {screen}: {target}")

    missing = [s for s in declared_screens if s not in mapped]
    lines.extend(["", "## Coverage Notes"])
    if missing:
        for s in missing:
            lines.append(f"- Missing router case: {s}")
    else:
        lines.append("- All declared screen IDs have router cases.")
    lines.append("")
    return "\n".join(lines)


def render_audit_checklist() -> str:
    screens = parse_screen_union()
    lines = [
        "# Visual Audit Checklist",
        "",
        "Use this checklist while reviewing the clean capture files.",
        "",
        "## Global Checks",
        "- Header alignment and title spacing are stable.",
        "- Footer is pinned and readable (shortcuts + version).",
        "- No overlapping text between body and footer.",
        "- No clipped status/info lines in the center body region.",
        "- Slash command menu is modal and does not leak background input.",
        "- Theme menu is modal and closes cleanly with Esc/Ctrl+C.",
        "",
        "## Screen-by-Screen Checks",
    ]
    for s in screens:
        lines.append(f"- [ ] {s}")
    lines.append("")
    return "\n".join(lines)


def set_nonblocking(fd: int) -> None:
    flags = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, flags | os.O_NONBLOCK)


def read_available(fd: int, timeout: float = 0.15) -> str:
    chunks: list[bytes] = []
    while True:
        rlist, _, _ = select.select([fd], [], [], timeout)
        if not rlist:
            break
        try:
            data = os.read(fd, 65536)
        except BlockingIOError:
            break
        if not data:
            break
        chunks.append(data)
        timeout = 0.03
    return b"".join(chunks).decode("utf-8", errors="replace")


def extract_last_frame(raw: str) -> str:
    # Split by clear-screen commands and keep the last repaint segment.
    parts = re.split(r"\x1b\[[0-9;?]*2J", raw)
    segment = parts[-1] if parts else raw
    # Keep from last "home" cursor if present.
    home_parts = re.split(r"\x1b\[[0-9;?]*H", segment)
    segment = home_parts[-1] if home_parts else segment
    return segment


def clean_frame(raw: str) -> str:
    text = ANSI_RE.sub("", raw)
    text = text.replace("\r", "")
    text = CTRL_RE.sub("", text)
    lines = text.splitlines()

    # Trim leading/trailing blank lines, keep body spacing in between.
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return "\n".join(lines) + ("\n" if lines else "")


def run_capture_sequence(raw_dir: pathlib.Path, clean_dir: pathlib.Path) -> list[str]:
    def send_keys(fd: int, keys: str, inter_key_delay_seconds: float) -> str:
        if not keys:
            return ""
        captured = ""
        if inter_key_delay_seconds <= 0:
            os.write(fd, keys.encode("utf-8", errors="ignore"))
            return read_available(fd, timeout=0.03)
        for ch in keys:
            os.write(fd, ch.encode("utf-8", errors="ignore"))
            time.sleep(inter_key_delay_seconds)
            captured += read_available(fd, timeout=0.03)
        return captured

    master_fd, slave_fd = pty.openpty()
    set_nonblocking(master_fd)

    proc = subprocess.Popen(
        ["node", "dist/index.js"],
        cwd=TUI_DIR,
        stdin=slave_fd,
        stdout=slave_fd,
        stderr=slave_fd,
        text=False,
        close_fds=True,
    )
    os.close(slave_fd)

    notes: list[str] = []
    all_output = ""
    time.sleep(0.8)
    all_output += read_available(master_fd, timeout=0.4)

    try:
        for idx, step in enumerate(CAPTURE_STEPS, start=1):
            if step.actions:
                for keys, delay_seconds in step.actions:
                    all_output += send_keys(
                        master_fd, keys, step.inter_key_delay_seconds
                    )
                    time.sleep(delay_seconds)
                    all_output += read_available(master_fd, timeout=0.05)
            elif step.keys:
                if step.enter_delay_seconds > 0 and "\r" in step.keys:
                    enter_idx = step.keys.find("\r")
                    typed = step.keys[:enter_idx]
                    submitted = step.keys[enter_idx:]
                    if typed:
                        all_output += send_keys(
                            master_fd, typed, step.inter_key_delay_seconds
                        )
                    time.sleep(step.enter_delay_seconds)
                    all_output += send_keys(
                        master_fd, submitted, step.inter_key_delay_seconds
                    )
                else:
                    all_output += send_keys(
                        master_fd, step.keys, step.inter_key_delay_seconds
                    )
            time.sleep(step.wait_seconds)
            all_output += read_available(master_fd, timeout=0.35)

            frame = extract_last_frame(all_output)
            cleaned = clean_frame(frame)

            name = f"{idx:02d}_{step.name}"
            raw_path = raw_dir / f"{name}.txt"
            clean_path = clean_dir / f"{name}.txt"
            write_text(raw_path, frame)
            write_text(clean_path, cleaned)
            notes.append(f"- {name}: captured")

            if proc.poll() is not None:
                notes.append(
                    f"- Process exited early at step {name} with code {proc.returncode}"
                )
                break
    finally:
        if proc.poll() is None:
            try:
                os.write(master_fd, b"\x03")
                time.sleep(0.2)
            except OSError:
                pass
            proc.terminate()
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                proc.kill()
        os.close(master_fd)

    return notes


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    ensure_dir(output_dir, replace=args.replace)

    inventory_dir = output_dir / "inventory"
    captures_raw_dir = output_dir / "captures" / "raw"
    captures_clean_dir = output_dir / "captures" / "clean"
    reports_dir = output_dir / "reports"
    for d in (inventory_dir, captures_raw_dir, captures_clean_dir, reports_dir):
        d.mkdir(parents=True, exist_ok=True)

    write_text(inventory_dir / "hierarchy.md", render_hierarchy())
    write_text(inventory_dir / "repo_tree.txt", build_repo_tree(ROOT, max_depth=5))
    write_text(inventory_dir / "screen_map.md", render_screen_map())
    write_text(inventory_dir / "audit_checklist.md", render_audit_checklist())

    capture_notes = run_capture_sequence(captures_raw_dir, captures_clean_dir)
    summary = [
        "# TUI Audit Bundle Summary",
        "",
        f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Output Structure",
        f"- inventory: {inventory_dir}",
        f"- captures/raw: {captures_raw_dir}",
        f"- captures/clean: {captures_clean_dir}",
        "",
        "## Capture Notes",
    ]
    summary.extend(capture_notes or ["- No captures were created."])
    summary.append("")
    write_text(reports_dir / "summary.md", "\n".join(summary))

    print(f"Audit bundle created at: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
