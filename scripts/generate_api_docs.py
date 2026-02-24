#!/usr/bin/env python3
"""Generate an API overview page for docs from runtime introspection."""

from __future__ import annotations

import argparse
import ast
import importlib
import inspect
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import click
from typer.main import get_command

from hedgehog.main import app as hedgehog_app
from hedgehog.tui_backend.server import JsonRpcServer

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "docs" / "content" / "api" / "index.mdx"
SRC_ROOT = REPO_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "hedgehog"

RPC_ERROR_CODES = [
    ("-32601", "Method not found", "Unknown RPC method name."),
    ("-32602", "Invalid params", "Validation/value error in method parameters."),
    ("-32000", "Server error", "Unhandled internal exception."),
    ("-32001", "File not found", "Requested file/path does not exist."),
    ("-32002", "Permission denied", "Insufficient permissions for operation."),
    ("-32003", "Not a directory", "Path expected to be a directory is not one."),
]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate docs/content/api/index.mdx from live API metadata."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output MDX file path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def _md_cell(value: Any) -> str:
    text = "-" if value is None else str(value)
    text = " ".join(text.split())
    return text.replace("|", "\\|")


def _first_line(text: str | None) -> str:
    if not text:
        return ""
    return text.strip().splitlines()[0].strip()


def _format_default(default: Any) -> str:
    if default is None:
        return "None"
    if default is ...:
        return "..."
    if callable(default):
        return "<callable>"

    if isinstance(default, Path):
        try:
            return str(default.resolve().relative_to(REPO_ROOT))
        except ValueError:
            return str(default)

    if isinstance(default, str):
        maybe_path = Path(default)
        if maybe_path.is_absolute():
            try:
                return str(maybe_path.relative_to(REPO_ROOT))
            except ValueError:
                return default
        return default

    if isinstance(default, (list, tuple, set)):
        return repr(list(default))

    return repr(default)


def _extract_cli() -> list[dict[str, Any]]:
    root_cmd = get_command(hedgehog_app)
    commands: list[dict[str, Any]] = []

    def walk(cmd: click.Command, parts: list[str]) -> None:
        params: list[dict[str, str]] = []
        for param in cmd.params:
            if getattr(param, "hidden", False):
                continue

            if isinstance(param, click.Option):
                flags = ", ".join([*param.opts, *param.secondary_opts]) or "-"
                params.append(
                    {
                        "name": param.name or "-",
                        "kind": "option",
                        "flags": flags,
                        "type": str(param.type),
                        "required": "yes" if param.required else "no",
                        "default": _format_default(param.default),
                        "description": _first_line(param.help),
                    }
                )
            elif isinstance(param, click.Argument):
                params.append(
                    {
                        "name": param.name or "-",
                        "kind": "argument",
                        "flags": "-",
                        "type": str(param.type),
                        "required": "yes" if param.required else "no",
                        "default": _format_default(param.default),
                        "description": _first_line(getattr(param, "help", "")),
                    }
                )

        commands.append(
            {
                "path": " ".join(parts),
                "kind": "group" if isinstance(cmd, click.Group) else "command",
                "description": _first_line(cmd.help),
                "params": params,
            }
        )

        if isinstance(cmd, click.Group):
            for sub_name in sorted(cmd.commands):
                walk(cmd.commands[sub_name], [*parts, sub_name])

    for name in sorted(root_cmd.commands):
        walk(root_cmd.commands[name], [name])

    return commands


def _extract_tui_rpc() -> list[dict[str, str]]:
    server = JsonRpcServer()
    methods: list[dict[str, str]] = []
    for method_name in sorted(server.handlers):
        handler = server.handlers[method_name]
        signature = str(inspect.signature(handler))
        qualname = getattr(handler, "__qualname__", getattr(handler, "__name__", ""))
        handler_ref = f"{handler.__module__}.{qualname}"
        methods.append(
            {
                "method": method_name,
                "signature": signature,
                "handler": handler_ref,
            }
        )
    return methods


def _module_name_from_init(init_file: Path) -> str:
    relative = init_file.relative_to(SRC_ROOT)
    module_path = relative.with_suffix("")
    parts = list(module_path.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _extract_all_symbols(init_file: Path) -> list[str]:
    source = init_file.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(init_file))

    for node in tree.body:
        if isinstance(node, ast.Assign):
            if not any(
                isinstance(target, ast.Name) and target.id == "__all__"
                for target in node.targets
            ):
                continue
            try:
                value = ast.literal_eval(node.value)
            except Exception:
                return []
            if isinstance(value, (list, tuple)):
                return [str(item) for item in value if isinstance(item, str)]

        if isinstance(node, ast.AnnAssign):
            if not (isinstance(node.target, ast.Name) and node.target.id == "__all__"):
                continue
            if node.value is None:
                return []
            try:
                value = ast.literal_eval(node.value)
            except Exception:
                return []
            if isinstance(value, (list, tuple)):
                return [str(item) for item in value if isinstance(item, str)]

    return []


def _classify_symbol(obj: Any) -> str:
    if inspect.isclass(obj):
        return "class"
    if inspect.isfunction(obj) or inspect.ismethod(obj) or inspect.isbuiltin(obj):
        return "function"
    return "other"


def _extract_python_public_api() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    init_files = sorted(PACKAGE_ROOT.rglob("__init__.py"))

    for init_file in init_files:
        module_name = _module_name_from_init(init_file)
        if module_name.startswith("hedgehog.vendor"):
            continue
        if module_name.startswith("hedgehog.workers"):
            continue

        exports = _extract_all_symbols(init_file)
        if not exports:
            continue

        try:
            module = importlib.import_module(module_name)
            import_error = ""
        except Exception as exc:  # noqa: BLE001
            module = None
            import_error = str(exc)

        for symbol in sorted(exports):
            if module is None:
                rows.append(
                    {
                        "module": module_name,
                        "symbol": symbol,
                        "kind": "unresolved",
                        "notes": f"import failed: {import_error}",
                    }
                )
                continue

            try:
                obj = getattr(module, symbol)
            except AttributeError:
                rows.append(
                    {
                        "module": module_name,
                        "symbol": symbol,
                        "kind": "missing",
                        "notes": "export not found on module at import time",
                    }
                )
                continue
            except Exception as exc:  # noqa: BLE001
                rows.append(
                    {
                        "module": module_name,
                        "symbol": symbol,
                        "kind": "unresolved",
                        "notes": f"symbol load failed: {exc}",
                    }
                )
                continue

            rows.append(
                {
                    "module": module_name,
                    "symbol": symbol,
                    "kind": _classify_symbol(obj),
                    "notes": "",
                }
            )

    rows.sort(key=lambda row: (row["module"], row["symbol"]))
    return rows


def _render_cli_section(cli_commands: list[dict[str, Any]]) -> list[str]:
    lines = [
        "## CLI",
        "",
        "The CLI metadata is collected from `hedgehog.main:app` via Typer/Click introspection.",
        "",
        "```bash",
        "uv run hedgehog --help",
        "```",
        "",
        "### Commands",
        "",
        "| Command | Kind | Description |",
        "| --- | --- | --- |",
    ]

    for command in cli_commands:
        lines.append(
            "| `{}` | {} | {} |".format(
                _md_cell(command["path"]),
                _md_cell(command["kind"]),
                _md_cell(command["description"] or "-"),
            )
        )

    for command in cli_commands:
        lines.extend(
            [
                "",
                f"### `hedgehog {command['path']}`",
                "",
                command["description"] or "No description.",
                "",
            ]
        )

        params = command["params"]
        if not params:
            lines.append("No parameters.")
            continue

        lines.extend(
            [
                "| Parameter | Kind | Flags | Type | Required | Default | Description |",
                "| --- | --- | --- | --- | --- | --- | --- |",
            ]
        )

        for param in params:
            lines.append(
                "| `{}` | {} | `{}` | `{}` | {} | `{}` | {} |".format(
                    _md_cell(param["name"]),
                    _md_cell(param["kind"]),
                    _md_cell(param["flags"]),
                    _md_cell(param["type"]),
                    _md_cell(param["required"]),
                    _md_cell(param["default"]),
                    _md_cell(param["description"] or "-"),
                )
            )

    return lines


def _render_tui_section(rpc_methods: list[dict[str, str]]) -> list[str]:
    lines = [
        "## TUI JSON-RPC",
        "",
        "The TUI backend API is read from `JsonRpcServer.handlers` in `hedgehog.tui_backend.server`.",
        "",
        "```json",
        '{"jsonrpc":"2.0","id":1,"method":"preflight_pipeline","params":{"stages":["descriptors"]}}',
        "```",
        "",
        "### Methods",
        "",
        "| Method | Signature | Handler |",
        "| --- | --- | --- |",
    ]

    for method in rpc_methods:
        lines.append(
            "| `{}` | `{}` | `{}` |".format(
                _md_cell(method["method"]),
                _md_cell(method["signature"]),
                _md_cell(method["handler"]),
            )
        )

    lines.extend(
        [
            "",
            "### Standard Error Codes",
            "",
            "| Code | Name | Meaning |",
            "| --- | --- | --- |",
        ]
    )

    for code, name, meaning in RPC_ERROR_CODES:
        lines.append(f"| `{_md_cell(code)}` | {_md_cell(name)} | {_md_cell(meaning)} |")

    return lines


def _render_python_section(python_api: list[dict[str, str]]) -> list[str]:
    lines = [
        "## Python Public API",
        "",
        "Public Python API entries are discovered from `__all__` in `src/hedgehog/**/__init__.py`.",
        "",
        "```python",
        "from hedgehog.reporting import ReportGenerator",
        "from hedgehog.setup import ensure_aizynthfinder",
        "```",
        "",
        "Modules under `hedgehog.vendor.*` and `hedgehog.workers.*` are intentionally excluded.",
        "",
    ]

    if not python_api:
        lines.append("No public exports were discovered.")
        return lines

    lines.extend(
        [
            "| Module | Symbol | Kind | Notes |",
            "| --- | --- | --- | --- |",
        ]
    )

    for row in python_api:
        lines.append(
            "| `{}` | `{}` | {} | {} |".format(
                _md_cell(row["module"]),
                _md_cell(row["symbol"]),
                _md_cell(row["kind"]),
                _md_cell(row["notes"] or "-"),
            )
        )

    return lines


def _render_document(
    cli_commands: list[dict[str, Any]],
    rpc_methods: list[dict[str, str]],
    python_api: list[dict[str, str]],
) -> str:
    now_utc = datetime.now(timezone.utc)
    generated_at = now_utc.strftime("%Y-%m-%d %H:%M:%S UTC")
    generated_timestamp_ms = int(now_utc.timestamp() * 1000)
    lines: list[str] = [
        "---",
        "title: API Overview",
        f"timestamp: {generated_timestamp_ms}",
        "---",
        "",
        "# API Overview",
        "",
        f"Generated automatically at build time by `scripts/generate_api_docs.py` ({generated_at}).",
        "",
        "This page is intentionally high-level and reflects the current runtime API surface.",
        "",
    ]

    lines.extend(_render_cli_section(cli_commands))
    lines.append("")
    lines.extend(_render_tui_section(rpc_methods))
    lines.append("")
    lines.extend(_render_python_section(python_api))
    lines.append("")

    return "\n".join(lines)


def main() -> int:
    args = _parse_args()
    output_path = args.output.resolve()

    cli_commands = _extract_cli()
    rpc_methods = _extract_tui_rpc()
    python_api = _extract_python_public_api()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        _render_document(cli_commands, rpc_methods, python_api),
        encoding="utf-8",
    )

    print(
        f"[api-docs] Wrote {output_path} "
        f"(cli={len(cli_commands)}, rpc={len(rpc_methods)}, py={len(python_api)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
