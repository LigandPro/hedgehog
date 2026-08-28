"""Launch AiZynthFinder with concise ONNX Runtime logging."""

from __future__ import annotations

import os


def _import_onnxruntime_quietly():
    """Import ONNX Runtime without its non-actionable Linux GPU probe warning."""
    saved_stderr = os.dup(2)
    try:
        with open(os.devnull, "w", encoding="utf-8") as sink:
            os.dup2(sink.fileno(), 2)
            import onnxruntime
    finally:
        os.dup2(saved_stderr, 2)
        os.close(saved_stderr)
    return onnxruntime


def main() -> int:
    """Run the upstream CLI while keeping real ONNX errors visible."""
    onnxruntime = _import_onnxruntime_quietly()
    onnxruntime.set_default_logger_severity(3)
    from aizynthfinder.interfaces.aizynthcli import main as aizynthfinder_main

    result = aizynthfinder_main()
    return int(result or 0)


if __name__ == "__main__":
    raise SystemExit(main())
