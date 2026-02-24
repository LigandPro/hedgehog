# Tools Index

This index lists helper tools that reduce repeated agent mistakes and improve verification speed.
When adding or changing a script, update this table in the same change.

| Tool | Path | Use When | Input | Output |
| --- | --- | --- | --- | --- |
| Unified check pipeline | `scripts/check_pipeline.py` | validating CLI + TUI health before merge | `--mode quick|ci|full` | smoke/full verification run |
| CLI smoke commands | `uv run hedgehog ...` | touching CLI command behavior | command-line flags | help/version contract validation |
| TUI build + smoke commands | `tui/package.json` scripts + `node tui/dist/index.js` | touching TUI screens/backend bridge | Node environment, built dist | TUI startup/build confidence |
| Agent mode guide | `docs/tools/agent-mode.md` | preparing low-noise agent execution | runtime flags + command choices | reproducible agent workflow |

## Trigger Rule

If an agent repeats the same class of error twice, add or improve a tool before adding more prose instructions.

## Related Docs

- [`../conventions.md`](../conventions.md)
- [`../quality.md`](../quality.md)
