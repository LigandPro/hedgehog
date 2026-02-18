"""Tests for plain-text log formatting."""

from __future__ import annotations

import logging

from hedgehog.configs.logger import _PlainTextFormatter


def _format_message(message: str) -> str:
    formatter = _PlainTextFormatter()
    record = logging.LogRecord("test", logging.INFO, "", 0, message, (), None)
    return formatter.format(record)


def test_plain_text_formatter_strips_rich_tags_only():
    formatted = _format_message("Status [bold]ok[/bold] [#B29EEE]done[/#B29EEE]")
    assert "Status ok done" in formatted


def test_plain_text_formatter_preserves_chemical_brackets():
    formatted = _format_message("Rejected molecule [NH4+] due to score")
    assert "[NH4+]" in formatted
