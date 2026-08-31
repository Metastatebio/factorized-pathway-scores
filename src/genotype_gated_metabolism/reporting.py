"""Dependency-free rendering helpers for reproducible text reports."""

from __future__ import annotations

import pandas as pd


def _markdown_cell(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ")


def dataframe_to_markdown(frame: pd.DataFrame) -> str:
    """Render a small DataFrame as a GitHub-flavoured Markdown table."""
    columns = [_markdown_cell(value) for value in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(_markdown_cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])
