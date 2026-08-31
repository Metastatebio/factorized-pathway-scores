from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.reporting import dataframe_to_markdown


def test_dataframe_to_markdown_has_no_optional_dependency() -> None:
    frame = pd.DataFrame({"name": ["A|B", "C"], "value": [1, None]})

    rendered = dataframe_to_markdown(frame)

    assert "| name | value |" in rendered
    assert "| A\\|B | 1.0 |" in rendered
    assert "| C |  |" in rendered
