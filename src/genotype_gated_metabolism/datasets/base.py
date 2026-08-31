"""Common data structures for cohort and model-system adapters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class CrossOmicsDataset:
    """Sample-aligned metadata and named feature blocks.

    Every table is indexed by the same stable sample identifier. Adapters may
    initially load different sample sets; ``aligned`` explicitly takes their
    intersection and preserves metadata order.
    """

    sample_metadata: pd.DataFrame
    blocks: Mapping[str, pd.DataFrame]
    provenance: Mapping[str, object]

    def __post_init__(self) -> None:
        if not self.sample_metadata.index.is_unique:
            raise ValueError("sample_metadata index must contain unique sample identifiers.")
        if self.sample_metadata.index.hasnans:
            raise ValueError("sample_metadata index cannot contain missing identifiers.")
        if not self.blocks:
            raise ValueError("At least one omics block is required.")
        for name, block in self.blocks.items():
            if not name:
                raise ValueError("Omics block names cannot be empty.")
            if not block.index.is_unique:
                raise ValueError(f"Block {name!r} has duplicate sample identifiers.")
            if not block.columns.is_unique:
                raise ValueError(f"Block {name!r} has duplicate feature identifiers.")

    @property
    def sample_ids(self) -> pd.Index:
        """Return the metadata-ordered sample identifiers."""
        return self.sample_metadata.index

    def aligned(self, required_blocks: list[str] | None = None) -> CrossOmicsDataset:
        """Return a sample-intersected copy across metadata and selected blocks."""
        selected = required_blocks or list(self.blocks)
        unknown = sorted(set(selected) - set(self.blocks))
        if unknown:
            raise KeyError(f"Unknown omics blocks: {', '.join(unknown)}")

        shared = self.sample_metadata.index
        for name in selected:
            shared = shared.intersection(self.blocks[name].index, sort=False)
        if shared.empty:
            raise ValueError("No samples are shared across the requested omics blocks.")

        metadata = self.sample_metadata.loc[shared].copy()
        blocks = {name: block.loc[shared].copy() for name, block in self.blocks.items()}
        return CrossOmicsDataset(metadata, blocks, dict(self.provenance))

    def require_features(self, block_name: str, feature_ids: list[str]) -> pd.DataFrame:
        """Return requested features in order, raising on absent identifiers."""
        if block_name not in self.blocks:
            raise KeyError(f"Unknown omics block: {block_name}")
        block = self.blocks[block_name]
        missing = [feature for feature in feature_ids if feature not in block.columns]
        if missing:
            raise KeyError(f"Missing {block_name} features: {', '.join(missing)}")
        return block.loc[:, feature_ids].copy()
