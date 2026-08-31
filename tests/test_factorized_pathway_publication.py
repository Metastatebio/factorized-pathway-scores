from __future__ import annotations

import pandas as pd

from genotype_gated_metabolism.pipelines.factorized_pathway_publication import (
    _claim_ledger,
)


def test_claim_ledger_separates_supported_and_failed_gates() -> None:
    st = {
        "gate_results": {
            "rmse_beats_random_grouping": False,
        }
    }
    ccle = {
        "gate_results": {
            "beats_random_factorization": True,
        }
    }
    structural = {"decision": "ADAPTIVE_STRUCTURE_RESOLUTION_SIGNAL"}
    replication = {"decision": "EXTERNAL_STRUCTURAL_REPLICATION"}
    ledger = _claim_ledger(st, structural, replication, ccle)
    decisions = dict(zip(ledger["claim"], ledger["decision"], strict=True))
    assert decisions["Broad lipid-family scores are pathway-specific"] == "NOT_SUPPORTED"
    assert (
        decisions[
            "Direct HumanGEM/GPR features outperform size-matched random features"
        ]
        == "SUPPORTED"
    )
    assert (
        decisions["Structure-resolved lipid scores beat a degree-preserving null"]
        == "SUPPORTED_ADAPTIVE_SAME_COHORT"
    )
    assert (
        decisions[
            "Structure-resolved lipid scores replicate in a separate human cohort"
        ]
        == "SUPPORTED_EXTERNAL_POPULATION_HELD_OUT"
    )
    assert isinstance(ledger, pd.DataFrame)
