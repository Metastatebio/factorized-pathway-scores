import numpy as np

from genotype_gated_metabolism.datasets.ccle import canonicalize_growth_medium


def test_canonicalize_growth_medium_orders_overlapping_names() -> None:
    assert canonicalize_growth_medium("DMEM:F12 (1:1) + 5% FBS") == "DMEM_F12"
    assert canonicalize_growth_medium("RPMI1640+10%FBS") == "RPMI"
    assert canonicalize_growth_medium("EMEM +10%FBS") == "EMEM"
    assert canonicalize_growth_medium("MEM + NEAA") == "MEM"
    assert canonicalize_growth_medium(np.nan) == "unknown"
