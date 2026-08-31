from pathlib import Path

import numpy as np

from genotype_gated_metabolism.datasets.mwtab import load_mwtab


def test_load_mwtab_builds_aligned_cross_omics_dataset(tmp_path: Path) -> None:
    source = tmp_path / "study.txt"
    source.write_text(
        "#SUBJECT_SAMPLE_FACTORS:\tSUBJECT\tSAMPLE\tFACTORS\tAdditional sample data\n"
        "SUBJECT_SAMPLE_FACTORS\tperson-1\ts1\tCL4:Healthy | Sex:F\t"
        "RandomID=person; AgeAtCollection(month)=100\n"
        "SUBJECT_SAMPLE_FACTORS\tperson-2\ts2\tCL4:Infection | Sex:F\t"
        "RandomID=person; AgeAtCollection(month)=101\n"
        "MS_METABOLITE_DATA_START\n"
        "Samples\ts1\ts2\n"
        "Factors\ta\tb\n"
        "lipid-a\t1.0\t2.0\n"
        "lipid-b\t3.0\t\n"
        "MS_METABOLITE_DATA_END\n"
    )

    dataset = load_mwtab(source)

    assert dataset.sample_ids.tolist() == ["s1", "s2"]
    assert dataset.sample_metadata.loc["s1", "additional_RandomID"] == "person"
    assert dataset.sample_metadata.loc["s2", "additional_AgeAtCollection(month)"] == 101
    assert dataset.blocks["metabolomics"].shape == (2, 2)
    assert np.isnan(dataset.blocks["metabolomics"].loc["s2", "lipid-b"])
