# CCLE reaction-cluster robustness

**Decision:** `CCLE_REACTION_CLUSTER_ROBUST`

## Cluster-resampled inference

The property-matched expected-null effect was 0.1374 SD across 60 targets assigned to 20 subsystems. The 95% cluster-bootstrap interval was 0.0857 to 0.1937.

## Subsystem summary

| subsystem | targets | mean_effect_sd |
| --- | --- | --- |
| Aminoacyl-tRNA biosynthesis | 1 | -0.082238 |
| Valine, leucine, and isoleucine metabolism | 2 | -0.054048 |
| Inositol phosphate metabolism | 1 | -0.017768 |
| Bile acid biosynthesis | 2 | 0.007057 |
| Glycerophospholipid metabolism | 3 | 0.031032 |
| Glycine, serine and threonine metabolism | 6 | 0.041045 |
| Isolated | 1 | 0.041174 |
| Arginine and proline metabolism | 6 | 0.0754 |
| Pyruvate metabolism | 1 | 0.084757 |
| Alanine, aspartate and glutamate metabolism | 5 | 0.095622 |
| Phenylalanine, tyrosine and tryptophan biosynthesis | 3 | 0.123846 |
| Beta-alanine metabolism | 3 | 0.134779 |
| Purine metabolism | 10 | 0.161448 |
| Nucleotide metabolism | 2 | 0.188783 |
| Nicotinate and nicotinamide metabolism | 1 | 0.209914 |
| Peptide metabolism | 1 | 0.231698 |
| Pyrimidine metabolism | 6 | 0.237941 |
| Histidine metabolism | 1 | 0.255128 |
| Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism | 4 | 0.385613 |
| Phenylalanine metabolism | 1 | 0.607512 |

## Adjudication

| clusters | minimum_clusters_pass | cluster_ci_lower | cluster_ci_pass | minimum_leave_one_cluster_effect_sd | leave_one_cluster_pass |
| --- | --- | --- | --- | --- | --- |
| 20 | True | 0.085702 | True | 0.119672 | True |

## Boundary

This post-result audit tests whether the property-matched CCLE topology result survives subsystem-cluster resampling and leave-one-subsystem-out analysis. It does not establish target independence, causal reaction direction, physiological flux, drug response or clinical utility.
