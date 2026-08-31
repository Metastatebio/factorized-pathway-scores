# Human graph-null mixing diagnostic

**Decision:** `HUMAN_GRAPH_NULL_MIXING_ADEQUATE`

## Dataset summary

| dataset | null_panel_rows | minimum_edge_replacement | maximum_observed_null_jaccard | minimum_unique_nulls | maximum_pairwise_null_jaccard |
| --- | --- | --- | --- | --- | --- |
| ST000818 | 100 | 0.864776 | 0.072515 | 20 | 0.076291 |
| ST002081 | 100 | 0.770569 | 0.12958 | 20 | 0.131141 |

## Mask-level null diversity

| dataset | mask | expected_nulls | unique_nulls | minimum_pairwise_jaccard | mean_pairwise_jaccard | maximum_pairwise_jaccard |
| --- | --- | --- | --- | --- | --- | --- |
| ST002081 | 0 | 20 | 20 | 0.112052 | 0.119344 | 0.129812 |
| ST002081 | 1 | 20 | 20 | 0.111067 | 0.121596 | 0.131141 |
| ST002081 | 2 | 20 | 20 | 0.107899 | 0.119488 | 0.128909 |
| ST002081 | 3 | 20 | 20 | 0.109443 | 0.118091 | 0.127711 |
| ST002081 | 4 | 20 | 20 | 0.104167 | 0.11693 | 0.127886 |
| ST000818 | 0 | 20 | 20 | 0.047181 | 0.060083 | 0.072481 |
| ST000818 | 1 | 20 | 20 | 0.048711 | 0.06107 | 0.074574 |
| ST000818 | 2 | 20 | 20 | 0.050401 | 0.062097 | 0.076291 |
| ST000818 | 3 | 20 | 20 | 0.047887 | 0.061805 | 0.073903 |
| ST000818 | 4 | 20 | 20 | 0.044053 | 0.060178 | 0.07605 |

## Adjudication

| all_degrees_preserved | minimum_edge_replacement_fraction | edge_replacement_pass | all_nulls_unique | maximum_pairwise_null_jaccard | pairwise_diversity_pass |
| --- | --- | --- | --- | --- | --- |
| True | 0.770569 | True | True | 0.131141 | True |

## Boundary

This post-result diagnostic verifies exact degree preservation, graph movement and null diversity for the reported human structural ensembles. It does not prove perfectly uniform graph sampling, independent biological replication, physiology, genetics or clinical utility.
