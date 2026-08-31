# Adaptive ST002081 structure-resolved pathway-score benchmark

**Decision:** `ADAPTIVE_STRUCTURE_RESOLUTION_SIGNAL`  
**Independent replication:** false

## Result

The adaptive score used 95 filtered structural descriptors across 493 lipids. Every null preserved the visible feature and descriptor degree sequences exactly.

| model | samples | subjects | prediction_pairs | row_weighted_rmse_sd | equal_subject_weighted_rmse_sd | pooled_r2 | precision_at_k | ndcg_at_k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| degree_preserving_random_structural_ridge | 1539 | 112 | 758727 | 0.4096 | 0.4193 | 0.8365 | 0.744 | 0.804 |
| population_mean | 1539 | 112 | 758727 | 1.0129 | 1.0652 | 0.0 | 0.1545 | 0.3694 |
| structural_descriptor_ridge | 1539 | 112 | 758727 | 0.3091 | 0.3168 | 0.9069 | 0.8877 | 0.8909 |

## Paired participant bootstrap

- `structural_descriptor_ridge` versus `population_mean`: RMSE improvement 0.7484 SD (0.6974 to 0.8013); precision improvement 0.7276 (0.7001 to 0.7538).
- `structural_descriptor_ridge` versus `degree_preserving_random_structural_ridge`: RMSE improvement 0.1025 SD (0.0971 to 0.1078); precision improvement 0.1368 (0.1221 to 0.1522).
- `structural_descriptor_ridge` versus `family_score_ridge`: RMSE improvement 0.3512 SD (0.3314 to 0.3709); precision improvement 0.4938 (0.4599 to 0.5272).
- `structural_descriptor_ridge` versus `all_visible_ridge`: RMSE improvement -0.1388 SD (-0.1447 to -0.1333); precision improvement -0.0891 (-0.0993 to -0.0794).

## Boundary

This adaptive same-cohort follow-up can test whether declared lipid structural descriptors outperform a degree-preserving random descriptor graph. It is not independent replication and does not establish temporal direction, genetic conditioning, flux or clinical utility.
