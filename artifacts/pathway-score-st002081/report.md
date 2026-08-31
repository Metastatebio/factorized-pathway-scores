# ST002081 factorized pathway missing-panel benchmark

**Decision:** `RESOURCE_ONLY`  
**Directional forecasting tested:** false  
**Participant genomics available:** false

## Result

The locked benchmark evaluated 493 complete lipid features from 9 biochemical families in 1539 samples from 112 people. Every feature was hidden exactly once across 5 disjoint masks.

| model | row_weighted_rmse_sd | equal_subject_weighted_rmse_sd | pooled_r2 | precision_at_k | ndcg_at_k |
| --- | --- | --- | --- | --- | --- |
| all_visible_ridge | 0.1747 | 0.178 | 0.9703 | 0.9725 | 0.9619 |
| family_score_ridge | 0.6503 | 0.668 | 0.5878 | 0.4029 | 0.6331 |
| population_mean | 1.0129 | 1.0652 | 0.0 | 0.1545 | 0.3694 |
| random_group_score_ridge | 0.6271 | 0.6417 | 0.6166 | 0.3691 | 0.6063 |

## Paired participant bootstrap

- `family_score_ridge` versus `population_mean`: RMSE improvement 0.3972 SD (95% interval 0.3506 to 0.4460); precision@k improvement 0.2338 (0.2067 to 0.2623).
- `family_score_ridge` versus `random_group_score_ridge`: RMSE improvement -0.0263 SD (95% interval -0.0429 to -0.0104); precision@k improvement 0.0166 (-0.0247 to 0.0591).
- `family_score_ridge` versus `all_visible_ridge`: RMSE improvement -0.4900 SD (95% interval -0.5107 to -0.4694); precision@k improvement -0.5829 (-0.6137 to -0.5518).

## Interpretation boundary

Participant-held-out missing-panel reconstruction can validate compact biochemical-family representation and marker-priority enrichment. It does not establish directional dynamics, genetic conditioning, physiological flux, diagnosis, treatment response or clinical utility.
