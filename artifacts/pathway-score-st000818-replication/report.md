# ST000818 external structural pathway-score replication

**Decision:** `EXTERNAL_STRUCTURAL_REPLICATION`

## Cohort

The locked analysis retained 255 complete lipids and 81 descriptors across 450 participants from 15 population categories. Every outer fold excluded complete categories.

## Model performance

| model | samples | subjects | prediction_pairs | row_weighted_rmse_sd | equal_subject_weighted_rmse_sd | pooled_r2 | precision_at_k | equal_subject_precision_at_k | ndcg_at_k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_visible_ridge | 450 | 15 | 114750 | 0.8398 | 0.8412 | 0.8236 | 0.6858 | 0.6858 | 0.8007 |
| family_score_ridge | 450 | 15 | 114750 | 1.8777 | 1.8931 | 0.1184 | 0.2997 | 0.2997 | 0.5091 |
| population_mean | 450 | 15 | 114750 | 1.9998 | 2.0145 | 0.0 | 0.1363 | 0.1363 | 0.3282 |
| random_group_score_ridge | 450 | 15 | 114750 | 1.9126 | 1.9283 | 0.0853 | 0.2113 | 0.2113 | 0.4241 |
| degree_preserving_random_structural_ridge | 450 | 15 | 114750 | 1.7855 | 1.8006 | 0.2016 | 0.4785 |  | 0.6435 |
| structural_descriptor_ridge | 450 | 15 | 114750 | 1.2374 | 1.2418 | 0.6165 | 0.5739 |  | 0.7163 |

## Primary paired comparisons

| reference | challenger | subjects | draws | rmse_improvement_sd | rmse_ci_lower | rmse_ci_upper | precision_improvement | precision_ci_lower | precision_ci_upper |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| population_mean | structural_descriptor_ridge | 15 | 5000 | 0.772 | 0.3555 | 1.2048 | 0.4387 | 0.3996 | 0.4861 |
| family_score_ridge | structural_descriptor_ridge | 15 | 5000 | 0.6513 | 0.145 | 1.1065 | 0.2742 | 0.2373 | 0.3298 |
| all_visible_ridge | structural_descriptor_ridge | 15 | 5000 | -0.4005 | -0.7073 | -0.0739 | -0.1119 | -0.1259 | -0.0988 |

## Repeated degree-preserving nulls

Across 20 graph realizations, random RMSE ranged from 1.7289 to 1.8592 SD. Structural RMSE improvement ranged from 0.5027 to 0.6278 SD.

| null_seed | random_rmse_sd | rmse_improvement_sd | rmse_ci_lower | random_precision_at_k | precision_improvement | precision_ci_lower |
| --- | --- | --- | --- | --- | --- | --- |
| 20261200 | 1.7855 | 0.5588 | 0.0617 | 0.4785 | 0.0954 | 0.081 |
| 20261201 | 1.7638 | 0.538 | 0.0587 | 0.4806 | 0.0933 | 0.075 |
| 20261202 | 1.7702 | 0.5442 | 0.066 | 0.4618 | 0.1121 | 0.0892 |
| 20261203 | 1.7763 | 0.5483 | 0.0583 | 0.4676 | 0.1064 | 0.0859 |
| 20261204 | 1.7941 | 0.5681 | 0.069 | 0.4593 | 0.1147 | 0.0935 |
| 20261205 | 1.8592 | 0.6278 | 0.0631 | 0.476 | 0.0979 | 0.0799 |
| 20261206 | 1.7851 | 0.5592 | 0.0623 | 0.476 | 0.0979 | 0.0767 |
| 20261207 | 1.7289 | 0.5027 | 0.0576 | 0.4686 | 0.1053 | 0.0799 |
| 20261208 | 1.7751 | 0.5489 | 0.0656 | 0.4661 | 0.1079 | 0.0839 |
| 20261209 | 1.7591 | 0.5339 | 0.0625 | 0.4773 | 0.0966 | 0.0756 |
| 20261210 | 1.7774 | 0.5507 | 0.0582 | 0.4735 | 0.1004 | 0.0849 |
| 20261211 | 1.7964 | 0.5702 | 0.0619 | 0.4729 | 0.101 | 0.0844 |
| 20261212 | 1.8129 | 0.5853 | 0.059 | 0.4766 | 0.0973 | 0.0797 |
| 20261213 | 1.7767 | 0.5498 | 0.0632 | 0.4794 | 0.0945 | 0.0721 |
| 20261214 | 1.7823 | 0.5562 | 0.0635 | 0.473 | 0.1009 | 0.081 |
| 20261215 | 1.7843 | 0.5583 | 0.0618 | 0.4801 | 0.0938 | 0.076 |
| 20261216 | 1.777 | 0.5505 | 0.0574 | 0.483 | 0.091 | 0.0735 |
| 20261217 | 1.7866 | 0.5606 | 0.0652 | 0.4717 | 0.1022 | 0.0819 |
| 20261218 | 1.782 | 0.5563 | 0.0604 | 0.4702 | 0.1037 | 0.0844 |
| 20261219 | 1.7888 | 0.5626 | 0.0638 | 0.476 | 0.0979 | 0.0834 |

## Boundary

This separate public-cohort benchmark can replicate structure-resolved lipid-profile reconstruction across held-out population categories. It does not test participant genomics or transcriptomics and does not establish temporal direction, physiological flux, diagnosis, treatment response or clinical utility.
