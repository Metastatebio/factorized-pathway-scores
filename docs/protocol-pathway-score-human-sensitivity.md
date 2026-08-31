# Post-result human pathway-score sensitivity protocol

**Version:** 1.0  
**Freeze date:** 30 August 2026, after the primary and replication results  
**Role:** robustness analysis, not untouched confirmation

## Objective

Test whether the structure-resolved advantage depends on the ridge penalty, descriptor-frequency
filter, maximum descriptor prevalence or graph-randomization intensity. The ST002081 participant-
held-out discovery cohort and ST000818 population-held-out replication cohort are rerun without
changing samples, feature eligibility, hidden panels or outer folds.

## Frozen settings

Nine settings are evaluated in each cohort:

1. primary: alpha 10, minimum descriptor frequency 5, maximum prevalence 0.90, 10 swaps per edge;
2. alpha 1;
3. alpha 100;
4. minimum descriptor frequency 3;
5. minimum descriptor frequency 10;
6. maximum descriptor prevalence 0.80;
7. maximum descriptor prevalence 0.95;
8. 5 attempted swaps per edge; and
9. 25 attempted swaps per edge.

Only the named parameter changes from the primary setting. Every structural representation is
compared with a degree-preserving graph realization generated from a setting-specific deterministic
seed. All graphs must preserve exact feature and descriptor degrees.

## Endpoints

For each setting and cohort we report structural-minus-null improvement in group-weighted RMSE and
precision@3, with 2,000 paired group-bootstrap draws. Robustness requires positive point estimates
for both endpoints in every setting and wholly positive 95% intervals in at least 80% of settings
in each cohort. Complete results are reported; no setting is selected or omitted.

This grid addresses model-definition sensitivity only. It does not turn a post-result analysis into
prospective confirmation and does not test genetics, flux or clinical utility.
