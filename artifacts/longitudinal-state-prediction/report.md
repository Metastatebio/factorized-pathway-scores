# Subject-held-out longitudinal lipid-state prediction

**Decision:** RESOURCE_ONLY  
**Genetic comparison supported:** false

## Result

The analysis retained 100 people, 1000 strictly forward transitions, and 512 transition-complete lipids.
The prior-profile autoregressive model achieved pooled out-of-subject R²=0.452 and RMSE=0.755 SD.

## Paired subject-cluster bootstrap

- Versus static_ridge: RMSE improvement 0.280 SD (95% bootstrap interval 0.234 to 0.330).
- Versus persistence: RMSE improvement 0.103 SD (95% bootstrap interval 0.075 to 0.132).
- Recorded context added to the prior state: RMSE improvement -0.009 SD (95% bootstrap interval -0.016 to -0.003); negative values favor the prior-state-only model.
- Precision@3 for membership in the observed top 10% of changes: 0.592 versus prevalence 0.420; bootstrap difference 0.159 (0.119 to 0.198). Exact top-k precision was 0.225 versus 0.082.
- Temporal falsification: forward R²=0.452, reverse-time R²=0.459, and within-person shuffled-time R²=0.415.

## Interpretation boundary

Subject-held-out reconstruction can quantify how much a prior lipid profile predicts another profile from the same person. It is evidence of directional dynamics only if forward-time performance beats reverse-time and within-person shuffled-time falsification tests. Precision at k measures enrichment for membership in a broad high-change set, not exact top-k recovery. This study has no participant genomics, so it cannot compare nurture with genetic effects, test genotype-by-context interactions, establish metabolic mobility, or establish clinical utility.
