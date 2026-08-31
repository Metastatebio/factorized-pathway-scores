# Post-result sensitivity protocol: repeated degree-preserving graph nulls

**Protocol version:** 1.0  
**Registration date:** 30 August 2026, after opening the first adaptive structural result  
**Role:** post-result robustness analysis; not an untouched confirmatory test

## Motivation

The adaptive ST002081 structural benchmark used one deterministic degree-preserving randomized
feature–descriptor graph. Its result could depend on that graph realization. This analysis repeats
the entire participant-isolated benchmark over 20 prespecified null seeds while holding the assay,
eligible features, structural descriptors, hidden panels, participant folds, model class and ridge
penalty fixed.

## Null ensemble

For every seed and hidden panel, double-edge swaps preserve exactly the feature and descriptor
degree sequences. The declared structural graph is unchanged. Each null is fitted separately; null
predictions are not selected, pooled during training or tuned against outcomes.

## Endpoints and adjudication

For each seed, the primary adaptive endpoints are repeated:

- participant-weighted RMSE improvement of declared structure over the randomized graph; and
- participant-weighted precision@3 improvement for top-decile hidden deviations.

Participant bootstrap intervals use 5,000 draws per seed. Robustness requires both lower 95%
interval bounds to exceed zero for every one of the 20 graph realizations, with exact degree
preservation in every hidden panel. The complete distribution is reported; no favorable seed may be
selected.

This sensitivity can show that the adaptive same-cohort result is not specific to one graph draw.
It remains post-result discovery evidence and is not external replication, temporal direction,
genetic conditioning, flux validation or clinical utility.
