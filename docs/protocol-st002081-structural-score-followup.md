# Adaptive protocol: structure-resolved ST002081 pathway scores

**Protocol version:** 1.0  
**Freeze date:** 30 August 2026, after opening the coarse-family benchmark  
**Role:** explicitly adaptive mechanistic-resolution follow-up

## Trigger

The frozen coarse biochemical-family score reconstructed held-out markers better than the
population mean but failed against size-matched random grouping. This follow-up asks whether the
failure reflects insufficient biochemical resolution rather than absence of structured signal.
It is not an independent replication and cannot rescue the original gate retrospectively.

## Structural representation

Every eligible lipid name is parsed without outcome access into the following descriptors:

- headgroup/class;
- exact total carbon count;
- exact total double-bond count;
- headgroup-by-total-unsaturation state;
- observed acyl-chain identities; and
- observed chain carbon and chain unsaturation states.

For multi-chain phospholipids, total carbon and unsaturation are summed across named chains. For
triacylglycerols, the declared total composition and explicitly reported fatty-acid component are
used. Descriptors occurring in fewer than five eligible markers or more than 90% of markers are
excluded.

Within each outer training fold, visible markers are standardized and every descriptor score is
the median of its visible member markers. The target marker never contributes to its own score.

## Degree-preserving null

For each hidden panel, the visible feature-by-descriptor bipartite graph is randomized by double
edge swaps. This preserves exactly:

- the number of descriptors assigned to every visible feature; and
- the number of visible features assigned to every descriptor.

The random score therefore has the same dimension and incidence degrees as the structural score.
At least ten attempted swaps per incidence edge are used with a deterministic seed.

## Validation and endpoints

The participant folds, five disjoint hidden panels, alpha, priority definition and participant
bootstrap are identical to the frozen coarse benchmark. Models are population mean,
degree-preserving random structural scores and declared structural scores. The already generated
coarse-family and all-visible predictions are imported only for paired descriptive comparison.

An adaptive structure-resolution signal requires both the RMSE and precision@3 participant-
bootstrap intervals for structural versus degree-preserving random scores to be wholly above zero.
If only RMSE passes, the result supports reconstruction but not priority ranking. If neither
passes, biochemical descriptors do not validate pathway-specific compression in this cohort.

No result is confirmatory evidence of temporal direction, genetics, flux or clinical utility.
