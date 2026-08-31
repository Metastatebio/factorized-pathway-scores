# Post-result CCLE reaction-cluster robustness protocol

**Version:** 1.0  
**Freeze date:** 30 August 2026, after the property-matched null ensemble  
**Role:** dependence robustness analysis, not untouched confirmation

## Objective

Test whether target-bootstrap inference is falsely precise because related CCLE metabolites share
reaction systems. The analysis uses the expected random-minus-topology target effects from the
20-seed property-matched null ensemble; no model is selected or refitted.

## Frozen clustering rule

Each target is assigned to the HumanGEM subsystem occurring most often among its eligible direct
candidate rows. Ties are resolved alphabetically. This gives one deterministic cluster per target
while retaining all targets. The complete assignment is exported.

## Inference and gate

Ten thousand bootstrap draws resample subsystem clusters with replacement and retain all target
effects within each sampled cluster. We also recompute the mean effect after leaving out each
subsystem in turn. Robustness requires at least ten clusters, a positive cluster-bootstrap lower
bound and a positive minimum leave-one-cluster-out effect.

Subsystem assignment is an imperfect dependence proxy because HumanGEM pathways overlap. This
analysis tests concentration of the result across broad reaction systems; it does not make targets
independent or establish causal reaction direction, flux, drug response or clinical utility.
