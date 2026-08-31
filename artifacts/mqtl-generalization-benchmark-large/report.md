# METSIM 71-state topology benchmark

**Status:** `71_state_positive_control_topology_benchmark_complete`  
**States:** 71  
**Direct topology recovered:** 26/71  
**Matched permutation p:** 1e-05  
**Global degree-matched permutation p:** 1e-05  
**Within-biochemical_class permutation p:** 1e-05  
**Independent causal-gene validation:** false  

## Null distributions

| null_id | permutations | mean_recoveries | median_recoveries | q95_recoveries | q99_recoveries | maximum_recoveries | draws_at_least_observed | plus_one_p_value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| within_panel | 100000 | 2.07096 | 2.0 | 4.0 | 6.0 | 10 | 0 | 9.99990000099999e-06 |
| global_degree_matched | 100000 | 1.36779 | 1.0 | 3.0 | 4.0 | 7 | 0 | 9.99990000099999e-06 |
| within_biochemical_class | 100000 | 5.7921 | 6.0 | 9.0 | 10.0 | 16 | 0 | 9.99990000099999e-06 |

## Biochemical-class recovery

| biochemical_class | states | unique_genes | direct_recoveries | recovery_fraction |
| --- | --- | --- | --- | --- |
| AA | 40 | 40 | 19 | 0.475 |
| CA | 3 | 3 | 0 | 0.0 |
| CV | 5 | 5 | 1 | 0.2 |
| LI | 9 | 9 | 0 | 0.0 |
| NU | 10 | 10 | 5 | 0.5 |
| PE | 2 | 2 | 0 | 0.0 |
| XE | 2 | 2 | 1 | 0.5 |

These 71 states exhaust the eligible significant METSIM signals under the frozen one-state-per-nominated-gene and exact Human-GEM metabolite-mapping rules. The within-panel label, biochemical-class-stratified, and gene-fixed global metabolite-degree-matched nulls test topology enrichment, but the causal-gene nominations remain knowledge based. This is a broad positive-control audit, not independent causal-gene validation, directional flux prediction, or target nomination.
