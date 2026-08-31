# Post-result human graph-null mixing diagnostic protocol

**Version:** 1.0  
**Freeze date:** 30 August 2026, after the graph-null and human sensitivity results  
**Role:** null-quality diagnostic, not independent biological confirmation

## Objective

Verify that the degree-preserving lipid feature–descriptor nulls move materially away from the
observed incidence graph and from one another. Exact degree preservation alone does not establish
adequate randomization.

## Frozen diagnostic

For ST002081 and ST000818, regenerate the 20 reported null seeds independently for each of five
hidden panels using the primary ten attempted swaps per edge. For every null, report exact row and
column degree preservation, the fraction of observed edges replaced and edge-set Jaccard similarity
to the observed graph. Within each panel, report all pairwise Jaccard similarities among the 20
nulls and the number of unique edge sets.

The diagnostic passes when every degree sequence is exact, all 20 nulls are unique within every
panel, every null replaces more than 50% of observed edges, and the maximum pairwise null Jaccard
similarity is below 0.80. Complete rows are retained.

This post-result audit tests graph movement and diversity. It does not prove uniform sampling from
all binary matrices with the degree sequence, convert the adaptive analysis into prospective
confirmation, or establish physiology, genetics or clinical utility.
