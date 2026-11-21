ELECTRONIC SUPPLEMENTARY TABLE S5: FILE GUIDE

This archive contains two tab-delimited UTF-8 text files retained unchanged
from the frozen analytical delivery. No spreadsheet conversion or value editing
has been performed.

1. S5_candidates39_complete.tsv
   The 39 source-linked, direction-opposed sequence units. The file includes
   sequence-unit identifiers, representative-gene and source records, source
   GO provenance, local alignment summaries, the P0, P1, interaction, D0 and
   D1 effect estimates and adjusted q values, cross-route support, and
   descriptive normalized counts for the three pooled biological replicates in
   each of Pa, PaPs, PaHp, and PaPsHp.

2. S5_source_GO_relations.tsv
   Direct GO relations for the source entries linked to the candidate table.

All column names and categorical values are already in English. The following
status values are used in S5_candidates39_complete.tsv:

SOURCE_CHAIN_READY
  A frozen source-chain record was available for the sequence unit.
AT_LEAST_ONE_SOURCE_ENTRY_RETURNED
  At least one required source entry returned an annotation record.
SUPPORTED_ORIGINAL_OPPOSITE
  The specified sensitivity route supported the original opposite-direction
  pattern.
ORIGINAL_PATTERN_WITHOUT_FULL_SUPPORT
  The original pattern was retained, but the specified route did not meet the
  full-support criterion.
TRUE / FALSE
  Boolean full-support flags for the corresponding sensitivity route.

GO_evidence follows the source GO evidence codes (for example, IBA and IEA).
The expression audit note identifies the group values as descriptive normalized
counts from three pooled biological replicates; group means are not model
contrasts.
