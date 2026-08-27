# 06 — Precise mode groups junctions by cluster

**What to build:** Precise mode stops inferring which junctions belong to the same element
by stripping a numeric copy suffix off the IS name, and groups on the cluster column
instead.

The suffix it strips today was invented by the parser in the first place, and it gets the
reference genome wrong in both directions: the two `IS17` fragments collapse together
despite not aligning to each other at all, while their actual parent `ISAba12_1` stays
separate. They are the two opposite ends of one element with the middle 819 bp missing
from the assembly.

A table without a cluster column is a hard error naming the migration subcommand. That is
acceptable only because the remedy exists — without one it would be user-hostile — so this
ticket must land after the fallback path is reachable in the same release.

**Blocked by:** 04 — Similarity clusters computed and written to the table.

**Status:** done

- [x] Junction pairing groups on the cluster column; the copy-suffix regex is gone
- [x] A table with no cluster column fails fast with a message naming the migration subcommand
- [x] The junction-pairs golden is re-pinned, with the diff showing the `IS17`/`ISAba12` regrouping (see the comment below: two rows also move columns, as a consequence of it)
- [x] Frequency estimates for elements whose grouping is unchanged are numerically unchanged

## Comments

- The error names `ijump migrate-is-table`, in `is_table.MIGRATE_SUBCOMMAND`. Ticket 08
  builds the subcommand under that name.
- **Landed before 08.** Until 08 ships, the error names a remedy that does not exist yet.
  The release constraint stands: 08 has to land before this reaches a user.
- The golden diff carries one consequence beyond the relabelling. `IS17`'s two `NODE_1`
  junctions were the only members of their group and were all right-junctions, which sent
  `find_pairs` down its orphan-only branch — and that branch writes a right-junction into
  `Position_l` with a zero count, so both rows reported frequency 0.0. Grouped with
  `ISAba12` the group now has junctions on both sides, the ordinary path runs, and the two
  rows report real frequencies (0.0096 and 0.0026). The orphan-only branch is still wrong
  for any cluster whose junctions are all right-handed — `ISAba12` and `ISAba53` on
  `NODE_2` still take it. Filed as
  `.scratch/junction-pairing-orphans/issues/01-right-only-junctions-written-into-the-left-column.md`.
