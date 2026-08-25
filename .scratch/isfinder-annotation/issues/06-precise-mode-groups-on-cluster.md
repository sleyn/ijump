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

**Status:** ready-for-agent

- [ ] Junction pairing groups on the cluster column; the copy-suffix regex is gone
- [ ] A table with no cluster column fails fast with a message naming the migration subcommand
- [ ] The junction-pairs golden is re-pinned, with the diff showing the `IS17`/`ISAba12` regrouping and nothing else
- [ ] Frequency estimates for elements whose grouping is unchanged are numerically unchanged
