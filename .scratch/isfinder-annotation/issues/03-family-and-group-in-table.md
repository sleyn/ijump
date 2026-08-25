# 03 — IS table carries family and group

**What to build:** The IS table gains the annotation that the pipeline already has in hand
and currently discards, and moves to a self-describing format.

An ISFinder BLAST subject id is `name_family_group` — `ISAba18_IS3_IS51`, `ISAba1_IS4_IS10`.
The parser today reduces it to the name by splitting on the first underscore, throwing away
family and group at zero saving and truncating the 11 database entries whose *name* itself
contains an underscore (`ISBj2_B_IS5_IS5` becomes `ISBj2`). Parsing from the right recovers
all three fields and fixes those names.

The table becomes a headered TSV so columns can be added without breaking readers. Existing
headerless 4-column tables are sniffed and still accepted, so nobody's working setup stops
running when they upgrade.

This is why the reader and writer land together: splitting them leaves a table format that
one half of the codebase writes and the other cannot read. The Circos writer's strict
three-value unpack of a table row breaks on any added column, so fixing it belongs to this
change, not to a follow-up.

**Blocked by:** 01 — Characterization goldens for the current IS-annotation path.

**Status:** ready-for-agent

- [ ] Subject ids are parsed from the right into name, family and group; the 11 underscore-bearing database names survive intact
- [ ] `isfinder-db-parse` emits a headered TSV carrying name, contig, coordinates, family, group and percent identity
- [ ] `ijump run` reads the headered TSV in both estimation modes
- [ ] Legacy headerless 4-column tables are detected and still accepted, with the added columns empty
- [ ] Circos input generation works against the widened table and produces output identical to the golden
- [ ] Goldens re-pinned in the same commit, with the diff showing only the intended format and annotation changes
