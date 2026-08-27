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

**Status:** done

- [x] Subject ids are parsed from the right into name, family and group; the 11 underscore-bearing database names survive intact
- [x] `isfinder-db-parse` emits a headered TSV carrying name, contig, coordinates, family, group and percent identity
- [x] `ijump run` reads the headered TSV in both estimation modes
- [x] Legacy headerless 4-column tables are detected and still accepted, with the added columns empty
- [x] Circos input generation works against the widened table and produces output identical to the golden
- [x] Goldens re-pinned in the same commit, with the diff showing only the intended format and annotation changes

## Comments

- The format lives in one new module, `src/ijump/is_table.py`: `parse_subject_id`,
  `write_is_table`, `read_is_table`. `isfinder_db_parcer.py` writes through it and
  `ISClipped.iscollect` reads through it, so neither end spells the format out.
- Columns are the spec's `is_name, contig, start, stop, family, group, pident`. A column
  the reader does not know about (`cluster`, from ticket 04) is passed through rather than
  dropped, and a headered table missing an annotation column reads it as empty.
- Legacy tables are recognised by the *absence* of the header row, matched on the whole
  leading run of column names — so a four-column table whose first element is called
  `is_name` is still read as data, not eaten as a header.
- Circos: `ISClipped.is_coords` is now derived from the table and deliberately kept to
  three fields, and `circos.py` slices the coordinate triple instead of unpacking the whole
  entry. Either alone would have held; both means a later widening of `is_coords` cannot
  re-break it.
- The parser golden is re-pinned in the commit: its diff is the header row plus the three
  new columns, with no name or coordinate changed. The fixture BLAST output happens to
  contain no underscore-bearing database name, so that fix is covered by unit tests rather
  than by the golden. The end-to-end goldens for both modes, Circos files included, were
  re-run and are byte-identical — no re-pin needed.
- `CONTEXT.md` gains **IS table** as the name of this file, and the README section is
  renamed to match; it had been "mobile elements coordinates file" while the code called it
  three other things.
