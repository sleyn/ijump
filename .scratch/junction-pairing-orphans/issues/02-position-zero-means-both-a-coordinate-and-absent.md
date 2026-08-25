# 02 — Position 0 means both "first base" and "no partner"

**What to build:** One rule for what an unpaired junction looks like, applied on both of
`find_pairs`' exits, so a junction at the first base of a contig is not confused with a
missing one.

`find_pairs` uses `0` in `Position_l`/`Position_r` to mean "this junction has no partner".
Positions are 0-based, so `0` is also a real coordinate — the first base of a contig. The
two meanings collide, and the two exits disagree about what to do:

- The main path ends with a filter dropping every row where both positions are zero. A
  genuine junction at position 0 with no partner is indistinguishable from a blank row in
  the pre-sized frame, so it is dropped.
- The early return does not filter at all — it pre-sizes the frame to exactly the number of
  junctions and fills every row, so there are no blanks to clear. A genuine junction at
  position 0 survives there.

Same junction, two answers, depending on whether the other side of the contig happened to
have any junctions. Neither exit is obviously the right one: the question underneath is
whether `0` should be carrying the sentinel at all, or whether absence wants its own
representation the coordinate space cannot collide with.

Reachability is unconfirmed. It was noticed by reading the branch, not from a failing run,
and no fixture produces it. Origin-spanning elements are where to look first: an element
the assembler cut at a contig boundary has a span starting at the contig's first base
(isfinder-annotation 05), which is where a junction at position 0 would come from. Settle
whether it is reachable before deciding how much to change — if it is not, saying so in a
comment on the sentinel is a legitimate outcome.

**Blocked by:** 01 — Right-only junctions are written into the left column. Both
tickets rewrite the same early-return arm, so they are sequenced rather than concurrent.

**Status:** done

- [x] Whether a junction at position 0 is reachable is established and recorded, with the evidence
- [x] Both exits of `find_pairs` agree on what an unpaired junction at position 0 produces
- [x] The chosen rule is stated where the sentinel is written, not left to be re-derived from the two branches
- [x] Test coverage pins the agreed behaviour on both exits
- [x] Any movement in the end-to-end goldens is re-pinned; an empty diff is the expected outcome if position 0 proves unreachable

## Comments

**Reachable — established, not assumed.** `_clboundaries` reads a read's junction off
`pysam`'s `get_reference_positions`, which is 0-based, so a left-clipped read whose aligned
part starts at the contig's first base reports its junction at 0. Pinned as
`test_a_read_clipped_at_the_first_base_of_a_contig_yields_position_zero`. That is the
origin-spanning shape the ticket pointed at, and the reference genome comes within one base
of producing it: `IS17_2` starts at base 2 of `NODE_2`.

**It was worse than the ticket described.** The one-sided exit does keep such a row, but
every consumer downstream tests the sentinel too, so the row survives with its evidence
zeroed and reports frequency 0.0 — the same failure mode as ticket 01, reached by a
different route. Agreeing the two exits without fixing that would have shipped a phantom
row in place of a dropped one.

**The rule.** Absence is `junction_pairing.NO_JUNCTION`, outside the 0-based coordinate
space. The written file is untouched: it is 1-based, where 0 is not a position, so
`convert_zero_one_base` spells absence 0 there and a 0-based position 0 becomes 1.
Consumers updated: `interpos_distance`, `keep_pair`, `count_depth_unclipped`,
`convert_zero_one_base`, and the six sentinel tests in `frequency_estimation`.

Two consequences beyond the pairing itself, both previously silent:

- `keep_pair` compares each side of a pair against the regions of interest, and regions can
  begin at 0 (`max(x - 5, 0)`). An absent side spelled 0 fell inside such a region and kept
  the pair on evidence that was not there.
- The main path's trailing "remove empty rows" filter selected by value — keep rows with a
  position above zero — which is what dropped a real junction at 0. The frame is filled
  from the front, so the used rows are a leading slice; it now takes that slice.

**Goldens unchanged, and the empty diff is checked rather than assumed.** The e2e tier
passes byte-identical against the existing goldens. The lowest junction position in that
sample is 319, so it has neither a junction at 0 nor a region beginning at 0 — the sample
does not exercise either path. Reachability is a property of the code, not of this sample.
