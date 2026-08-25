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

**Status:** ready-for-agent

- [ ] Whether a junction at position 0 is reachable is established and recorded, with the evidence
- [ ] Both exits of `find_pairs` agree on what an unpaired junction at position 0 produces
- [ ] The chosen rule is stated where the sentinel is written, not left to be re-derived from the two branches
- [ ] Test coverage pins the agreed behaviour on both exits
- [ ] Any movement in the end-to-end goldens is re-pinned; an empty diff is the expected outcome if position 0 proves unreachable
