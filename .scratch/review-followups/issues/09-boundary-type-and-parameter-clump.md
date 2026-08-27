# 09 — Proposal: give the boundary clump a type and shrink a 13-parameter function

**What to build:** *Proposal, not yet approved work.* The ungapped clipped-read
table builder takes thirteen parameters, mutates several of them in place, and
returns a three-tuple; the boundaries it works with are bare five-element lists
unpacked by position. Giving that clump a name would make the function's contract
readable and let the positional unpacking go away.

**This ticket needs a scope decision before anyone starts.** See *Open questions*.

## Why

The review flagged this as the largest **Data Clumps** / **Primitive Obsession**
instance in the refactor, and the only finding that is a genuine design change
rather than a defect.

`_crtable_ungapped` in `clipped_read_search.py` takes 13 parameters. Some are
mutated in place; others are returned in a 3-tuple — so a caller cannot tell
from the signature what is input, what is output, and what is both. Its
`boundaries` entries are five-element lists whose fields are recovered by
positional unpacking at the consuming site, which means the meaning of each slot
lives only in whoever wrote both ends.

A small `Boundary` type would fix the positional unpacking and shorten the
signature. The module **already has a `SearchResult` dataclass**, so this is
applying an established in-repo pattern rather than importing a new idea.

The clump is not isolated. The same shape recurs at the Circos writer (8
parameters), the frequency estimator (7), the average-mode report (7), and the
junction pairer (7) — which is what makes the scope question real.

## Why this is `needs-triage` and not `ready-for-agent`

Two reasons, and they are the point of this ticket:

1. **No test safety net.** `.scratch/isclipped-refactor/notes/no-test-safety-net.md`
   documents that this repo lacks the coverage to make a change of this shape
   safely. The characterization tests added during the refactor cover junction
   pairing and region summary — not the ungapped clipped-read path. Refactoring a
   13-parameter function with in-place mutation and no tests over it is how silent
   behaviour changes get made.
2. **It is a judgement call.** Every other ticket in this batch fixes something
   demonstrably wrong. This one changes a design that currently works. That is a
   maintainer's call about whether the readability is worth the risk, not an
   agent's.

## Scope

Provisional — the final scope is the triage decision below, so treat this as the
candidate rather than the commitment. Smallest defensible version:

- Introduce a `Boundary` type for the five-field boundary entries, following the
  module's existing `SearchResult` dataclass.
- Replace the positional unpacking at the consuming site with attribute access.
- Thread the type through `_crtable_ungapped`'s signature so the parameter count
  drops, without changing which values the function computes.

Anything beyond that — the other four clump sites, the in-place mutation, the
3-tuple return — is explicitly deferred to triage.

## Open questions for triage

- **Scope**: this one function only, or the recurring clump across all five sites?
  One function is a contained change; all five is a wide refactor and should be
  sequenced expand–contract, with the migration batched per call site so each
  batch stays green.
- **Test-first?** Given the note above, the defensible order is: characterization
  tests over the ungapped path first (as its own ticket, mirroring how
  isclipped-refactor tickets 01 and 06 preceded their extractions), then the
  type. Confirm whether that prerequisite ticket should be filed.
- **Does `Boundary` earn its keep on its own**, independent of the parameter
  count? The positional five-list unpacking is arguably the worse of the two
  problems and is fixable without touching the signature at all — a smaller,
  safer first step if the full change is judged too risky.
- **Mutation**: should the function stop mutating its arguments in place, or is
  that a separate change? Doing both at once is what makes this risky; doing only
  the type is much safer.

## Out of scope regardless of how triage goes

- Changing what the clipped-read search computes. Any version of this is strictly
  behaviour-preserving.
- The Circos writer's signature, if ticket 07 has not landed yet — coordinate,
  since 07 introduces a method that wraps that call.

## Verification (whatever scope is chosen)

- Characterization coverage over the affected path exists **before** the
  refactor, and passes unchanged after. If triage decides to proceed without it,
  that decision should be written down here explicitly.
- Output byte-identical on a real sample.
- No behaviour change in the in-place-mutation semantics unless that was
  explicitly chosen in triage.
- `ruff`/`mypy` clean on `src/ijump/`.

**Blocked by:** 05 and 07.

- **05** touches `clipped_read_search.py`'s module-level constants; this ticket
  restructures the same module, so the small change should land first.
- **07** introduces the wrapper around the Circos writer, whose signature is one
  of the clump sites this ticket might extend to.

Both are cheap. Neither should be a reason to delay the *triage decision*, only
the implementation.

**Status:** needs-triage

- [ ] Scope decided: this function only, or the recurring clump across all sites.
- [ ] Decided whether characterization tests are a prerequisite ticket.
- [ ] Decided whether `Boundary` lands independently of the signature change.
- [ ] Decided whether in-place mutation is in or out of scope.
- [ ] Ticket re-labelled `ready-for-agent` or `wontfix` once the above are answered.

## Comments
