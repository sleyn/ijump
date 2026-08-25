# 11 — `average_depth` caches on the instance

**What to build:** The coverage cache stops holding its pipeline object alive, and the
`# noqa: B019` suppression standing in for that goes away.

`ISClipped.average_depth` is decorated with `@lru_cache`. Caching on an instance method
keeps a reference to `self` for as long as the cache lives, so every `ISClipped` ever built
stays reachable — with its alignment handle, its junction tables and its depth
dictionaries. One run does not care; anything that builds many instances leaks all of them.

Ruff's bugbear rule catches this (`B019`). It was found during the ruff/mypy configuration
work and deliberately left alone there, on the grounds that fixing it properly is a design
change rather than a lint fix, and that ticket's scope was tooling configuration only. It
was suppressed with a `noqa` and an inline comment saying so. No ticket followed, so the
suppression has been carrying the decision ever since.

The caching itself is worth keeping — depth over a region is recomputed often enough that
the cache is doing real work — so this is about where the cache lives, not whether it
exists. Whatever replaces it has to keep `average_depth`'s current results exactly: the
end-to-end goldens pin depth columns in both estimation modes, so a change in what is cached
must not be a change in what is computed.

**Blocked by:** None — can start immediately.

**Status:** done

- [x] Repeated calls for one region are still served from a cache, not recomputed
- [x] A discarded `ISClipped` becomes collectable, with its alignment handle and tables
- [x] The `# noqa: B019` and the note explaining it are gone, and ruff passes without them
- [x] The end-to-end goldens are unchanged, in both estimation modes

## Comments

The cache is a plain dict on the instance now, keyed by `(chrom, start, stop)`, and the
`# noqa: B019` is gone. `average_depth` looks the answer up or delegates to
`_measure_average_depth`, which is the old body unchanged.

**A dict rather than a per-instance `lru_cache`.** The obvious fix — building an
`lru_cache` around the bound method in `__init__` — would stop the *class* holding every
pipeline, but the cache would then hold the bound method, which holds `self`: a reference
cycle, collectable only when the garbage collector next runs rather than when the last real
reference goes. A dict of coordinates to floats refers to nothing, so there is no cycle at
all. The test pins that distinction: it disables the collector before dropping the pipeline,
so it fails for a cycle exactly as it failed for the class-level cache.

**Unbounded, deliberately.** The `lru_cache` held 128 entries, which turned out to be far
too few to be doing the job it was there for: Circos measures every annotated region — 7670
of them on the test genome — so an entry was evicted long before a second caller asked for it
again. One float per region is nothing beside the alignment already open.

**Four tests**, covering both halves of the ticket: a region is measured once, different
regions are measured separately, two pipelines do not answer from each other's cache, and a
discarded pipeline is collectable without a collector pass. The third passed before this
change — `lru_cache` keyed on `self`, so instances never got each other's answers — and is
kept as a guard on the replacement, which could plausibly have got that wrong.

**Verification.** `ruff` and `mypy` clean; full suite 202 passed, 4 skipped; the end-to-end
goldens are byte-identical in both estimation modes, so the depth columns did not move.
