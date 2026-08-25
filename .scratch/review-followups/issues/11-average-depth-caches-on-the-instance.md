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

**Status:** ready-for-agent

- [ ] Repeated calls for one region are still served from a cache, not recomputed
- [ ] A discarded `ISClipped` becomes collectable, with its alignment handle and tables
- [ ] The `# noqa: B019` and the note explaining it are gone, and ruff passes without them
- [ ] The end-to-end goldens are unchanged, in both estimation modes
