# Triage Labels

The skills speak in terms of five canonical triage roles. This file maps those roles to the actual label strings used in this repo's issue tracker.

Because this repo uses a local-markdown tracker, a "label" is the value of the `Status:` line near the top of an issue file.

| Label in mattpocock/skills | Label in our tracker | Meaning                                  |
| -------------------------- | -------------------- | ---------------------------------------- |
| `needs-triage`             | `needs-triage`       | Maintainer needs to evaluate this issue  |
| `needs-info`               | `needs-info`         | Waiting on reporter for more information |
| `ready-for-agent`          | `ready-for-agent`    | Fully specified, ready for an AFK agent  |
| `ready-for-human`          | `ready-for-human`    | Requires human implementation            |
| `wontfix`                  | `wontfix`            | Will not be actioned                     |
| *(none — local addition)*  | `done`                | Implemented, verified, and merged; not one of the five upstream roles, but the tracker's actual closed state | 

`done` may carry a trailing qualifier (`done, partially`, `done, mostly`, or a note like `done — superseded by 12, 13, 14`) when the ticket wasn't closed cleanly; the qualifier is free text, not a separate label. `resolved` (used once, in `isclipped-refactor/15`) is a synonym for `done` from before this convention settled — new tickets should use `done`.

When a skill mentions a role (e.g. "apply the AFK-ready triage label"), use the corresponding label string from this table.

Edit the right-hand column to match whatever vocabulary you actually use.
