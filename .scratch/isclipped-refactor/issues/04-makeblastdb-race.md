# 04 — makeblastdb is never waited on

Status: ready-for-agent
Blocked by: —

## The bug

`ijump.py:50-55`:

```python
def check_blast_ref(ref_name, ref_file):
    if os.path.isfile(ref_name + '.nsq'):
        pass
    else:
        makeblastdb_command = f'makeblastdb -in {ref_file} -dbtype nucl -out {ref_name}'
        makeblastdb = subprocess.Popen(makeblastdb_command.split(), stdout=subprocess.PIPE)
```

`Popen` is never awaited — no `.wait()`, no `.communicate()`. Execution continues
immediately while `makeblastdb` is still writing the database.

`runblast` (isclipped.py:458, called at ijump.py:194) then invokes `blastn` against that
database. On a first run against a fresh reference, `blastn` can start before
`makeblastdb` has finished. The larger the reference, the wider the window.

Also note: `stdout=subprocess.PIPE` with nothing ever reading the pipe. If `makeblastdb`
writes enough to fill the OS pipe buffer it will block indefinitely.

## Why it is rarely seen

Only fires when `<ref>.nsq` is absent — the first run against a given reference. On repeat
runs line 51 short-circuits. That is also why it survives in a codebase people use daily.

## Scope

- Wait for completion before returning. `subprocess.run(..., check=True)` handles the wait,
  the pipe, and the exit status in one call.
- Fail loudly if `makeblastdb` returns non-zero or is not installed
  (`FileNotFoundError`) — currently a missing BLAST+ surfaces as a confusing `blastn`
  failure much later.
- Consider verifying the database exists after the call rather than trusting the exit code.

## Verification

Awkward to test directly without BLAST+ installed, and ticket 01 deliberately commits
`tiny.nsq` so the suite needs no external tooling. Options, in preference order:

1. Extract the command construction from its execution and unit-test that the runner is
   invoked with a waiting call — no BLAST+ needed.
2. Integration test marked `skipif(shutil.which('makeblastdb') is None)`, pointed at a
   fixture directory with the `.nsq` deliberately removed.

Option 1 is the smaller change and keeps the suite tool-free.

## Comments

Implemented in bd7bfae (`fix: Wait for makeblastdb to finish before returning`): the
`makeblastdb` call now waits on completion before `check_blast_ref` returns.