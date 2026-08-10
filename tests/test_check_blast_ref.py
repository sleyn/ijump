"""Unit tests for check_blast_ref (ticket 04).

makeblastdb was previously launched with Popen and never waited on, so
blastn could race a still-writing database, and a filled stdout pipe could
hang forever. These tests exercise the fix without needing BLAST+
installed: makeblastdb_command is tested directly for correct argv
construction, and check_blast_ref's calls to subprocess.run are verified
with a stub instead of shelling out.
"""
import subprocess

import pytest

import ijump


def test_makeblastdb_command_builds_expected_argv():
    argv = ijump.makeblastdb_command('/path/to/ref', '/path/to/ref.fna')
    assert argv == ['makeblastdb', '-in', '/path/to/ref.fna', '-dbtype', 'nucl', '-out', '/path/to/ref']


def test_check_blast_ref_skips_makeblastdb_when_nsq_exists(tmp_path, monkeypatch):
    ref_name = str(tmp_path / "ref")
    (tmp_path / "ref.nsq").touch()

    called = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: called.append((a, k)))

    ijump.check_blast_ref(ref_name, "ref.fna")

    assert called == []


def test_check_blast_ref_waits_on_makeblastdb(tmp_path, monkeypatch):
    ref_name = str(tmp_path / "ref")

    calls = []

    def fake_run(argv, check, stdout, stderr):
        calls.append((argv, check))
        (tmp_path / "ref.nsq").touch()
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    ijump.check_blast_ref(ref_name, "ref.fna")

    assert len(calls) == 1
    argv, check = calls[0]
    assert argv == ijump.makeblastdb_command(ref_name, "ref.fna")
    assert check is True


def test_check_blast_ref_raises_on_missing_makeblastdb(tmp_path, monkeypatch):
    ref_name = str(tmp_path / "ref")

    def fake_run(*a, **k):
        raise FileNotFoundError("makeblastdb not found")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(FileNotFoundError):
        ijump.check_blast_ref(ref_name, "ref.fna")


def test_check_blast_ref_raises_on_nonzero_exit(tmp_path, monkeypatch):
    ref_name = str(tmp_path / "ref")

    def fake_run(argv, check, stdout, stderr):
        raise subprocess.CalledProcessError(1, argv, output=b"", stderr=b"bad input file")

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(subprocess.CalledProcessError):
        ijump.check_blast_ref(ref_name, "ref.fna")


def test_check_blast_ref_raises_if_nsq_missing_after_success(tmp_path, monkeypatch):
    ref_name = str(tmp_path / "ref")

    def fake_run(argv, check, stdout, stderr):
        # Simulate a makeblastdb that reports success but didn't write the .nsq.
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)

    with pytest.raises(RuntimeError):
        ijump.check_blast_ref(ref_name, "ref.fna")
