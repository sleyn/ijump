"""Regression test for ticket 03: 'IS pos' must land on the same coordinate
base as 'Position' in precise-mode ijump_junctions.txt output.

check_junctions_presence() (now isclipped.py, moved there by ticket 07) used
to compare est_mode against the misspelled literal 'presice', so the
'IS pos' -> 1-base conversion never ran. Drives the function directly
against a constructed junctions DataFrame, since the tiny.bam fixture used
elsewhere in this suite produces no clipped reads and never reaches
junction output (see ticket 01/02 fixture notes).
"""
import pandas as pd
import pytest

from ijump.isclipped import EstimationMode, check_junctions_presence


def _junctions_df():
    return pd.DataFrame({
        'IS pos': [10],
        'Position': [20],
    })


def test_precise_mode_converts_is_pos_to_1_base(tmp_path):
    check_junctions_presence(_junctions_df(), str(tmp_path), EstimationMode.PRECISE)

    written = pd.read_csv(tmp_path / 'ijump_junctions.txt', sep='\t')
    assert written.loc[0, 'IS pos'] == 11
    assert written.loc[0, 'Position'] == 21


def test_precise_mode_keeps_is_pos_and_position_on_the_same_base(tmp_path):
    junc_tbl = _junctions_df()
    original_gap = junc_tbl.loc[0, 'Position'] - junc_tbl.loc[0, 'IS pos']

    check_junctions_presence(junc_tbl, str(tmp_path), EstimationMode.PRECISE)

    written = pd.read_csv(tmp_path / 'ijump_junctions.txt', sep='\t')
    assert written.loc[0, 'Position'] - written.loc[0, 'IS pos'] == original_gap
