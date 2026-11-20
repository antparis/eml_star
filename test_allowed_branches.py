"""Regression tests of the branch mutant (review of 2026-11-19). Run: python3 -m pytest test_allowed_branches.py -q"""
import mpmath as mp
import allowed_branches as ba

TURN = 2j * mp.pi
SMALL = mp.mpf(10) ** -40


def test_twenty_logarithms_eight_of_them_with_a_variable_argument():
    _, logs = ba.evaluate(ba.E_mul(ba.X, ba.Y), mp.mpc(0.7, 0.3), mp.mpc(-1.3, 0.4))
    assert len(logs) == 20 and [n for n, variable, _ in logs if variable] == [0, 2, 3, 5, 12, 14, 16, 18]


def test_the_deepest_log_1_of_E_zero_changes_nothing():
    assert abs(ba.effect_of_a_shift(6)) < SMALL             # no. 6: e − 2πi, exp(e − 2πi) = exp(e), E_zero stays 0


def test_the_log_1_of_the_root_shifts_the_result_by_one_turn():
    assert abs(ba.effect_of_a_shift(19) + TURN) < SMALL     # no. 19: exactly −2πi


def test_the_eight_variable_logs_one_by_one():
    for n in (0, 2, 3, 5, 12, 14, 16, 18):
        assert abs(ba.effect_of_a_shift(n)) < SMALL, n


def test_some_logs_of_constants_break_the_result():
    assert abs(ba.effect_of_a_shift(1)) > 1e-3 and abs(ba.effect_of_a_shift(7)) > 1e-3
