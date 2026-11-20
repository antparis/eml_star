"""Tests that REALLY test eml_toolkit/core.py (direct import: if a name is missing, the test fails loudly).
Run: python3 -m pytest eml_toolkit/test_core.py -q"""
import random

import mpmath as mp
import pytest

from eml_toolkit import core as c

TOL = mp.mpf(10) ** -40
rng = random.Random(2026)
STRIP = [mp.mpc(rng.uniform(-5, 5), rng.uniform(-3.1, 3.1)) for _ in range(200)] + [mp.mpc(0.5, -mp.pi), mp.mpc(3, 4) / 2]
PAIRS = [(mp.mpc(rng.uniform(-3, 3), rng.uniform(-3, 3)), mp.mpc(rng.uniform(-3, 3), rng.uniform(-3, 3))) for _ in range(200)]


def test_expected_names_exist():
    for name in ("eml", "eml_star", "eml_exp", "eml_ln", "eml_zero", "eml_neg", "eml_sub", "eml_add", "eml_inv", "eml_mul",
                 "conjugate_formula", "real_part", "imag_part", "modulus_squared", "modulus", "alt_conjugate",
                 "alt_modulus_squared", "E_add", "E_mul", "E_conj", "fold_to_strip", "ONE"):
        assert hasattr(c, name), name


def test_arithmetic():
    assert c.eml_zero() == 0
    for a, b in PAIRS:
        assert abs(c.eml_mul(a, b) - a * b) < TOL * max(1, abs(a * b))          # everywhere except 0 (Theorem B)
        assert abs(c.eml_exp(a) - mp.exp(a)) < TOL * max(1, abs(mp.exp(a)))
        assert abs(c.eml_inv(a) - 1 / a) < TOL * max(1, abs(1 / a))
        assert abs(c.eml_neg(b) + b) < TOL and abs(c.eml_sub(a, b) - (a - b)) < TOL   # |Im b| ≤ 3 < π
        assert abs(c.eml_add(a, b) - (a + b)) < TOL
        assert abs(c.E_add(a, b) - (a + b)) < TOL and abs(c.E_mul(a, b) - a * b) < TOL * max(1, abs(a * b))


def test_addition_outside_strip_is_wrapped():
    assert abs(c.eml_add(mp.mpc(1), mp.mpc(0, 5)) - (1 + (5 - 2 * mp.pi) * 1j)) < TOL    # Theorem C


def test_mirror_in_the_strip():
    for z in STRIP:
        assert abs(c.conjugate_formula(z) - mp.conj(z)) < TOL and abs(c.E_conj(z) - mp.conj(z)) < TOL
        assert abs(c.real_part(z) - z.real) < TOL and abs(c.imag_part(z) - z.imag) < TOL
        assert abs(c.modulus_squared(z) - abs(z) ** 2) < TOL * max(1, abs(z) ** 2)
        assert abs(c.modulus(z) - abs(z)) < TOL * max(1, abs(z))


def test_the_case_that_was_wrong():
    assert abs(c.modulus_squared(mp.mpc(3, 4) / 2) - mp.mpf(25) / 4) < TOL
    assert abs(c.modulus(mp.mpc(3, 2)) - mp.sqrt(13)) < TOL


def test_mirror_outside_strip_is_wrapped():
    z = mp.mpc(1, 5)                                   # outside [−π, π): z̄ plus one turn, not z̄
    assert abs(c.conjugate_formula(z) - (1 + (2 * mp.pi - 5) * 1j)) < TOL
    assert abs(c.real_part(z).real - 1) < TOL          # the real part, for its part, stays correct everywhere


def test_alternative_basis_unconditional():
    for z in STRIP + [mp.mpc(1, 5), mp.mpc(0, 100), mp.mpc(2, -4)]:
        assert abs(c.alt_conjugate(z) - mp.conj(z)) < TOL and abs(c.alt_modulus_squared(z) - abs(z) ** 2) < TOL * max(1, abs(z) ** 2)


def test_fold_to_strip():
    assert abs(c.fold_to_strip(mp.mpc(1, 5)) - mp.mpc(1, 5 - 2 * mp.pi)) < TOL
    assert c.fold_to_strip(mp.mpc(1, 0.3)) == mp.mpc(1, 0.3)


# ---------- edges of the mirror strip (review of 2026-11-19) ----------

def test_explicit_folding_at_the_edges():
    pi = mp.pi
    for re in (mp.mpf("0.5"), mp.mpf(-2)):
        assert c.conjugate_folded(mp.mpc(re, pi)) == mp.mpc(re, pi)          # Im z = +π: Re z + iπ (this is NOT conj z)
        assert c.conjugate_folded(mp.mpc(re, -pi)) == mp.mpc(re, pi)         # Im z = −π: Re z + iπ = conj z
    d = mp.mpf(10) ** -30
    assert abs(c.conjugate_folded(mp.mpc(1, pi - d)) - mp.mpc(1, -(pi - d))) < TOL       # just inside: conj z
    assert abs(c.conjugate_folded(mp.mpc(1, pi + d)) - mp.mpc(1, pi - d)) < TOL          # just above: one turn added
    assert abs(c.conjugate_folded(mp.mpc(1, -pi + d)) - mp.mpc(1, pi - d)) < TOL         # just inside: conj z
    assert abs(c.conjugate_folded(mp.mpc(1, -pi - d)) - mp.mpc(1, -(pi - d))) < TOL      # just below: one turn removed


def test_raw_tree_and_folding_coincide_away_from_the_edges():
    for z in STRIP + [mp.mpc(1, 5), mp.mpc(-2, 9.3), mp.mpc(0.3, -14), mp.mpc(2, -4)]:
        if abs(abs(z.imag) - mp.pi) < mp.mpf(10) ** -20:
            continue
        assert abs(c.conjugate_formula(z) - c.conjugate_folded(z)) < TOL, z


def test_in_the_strip_the_folding_is_the_conjugate():
    for z in STRIP:
        assert abs(c.conjugate_folded(z) - mp.conj(z)) < TOL
