"""Essais qui testent VRAIMENT eml_toolkit/core.py (import direct : s'il manque un nom, l'essai échoue bruyamment).
Lancer : python3 -m pytest eml_toolkit/test_core.py -q"""
import random

import mpmath as mp
import pytest

from eml_toolkit import core as c

TOL = mp.mpf(10) ** -40
rng = random.Random(2026)
BANDE = [mp.mpc(rng.uniform(-5, 5), rng.uniform(-3.1, 3.1)) for _ in range(200)] + [mp.mpc(0.5, -mp.pi), mp.mpc(3, 4) / 2]
PAIRES = [(mp.mpc(rng.uniform(-3, 3), rng.uniform(-3, 3)), mp.mpc(rng.uniform(-3, 3), rng.uniform(-3, 3))) for _ in range(200)]


def test_les_noms_attendus_existent():
    for nom in ("eml", "eml_star", "eml_exp", "eml_ln", "eml_zero", "eml_neg", "eml_sub", "eml_add", "eml_inv", "eml_mul",
                "conjugate_formula", "real_part", "imag_part", "modulus_squared", "modulus", "alt_conjugate",
                "alt_modulus_squared", "E_add", "E_mul", "E_conj", "fold_to_strip", "ONE"):
        assert hasattr(c, nom), nom


def test_arithmetique():
    assert c.eml_zero() == 0
    for a, b in PAIRES:
        assert abs(c.eml_mul(a, b) - a * b) < TOL * max(1, abs(a * b))          # partout sauf 0 (théorème B)
        assert abs(c.eml_exp(a) - mp.exp(a)) < TOL * max(1, abs(mp.exp(a)))
        assert abs(c.eml_inv(a) - 1 / a) < TOL * max(1, abs(1 / a))
        assert abs(c.eml_neg(b) + b) < TOL and abs(c.eml_sub(a, b) - (a - b)) < TOL   # |Im b| ≤ 3 < π
        assert abs(c.eml_add(a, b) - (a + b)) < TOL
        assert abs(c.E_add(a, b) - (a + b)) < TOL and abs(c.E_mul(a, b) - a * b) < TOL * max(1, abs(a * b))


def test_addition_hors_bande_est_enroulee():
    assert abs(c.eml_add(mp.mpc(1), mp.mpc(0, 5)) - (1 + (5 - 2 * mp.pi) * 1j)) < TOL    # théorème C


def test_miroir_dans_la_bande():
    for z in BANDE:
        assert abs(c.conjugate_formula(z) - mp.conj(z)) < TOL and abs(c.E_conj(z) - mp.conj(z)) < TOL
        assert abs(c.real_part(z) - z.real) < TOL and abs(c.imag_part(z) - z.imag) < TOL
        assert abs(c.modulus_squared(z) - abs(z) ** 2) < TOL * max(1, abs(z) ** 2)
        assert abs(c.modulus(z) - abs(z)) < TOL * max(1, abs(z))


def test_le_cas_qui_etait_faux():
    assert abs(c.modulus_squared(mp.mpc(3, 4) / 2) - mp.mpf(25) / 4) < TOL
    assert abs(c.modulus(mp.mpc(3, 2)) - mp.sqrt(13)) < TOL


def test_miroir_hors_bande_est_enroule():
    z = mp.mpc(1, 5)                                   # hors de [−π, π) : z̄ plus un tour, pas z̄
    assert abs(c.conjugate_formula(z) - (1 + (2 * mp.pi - 5) * 1j)) < TOL
    assert abs(c.real_part(z).real - 1) < TOL          # la partie réelle, elle, reste juste partout


def test_base_alternative_sans_condition():
    for z in BANDE + [mp.mpc(1, 5), mp.mpc(0, 100), mp.mpc(2, -4)]:
        assert abs(c.alt_conjugate(z) - mp.conj(z)) < TOL and abs(c.alt_modulus_squared(z) - abs(z) ** 2) < TOL * max(1, abs(z) ** 2)


def test_pliage():
    assert abs(c.fold_to_strip(mp.mpc(1, 5)) - mp.mpc(1, 5 - 2 * mp.pi)) < TOL
    assert c.fold_to_strip(mp.mpc(1, 0.3)) == mp.mpc(1, 0.3)


# ---------- bords de la bande du miroir (relecture du 2026-11-19) ----------

def test_repliement_explicite_aux_bords():
    pi = mp.pi
    for re in (mp.mpf("0.5"), mp.mpf(-2)):
        assert c.conjugate_folded(mp.mpc(re, pi)) == mp.mpc(re, pi)          # Im z = +π : Re z + iπ (ce n'est PAS conj z)
        assert c.conjugate_folded(mp.mpc(re, -pi)) == mp.mpc(re, pi)         # Im z = −π : Re z + iπ = conj z
    d = mp.mpf(10) ** -30
    assert abs(c.conjugate_folded(mp.mpc(1, pi - d)) - mp.mpc(1, -(pi - d))) < TOL       # juste dedans : conj z
    assert abs(c.conjugate_folded(mp.mpc(1, pi + d)) - mp.mpc(1, pi - d)) < TOL          # juste au-dessus : un tour ajouté
    assert abs(c.conjugate_folded(mp.mpc(1, -pi + d)) - mp.mpc(1, pi - d)) < TOL         # juste dedans : conj z
    assert abs(c.conjugate_folded(mp.mpc(1, -pi - d)) - mp.mpc(1, -(pi - d))) < TOL      # juste en dessous : un tour ôté


def test_l_arbre_brut_et_le_repliement_coincident_hors_des_bords():
    for z in BANDE + [mp.mpc(1, 5), mp.mpc(-2, 9.3), mp.mpc(0.3, -14), mp.mpc(2, -4)]:
        if abs(abs(z.imag) - mp.pi) < mp.mpf(10) ** -20:
            continue
        assert abs(c.conjugate_formula(z) - c.conjugate_folded(z)) < TOL, z


def test_dans_la_bande_le_repliement_est_le_conjugue():
    for z in BANDE:
        assert abs(c.conjugate_folded(z) - mp.conj(z)) < TOL
