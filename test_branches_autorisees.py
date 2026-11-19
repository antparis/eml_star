"""Régressions du mutant de branche (relecture du 2026-11-19). Lancer : python3 -m pytest test_branches_autorisees.py -q"""
import mpmath as mp
import branches_autorisees as ba

TOUR = 2j * mp.pi
PETIT = mp.mpf(10) ** -40


def test_vingt_logarithmes_dont_huit_a_argument_variable():
    _, logs = ba.evalue(ba.E_mul(ba.X, ba.Y), mp.mpc(0.7, 0.3), mp.mpc(-1.3, 0.4))
    assert len(logs) == 20 and [n for n, variable, _ in logs if variable] == [0, 2, 3, 5, 12, 14, 16, 18]


def test_le_log_1_le_plus_profond_de_E_zero_ne_change_rien():
    assert abs(ba.effet_d_un_decalage(6)) < PETIT           # n° 6 : e − 2πi, exp(e − 2πi) = exp(e), E_zero reste 0


def test_le_log_1_de_la_racine_decale_le_resultat_d_un_tour():
    assert abs(ba.effet_d_un_decalage(19) + TOUR) < PETIT   # n° 19 : exactement −2πi


def test_les_huit_logs_variables_un_par_un():
    for n in (0, 2, 3, 5, 12, 14, 16, 18):
        assert abs(ba.effet_d_un_decalage(n)) < PETIT, n


def test_certains_logs_de_constantes_cassent_le_resultat():
    assert abs(ba.effet_d_un_decalage(1)) > 1e-3 and abs(ba.effet_d_un_decalage(7)) > 1e-3
