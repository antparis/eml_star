"""Essais du certificat de branch-safety (certify_branch_safety.py).
Lancer : python3 -m pytest test_certify_branch_safety.py -q   (ou l'environnement de ~/cutlattice : .venv/bin/python -m pytest)"""
import mpmath as mp
import pytest
import certify_branch_safety as cb

X, Y, ONE = cb.X, cb.Y, cb.ONE
eml = cb.eml


def B(a, b, c, d):
    return cb.Boite(a, b, c, d)


# ---------- constantes (aucune variable) : valeurs exactes sur les réels étendus ----------

def test_zero_est_exact():
    assert cb.constante(cb.E_zero()) == 0


def test_la_chaine_de_E_neg_passe_par_plus_et_moins_l_infini():
    a = eml(ONE, cb.E_zero())                    # e − ln 0
    assert cb.constante(a) == cb.PLUS_INF
    assert cb.constante(eml(a, ONE)) == cb.PLUS_INF      # exp(+∞) − ln 1
    assert cb.constante(eml(ONE, eml(a, ONE))) == cb.MOINS_INF  # e − ln(+∞)


def test_constante_refuse_l_indefini():
    a = eml(ONE, cb.E_zero())                    # +∞
    with pytest.raises(cb.Indefini):
        cb.constante(eml(a, a))                  # exp(+∞) − ln(+∞) = ∞ − ∞


# ---------- les arbres calculent bien ce qu'ils disent ----------

@pytest.mark.parametrize("a,b", [(mp.mpc(1.3, 0.1), mp.mpc(0.7, 0.2)), (mp.mpc(0.4, -0.5), mp.mpc(2.5, 0.3)), (mp.mpc(0.7, 1e-12), mp.mpc(1.3))])
def test_arbres_justes_en_un_point(a, b):
    assert abs(cb.valeur_point(cb.E_add(X, Y), a, b) - (a + b)) < mp.mpf(10) ** -30
    assert abs(cb.valeur_point(cb.E_mul(X, Y), a, b) - a * b) < mp.mpf(10) ** -30


# ---------- la règle des deux portes ----------

GAUCHE = B(-2, -1, -0.1, 0.1)                  # à cheval sur l'axe réel négatif, loin de 0


def test_traversee_en_porte_x_est_absorbee():
    t = eml(eml(ONE, X), ONE)                  # le nœud eml(1, x) entre par la porte x : exp(e − ln x) = e^e / x
    r = cb.certifier(t, {"x": GAUCHE}, profondeur=4)
    assert r["certifie"] and r["absorbees"] >= 1
    v = cb.enclos(t, {"x": GAUCHE.point(-1.5, 0.05)})
    attendu = mp.exp(mp.e) / mp.mpc(-1.5, 0.05)
    assert v.real.a <= attendu.real <= v.real.b and v.imag.a <= attendu.imag <= v.imag.b


def test_traversee_en_porte_y_est_refusee():
    t = eml(ONE, eml(ONE, X))                  # le nœud eml(1, x) entre par la porte y : son log est repris, le saut reste
    r = cb.certifier(t, {"x": GAUCHE}, profondeur=4)
    assert not r["certifie"] and "porte y" in r["raison"]


def test_traversee_a_la_racine_est_refusee():
    r = cb.certifier(eml(ONE, X), {"x": GAUCHE}, profondeur=4)
    assert not r["certifie"] and "racine" in r["raison"]


def test_log_de_zero_variable_est_refuse():
    r = cb.certifier(eml(eml(ONE, X), ONE), {"x": B(-0.1, 0.1, -0.1, 0.1)}, profondeur=3)
    assert not r["certifie"] and "zéro" in r["raison"]


def test_rotation_garde_le_demi_tour():
    """ln(−y) + iπ : écart 0 ou 2πi avec le log principal ; sans le + iπ le signe serait faux."""
    for y in (mp.mpc(-1.5, 0.05), mp.mpc(-1.5, -0.05)):
        k = (cb.log_tourne_point(y) - mp.log(y)) / (2j * mp.pi)
        assert abs(k - mp.nint(k.real)) < 1e-40 and abs(mp.exp(cb.log_tourne_point(y)) - y) < 1e-40


# ---------- les vrais arbres ----------

def test_addition_certifiee():
    r = cb.certifier(cb.E_add(X, Y), {"x": B(0.5, 2, -0.3, 0.3), "y": B(0.5, 2, -0.3, 0.3)}, profondeur=6)
    assert r["certifie"], r


def test_multiplication_echoue_pres_de_x_egal_1_sans_disque():
    r = cb.certifier(cb.E_mul(X, Y), {"x": B(0.9, 1.1, -0.1, 0.1), "y": B(1.1, 1.5, -0.1, 0.1)}, profondeur=5, budget=60)
    assert not r["certifie"]


def test_multiplication_certifiee_hors_du_disque_autour_de_1():
    r = cb.certifier(cb.E_mul(X, Y), {"x": B(0.5, 1.5, -0.3, 0.3), "y": B(1.1, 1.5, -0.1, 0.1)}, profondeur=8, budget=240,
                     disque_x=(1, 0.1))
    assert r["certifie"], r
    assert r["exclues"] >= 1 and r["absorbees"] >= 1        # on a bien traversé (x réel < 1) et bien exclu le disque


def test_le_disque_ne_cache_pas_tout():
    with pytest.raises(ValueError):
        cb.certifier(cb.E_mul(X, Y), {"x": B(0.95, 1.05, -0.05, 0.05), "y": B(1.1, 1.5, -0.1, 0.1)}, disque_x=(1, 0.2))


# ---------- le disque autour de x = 1 : lemme exp(E_ln(w)) = w, puis le programme vérifie tout le reste ----------

@pytest.mark.parametrize("w", [mp.mpc(-0.3, 1e-12), mp.mpc(-0.3, -1e-12), mp.mpc(2, 0.5), mp.mpc(-4, 0), mp.mpc(1e-8, -1e-9)])
def test_lemme_exp_E_ln_est_l_identite(w):
    assert abs(mp.exp(cb.valeur_point(cb.E_ln(X), w)) - w) < mp.mpf(10) ** -35


@pytest.mark.parametrize("a,b", [(mp.mpc(1 + 1e-6, 1e-7), mp.mpc(1.3, 0.1)), (mp.mpc(0.97, -0.02), mp.mpc(1.2, 0)), (mp.mpc(1.03, 0.04), mp.mpc(1.4, -0.05))])
def test_la_version_lemme_calcule_le_meme_produit(a, b):
    assert abs(cb.valeur_point(cb.E_mul_lemme(X, Y), a, b) - a * b) < mp.mpf(10) ** -30
    assert abs(cb.valeur_point(cb.E_mul_lemme(X, Y), a, b) - cb.valeur_point(cb.E_mul(X, Y), a, b)) < mp.mpf(10) ** -30


def test_le_disque_est_certifie_avec_le_lemme():
    r = cb.certifier(cb.E_mul_lemme(X, Y), {"x": B(0.9, 1.1, -0.1, 0.1), "y": B(1.1, 1.5, -0.1, 0.1)}, profondeur=6, budget=120)
    assert r["certifie"], r


def test_moins_ln_une_valeur_qui_saute_n_est_permise_qu_en_porte_x():
    t = cb.moins_ln(cb.eml(ONE, X), ONE)        # a = e − ln x saute de 2πi pour x réel négatif
    r = cb.certifier(t, {"x": GAUCHE}, profondeur=3)
    assert not r["certifie"] and "racine" in r["raison"]                          # à la racine : refusé
    r = cb.certifier(cb.E_exp(t), {"x": GAUCHE}, profondeur=3)
    assert r["certifie"] and r["absorbees"] >= 1                                  # sous une exponentielle : effacé
    r = cb.certifier(cb.eml(ONE, t), {"x": GAUCHE}, profondeur=3)
    assert not r["certifie"]                                                      # sous un log : refusé


def test_produit_certifie_pour_x_negatif_avec_le_lemme():
    r = cb.certifier(cb.E_mul_lemme(X, Y), {"x": B(-3, -0.1, -1, 1), "y": B(0.1, 3, -1, 1)}, profondeur=7, budget=200)
    assert r["certifie"] and r["absorbees"] >= 1, r


# ---------- rigueur des bornes et des constantes (relecture du 2026-11-19) ----------

def test_une_borne_decimale_en_texte_est_arrondie_vers_l_exterieur():
    b = B("0.1", "0.3", "-0.1", "0.7")
    dixieme, trois, sept = mp.mpf(1) / 10, mp.mpf(3) / 10, mp.mpf(7) / 10
    with mp.workdps(80):
        assert b.r[0] <= mp.mpf(1) / 10 <= mp.mpf(3) / 10 <= b.r[1] and b.i[0] <= -mp.mpf(1) / 10 and mp.mpf(7) / 10 <= b.i[1]
    assert B(0.5, 2, -1, 1).r == (mp.mpf(0.5), mp.mpf(2))            # un nombre Python est pris tel quel


def test_le_disque_exclu_est_garanti_dans_le_disque_annonce():
    assert B(0.97, 1.03, -0.03, 0.03).dans_disque(1, "0.05")
    assert not B(0.96, 1.04, -0.03, 0.03).dans_disque(1, "0.05")     # coin à distance exactement 0,05 : pas STRICTEMENT dedans
    assert not B(0.97, 1.03, -0.045, 0.045).dans_disque(1, "0.05")


def test_l_intervalle_d_une_constante_contient_sa_valeur():
    t = cb.eml(cb.eml(ONE, ONE), ONE)                                 # e^e ≈ 15,15 : fini, non nul
    z = cb._intervalle_constant(t)
    with mp.workdps(90):
        v = mp.exp(mp.e)
        assert mp.mpf(z.real._mpi_[0]) < v < mp.mpf(z.real._mpi_[1])
    assert mp.mpf(z.real._mpi_[1]) - mp.mpf(z.real._mpi_[0]) > 0       # largeur non nulle


def test_les_arbres_plus_et_fois_n_utilisent_aucune_constante_finie_non_nulle():
    cb.CONSTANTES_FINIES[0] = 0
    d = {"x": B(1.2, 1.4, 0.1, 0.2), "y": B(0.6, 0.8, 0.1, 0.2)}
    for t in (cb.E_add(X, Y), cb.E_mul(X, Y), cb.E_mul_lemme(X, Y), cb.E_add_lemme(X, Y)):
        assert cb.certifier(t, d, profondeur=4)["certifie"]
    assert cb.CONSTANTES_FINIES[0] == 0


def test_tout_le_domaine_dans_le_disque_est_refuse():
    with pytest.raises(ValueError):
        cb.certifier(cb.E_exp(X), {"x": B(0.99, 1.01, -0.01, 0.01)}, disque_x=(1, "0.5"))
