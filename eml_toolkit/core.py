"""eml★ toolkit — fonctions de base (Monnerot 2026), à partir de eml (Odrzywołek 2026) et de eml★.

Corrigé le 2026-11-19 (original : core.py.bak_2026-11-19). L'ancienne version donnait des résultats faux pour
E_sub, E_add, E_mul et modulus_squared (modulus_squared(z) valait exp(z) − Log(exp(z)), soit exp(z) − z dans la
bande |Im z| < π), et la suite d'essais ne
l'importait pas réellement. Essais : eml_toolkit/test_core.py.

Domaines de validité (mesurés et certifiés, voir FINDINGS_20261119_branch_safety_certificate.md) :
- eml_mul(x, y) = x·y pour tous x, y non nuls, sans condition sur les parties imaginaires ;
- eml_add(x, y) = x + y si |Im y| < π (sinon x + y à un nombre entier de tours 2πi près) ; aucune condition sur x ;
- conjugate_formula(z) = conj(z) si −π ≤ Im z < π (Théorème 3.1) ; hors de cette bande, conj(z) à des tours 2πi près ;
  real_part, imag_part, modulus_squared, modulus (et l'ancien nom E_conj) appellent cet arbre BRUT : ils héritent de
  cette bande et de son instabilité numérique exactement sur les bords Im z = ±π. La seule fonction à convention
  explicite et garantie aux bords est conjugate_folded.
Les constantes passent par les réels étendus : ln 0 = −∞, exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞."""
import mpmath as mp

mp.mp.dps = 60

ONE = mp.mpf(1)


def eml(x, y):
    """eml(x, y) = exp(x) - ln(y)"""
    return mp.exp(x) - mp.log(y)


def eml_star(x, y):
    """eml★(x, y) = exp(x) - ln(conj(y))"""
    return mp.exp(x) - mp.log(mp.conj(y))


def fold_to_strip(z):
    """Ramène Im(z) dans [-π, π) par des tours entiers de 2π."""
    re = mp.re(z)
    im = mp.im(z)
    folded_im = mp.fmod(im + mp.pi, mp.mpf(2) * mp.pi) - mp.pi
    return mp.mpc(re, folded_im)


# ── Arbres holomorphes (eml seul) ─────────────────────────────────────────────

def eml_exp(z):
    return eml(z, ONE)


def eml_ln(z):
    return eml(ONE, eml(eml(ONE, z), ONE))


def eml_zero():
    """eml(1, eml(eml(1, 1), 1)) = e − ln(e^e) = 0, exactement (en virgule flottante on obtiendrait ~1e-60)."""
    return mp.mpf(0)


def eml_neg(z):
    """−z = eml(C, eml(z, 1)) avec C = eml(1, eml(eml(1, 0), 1)) = −∞ : exp(−∞) − ln(exp z). Valide pour |Im z| < π."""
    return eml(eml(ONE, eml(eml(ONE, eml_zero()), ONE)), eml(z, ONE))


def eml_sub(a, b):
    return eml(eml_ln(a), eml_exp(b))


def eml_add(a, b):
    return eml_sub(a, eml_neg(b))


def eml_inv(z):
    return eml_exp(eml_neg(eml_ln(z)))


def eml_mul(a, b):
    return eml_exp(eml_add(eml_ln(a), eml_ln(b)))


# ── Partie anti-holomorphe (eml★) ─────────────────────────────────────────────

def conjugate_formula(z):
    """Évaluation BRUTE de l'arbre 1 − eml★(0, eml(z, 1)) (Théorème 3.1).

    Valeur exacte : Log(exp(conj z)), c'est-à-dire conj(z) ramené modulo 2πi dans −π < Im ≤ π ; égale à conj(z)
    si et seulement si −π ≤ Im z < π. En virgule flottante, quand Im z vaut ±π (à l'arrondi près), exp(z) tombe à
    un cheveu de l'axe réel négatif et le logarithme peut choisir l'un ou l'autre côté : résultat non fiable SUR le
    bord. Pour une convention garantie aux bords, utiliser conjugate_folded."""
    return ONE - eml_star(mp.mpc(0), eml(z, ONE))


def conjugate_folded(z):
    """Valeur exacte de la formule du miroir, calculée par repliement explicite (sans exp ni log) :
    conj(z) ramené modulo 2πi dans −π < Im ≤ π. Convention aux bords : Im z = +π et Im z = −π donnent tous deux
    Re z + iπ. Coïncide avec conj(z) pour −π ≤ Im z < π, et avec conjugate_formula(z) hors des bords."""
    im = -mp.im(z)
    k = mp.ceil((im - mp.pi) / (2 * mp.pi))
    return mp.mpc(mp.re(z), im - 2 * mp.pi * k)


def real_part(z):
    """Re(z) = (z + conj(z)) / 2, avec conj par eml★ (bande du Théorème 3.1)."""
    return (z + conjugate_formula(z)) / 2


def imag_part(z):
    """Im(z) = (z − conj(z)) / (2i), avec conj par eml★ (bande du Théorème 3.1)."""
    return (z - conjugate_formula(z)) / (2 * mp.j)


def modulus_squared(z):
    """|z|² = z · conj(z), avec conj par eml★ (bande du Théorème 3.1)."""
    return z * conjugate_formula(z)


def modulus(z):
    return mp.sqrt(modulus_squared(z))


# ── Base alternative {eml, Re, 1} : sans condition de bande ───────────────────

def alt_conjugate(z):
    return 2 * mp.re(z) - z


def alt_modulus_squared(z):
    return mp.re(z) ** 2 + mp.im(z) ** 2


# Anciens noms, gardés pour compatibilité
E_neg, E_sub, E_add, E_mul, E_conj = eml_neg, eml_sub, eml_add, eml_mul, conjugate_formula

__all__ = ['eml', 'eml_star', 'fold_to_strip', 'eml_exp', 'eml_ln', 'eml_zero', 'eml_neg', 'eml_sub', 'eml_add',
           'eml_inv', 'eml_mul', 'conjugate_formula', 'conjugate_folded', 'real_part', 'imag_part', 'modulus_squared', 'modulus',
           'alt_conjugate', 'alt_modulus_squared', 'E_neg', 'E_sub', 'E_add', 'E_mul', 'E_conj', 'ONE']
