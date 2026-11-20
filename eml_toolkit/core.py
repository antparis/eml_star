"""eml★ toolkit — basic functions (Monnerot 2026), built from eml (Odrzywołek 2026) and eml★.

Corrected on 2026-11-19 (original: core.py.bak_2026-11-19). The old version gave wrong results for
E_sub, E_add, E_mul and modulus_squared (modulus_squared(z) was equal to exp(z) − Log(exp(z)), i.e. exp(z) − z in the
strip |Im z| < π), and the test suite did not actually import it. Tests: eml_toolkit/test_core.py.

Domains of validity (measured and certified, see FINDINGS_20261119_branch_safety_certificate.md):
- eml_mul(x, y) = x·y for all nonzero x, y, with no condition on the imaginary parts;
- eml_add(x, y) = x + y if |Im y| < π (otherwise x + y up to an integer number of turns of 2πi); no condition on x;
- conjugate_formula(z) = conj(z) if −π ≤ Im z < π (Theorem 3.1); outside this strip, conj(z) up to turns of 2πi;
  real_part, imag_part, modulus_squared, modulus (and the old name E_conj) call this RAW tree: they inherit this
  strip and its numerical instability exactly on the edges Im z = ±π. The only function with an explicit
  convention that is guaranteed at the edges is conjugate_folded.
The constants go through the extended reals: ln 0 = −∞, exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞."""
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
    """Brings Im(z) back into [-π, π) by whole turns of 2π."""
    re = mp.re(z)
    im = mp.im(z)
    folded_im = mp.fmod(im + mp.pi, mp.mpf(2) * mp.pi) - mp.pi
    return mp.mpc(re, folded_im)


# ── Holomorphic trees (eml alone) ─────────────────────────────────────────────

def eml_exp(z):
    return eml(z, ONE)


def eml_ln(z):
    return eml(ONE, eml(eml(ONE, z), ONE))


def eml_zero():
    """eml(1, eml(eml(1, 1), 1)) = e − ln(e^e) = 0, exactly (in floating point one would get ~1e-60)."""
    return mp.mpf(0)


def eml_neg(z):
    """−z = eml(C, eml(z, 1)) with C = eml(1, eml(eml(1, 0), 1)) = −∞: exp(−∞) − ln(exp z). Valid for |Im z| < π."""
    return eml(eml(ONE, eml(eml(ONE, eml_zero()), ONE)), eml(z, ONE))


def eml_sub(a, b):
    return eml(eml_ln(a), eml_exp(b))


def eml_add(a, b):
    return eml_sub(a, eml_neg(b))


def eml_inv(z):
    return eml_exp(eml_neg(eml_ln(z)))


def eml_mul(a, b):
    return eml_exp(eml_add(eml_ln(a), eml_ln(b)))


# ── Anti-holomorphic part (eml★) ──────────────────────────────────────────────

def conjugate_formula(z):
    """RAW evaluation of the tree 1 − eml★(0, eml(z, 1)) (Theorem 3.1).

    Exact value: Log(exp(conj z)), that is, conj(z) brought back modulo 2πi into −π < Im ≤ π; equal to conj(z)
    if and only if −π ≤ Im z < π. In floating point, when Im z is equal to ±π (up to rounding), exp(z) lands
    a hair's breadth from the negative real axis and the logarithm may choose either side: unreliable result ON the
    edge. For a convention that is guaranteed at the edges, use conjugate_folded."""
    return ONE - eml_star(mp.mpc(0), eml(z, ONE))


def conjugate_folded(z):
    """Exact value of the mirror formula, computed by explicit folding (without exp or log):
    conj(z) brought back modulo 2πi into −π < Im ≤ π. Convention at the edges: Im z = +π and Im z = −π both give
    Re z + iπ. Coincides with conj(z) for −π ≤ Im z < π, and with conjugate_formula(z) away from the edges."""
    im = -mp.im(z)
    k = mp.ceil((im - mp.pi) / (2 * mp.pi))
    return mp.mpc(mp.re(z), im - 2 * mp.pi * k)


def real_part(z):
    """Re(z) = (z + conj(z)) / 2, with conj via eml★ (strip of Theorem 3.1)."""
    return (z + conjugate_formula(z)) / 2


def imag_part(z):
    """Im(z) = (z − conj(z)) / (2i), with conj via eml★ (strip of Theorem 3.1)."""
    return (z - conjugate_formula(z)) / (2 * mp.j)


def modulus_squared(z):
    """|z|² = z · conj(z), with conj via eml★ (strip of Theorem 3.1)."""
    return z * conjugate_formula(z)


def modulus(z):
    return mp.sqrt(modulus_squared(z))


# ── Alternative basis {eml, Re, 1}: no strip condition ────────────────────────

def alt_conjugate(z):
    return 2 * mp.re(z) - z


def alt_modulus_squared(z):
    return mp.re(z) ** 2 + mp.im(z) ** 2


# Old names, kept for compatibility
E_neg, E_sub, E_add, E_mul, E_conj = eml_neg, eml_sub, eml_add, eml_mul, conjugate_formula

__all__ = ['eml', 'eml_star', 'fold_to_strip', 'eml_exp', 'eml_ln', 'eml_zero', 'eml_neg', 'eml_sub', 'eml_add',
           'eml_inv', 'eml_mul', 'conjugate_formula', 'conjugate_folded', 'real_part', 'imag_part', 'modulus_squared', 'modulus',
           'alt_conjugate', 'alt_modulus_squared', 'E_neg', 'E_sub', 'E_add', 'E_mul', 'E_conj', 'ONE']
