"""
eml_toolkit/core.py
Core operators for the eml★ system.

eml(x, y)  = exp(x) - ln(y)        [Odrzywołek 2026, arXiv:2603.21852v2]
eml★(x, y) = exp(x) - ln(conj(y))  [Monnerot 2026]

Reference: Zenodo DOI:10.5281/zenodo.19183008
"""
import mpmath

# Default precision: 50 decimal places
mpmath.mp.dps = 50


# ── Core operators ─────────────────────────────────────────────────────────────

def eml(x, y):
    """
    EML Sheffer operator: eml(x, y) = exp(x) - ln(y).
    Generates all classical elementary functions when composed with constant 1.
    Holomorphic by construction.
    """
    return mpmath.exp(x) - mpmath.log(y)


def eml_star(x, y):
    """
    eml★ anti-holomorphic extension: eml★(x, y) = exp(x) - ln(conj(y)).
    Injects anti-holomorphic structure into the Sheffer basis.
    Minimal change: single conjugate on the second argument.
    """
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))


# ── EML arithmetic trees (source: eml_compiler_v4.py, Zenodo) ─────────────────

ONE = mpmath.mpc(1)


def eml_exp(z):
    """exp(z) via EML tree, K=1."""
    return eml(z, ONE)


def eml_ln(z):
    """ln(z) via EML tree, K=3."""
    return eml(ONE, eml(eml(ONE, z), ONE))


def eml_zero():
    """Constant 0 via EML tree, K=3."""
    return eml(ONE, eml(eml(ONE, ONE), ONE))


def eml_neg(z):
    """Negation -z via EML tree, K=8."""
    return eml(eml(ONE, eml(eml(ONE, eml_zero()), ONE)), eml(z, ONE))


def eml_sub(a, b):
    """Subtraction a - b via EML tree, K=5."""
    return eml(eml_ln(a), eml_exp(b))


def eml_add(a, b):
    """Addition a + b via EML tree, K=13."""
    return eml_sub(a, eml_neg(b))


def eml_inv(z):
    """Reciprocal 1/z via EML tree, K=12."""
    return eml_exp(eml_neg(eml_ln(z)))


def eml_mul(a, b):
    """Multiplication a * b via EML tree, K=20."""
    return eml_exp(eml_add(eml_ln(a), eml_ln(b)))


def eml_div(a, b):
    """Division a / b via EML tree."""
    return eml_mul(a, eml_inv(b))


def eml_pow(a, b):
    """Power a^b via EML tree."""
    return eml_exp(eml_mul(b, eml_ln(a)))


# ── eml★ formulas (Monnerot 2026) ─────────────────────────────────────────────

def conjugate_formula(z):
    """
    Complex conjugate via eml★ (Theorem 3.1).

    z_bar = 1 - eml★(0, eml(z, 1))

    Valid for Im(z) in (-pi, pi).
    Returns z_bar with error < 10^-50 for Im(z) strictly inside the strip.
    Outside the strip: error = n * 2*pi*i for integer n (branch jump).
    """
    return 1 - eml_star(mpmath.mpc(0), eml(z, ONE))


def real_part(z):
    """
    Real part Re(z) via eml★, relative depth 3.
    Re(z) = (z + conj(z)) / 2
    Valid for Im(z) in (-pi, pi).
    """
    return (z + conjugate_formula(z)) / 2


def imag_part(z):
    """
    Imaginary part Im(z) via eml★, relative depth 3.
    Im(z) = (z - conj(z)) / (2i)
    Valid for Im(z) in (-pi, pi).
    """
    return (z - conjugate_formula(z)) / (2 * mpmath.j)


def modulus_squared(z):
    """
    |z|^2 via eml★, relative depth 3.
    |z|^2 = z * conj(z)
    Valid for Im(z) in (-pi, pi).
    """
    return z * conjugate_formula(z)


def modulus(z):
    """
    |z| via eml★, relative depth 4.
    |z| = sqrt(z * conj(z))
    Valid for Im(z) in (-pi, pi).
    """
    return mpmath.sqrt(modulus_squared(z))


# ── Alternative: {eml, Re, 1} basis (Grok's proposal) ────────────────────────

def alt_conjugate(z):
    """
    Complex conjugate via {eml, Re, 1} basis.
    z_bar = 2 * Re(z) - z
    Unconditional (no strip restriction).
    """
    return 2 * mpmath.re(z) - z


def alt_modulus_squared(z):
    """
    |z|^2 via {eml, Re, 1} basis.
    |z|^2 = Re(z)^2 + Im(z)^2
    Unconditional.
    """
    return mpmath.re(z) ** 2 + mpmath.im(z) ** 2


# ── Utility ────────────────────────────────────────────────────────────────────

def is_in_strip(z, margin=1e-10):
    """Check whether Im(z) is strictly inside (-pi, pi)."""
    return abs(float(mpmath.im(z))) < float(mpmath.pi) - margin


if __name__ == "__main__":
    print("eml_toolkit/core.py — sanity check")
    z = mpmath.mpc(1.5, 0.8)
    print(f"  z          = {z}")
    print(f"  conj(z)    = {conjugate_formula(z)}")
    print(f"  Re(z)      = {real_part(z)}")
    print(f"  Im(z)      = {imag_part(z)}")
    print(f"  |z|^2      = {modulus_squared(z)}")
    print(f"  |z|        = {modulus(z)}")
    print(f"  Expected:    conj={mpmath.conj(z)}, Re={mpmath.re(z)}, Im={mpmath.im(z)}")
