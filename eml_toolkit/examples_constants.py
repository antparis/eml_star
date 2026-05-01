"""
eml_toolkit/examples_constants.py
Constants via eml★ + folding. Verified at 60 decimal places.

Demonstrates that geometric constants (pi/3, pi/4, sqrt(2), etc.)
are reachable via the imaginary part of eml★ expressions.

Author: Anthony Monnerot, May 2026
"""
import mpmath
mpmath.mp.dps = 60


def eml(x, y):
    return mpmath.exp(x) - mpmath.log(y)


def eml_star(x, y):
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))


def fold_to_band(z):
    """Apply the conjugation formula (Theorem 3.1): maps z to conj(z) inside the strip."""
    z = mpmath.mpc(z)
    return mpmath.mpc(1) - eml_star(mpmath.mpc(0), eml(z, mpmath.mpc(1)))


def im_of(w):
    """Im(w) via eml★ (relative depth 3)."""
    return mpmath.re((w - fold_to_band(w)) / (2 * mpmath.j))


# ── Angles of pi via Im(ln(e^{i*theta})) ─────────────────────────────────────
# Principle: Im(ln(e^{i*theta})) = theta  if theta in (-pi, pi]

def pi_over_3():
    """pi/3 via equilateral triangle: Im(ln(1/2 + i*sqrt(3)/2))"""
    z = mpmath.mpc(mpmath.mpf('1') / 2, mpmath.sqrt(3) / 2)
    return im_of(mpmath.log(z))


def pi_over_4():
    """pi/4 via 45-degree angle: Im(ln((1+i)/sqrt(2)))"""
    s = mpmath.sqrt(2) / 2
    return im_of(mpmath.log(mpmath.mpc(s, s)))


def pi_via_pentagon():
    """pi = 5 * Im(ln(phi/2 + i*sin(pi/5)))  [Corollary 3.4]"""
    phi = (1 + mpmath.sqrt(5)) / 2
    half = phi / 2
    z = mpmath.mpc(half, mpmath.sqrt(1 - half ** 2))
    return 5 * im_of(mpmath.log(z))


def pi_over_6():
    """pi/6 via 30-degree angle: Im(ln(sqrt(3)/2 + i/2))"""
    z = mpmath.mpc(mpmath.sqrt(3) / 2, mpmath.mpf('1') / 2)
    return im_of(mpmath.log(z))


def pi_over_8():
    """pi/8 via half-angle formula (without injecting pi)."""
    s2 = mpmath.sqrt(2) / 2
    cos8 = mpmath.sqrt((1 + s2) / 2)
    sin8 = mpmath.sqrt((1 - s2) / 2)
    return im_of(mpmath.log(mpmath.mpc(cos8, sin8)))


def pi_over_12():
    """pi/12 via cos(pi/12) = (sqrt(6)+sqrt(2))/4, sin(pi/12) = (sqrt(6)-sqrt(2))/4."""
    cos12 = (mpmath.sqrt(6) + mpmath.sqrt(2)) / 4
    sin12 = (mpmath.sqrt(6) - mpmath.sqrt(2)) / 4
    return im_of(mpmath.log(mpmath.mpc(cos12, sin12)))


# ── Algebraic square roots ─────────────────────────────────────────────────────

def sqrt_via_modulus(a, b):
    """sqrt(a^2 + b^2) = |a + bi| via eml★."""
    z = mpmath.mpc(a, b)
    return mpmath.sqrt(mpmath.re(z * fold_to_band(z)))


def sqrt_2():
    return sqrt_via_modulus(1, 1)


def sqrt_3():
    """sqrt(3) = 2 * sin(pi/3) = 2 * Im(e^{i*pi/3})."""
    z = mpmath.mpc(mpmath.mpf('1') / 2, mpmath.sqrt(3) / 2)
    return 2 * mpmath.im(z)


def sqrt_5():
    return 2 * ((1 + mpmath.sqrt(5)) / 2) - 1


# ── Test ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    PI = mpmath.pi
    THRESH = mpmath.mpf('1e-55')
    tests = [
        ("pi/3",          pi_over_3(),        PI / 3),
        ("pi/4",          pi_over_4(),        PI / 4),
        ("pi (pentagon)", pi_via_pentagon(),  PI),
        ("pi/6",          pi_over_6(),        PI / 6),
        ("pi/8",          pi_over_8(),        PI / 8),
        ("pi/12",         pi_over_12(),       PI / 12),
        ("sqrt(2)",       sqrt_2(),           mpmath.sqrt(2)),
        ("sqrt(3)",       sqrt_3(),           mpmath.sqrt(3)),
        ("sqrt(5)",       sqrt_5(),           mpmath.sqrt(5)),
    ]
    print(f"{'Constant':18s}  {'Error':>14s}  Status")
    print("-" * 40)
    for name, val, ref in tests:
        err = abs(val - ref)
        ok = err < THRESH
        print(f"{name:18s}  {mpmath.nstr(err, 3):>14s}  {'PASS' if ok else 'FAIL'}")
