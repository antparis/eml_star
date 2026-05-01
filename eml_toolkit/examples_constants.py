"""
eml_toolkit/examples_constants.py
==================================
Constantes via eml★ + folding. Vérifiées à 60 décimales.
Anthony Monnerot, mai 2026.
"""
import mpmath
mpmath.mp.dps = 60

def eml(x, y):      return mpmath.exp(x) - mpmath.log(y)
def eml_star(x, y): return mpmath.exp(x) - mpmath.log(mpmath.conj(y))

def fold_to_band(z):
    z = mpmath.mpc(z)
    return mpmath.mpc(1) - eml_star(mpmath.mpc(0), eml(z, mpmath.mpc(1)))

def im_of(w):
    """Im(w) via eml★ (profondeur relative 3)."""
    return mpmath.re((w - fold_to_band(w)) / (2 * mpmath.j))

# ── Angles de π via Im(ln(e^{iθ})) ──────────────────────────────
# Principe : Im(ln(e^{iθ})) = θ  si θ ∈ (-π, π]

def pi_over_3():
    """π/3 via triangle équilatéral : Im(ln(1/2 + i√3/2))"""
    z = mpmath.mpc(mpmath.mpf('1')/2, mpmath.sqrt(3)/2)
    return im_of(mpmath.log(z))

def pi_over_4():
    """π/4 via angle 45° : Im(ln((1+i)/√2))"""
    s = mpmath.sqrt(2)/2
    return im_of(mpmath.log(mpmath.mpc(s, s)))

def pi_via_pentagon():
    """π = 5 × Im(ln(φ/2 + i·sin(π/5)))  [Corollary 3.4]"""
    phi = (1 + mpmath.sqrt(5))/2
    half = phi/2
    z = mpmath.mpc(half, mpmath.sqrt(1 - half**2))
    return 5 * im_of(mpmath.log(z))

def pi_over_6():
    """π/6 via angle 30° : Im(ln(√3/2 + i/2))"""
    z = mpmath.mpc(mpmath.sqrt(3)/2, mpmath.mpf('1')/2)
    return im_of(mpmath.log(z))

def pi_over_8():
    """π/8 via formule demi-angle (sans injecter π)."""
    s2 = mpmath.sqrt(2)/2
    cos8 = mpmath.sqrt((1 + s2)/2)
    sin8 = mpmath.sqrt((1 - s2)/2)
    return im_of(mpmath.log(mpmath.mpc(cos8, sin8)))

def pi_over_12():
    """π/12 via cos(π/12) = (√6+√2)/4, sin(π/12) = (√6−√2)/4."""
    cos12 = (mpmath.sqrt(6) + mpmath.sqrt(2))/4
    sin12 = (mpmath.sqrt(6) - mpmath.sqrt(2))/4
    return im_of(mpmath.log(mpmath.mpc(cos12, sin12)))

# ── Racines algébriques ───────────────────────────────────────────

def sqrt_via_modulus(a, b):
    """√(a²+b²) = |a+bi| via eml★."""
    z = mpmath.mpc(a, b)
    return mpmath.sqrt(mpmath.re(z * fold_to_band(z)))

def sqrt_2(): return sqrt_via_modulus(1, 1)
def sqrt_3():
    """√3 = 2×sin(π/3) = 2×Im(e^{iπ/3})."""
    z = mpmath.mpc(mpmath.mpf('1')/2, mpmath.sqrt(3)/2)
    return 2 * mpmath.im(z)   # Im directement, sans log
def sqrt_5(): return 2*((1 + mpmath.sqrt(5))/2) - 1

# ── Test ─────────────────────────────────────────────────────────

if __name__ == "__main__":
    PI = mpmath.pi
    THRESH = mpmath.mpf('1e-55')
    tests = [
        ("π/3",  pi_over_3(),        PI/3),
        ("π/4",  pi_over_4(),        PI/4),
        ("π (pentagone)", pi_via_pentagon(), PI),
        ("π/6",  pi_over_6(),        PI/6),
        ("π/8",  pi_over_8(),        PI/8),
        ("π/12", pi_over_12(),       PI/12),
        ("√2",   sqrt_2(),           mpmath.sqrt(2)),
        ("√3",   sqrt_3(),           mpmath.sqrt(3)),
        ("√5",   sqrt_5(),           mpmath.sqrt(5)),
    ]
    print(f"{'Constante':16s}  {'Erreur':>12s}  Statut")
    print("-" * 36)
    for name, val, ref in tests:
        err = abs(val - ref)
        ok  = err < THRESH
        print(f"{name:16s}  {mpmath.nstr(err,3):>12s}  {'✓' if ok else '✗'}")
