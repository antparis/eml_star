"""
verify_theorem4.py
==================
Ancillary file — Monnerot (2026), arXiv comment on 2603.21852v2.

Theorem 3.2 (half-open domain [−π, π)):
  - Interior (−π, π)  : verified NUMERICALLY to 50 decimal digits.
  - Boundary Im(z)=−π : verified SYMBOLICALLY (exact algebra). Numerically
    unstable: exp(−iπ) has tiny positive imaginary part in floating-point,
    causing conj(exp(z)) to land below the branch cut → spurious 2πi jump.
  - Boundary Im(z)=+π : fails mathematically and numerically.
  - Exterior           : fails with jump = 2πi.

Requirements: mpmath >= 1.3.0, sympy >= 1.12
"""

import sys
try:
    import mpmath
    import sympy as sp
except ImportError as e:
    sys.exit(f"Missing: {e}. pip install mpmath sympy")

mpmath.mp.dps = 50
THRESHOLD = mpmath.mpf("1e-40")
TWO_PI = 2 * mpmath.pi

def eml(x, y):      return mpmath.exp(x) - mpmath.log(y)
def eml_star(x, y): return mpmath.exp(x) - mpmath.log(mpmath.conj(y))
def formula(z):     return 1 - eml_star(0, eml(z, 1))
def err(z):         return abs(formula(z) - mpmath.conj(z))

def symbolic_minus_pi():
    """Exact verification at Im(z)=-π using SymPy (no floating point)."""
    x = sp.Symbol('x', real=True)
    # ln(-exp(x)) = x + iπ under principal branch (since -exp(x) < 0)
    result = sp.simplify(1 - (1 - sp.log(-sp.exp(x))))
    expected = x + sp.I * sp.pi  # = conj(x - iπ)
    return sp.simplify(result - expected) == 0

def run():
    print("="*65)
    print(f"verify_theorem4.py  —  {mpmath.mp.dps} decimal digits")
    print("="*65)

    p = t = 0

    # [1] Interior
    print("\n[1] Interior Im(z) ∈ (−π, π) — numerical — expect PASS\n")
    for label, z in [
        ("z = 1+0i",            mpmath.mpc(1, 0)),
        ("z = 0.5+1i",          mpmath.mpc("0.5", "1")),
        ("z = -2+2.5i",         mpmath.mpc("-2", "2.5")),
        ("z = 3-2.9i",          mpmath.mpc("3", "-2.9")),
        ("z = 1e6+3.14i",       mpmath.mpc("1e6", "3.14")),
        ("Im(z) = π/2",         mpmath.mpc(1, mpmath.pi/2)),
        ("Im(z) = -π/2",        mpmath.mpc(1, -mpmath.pi/2)),
        ("Im(z) = π - 1e-30",   mpmath.mpc(1, mpmath.pi - mpmath.mpf("1e-30"))),
        ("Im(z) = -π + 1e-30",  mpmath.mpc(1, -mpmath.pi + mpmath.mpf("1e-30"))),
    ]:
        e = err(z); ok = e < THRESHOLD; t += 1; p += ok
        print(f"  {'✓' if ok else '✗'}  {label:38s}  err={float(e):.2e}  {'PASS' if ok else 'FAIL'}")

    # [2] Boundary Im(z) = -π : symbolic
    print("\n[2] Boundary Im(z) = −π — SYMBOLIC (exact, no floating point)\n")
    ok_sym = symbolic_minus_pi(); t += 1; p += ok_sym
    print(f"  {'✓' if ok_sym else '✗'}  Symbolic proof: formula(x−iπ) = x+iπ = z̄   {'PASS' if ok_sym else 'FAIL'}")
    fp_err = err(mpmath.mpc(1, -mpmath.pi))
    print(f"\n  ⚠  Floating-point result at Im(z)=−π: err = {float(fp_err):.4f} ≈ 2π")
    print(f"     Cause: exp(−iπ) has Im ≈ +{float(mpmath.im(mpmath.exp(-mpmath.pi*1j))):.2e} (rounding).")
    print(f"     This pushes conj(exp(z)) below the branch cut → spurious jump.")
    print(f"     Mathematical domain is [−π, π); numerical domain is (−π, π).")
    print(f"     The boundary −π is in the domain but requires symbolic/exact arithmetic.")

    # [3] Boundary Im(z) = +π : fails
    print("\n[3] Boundary Im(z) = +π — expect FAIL\n")
    e = err(mpmath.mpc(1, mpmath.pi)); fails = e > THRESHOLD; t += 1; p += fails
    print(f"  {'✓' if fails else '✗'}  err = {float(e):.4f} ≈ {float(e/TWO_PI):.4f}×2π   {'CONFIRMED' if fails else 'UNEXPECTED'}")

    # [4] Exterior
    print("\n[4] Exterior |Im(z)| > π — expect FAIL (2π jump)\n")
    for label, z in [
        ("Im(z) = π+1e-10",   mpmath.mpc(1, mpmath.pi + mpmath.mpf("1e-10"))),
        ("Im(z) = π+0.1",     mpmath.mpc(1, mpmath.pi + mpmath.mpf("0.1"))),
        ("Im(z) = 2π",        mpmath.mpc(1, 2*mpmath.pi)),
        ("Im(z) = −π−0.1",    mpmath.mpc(1, -mpmath.pi - mpmath.mpf("0.1"))),
    ]:
        e = err(z); ok = e > THRESHOLD; t += 1; p += ok
        print(f"  {'✓' if ok else '✗'}  {label:25s}  err={float(e):.4f} = {float(e/TWO_PI):.4f}×2π")

    # Summary
    print("\n"+"="*65)
    print("SUMMARY")
    print("="*65)
    print(f"  Outcomes as expected : {p}/{t}")
    print()
    print("  MATHEMATICAL DOMAIN : Im(z) ∈ [−π, π)  (half-open)")
    print("  NUMERICAL DOMAIN    : Im(z) ∈ (−π, π)  (open — boundary −π unstable)")
    print()
    print("  Interior (−π, π)   ✓ numerical,  50 decimal digits, err < 1e-40")
    print("  Boundary Im(z)=−π  ✓ symbolic (exact) / ✗ numerical (unstable)")
    print("  Boundary Im(z)=+π  ✗ fails (branch cut, jump = 2πi)")
    print("  Exterior           ✗ fails (jump = 2πi)")
    return p == t

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
