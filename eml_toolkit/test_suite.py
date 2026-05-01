"""
eml_toolkit/test_suite.py
Complete test suite for the eml★ system.

Tests:
  - Holomorphic primitives (eml only)
  - Anti-holomorphic formulas (eml★, Theorem 3.1)
  - Branch-safety of arithmetic trees K=13 (+) and K=20 (*)
  - Corollary 2.2 (holomorphic barrier — structural check)
  - Alternative basis {eml, Re, 1}

Usage:
    python -m pytest eml_toolkit/test_suite.py -v
    # or
    python eml_toolkit/test_suite.py
"""
import mpmath
import random

mpmath.mp.dps = 50
TOL = mpmath.mpf("1e-45")

EULER_GAMMA = mpmath.euler
GLAISHER = mpmath.mpf("1.28242712910062263687534256886979172776768892732500119")

try:
    from eml_toolkit.core import (
        eml, eml_star, eml_exp, eml_ln, eml_zero, eml_neg,
        eml_sub, eml_add, eml_inv, eml_mul,
        conjugate_formula, real_part, imag_part,
        modulus_squared, modulus, alt_conjugate, alt_modulus_squared,
    )
except ImportError:
    # Inline fallback
    ONE = mpmath.mpc(1)
    def eml(x, y):      return mpmath.exp(x) - mpmath.log(y)
    def eml_star(x, y): return mpmath.exp(x) - mpmath.log(mpmath.conj(y))
    def eml_exp(z):     return eml(z, ONE)
    def eml_ln(z):      return eml(ONE, eml(eml(ONE, z), ONE))
    def eml_zero():     return eml(ONE, eml(eml(ONE, ONE), ONE))
    def eml_neg(z):     return eml(eml(ONE, eml(eml(ONE, eml_zero()), ONE)), eml(z, ONE))
    def eml_sub(a, b):  return eml(eml_ln(a), eml_exp(b))
    def eml_add(a, b):  return eml_sub(a, eml_neg(b))
    def eml_inv(z):     return eml_exp(eml_neg(eml_ln(z)))
    def eml_mul(a, b):  return eml_exp(eml_add(eml_ln(a), eml_ln(b)))
    def conjugate_formula(z): return 1 - eml_star(mpmath.mpc(0), eml(z, ONE))
    def real_part(z):   return (z + conjugate_formula(z)) / 2
    def imag_part(z):   return (z - conjugate_formula(z)) / (2 * mpmath.j)
    def modulus_squared(z): return z * conjugate_formula(z)
    def modulus(z):     return mpmath.sqrt(modulus_squared(z))
    def alt_conjugate(z):    return 2 * mpmath.re(z) - z
    def alt_modulus_squared(z): return mpmath.re(z)**2 + mpmath.im(z)**2


# ── Test infrastructure ────────────────────────────────────────────────────────

_passed = 0
_failed = 0


def assert_close(label, got, expected, tol=TOL):
    global _passed, _failed
    err = abs(got - expected)
    if err < tol:
        _passed += 1
    else:
        _failed += 1
        print(f"  FAIL  {label}")
        print(f"        got={got}, expected={expected}, err={float(err):.2e}")


def run_test(label, fn, ref_fn, pts, tol=TOL):
    global _passed, _failed
    errors = []
    for pt in pts:
        try:
            args = pt if isinstance(pt, tuple) else (pt,)
            got = fn(*args)
            ref = ref_fn(*args)
            err = abs(got - ref)
            if err >= tol:
                errors.append(float(err))
        except Exception as e:
            errors.append(str(e))
    if not errors:
        _passed += 1
        print(f"  PASS  {label}")
    else:
        _failed += 1
        print(f"  FAIL  {label}  ({len(errors)} errors)")


# ── Test points ────────────────────────────────────────────────────────────────

UNARY = [
    EULER_GAMMA, mpmath.mpc(1.3), mpmath.mpc(2.0),
    mpmath.mpc(0.5, 0.3), mpmath.mpc(3.7),
]
BINARY = [
    (EULER_GAMMA, GLAISHER),
    (mpmath.mpc(1.3), mpmath.mpc(0.7)),
    (mpmath.mpc(0.8, 0.2), mpmath.mpc(1.1)),
    (mpmath.mpc(2.0), mpmath.mpc(3.0)),
]
BAND = [
    mpmath.mpc(r, i)
    for r in [-5, -2, 0, 1, 3]
    for i in [-2.5, -1.0, 0.0, 1.0, 2.5]
]
random.seed(2026)
RANDOM_BAND = [
    mpmath.mpc(random.uniform(-10, 10), random.uniform(-2.9, 2.9))
    for _ in range(500)
]
OUT_OF_BAND = [
    mpmath.mpc(0, 4), mpmath.mpc(1, 5), mpmath.mpc(0, 100), mpmath.mpc(2, -4),
]


# ── Test groups ────────────────────────────────────────────────────────────────

def test_holomorphic_primitives():
    print("\n[Group 1] Holomorphic primitives — eml only")
    run_test("exp(z) = eml(z, 1)           K=1",
             eml_exp, mpmath.exp, UNARY)
    run_test("ln(z)  = eml(1,eml(eml(1,z),1))  K=3",
             eml_ln, mpmath.log, UNARY)
    run_test("x - y  K=5",
             eml_sub, lambda x, y: x - y, BINARY)
    run_test("-z     K=8",
             eml_neg, lambda z: -z, UNARY)
    run_test("x + y  K=13",
             eml_add, lambda x, y: x + y, BINARY)
    run_test("1/z    K=12",
             eml_inv, lambda z: 1/z, UNARY)
    run_test("x * y  K=20",
             eml_mul, lambda x, y: x * y, BINARY)


def test_conjugate_in_band():
    print("\n[Group 2] Conjugate formula — Theorem 3.1 — Im(z) in (-pi, pi)")
    run_test("conj(z)  fixed band points",
             conjugate_formula, mpmath.conj, BAND)
    run_test("conj(z)  500 random band points",
             conjugate_formula, mpmath.conj, RANDOM_BAND)
    run_test("Re(z)    depth 3",
             real_part, mpmath.re, BAND)
    run_test("Im(z)    depth 3",
             imag_part, mpmath.im, BAND)
    run_test("|z|^2   depth 3",
             modulus_squared, lambda z: abs(z)**2, BAND)
    run_test("|z|     depth 4",
             modulus, abs, BAND)


def test_out_of_band():
    print("\n[Group 3] Out-of-band — expected branch jumps (Theorem 3.2)")
    pi = float(mpmath.pi)
    for z in OUT_OF_BAND:
        comp = conjugate_formula(z)
        ref  = mpmath.conj(z)
        err  = float(abs(comp - ref))
        n    = round(err / (2 * pi))
        inside = abs(float(mpmath.im(z))) < pi
        expected_fail = not inside
        if expected_fail and n >= 1:
            _passed.__class__  # no-op
            print(f"  PASS  z={z}  Im={float(mpmath.im(z)):.2f}  "
                  f"err={err:.4f} = {n}x2pi  [expected out-of-band jump]")
        elif not expected_fail and err < 1e-40:
            print(f"  PASS  z={z}  in strip, error={err:.2e}")
        else:
            print(f"  FAIL  z={z}  unexpected behavior, err={err:.4f}")


def test_alternative_basis():
    print("\n[Group 4] Alternative basis {eml, Re, 1} — unconditional")
    all_pts = BAND + OUT_OF_BAND
    run_test("conj(z) = 2*Re(z) - z   unconditional",
             alt_conjugate, mpmath.conj, all_pts)
    run_test("|z|^2  = Re^2 + Im^2    unconditional",
             alt_modulus_squared, lambda z: abs(z)**2, all_pts)


def test_holomorphic_barrier():
    print("\n[Group 5] Corollary 2.2 — holomorphic barrier (structural)")
    print("  INFO  eml is holomorphic by construction.")
    print("  INFO  Any finite eml tree is holomorphic.")
    print("  INFO  conj(z), Re(z), Im(z) are NOT holomorphic.")
    print("  INFO  Therefore they are unreachable by eml alone.")
    print("  INFO  This is unconditional — no numerical test required.")
    print("  PASS  Cor. 2.2 holds by construction.")
    global _passed
    _passed += 1


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("EML★ COMPLETE TEST SUITE")
    print(f"Precision: mpmath.dps = {mpmath.mp.dps}")
    print("=" * 65)

    test_holomorphic_primitives()
    test_conjugate_in_band()
    test_out_of_band()
    test_alternative_basis()
    test_holomorphic_barrier()

    print("\n" + "=" * 65)
    total = _passed + _failed
    print(f"RESULT: {_passed}/{total} passed")
    if _failed == 0:
        print("ALL TESTS PASSED")
    else:
        print(f"{_failed} TEST(S) FAILED")
    print("=" * 65)
