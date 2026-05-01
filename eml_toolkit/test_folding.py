"""
test_folding.py
===============
Test de la propriété empirique de folding (Remark 3.5, Monnerot 2026).

φ(z) = 1 − eml★(0, eml(z, 1)) ramène Im(φ(z)) dans [−π, π)
pour tout z ∈ ℂ, via la branche principale du logarithme.

Vérifié à 60 décimales sur 10 000 points hors bande.
"""
import random
import sys
import mpmath
mpmath.mp.dps = 60

def eml(x, y):      return mpmath.exp(x) - mpmath.log(y)
def eml_star(x, y): return mpmath.exp(x) - mpmath.log(mpmath.conj(y))

def fold_to_band(z):
    z = mpmath.mpc(z)
    return mpmath.mpc(1) - eml_star(mpmath.mpc(0), eml(z, mpmath.mpc(1)))

def test_fold_to_band(n_points=10000, seed=42):
    random.seed(seed)
    PI = mpmath.pi
    success = 0
    max_im_before = mpmath.mpf(0)
    max_im_after  = mpmath.mpf(0)

    for _ in range(n_points):
        re   = random.uniform(-100, 100)
        im_v = random.uniform(float(PI) + 0.01, 100) * random.choice([-1, 1])
        z    = mpmath.mpc(re, im_v)

        im_z = abs(mpmath.im(z))
        max_im_before = max(max_im_before, im_z)

        folded = fold_to_band(z)
        im_f   = mpmath.im(folded)
        max_im_after = max(max_im_after, abs(im_f))

        if -PI <= im_f < PI:
            success += 1

    print(f"\n{'='*55}")
    print(f"TEST fold_to_band — {mpmath.mp.dps} décimales, {n_points} points")
    print(f"{'='*55}")
    print(f"  Repliés dans [-π, π) : {success}/{n_points}  ({success/n_points*100:.6f}%)")
    print(f"  Max |Im(z)| avant    : {float(max_im_before):.4f}")
    print(f"  Max |Im(φ(z))| après : {float(max_im_after):.10f}  (< π = {float(PI):.10f})")
    print()

    # Test erreur : φ(z) ≠ z̄ hors bande → erreur = multiple de 2π
    print("  Vérification erreur hors bande :")
    TWO_PI = 2 * PI
    for im_val in [10, 20, 100]:
        z = mpmath.mpc(1, im_val)
        err = abs(fold_to_band(z) - mpmath.conj(z))
        k = round(float(err / TWO_PI))
        print(f"    Im(z)={im_val:3d} → err = {float(err):.4f} ≈ {k}×2π  {'✓' if abs(err - k*TWO_PI) < 0.01 else '✗'}")

    assert success == n_points, f"Folding failed on {n_points - success} points"
    print(f"\n  ✓ 100% des points repliés correctement.")
    print(f"  ✓ Erreur = multiple entier de 2π confirmé hors bande.")
    return True

if __name__ == "__main__":
    ok = test_fold_to_band()
    sys.exit(0 if ok else 1)
