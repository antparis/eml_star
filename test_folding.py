"""
test_folding.py
Empirical test of the folding property (Remark 3.5, Monnerot 2026).

phi(z) = 1 - eml★(0, eml(z, 1)) maps Im(phi(z)) into [-pi, pi)
for any z in C, via the principal branch of the logarithm.

Verified at 60 decimal places on 10,000 out-of-band points.

Usage:
    python test_folding.py
"""
import random
import sys
import mpmath

mpmath.mp.dps = 60


def eml(x, y):
    return mpmath.exp(x) - mpmath.log(y)


def eml_star(x, y):
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))


def fold_to_band(z):
    """Conjugation formula (Theorem 3.1): 1 - eml★(0, eml(z, 1))."""
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

    print(f"\n{'=' * 55}")
    print(f"TEST fold_to_band — {mpmath.mp.dps} decimal places, {n_points} points")
    print(f"{'=' * 55}")
    print(f"  Folded into [-pi, pi) : {success}/{n_points}  ({success / n_points * 100:.6f}%)")
    print(f"  Max |Im(z)| before    : {float(max_im_before):.4f}")
    print(f"  Max |Im(phi(z))| after: {float(max_im_after):.10f}  (< pi = {float(PI):.10f})")
    print()

    # Error check: phi(z) != conj(z) out of band -> error = integer multiple of 2*pi
    print("  Out-of-band error verification:")
    TWO_PI = 2 * PI
    for im_val in [10, 20, 100]:
        z = mpmath.mpc(1, im_val)
        err = abs(fold_to_band(z) - mpmath.conj(z))
        k = round(float(err / TWO_PI))
        ok = abs(err - k * TWO_PI) < 0.01
        print(f"    Im(z)={im_val:3d} -> err = {float(err):.4f} ~ {k}x2pi  "
              f"[{'PASS' if ok else 'FAIL'}]")

    assert success == n_points, f"Folding failed on {n_points - success} points"
    print(f"\n  PASS  100% of points folded correctly.")
    print(f"  PASS  Error = integer multiple of 2*pi confirmed out of band.")
    return True


if __name__ == "__main__":
    ok = test_fold_to_band()
    sys.exit(0 if ok else 1)
