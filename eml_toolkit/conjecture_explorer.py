"""
eml_toolkit/conjecture_explorer.py
Heuristic conjecture explorer for {eml, eml★, 1}.

Generates random candidate expressions, evaluates them numerically at high
precision, and checks whether they match a given target.

This is a parallel random search over the space of eml/eml★ expression trees.
Useful for:
  - Discovering unknown identities in the eml / eml★ system
  - Numerically validating conjectures at 50+ decimal places
  - Exploring the expression space up to a given depth

Author: Anthony Monnerot (2026)
Extension of: Odrzywołek (arXiv:2603.21852v2)
"""

import random
import mpmath
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from eml_toolkit.core import (
        eml, eml_star, eml_exp, eml_ln, eml_zero, eml_neg,
        eml_sub, eml_add, eml_mul, conjugate_formula,
    )
except ImportError:
    ONE = mpmath.mpc(1)
    def eml(x, y):       return mpmath.exp(x) - mpmath.log(y)
    def eml_star(x, y):  return mpmath.exp(x) - mpmath.log(mpmath.conj(y))
    def eml_exp(z):      return eml(z, ONE)
    def eml_ln(z):       return eml(ONE, eml(eml(ONE, z), ONE))
    def eml_zero():      return eml(ONE, eml(eml(ONE, ONE), ONE))
    def eml_neg(z):      return eml(eml(ONE, eml(eml(ONE, eml_zero()), ONE)), eml(z, ONE))
    def eml_sub(a, b):   return eml(eml_ln(a), eml_exp(b))
    def eml_add(a, b):   return eml_sub(a, eml_neg(b))
    def eml_mul(a, b):   return eml_exp(eml_add(eml_ln(a), eml_ln(b)))
    def conjugate_formula(z): return 1 - eml_star(mpmath.mpc(0), eml(z, mpmath.mpc(1)))

_EPS = mpmath.mpf("1e-40")


# ── Random expression generator ────────────────────────────────────────────────

def _random_leaf(z, rng):
    """Return a random terminal: z, 1, 0, or a small integer."""
    choice = rng.randint(0, 3)
    if choice == 0:
        return z
    elif choice == 1:
        return mpmath.mpc(1)
    elif choice == 2:
        return mpmath.mpc(0)
    else:
        return mpmath.mpc(rng.choice([2, -1, 3]))


def _random_expr(z, depth, max_depth, rng):
    """
    Recursively build a random eml / eml★ expression tree.

    Parameters
    ----------
    z         : the input value (mpmath complex number)
    depth     : current recursion depth
    max_depth : maximum allowed depth
    rng       : random.Random instance for reproducibility
    """
    if depth >= max_depth:
        return _random_leaf(z, rng)

    op = rng.choice(["eml", "eml_star", "leaf"])
    if op == "leaf":
        return _random_leaf(z, rng)
    elif op == "eml":
        left  = _random_expr(z, depth + 1, max_depth, rng)
        right = _random_expr(z, depth + 1, max_depth, rng)
        try:
            return eml(left, right)
        except Exception:
            return _random_leaf(z, rng)
    else:  # eml_star
        left  = _random_expr(z, depth + 1, max_depth, rng)
        right = _random_expr(z, depth + 1, max_depth, rng)
        try:
            return eml_star(left, right)
        except Exception:
            return _random_leaf(z, rng)


# ── Single-candidate evaluator (runs in worker process) ───────────────────────

def _evaluate_one(args):
    """
    Evaluate one random candidate against the target.

    Returns (candidate_value, error) if error < EPS, else None.
    """
    seed, z_re, z_im, target_re, target_im, max_depth, dps = args
    mpmath.mp.dps = dps
    rng = random.Random(seed)
    z      = mpmath.mpc(z_re, z_im)
    target = mpmath.mpc(target_re, target_im)
    try:
        val = _random_expr(z, 0, max_depth, rng)
        err = abs(val - target)
        if err < _EPS:
            return (val, float(err))
    except Exception:
        pass
    return None


# ── ConjectureExplorer class ───────────────────────────────────────────────────

class ConjectureExplorer:
    """
    Parallel heuristic explorer for eml / eml★ conjectures.

    Generates random candidate expressions and checks whether any of them
    numerically matches a given target function at a test point.

    Parameters
    ----------
    max_depth    : maximum depth of the random expression trees
    num_rollouts : total number of random candidates to evaluate
    num_workers  : number of parallel worker processes
    dps          : mpmath precision in decimal places
    """

    def __init__(
        self,
        max_depth: int = 4,
        num_rollouts: int = 200,
        num_workers: int = 4,
        dps: int = 50,
    ):
        self.max_depth    = max_depth
        self.num_rollouts = num_rollouts
        self.num_workers  = num_workers
        self.dps          = dps
        self.results      = []

    def explore(self, target_fn, test_point=None):
        """
        Launch parallel exploration for a target function.

        Parameters
        ----------
        target_fn  : callable z -> mpmath complex, the function to match
        test_point : mpmath complex test point (default: Euler-Mascheroni + 0.3i)

        Returns
        -------
        list of (value, error) sorted by increasing error
        """
        mpmath.mp.dps = self.dps

        if test_point is None:
            test_point = mpmath.mpc(float(mpmath.euler), 0.3)

        target_val = target_fn(test_point)
        z_re       = float(mpmath.re(test_point))
        z_im       = float(mpmath.im(test_point))
        t_re       = float(mpmath.re(target_val))
        t_im       = float(mpmath.im(target_val))

        args_list = [
            (seed, z_re, z_im, t_re, t_im, self.max_depth, self.dps)
            for seed in range(self.num_rollouts)
        ]

        self.results = []
        with ProcessPoolExecutor(max_workers=self.num_workers) as pool:
            futures = [pool.submit(_evaluate_one, a) for a in args_list]
            for f in as_completed(futures):
                res = f.result()
                if res is not None:
                    self.results.append(res)

        self.results.sort(key=lambda x: x[1])
        return self.results[:10]

    def report(self):
        """Print the top exploration results."""
        if not self.results:
            print("No matching candidates found.")
            return
        print(f"{'Rank':>4}  {'Error':>14}  Value")
        print("-" * 55)
        for i, (val, err) in enumerate(self.results[:10], 1):
            print(f"{i:>4}  {err:>14.3e}  {val}")


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("ConjectureExplorer — demo: searching for conj(z)")
    print()

    explorer = ConjectureExplorer(
        max_depth=3,
        num_rollouts=500,
        num_workers=4,
        dps=50,
    )

    results = explorer.explore(
        target_fn=mpmath.conj,
        test_point=mpmath.mpc(1.2, 0.8),
    )

    if results:
        print(f"Found {len(results)} matching candidate(s):")
        explorer.report()
    else:
        print("No exact match found in this run.")
        print("(Expected: the conjugate formula 1 - eml★(0, eml(z,1)) may appear)")
