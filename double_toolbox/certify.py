"""Interval certificate of the "wall-free" recipes, with the already reviewed certifier of the repository
(certify_branch_safety.py: splitting into boxes, logarithms checked one by one, gate rule).

What the certificate establishes: on the whole domain, NO logarithm of the tree crosses its cut without the
jump being erased by an exponential (gate), and none sees zero. The tree therefore makes no branch jump there.
What it does not establish on its own: equality with the target function. That comes from the formal identity, plus the
witness points compared here, plus the continuation argument (to be reviewed, as in the report).

Two preparations, both already admitted in the report:
  1. the absorption lemma: eml(E_ln(a), c) is rewritten a − ln c (proof by hand, report § 4);
  2. the constants (i, π, ½...): EXACT value by SymPy, then enclosure by interval arithmetic on the
     exact form (integers, rationals, π, e, i). The original certifier refuses i because its construction goes
     deliberately ON the cut (ln(−1)); the exactness of SymPy is already the admitted basis for 0 and ±∞.
"""
import os
import sys
import time

import mpmath as mp
import sympy as sp
from mpmath import iv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import certify_branch_safety as cb  # noqa: E402

from . import eml as E  # noqa: E402


# ── 1. Translation of the trees, with the absorption lemma ────────────────────
def _is_ln(n):
    """n = ln_(a) = eml(1, eml(eml(1, a), 1))? Returns a, or None."""
    if isinstance(n, str) or n.op != "eml" or n.left != E.ONE:
        return None
    m = n.right
    if isinstance(m, str) or m.op != "eml" or m.right != E.ONE:
        return None
    k = m.left
    if isinstance(k, str) or k.op != "eml" or k.left != E.ONE:
        return None
    return k.right


def translate(t, _memo=None):
    memo = {} if _memo is None else _memo
    if isinstance(t, str):
        return cb.ONE if t == E.ONE else (t,)
    if id(t) in memo:
        return memo[id(t)]
    a = _is_ln(t.left) if t.op == "eml" else None
    if a is not None and E.variables(t):
        r = cb.minus_ln(translate(a, memo), translate(t.right, memo))
    else:
        r = (t.op, translate(t.left, memo), translate(t.right, memo))
    memo[id(t)] = r
    return r


# ── 2. Exact constants → intervals ────────────────────────────────────────────
def _iv_exact(e):
    """Guaranteed enclosure of an exact SymPy expression made of integers, rationals, π, e, i, +, ×, integer
    powers, exp, and log of a strictly positive real. Any other form is refused."""
    if e.is_Integer:
        return iv.mpc(int(e), 0)
    if e.is_Rational:
        return iv.mpc(iv.mpf(int(e.p)) / iv.mpf(int(e.q)), 0)
    if e is sp.pi:
        return iv.mpc(iv.pi, 0)
    if e is sp.E:
        return iv.mpc(iv.e, 0)
    if e is sp.I:
        return iv.mpc(0, 1)
    if e.is_Add or e.is_Mul:
        pieces = [_iv_exact(a) for a in e.args]
        r = pieces[0]
        for m in pieces[1:]:
            r = r + m if e.is_Add else r * m
        return r
    if e.is_Pow and e.exp.is_Integer:
        base, n = _iv_exact(e.base), int(e.exp)
        r = iv.mpc(1, 0)
        for _ in range(abs(n)):
            r = r * base
        return r if n >= 0 else iv.mpc(1, 0) / r
    if isinstance(e, sp.exp):
        return iv.exp(_iv_exact(e.args[0]))
    if isinstance(e, sp.log):
        a = _iv_exact(e.args[0])
        if a.real.a > 0 and a.imag.a == 0 == a.imag.b:
            return iv.mpc(iv.log(a.real), 0)
    raise cb.Failure(f"exact constant of an unforeseen form: {e}")


_original = cb._constant_interval


def _constant_interval(t):
    try:
        return _original(t)
    except cb.Failure:
        return _iv_exact(cb.constant(t))


cb._constant_interval = _constant_interval


# ── 3. The certificates ───────────────────────────────────────────────────────
def certify(name, tree, reference, domains, witnesses, depth=7, budget=240, disc=None):
    """Certifies the tree on each product of boxes of `domains`, then compares with the witness points.
    disc = (center, radius): the boxes of x entirely inside this open disc are set aside (counted)."""
    t = translate(tree)
    summary = {"name": name, "symbols": E.size(tree), "certified": True, "boxes": 0, "absorbed": 0, "excluded": 0, "reason": ""}
    t0 = time.time()
    for d in domains:
        r = cb.certify(t, d, depth=depth, budget=budget, disc_x=disc)
        summary["boxes"] += r["boxes"]
        summary["excluded"] += r["excluded"]
        summary["absorbed"] += r["absorbed"]
        if not r["certified"]:
            summary.update(certified=False, reason=r["reason"][:300])
            break
    worst = mp.mpf(0)
    for p in witnesses:
        worst = max(worst, abs(E.evaluate(tree, p) - reference(**{k: mp.mpmathify(v) for k, v in p.items()})))
    summary.update(worst_witness_gap=worst, witnesses=len(witnesses), seconds=round(time.time() - t0, 1))
    return summary
