"""Strict mirror in eml★ alone, 101 symbols — standalone verification script (SymPy only).

Discovery: ChatGPT, 20 November 2026 (see FINDINGS_20261119_branch_safety_certificate.md, § 17).
This file rebuilds the tree without depending on any other file of the repository, writes it in postfix notation, then judges it.

    S(a, b) = eml★(a, b) = exp(a) − Log(conj b),   principal log, Arg in (−π, π].
    Leaves: 1 and z. No other gate, no constant supplied.

STRICT evaluation: every node must have a finite value; Log(0) is refused (Undefined error).

What the script establishes: the grammar, the length, the identity J(z) = conj(z) in exact computation on a nonzero
complex symbol and at exact points, the refusal at zero, and eml★(x, J(y)) = eml(x, y).
What it does NOT establish: that no intermediate node vanishes for EVERY nonzero z. That comes from the domain
argument of § 17 (the imaginary part of S(b, z) is in [−3π, −π)), to be reviewed. Nor the minimality of 101.
Exit code 0 only if all the checks pass.   Usage: python3 star_only_mirror_101.py
"""
import sys

import sympy as sp


class Undefined(ArithmeticError):
    pass


def S(a, b):
    return ("S", a, b)


def X(t):                      # exp(t)
    return S(t, "1")


def Q(t):                      # e^e / conj(t)
    return X(S("1", t))


def L(t):                      # a logarithm of t: exp(L(t)) = t
    return S("1", Q(t))


def d(t, r):                   # t − r, r real
    return S(L(t), X(r))


E_ = X("1")                                    # e
O = S("1", X(E_))                              # 0
MINUS_ONE = d(S("1", E_), E_)                  # −1
K = d(S(O, MINUS_ONE), "1")                    # −iπ
H = S(L(K), MINUS_ONE)                         # −2πi: one full turn
B = L(H)                                       # exp(B) = −2πi
A = L(S("1", MINUS_ONE))                       # exp(exp(A)) = −e^e


def P(t):                      # −t for nonzero t
    return X(S(A, Q(t)))


def J(t):                      # conj(t) for nonzero t
    return X(P(S(B, t)))


def postfix(t):
    return [t] if isinstance(t, str) else postfix(t[1]) + postfix(t[2]) + ["eml★"]


def evaluate(t, values, star=True, _memo=None):
    """Exact value (SymPy), strict. star=False replaces every eml★ by eml (check)."""
    memo = {} if _memo is None else _memo
    if isinstance(t, str):
        return sp.Integer(1) if t == "1" else values[t]
    if t not in memo:
        a, b = evaluate(t[1], values, star, memo), evaluate(t[2], values, star, memo)
        if b == 0:
            raise Undefined("Log(0)")
        v = sp.simplify(sp.exp(a) - sp.log(sp.conjugate(b) if star else b))
        if not v.is_finite and v.is_finite is not None:
            raise Undefined("infinite value")
        memo[t] = v
    return memo[t]


def main():
    ok = []

    def check(name, condition):
        ok.append(bool(condition))
        print(f"[{'OK' if condition else 'FAIL'}] {name}")

    tree = J("z")
    tokens = postfix(tree)
    print("Tree in postfix notation:\n" + " ".join(tokens) + "\n")
    check(f"grammar: only 1, z and eml★ (found: {sorted(set(tokens))})", set(tokens) == {"1", "z", "eml★"})
    check(f"unfolded length = {len(tokens)} symbols (expected 101)", len(tokens) == 101)

    for name, c, expected in (("O", O, 0), ("−1", MINUS_ONE, -1), ("k", K, -sp.I * sp.pi), ("h", H, -2 * sp.I * sp.pi)):
        check(f"constant {name} = {expected}", sp.simplify(evaluate(c, {}) - expected) == 0)
    check("constant b: exp(b) = −2πi", sp.simplify(sp.exp(evaluate(B, {})) + 2 * sp.I * sp.pi) == 0)
    check("constant a: exp(exp(a)) = −e^e", sp.simplify(sp.exp(sp.exp(evaluate(A, {}))) + sp.exp(sp.E)) == 0)

    z = sp.Symbol("z", complex=True, zero=False)
    check("J(z) = conj(z) on a nonzero complex symbol", sp.simplify(evaluate(tree, {"z": z}) - sp.conjugate(z)) == 0)
    points = (sp.Integer(1), sp.Integer(-1), sp.I, sp.exp(sp.E), sp.Rational(-1, 3), sp.Integer(-50),
              sp.Rational(3, 2) + 2 * sp.I, -5 * sp.I, sp.Rational(1, 7) - 40 * sp.I)
    for q in points:
        check(f"J({q}) = conj({q})", sp.simplify(evaluate(tree, {"z": q}) - sp.conjugate(q)) == 0)

    try:
        evaluate(tree, {"z": sp.Integer(0)})
        check("strict refusal at z = 0", False)
    except Undefined:
        check("strict refusal at z = 0", True)

    x, y = sp.Symbol("x", complex=True), sp.Symbol("y", complex=True, zero=False)
    gate = S("x", J("y"))
    check(f"eml★(x, J(y)) = eml(x, y) on symbols ({len(postfix(gate))} symbols)",
          sp.simplify(evaluate(gate, {"x": x, "y": y}) - (sp.exp(x) - sp.log(y))) == 0)
    for q in (sp.Integer(-2), sp.Integer(-50), sp.Rational(1, 2) + sp.I):
        check(f"eml★(1/3, J({q})) = eml(1/3, {q})",
              sp.simplify(evaluate(gate, {"x": sp.Rational(1, 3), "y": q}) - (sp.exp(sp.Rational(1, 3)) - sp.log(q))) == 0)

    try:                                         # negative check: without the star, the tree is NOT a mirror
        without_star = evaluate(tree, {"z": sp.I}, star=False)
        check("check: the same tree without the star does not return conj(i)", sp.simplify(without_star + sp.I) != 0)
    except Undefined:
        check("check: the same tree without the star does not return conj(i) (it is undefined)", True)

    print(f"\nSUMMARY: {sum(ok)}/{len(ok)} checks passed. Python {sys.version.split()[0]}, SymPy {sp.__version__}.")
    print("Not established by this script: absence of an intermediate zero for EVERY nonzero z (domain argument, § 17); minimality.")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
