"""Guaranteed-rounding certificate of the minimality of the mirror (2026-11-20, requested in review by ChatGPT).

Exact RATIONAL witness point z0 = (37 + 61 i)/100. For each of the 714 trees {eml, eml★} with leaves {1, z} of RPN
length ≤ 7, T(z0) − conj(z0) is enclosed in interval arithmetic (mpmath.iv, outward rounding) and the tree
is ELIMINATED as soon as the interval of the real part or that of the imaginary part excludes 0. A tree equal to conj
everywhere would be equal to it at z0: 714 eliminations = no tree of length ≤ 7 realizes the mirror.

Convention (that of the manuscript): principal log, Arg ∈ (−π, π]; extended evaluation node by node: ln 0 = −∞,
exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞; ∞ − ∞ and ln(−∞) undefined.
- CONSTANT subtrees (without z): decided exactly with SymPy (equality to 0, infinite, negative real). An EXACT
  negative real has principal log ln|b| + iπ: enclosed that way, without depending on a floating-point branch cut.
- Subtrees depending on z: if the argument of a log is a rectangle that contains 0 or that straddles the negative
  real axis, the tree is declared UNDECIDED (nothing is guessed).
- A tree that is undefined or infinite at z0 cannot be equal to conj(z0) (finite): it is eliminated, with that reason."""
import sympy as sp
from mpmath import iv

iv.dps = 60
Z0_EXACT = sp.Rational(37, 100) + sp.I * sp.Rational(61, 100)
Z0 = iv.mpc(iv.mpf("0.37"), iv.mpf("0.61"))
CONJ = iv.mpc(iv.mpf("0.37"), -iv.mpf("0.61"))
PI_I = iv.mpc(0, iv.pi)


class Undecided(Exception):
    pass


class Undefined(Exception):
    pass


def exact_const(op, a, b):
    """Exact value (SymPy, with sp.oo / -sp.oo) of a node whose two children are constant."""
    if op == "eml★" and b not in (sp.oo, -sp.oo):
        b = sp.conjugate(b)
    if b == -sp.oo:
        raise Undefined("ln(−∞)")
    lb = -sp.oo if sp.simplify(b) == 0 else (sp.oo if b == sp.oo else sp.log(b))
    ea = 0 if a == -sp.oo else (sp.oo if a == sp.oo else sp.exp(a))
    if ea == sp.oo and lb == sp.oo:
        raise Undefined("∞ − ∞")
    return sp.oo if ea == sp.oo else (-sp.oo if lb == sp.oo else (sp.oo if lb == -sp.oo else sp.simplify(ea - lb)))


def _is_negative_real(c):
    """Exact decision "c is a strictly negative real". If SymPy cannot decide, the tree is UNDECIDED
    (we do not guess, and we do not crash)."""
    im_zero = sp.im(c) == 0
    if im_zero is not True and im_zero is not False:
        raise Undecided(f"imaginary part of a constant not decided: {c}")
    if not im_zero:
        return False
    negative = sp.re(c).is_negative
    if negative is None:
        raise Undecided(f"sign of a constant not decided: {c}")
    return negative


def ev(t):
    """Returns (enclosure, exact value or None). The enclosure is an iv.mpc, or '+inf' / '-inf'."""
    if t == "1":
        return iv.mpc(1, 0), sp.Integer(1)
    if t == "z":
        return Z0, None
    op, ta, tb = t
    A, sa = ev(ta)
    B, sb = ev(tb)
    if sa is not None and sb is not None:                      # constant node: exact first
        s = exact_const(op, sa, sb)
        if s in (sp.oo, -sp.oo):
            return ("+inf" if s == sp.oo else "-inf"), s
    else:
        s = None
    # logarithm of the second argument
    if isinstance(B, str):
        if B == "-inf":
            raise Undefined("ln(−∞)")
        LB = "+inf"
    else:
        Bc = iv.mpc(B.real, -B.imag) if op == "eml★" else B
        if sb is not None:                                      # constant argument: exact decision
            bc = sp.conjugate(sb) if op == "eml★" else sb
            if sp.simplify(bc) == 0:
                LB = "-inf"
            elif _is_negative_real(bc):                           # EXACT negative real: Log = ln|b| + iπ
                LB = iv.log(iv.mpc(-Bc.real, 0)) + PI_I
            else:
                if Bc.real.a <= 0 and Bc.imag.a <= 0 <= Bc.imag.b:
                    raise Undecided(f"constant not decided near the branch cut: {bc}")
                LB = iv.log(Bc)
        else:
            contains0 = Bc.real.a <= 0 <= Bc.real.b and Bc.imag.a <= 0 <= Bc.imag.b
            if contains0:
                raise Undecided("log of a rectangle that contains 0")
            if not (Bc.real.a > 0 or Bc.imag.a > 0 or Bc.imag.b < 0):
                raise Undecided("log of a rectangle straddling the negative real axis")
            LB = iv.log(Bc)
    # exponential of the first argument
    if isinstance(A, str):
        EA = iv.mpc(0, 0) if A == "-inf" else "+inf"
    else:
        EA = iv.exp(A)
    if EA == "+inf":
        if LB == "+inf":
            raise Undefined("∞ − ∞")
        return "+inf", s
    if LB == "+inf":
        return "-inf", s
    if LB == "-inf":
        return "+inf", s
    return EA - LB, s


def trees(n):
    if n == 1:
        yield "1"; yield "z"; return
    for k in range(1, n):
        for left in trees(k):
            for right in trees(n - k):
                for op in ("eml", "eml★"):
                    yield (op, left, right)


def separated_from_zero(d):
    """Elimination decision, USED BY judge AND BY CHECK 3: returns "real" or "imaginary" if the corresponding
    part of the interval d strictly excludes 0 (the real part is examined first), otherwise None."""
    if d.real.a > 0 or d.real.b < 0:
        return "real"
    if d.imag.a > 0 or d.imag.b < 0:
        return "imaginary"
    return None


REASON_CONTAINS_ZERO = "UNDECIDED (both intervals contain 0)"


def judge(t):
    try:
        v, _ = ev(t)
    except Undefined as e:
        return "eliminated (undefined at z0: " + str(e) + ")"
    except Undecided as e:
        return "UNDECIDED (" + str(e) + ")"
    if isinstance(v, str):
        return "eliminated (infinite value at z0)"
    part = separated_from_zero(v - CONJ)
    if part is not None:
        return f"eliminated ({part} part of T(z0) − conj(z0) separated from 0)"
    return REASON_CONTAINS_ZERO


def checks():
    """Checks of the code, requested in review: coverage of the generator, enclosure of the exact point, elimination rule."""
    import sympy as sp
    all_trees = [t for n in range(1, 5) for t in trees(n)]
    distinct = len(set(all_trees))
    from math import comb
    catalan = lambda m: comb(2 * m, m) // (m + 1)
    expected = sum(catalan(n - 1) * 2 ** n * 2 ** (n - 1) for n in range(1, 5))     # shapes × leaves × operators
    print(f"Check 1, generator: {len(all_trees)} trees produced, {distinct} distinct, counting formula {expected}")
    assert len(all_trees) == distinct == expected == 714
    def exact_rational(raw):
        """EXACT binary value of an interval bound, read directly from its internal representation
        (sign, mantissa, exponent), as a SymPy rational. (Do not go through mpmath.mpf: it rounds to the current
        precision, 53 bits by default, and would wrongly suggest that the interval is reduced to a point.)"""
        sign, mantissa, exponent, _ = raw
        return (-1) ** sign * sp.Integer(mantissa) * sp.Integer(2) ** exponent
    for v, exact in ((Z0.real, sp.Rational(37, 100)), (Z0.imag, sp.Rational(61, 100))):
        low, high = exact_rational(v._mpi_[0]), exact_rational(v._mpi_[1])
        assert low < exact < high, "z0 outside its interval"
    print("Check 2, enclosure: the interval of z0 strictly contains 37/100 and 61/100 (exact comparison in rationals)")
    zero_real = iv.mpc(iv.mpf([-1e-70, 1e-70]), iv.mpf([1, 2]))       # real part contains 0, imaginary part does not
    zero_imag = iv.mpc(iv.mpf([1, 2]), iv.mpf([-1e-70, 1e-70]))       # the reverse
    zero_both = iv.mpc(iv.mpf([-1e-70, 1e-70]), iv.mpf([-1e-70, 1e-70]))
    edge = iv.mpc(iv.mpf([0, 1]), iv.mpf([0, 1]))                      # 0 on the edge of both intervals: not separated
    assert separated_from_zero(zero_real) == "imaginary" and separated_from_zero(zero_imag) == "real"
    assert separated_from_zero(zero_both) is None and separated_from_zero(edge) is None
    print("Check 3, rule (the function separated_from_zero, the one judge uses): a deviation whose two intervals both contain 0, even at the edge, is NOT eliminated; a single interval strictly separated from 0 suffices")
    return all_trees


if __name__ == "__main__":
    import collections, sys
    checks()
    tally, by_length, total = collections.Counter(), collections.Counter(), 0
    undecided, by_imag = [], []
    for n in range(1, 5):
        for t in trees(n):
            total += 1
            by_length[2 * n - 1] += 1
            r = judge(t)
            tally[r] += 1
            if "imaginary part" in r:
                by_imag.append(t)
            if r.startswith("UNDECIDED"):
                undecided.append((t, r))
    print(f"z0 = (37 + 61 i)/100 exact; trees enumerated: {total}; by RPN length: {dict(sorted(by_length.items()))}")
    for k, v in sorted(tally.items()):
        print(f"   {v:4d} × {k}")
    print(f"UNDECIDED: {len(undecided)}")
    print(f"Tree(s) eliminated by the imaginary part: {by_imag} (predicted in review: the leaf z alone, since z0 − conj(z0) = 61i/50)")
    for t, r in undecided[:10]:
        print("   ", t, r)
    F = ("eml", ("eml", "1", ("eml", ("eml★", "1", "z"), "1")), "1")
    verdict_F = judge(F)
    print("Positive check, the 9-symbol tree F:", verdict_F, "(required: \"both intervals contain 0\", since it is exactly equal to conj(z0))")
    n_real = sum(v for k, v in tally.items() if "real part" in k)
    expectations = {
        "714 trees enumerated": total == 714,
        "0 undecided": not undecided,
        "713 eliminated by the real part": n_real == 713,
        "only the leaf z eliminated by the imaginary part": by_imag == ["z"],
        "positive check F not eliminated, exact reason": verdict_F == REASON_CONTAINS_ZERO,
    }
    for name, passed in expectations.items():
        print(f"   [{'OK' if passed else 'FAIL'}] {name}")
    ok = all(expectations.values())
    print("CONCLUSION:", "714 trees eliminated, 0 undecided: no tree of RPN length ≤ 7 realizes the mirror; 9 is minimal." if ok else "INCONCLUSIVE")
    sys.exit(0 if ok else 1)
