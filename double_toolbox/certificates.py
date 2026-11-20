"""Runs all the certificates of the wall-free recipes:  python -m double_toolbox.certificates
Exit code 0 only if all the expected verdicts are obtained (negative controls included)."""
import sys

import mpmath as mp
import sympy as sp

import certify_branch_safety as cb

from . import certify as C
from . import eml as E

B = cb.Box
EXPECTED = []


def report(summary, expected=True):
    ok = summary["certified"] == expected and (not expected or summary["worst_witness_gap"] < mp.mpf(10) ** -40)
    EXPECTED.append(ok)
    print(f"[{'OK' if ok else 'FAILED'}] {summary['name']} ({summary['symbols']} symbols): "
          f"{'CERTIFIED with no jump' if summary['certified'] else 'not certified'} — {summary['boxes']} boxes, "
          f"{summary['absorbed']} absorbed crossings, {summary['excluded']} boxes excluded (disc), {summary['witnesses']} witnesses (worst gap "
          f"{mp.nstr(summary['worst_witness_gap'], 3)}), {summary['seconds']} s"
          + (f"\n       reason: {summary['reason']}" if summary["reason"] else ""), flush=True)


def grid(box, n=3):
    (a, b), (c, d) = box.r, box.i
    return [mp.mpc(a + (b - a) * (k + mp.mpf("0.37")) / n, c + (d - c) * (j + mp.mpf("0.61")) / n) for k in range(n) for j in range(n)]


def witnesses(domains):
    pts = []
    for d in domains:
        if "y" in d:
            pts += [{"x": x, "y": y} for x in grid(d["x"], 2) for y in grid(d["y"], 3)]
        else:
            pts += [{"x": x} for x in grid(d["x"], 4)]
    return pts


def exact_checks(name, tree, points, expected):
    """Exact judge (SymPy) at the points where intervals can say nothing."""
    good = []
    for p in points:
        try:
            good.append(sp.simplify(E.exact(tree, p) - expected(**p)) == 0)
        except Exception as e:                       # noqa: BLE001 — we report, we do not hide
            print(f"       at {p}: {type(e).__name__} {e}", flush=True)
            good.append(False)
    EXPECTED.append(all(good))
    print(f"[{'OK' if all(good) else 'FAILED'}] {name}, exact computation (SymPy): {sum(good)}/{len(good)} points correct", flush=True)


def main(budget=600):
    ln_ref, add_ref = (lambda x: mp.log(x)), (lambda x, y: x + y)

    print("\n— Wall-free logarithm: ln_exact(z) = mirror(ln(mirror z)) —", flush=True)
    d = [{"x": B(-4, 4, "0.001", 8)}, {"x": B(-4, 4, -8, "-0.001")}, {"x": B("0.001", 4, -8, 8)}]
    t = [p for p in witnesses(d) if abs(p["x"] - 1) > mp.mpf("0.05")]
    report(C.certify("ln_exact above, below and to the right of the negative axis, outside the disc |z − 1| < 0.05",
                     E.ln_exact("x"), ln_ref, d, t, 14, budget, disc=(1, "0.05")))
    axis = [{"x": q} for q in (sp.Rational(-1, 3), sp.Integer(-2), sp.Integer(-7), sp.Rational(-22, 7))]
    exact_checks("ln_exact ON the negative axis (must equal ln|z| + iπ)", E.ln_exact("x"), axis, lambda x: sp.log(-x) + sp.I * sp.pi)
    exact_checks("ln_exact at z = 1 (center of the excluded disc; passes through ln 0 = −∞)", E.ln_exact("x"), [{"x": sp.Integer(1)}], lambda x: sp.Integer(0))
    wrong = [sp.simplify(E.exact(E.ln_("x"), p) - (sp.log(-p["x"]) + sp.I * sp.pi)) == 0 for p in axis]
    EXPECTED.append(not any(wrong))
    print(f"[{'OK' if not any(wrong) else 'FAILED'}] control: the short ln, for its part, is wrong on the axis at {len(wrong) - sum(wrong)}/{len(wrong)} points", flush=True)

    print("\n— Wall-free addition: add_safe(x, y) = y × (x/y + 1), on domains where x, y and x + y stay far from 0 —", flush=True)
    X = B(1, 3, -3, 3)
    d = [{"x": X, "y": B("0.5", 3, -8, 8)}, {"x": X, "y": B("-0.5", "0.5", "0.5", 8)},
         {"x": X, "y": B("-0.5", "0.5", -8, "-0.5")}, {"x": X, "y": B(-8, -4, -8, 8)}]
    report(C.certify("add_safe, y up to ±8i (well beyond the wall ±πi)", E.add_safe("x", "y"), add_ref, d, witnesses(d), 7, budget))

    print("\n— Negative controls (they MUST fail) —", flush=True)
    d = [{"x": X, "y": B("0.5", 3, 2, 4)}]
    report(C.certify("SHORT addition straddling the wall Im y = π", E.add("x", "y"), add_ref, d, [], 4, 60), expected=False)
    d = [{"x": X, "y": B("0.5", 3, 4, 6)}]
    b = C.certify("SHORT addition beyond the wall (no jump, but wrong value: it is the witness that unmasks it)",
                  E.add("x", "y"), add_ref, d, witnesses(d), 4, 60)
    ok = b["certified"] and abs(b["worst_witness_gap"] - 2 * mp.pi) < mp.mpf(10) ** -40
    EXPECTED.append(ok)
    print(f"[{'OK' if ok else 'FAILED'}] {b['name']}: no jump = {b['certified']}, gap at the witnesses = {mp.nstr(b['worst_witness_gap'], 8)} (= 2π)", flush=True)

    print("\n— Total addition (orthogonal axis): add_total(x, y) = i(−ix + Im y) + Re y —", flush=True)
    q = [{"x": B("0.25", 3, "0.25", 8)}, {"x": B("0.25", 3, -8, "-0.25")}, {"x": B(-3, "-0.25", "0.25", 8)}, {"x": B(-3, "-0.25", -8, "-0.25")}]
    report(C.certify("piece 1: wall-free real part, four quadrants up to ±8i", E.real_part_safe("x"), lambda x: mp.re(x), q, witnesses(q), 12, budget))
    report(C.certify("piece 2: imaginary part, four quadrants up to ±8i", E.imag_part("x"), lambda x: mp.im(x), q, witnesses(q), 12, budget))
    d = [{"x": B(1, "1.125", 0, "0.125"), "y": B(1, "1.125", 5, "5.125")}]
    report(C.certify("whole tree, small domain beyond the wall (Im y between 5 and 5.125)", E.add_total("x", "y"), add_ref, d, witnesses(d), 5, budget))
    zeros = [(0, 0), (0, sp.Rational(3, 2)), (sp.Rational(-5, 2), 0), (2, -2), (0, 5 * sp.I), (3 * sp.I, 0), (sp.Rational(1, 2), 7 * sp.I)]
    exact_checks("add_total at the points at zero or on the axes (where the intervals see nothing)", E.add_total("x", "y"),
                 [{"x": sp.sympify(a), "y": sp.sympify(c)} for a, c in zeros], lambda x, y: x + y)

    print(f"\nSUMMARY: {sum(EXPECTED)}/{len(EXPECTED)} expected verdicts obtained.")
    return 0 if all(EXPECTED) else 1


if __name__ == "__main__":
    sys.exit(main())
