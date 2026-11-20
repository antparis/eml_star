"""Which logarithms of the × tree can be moved without changing the result? (review objection of 2026-11-19:
"for every choice of branch at each logarithm" is false; counterexample: the log 1 of the last node.)
We number the logarithms of the expanded tree in evaluation order, we classify each one according to whether its argument
depends on the variables or not, and we shift ONE SINGLE logarithm by 2πi at a time."""
import mpmath as mp

mp.mp.dps = 60
TURN = 2j * mp.pi
X, Y, ONE = ("x",), ("y",), ("1",)
eml = lambda a, b: ("eml", a, b)
E_exp = lambda z: eml(z, ONE)
E_ln = lambda z: eml(ONE, eml(eml(ONE, z), ONE))
E_zero = lambda: eml(ONE, eml(eml(ONE, ONE), ONE))
E_neg = lambda z: eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))
E_sub = lambda a, b: eml(E_ln(a), E_exp(b))
E_add = lambda a, b: E_sub(a, E_neg(b))
E_mul = lambda a, b: E_exp(E_add(E_ln(a), E_ln(b)))


def has_variable(t):
    return t in (X, Y) or (len(t) == 3 and (has_variable(t[1]) or has_variable(t[2])))


def evaluate(t, x, y, shift=None, turns=None):
    """Evaluates the tree (principal log, extended reals). `shift` = number of the logarithm that receives + 2πi.
    Returns (value, list of the logs: (number, does the argument depend on the variables?, is it the log of the leaf 1?))."""
    logs = []

    def ev(t):
        if t == ONE:
            return mp.mpf(1)
        if t == X:
            return x
        if t == Y:
            return y
        a, b = ev(t[1]), ev(t[2])
        n = len(logs)
        logs.append((n, has_variable(t[2]), t[2] == ONE))
        if b == 0:
            lb = -mp.inf
        elif b == mp.inf:
            lb = mp.inf
        else:
            lb = mp.log(b)
        if shift == n and mp.isfinite(lb):
            lb = lb + TURN
        ea = mp.mpf(0) if a == -mp.inf else (mp.inf if a == mp.inf else mp.exp(a))
        if abs(mp.re(ea)) < mp.mpf(10) ** -55 and mp.im(ea) == 0 and not has_variable(t[1]) and t[1] != ONE:
            pass
        return ea - lb

    # E_zero equals 0 exactly: we force the exact value of this constant subtree (in floating point: ~1e-60)
    zero = E_zero()

    def ev2(t):
        if t == zero:
            n0 = len(logs)
            for k in range(3):                                   # its three logs, numbered so that they can be shifted
                logs.append((n0 + k, False, k != 2))
            # E_zero = eml(1, eml(eml(1, 1), 1)). EXACT value depending on the shifted logarithm (computed by hand):
            if shift == n0:           # ln 1 of the deepest node: e − 2πi, then exp(e − 2πi) − 0 = exp(e), then e − ln(exp e) = 0
                return mp.mpf(0)
            if shift == n0 + 1:       # ln 1 of the middle node: exp(e) − 2πi, then e − ln(exp(e) − 2πi): finite, nonzero
                return mp.e - mp.log(mp.exp(mp.e) - TURN)
            if shift == n0 + 2:       # ln(exp e) of the root of E_zero: e − (e + 2πi) = −2πi
                return mp.mpc(0) - TURN
            return mp.mpf(0)
        if t in (ONE, X, Y):
            return ev(t)
        a, b = ev2(t[1]), ev2(t[2])
        if a is None or b is None:
            return None
        n = len(logs)
        logs.append((n, has_variable(t[2]), t[2] == ONE))
        lb = -mp.inf if b == 0 else (mp.inf if b == mp.inf else mp.log(b))
        if shift == n and mp.isfinite(lb):
            lb = lb + TURN
        if turns and mp.isfinite(lb):
            lb = lb + TURN * turns.get(n, 0)               # several logarithms shifted AT THE SAME TIME
        ea = mp.mpf(0) if a == -mp.inf else (mp.inf if a == mp.inf else mp.exp(a))
        positions.append((n, gate.get(id(t), "root")))
        return ea - lb

    positions, gate = [], {}

    def mark(t):
        if len(t) == 3:
            gate[id(t[1])], gate[id(t[2])] = "x (exponentiated by the parent)", "y (its log is taken by the parent)"
            mark(t[1]); mark(t[2])
    mark(t)
    v = ev2(t)
    return (v, logs, dict(positions)) if turns is not None else (v, logs)


def effect_of_a_shift(n, x=mp.mpc(0.7, 0.3), y=mp.mpc(-1.3, 0.4)):
    """(result − x·y) when the single logarithm no. n of the × tree receives + 2πi."""
    v, _ = evaluate(E_mul(X, Y), x, y, shift=n)
    return v - x * y


if __name__ == "__main__":
    x, y = mp.mpc(0.7, 0.3), mp.mpc(-1.3, 0.4)
    t = E_mul(X, Y)
    base, logs = evaluate(t, x, y)
    print(f"× tree: {len(logs)} logarithms. Without shift: difference from x·y = {mp.nstr(abs(base - x * y), 3)}")
    summary = {}
    for n, variable, leaf1 in logs:
        v, _ = evaluate(t, x, y, shift=n)
        if not mp.isfinite(abs(v)):
            effect = "infinite result"
        else:
            d = (v - x * y)
            effect = "no effect" if abs(d) < mp.mpf(10) ** -40 else (f"result shifted by {mp.nstr(d / TURN, 4)} turn(s)" if abs(d / TURN - mp.nint((d / TURN).real)) < 1e-30 else "wrong result")
        kind = "argument is variable" if variable else ("log of the leaf 1" if leaf1 else "argument is constant")
        summary.setdefault((kind, effect), []).append(n)
    for (kind, effect), ns in sorted(summary.items()):
        print(f"   {kind:22s}: {effect:32s} — logs no. {ns}")
    v, _ = evaluate(t, x, y, shift=len(logs) - 1)
    print(f"Last node alone (log 1 of the root shifted by +2πi): result − x·y = {mp.nstr(v - x * y, 12)} (expected −2πi = {mp.nstr(-TURN, 12)})")

    # ---- 2) where do the 8 free logarithms go? and what happens when they are ALL moved together?
    import random
    _, logs, pos = evaluate(t, x, y, turns={})
    free = [n for n, variable, _ in logs if variable]
    print("Position, in its parent, of the node carrying each logarithm with a variable argument:")
    for n in free:
        print(f"   log no. {n:2d}: gate {pos[n]}")
    rng = random.Random(1119)
    worst, n_pts = mp.mpf(0), 0
    for _ in range(3000):
        xx = mp.mpc(rng.uniform(-3, 3), rng.uniform(-6, 6)); yy = mp.mpc(rng.uniform(-3, 3), rng.uniform(-6, 6))
        k = {n: rng.randint(-3, 3) for n in free}
        v, _, _ = evaluate(t, xx, yy, turns=k)
        worst = max(worst, abs(v - xx * yy) / abs(xx * yy)); n_pts += 1
    print(f"Random SIMULTANEOUS shifts (−3…3 turns) on the {len(free)} free logs, {n_pts} points: worst relative difference from x·y = {mp.nstr(worst, 3)}")
    for kk in (0, 1, -2):
        v, _, _ = evaluate(t, mp.mpf(1), mp.mpc(-1.3, 0.4), turns={n: kk for n in free})
        print(f"   at x = 1 exactly, all the free logs shifted by {kk:+d} turn(s): difference from x·y = {mp.nstr(abs(v - mp.mpc(-1.3, 0.4)), 3)}")
