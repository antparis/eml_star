"""Branch-safety certificate by interval arithmetic (mpmath.iv) for eml / eml★ trees.

What a successful certificate proves, on the whole given domain (a product of complex boxes for x and y, possibly minus
a disc around a point of x): at every point, the tree can be evaluated and the computed function is
CONTINUOUS. To that end, every logarithm of the tree is checked:

- its argument avoids the branch cut (−∞, 0]: nothing to say;
- its argument crosses the negative real axis (away from 0) AND its node enters the parent through the X GATE (the first
  argument, the one that gets exponentiated): the 2πi jump is erased by the parent's exponential, because
  exp(c − log y) = e^c / y whatever the branch of the log. We then compute with a rotated branch,
  ln(−y) + iπ, whose cut avoids the box (difference 0 or 2πi from the principal log);
- any other case (crossing at the y gate, crossing at the root, 0 in the box): failure, we subdivide; definitive
  failure if the maximum depth is reached.

The CONSTANT subtrees (no variable) are computed separately, exactly, on the extended reals [−∞, +∞]:
ln 0 = −∞, exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞; ∞ − ∞ and ln(−∞) are refused. Equality to 0 is decided by
SymPy (e.g. E_zero = e − ln(e^e) = 0). A finite nonzero constant is recomputed node by node in interval arithmetic
(directed rounding); if one of its logarithms is not strictly off the branch cut, it is refused rather than guessed.
The + and × trees use no such constant (counter FINITE_CONSTANTS).

Trees: exact copies of branch_safety_final.py (E_exp, E_ln, E_zero, E_neg, E_sub, E_add, E_mul), in tree form.
Origin of the work: CutLattice handover docs of 2026-11-19 (discovery of the crossing of ln(ln x) in E_mul)."""
import time

import mpmath as mp
import sympy as sp
from mpmath import iv

iv.dps = 50
mp.mp.dps = 50

PLUS_INF, MINUS_INF = "+inf", "-inf"


class Undefined(ValueError):
    pass


class Failure(Exception):
    """A logarithm could not be checked on this box."""


# ---------- trees ----------

ONE, X, Y = ("1",), ("x",), ("y",)


def eml(a, b):
    return ("eml", a, b)


def emls(a, b):
    return ("eml*", a, b)


def E_exp(z):    return eml(z, ONE)
def E_ln(z):     return eml(ONE, eml(eml(ONE, z), ONE))
def E_zero():    return eml(ONE, eml(eml(ONE, ONE), ONE))
def E_neg(z):    return eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))
def E_sub(a, b): return eml(E_ln(a), E_exp(b))
def E_add(a, b): return E_sub(a, E_neg(b))
def E_mul(a, b): return E_exp(E_add(E_ln(a), E_ln(b)))


def minus_ln(a, c):
    """Node "a − ln c". Used ONLY to replace eml(E_ln(a), c), by virtue of the lemma (proved by hand, see
    FINDINGS_20261119_branch_safety_certificate.md): exp(E_ln(w)) = w for every w ≠ 0, whatever the branch, and
    = 0 = w at w = 0 with the conventions of the extended reals."""
    return ("-ln", a, c)


def E_mul_lemma(a, b):
    """E_mul(a, b) where the only passage through ln(ln a) is replaced thanks to the lemma: same value, without that passage."""
    return E_exp(minus_ln(E_ln(a), E_exp(E_neg(E_ln(b)))))


def variables(t):
    if t[0] in ("x", "y"):
        return {t[0]}
    if t[0] == "1":
        return set()
    return variables(t[1]) | variables(t[2])


def size(t):
    """Number of symbols in reverse Polish notation (leaves + operators)."""
    return 1 if len(t) == 1 else 1 + size(t[1]) + size(t[2])


def text(t):
    if len(t) == 1:
        return t[0]
    return f"{t[0]}({text(t[1])}, {text(t[2])})"


# ---------- exact constants on the extended reals ----------

_cache = {}


def constant(t):
    """Exact value of a subtree without variables: finite SymPy expression, PLUS_INF or MINUS_INF."""
    if t in _cache:
        return _cache[t]
    if t == ONE:
        v = sp.Integer(1)
    else:
        a, b = constant(t[1]), constant(t[2])
        if b == MINUS_INF:
            raise Undefined(f"ln(−∞) in {text(t)}")
        if isinstance(b, str):                      # b = +∞: ln b = +∞
            if a == PLUS_INF:
                raise Undefined(f"∞ − ∞ in {text(t)}")
            v = MINUS_INF
        elif b == 0:                                # ln 0 = −∞
            v = PLUS_INF
        elif a == PLUS_INF:
            v = PLUS_INF
        else:
            ea = sp.Integer(0) if a == MINUS_INF else sp.exp(a)
            arg = sp.conjugate(b) if t[0] == "eml*" else b
            v = sp.simplify(ea - sp.log(arg))
    _cache[t] = v
    return v


FINITE_CONSTANTS = [0]      # counter: how many nonzero finite constants have been enclosed (none for + and ×)


def _constant_interval(t):
    """GUARANTEED enclosure of a constant subtree. ±∞ and 0 come from the exact computation (SymPy). A nonzero finite
    constant is recomputed node by node in interval arithmetic (directed rounding); if one of its logarithms
    is not strictly off the branch cut, we refuse (Failure) rather than guess."""
    v = constant(t)
    if isinstance(v, str):
        return v
    if v == 0:
        return iv.mpc(0, 0)
    if t == ONE:
        return iv.mpc(1, 0)
    FINITE_CONSTANTS[0] += 1
    a, b = _constant_interval(t[1]), _constant_interval(t[2])
    if isinstance(a, str) or isinstance(b, str):
        raise Failure(f"finite constant going through infinity, not handled: {text(t)[:80]}")
    arg = iv.mpc(b.real, -b.imag) if t[0] == "eml*" else b
    if _contains_zero(arg) or not _off_cut(arg):
        raise Failure(f"constant with a logarithm touching the branch cut, not handled: {text(t)[:80]}")
    return iv.exp(a) - iv.log(arg)


# ---------- boxes ----------

def _bound(v, downward):
    """A decimal bound given as TEXT ("0.1") is rounded towards the outside of the box, so that the computed box
    contains the stated box. A Python number (float, int) is taken as is: its exact binary value."""
    if isinstance(v, str):
        return mp.mpf(iv.mpf(v)._mpi_[0 if downward else 1])
    return mp.mpf(v)


class Box:
    def __init__(self, re_min, re_max, im_min, im_max):
        self.r = (_bound(re_min, True), _bound(re_max, False))
        self.i = (_bound(im_min, True), _bound(im_max, False))

    def iv(self):
        return iv.mpc(iv.mpf([self.r[0], self.r[1]]), iv.mpf([self.i[0], self.i[1]]))

    def point(self, re, im):
        p = Box(re, re, im, im)
        return p

    def split(self):
        mr, mi = (self.r[0] + self.r[1]) / 2, (self.i[0] + self.i[1]) / 2
        return [Box(self.r[0], mr, self.i[0], mi), Box(mr, self.r[1], self.i[0], mi),
                Box(self.r[0], mr, mi, self.i[1]), Box(mr, self.r[1], mi, self.i[1])]

    def in_disc(self, center, radius):
        """Is the whole box, in a GUARANTEED way, inside the open disc? (a disc is convex: the four corners
        suffice.) Comparison of the squares in interval arithmetic: upper bound of distance² < lower bound of
        radius². A decimal radius is given as text ("0.05")."""
        c = mp.mpc(center)
        r = iv.mpf(radius) if isinstance(radius, str) else iv.mpf(mp.mpf(radius))
        r2 = r * r
        for a in self.r:
            for b in self.i:
                d2 = (iv.mpf(a) - iv.mpf(c.real)) ** 2 + (iv.mpf(b) - iv.mpf(c.imag)) ** 2
                if not d2.b < r2.a:
                    return False
        return True

    def __repr__(self):
        f = lambda v: mp.nstr(v, 6)
        return f"[{f(self.r[0])}, {f(self.r[1])}] + i[{f(self.i[0])}, {f(self.i[1])}]"


# ---------- interval evaluation ----------

def _contains_zero(z):
    return z.real.a <= 0 <= z.real.b and z.imag.a <= 0 <= z.imag.b


def _off_cut(z):
    return z.real.a > 0 or z.imag.a > 0 or z.imag.b < 0


def rotated_log_point(y):
    """Rotated branch of the log, cut on the positive real axis: ln(−y) + iπ."""
    return mp.log(-y) + 1j * mp.pi


def _ev(t, env, gate, count):
    if t == ONE:
        return iv.mpc(1, 0)
    if t in (X, Y):
        return env[t[0]]
    if not variables(t):
        return _constant_interval(t)
    if t[0] == "-ln":                              # a − ln c: a 2πik jump in a becomes a jump of the node's value,
        a = _ev(t[1], env, gate, count)            # allowed exactly when the node itself is at the x gate
        b = _ev(t[2], env, "y", count)
        if isinstance(a, str) or isinstance(b, str):
            raise Failure(f"infinite value in {text(t)[:80]}")
        if _contains_zero(b):
            raise Failure(f"log of a box that contains zero (node {text(t)[:80]})")
        if _off_cut(b):
            return a - iv.log(b)
        if gate == "x":
            count["absorbed"] += 1
            return a - (iv.log(iv.mpc(-b.real, -b.imag)) + iv.mpc(0, iv.pi))
        raise Failure(f"branch-cut crossing in a − ln c (node {text(t)[:80]})")
    a = _ev(t[1], env, "x", count)
    b = _ev(t[2], env, "y", count)
    if isinstance(b, str):
        if b == MINUS_INF:
            raise Failure(f"ln(−∞) in {text(t)[:80]}")
        if a == PLUS_INF:
            raise Failure(f"∞ − ∞ in {text(t)[:80]}")
        return MINUS_INF
    if a == PLUS_INF:
        return PLUS_INF
    ea = iv.mpc(0, 0) if a == MINUS_INF else iv.exp(a)
    arg = iv.mpc(b.real, -b.imag) if t[0] == "eml*" else b
    if _contains_zero(arg):
        raise Failure(f"log of a box that contains zero (node {text(t)[:80]})")
    if _off_cut(arg):
        return ea - iv.log(arg)
    if gate == "x":                                  # crossing absorbed by the parent's exp
        count["absorbed"] += 1
        return ea - (iv.log(iv.mpc(-arg.real, -arg.imag)) + iv.mpc(0, iv.pi))
    where = {"root": "at the root", "direct": "in a value used as is"}.get(gate, "at the y gate")
    raise Failure(f"branch-cut crossing {where} (node {text(t)[:80]})")


def enclose(t, boxes):
    """Enclosure of the value of the tree on the given boxes (raises Failure if a log is not checked)."""
    r = _ev(t, {k: v.iv() for k, v in boxes.items()}, "root", {"absorbed": 0})
    if isinstance(r, str):
        raise Failure(f"infinite value at the root ({r})")
    return r


def point_value(t, x, y=None):
    """Ordinary value (mpmath, principal log) at a point: to verify that the tree does compute its function."""
    def ev(t):
        if t == ONE:
            return mp.mpc(1)
        if t == X:
            return mp.mpc(x)
        if t == Y:
            return mp.mpc(y)
        if not variables(t):
            v = constant(t)
            return {PLUS_INF: mp.inf, MINUS_INF: -mp.inf}.get(v, v) if isinstance(v, str) else mp.mpc(complex(sp.N(v, 60)))
        a, b = ev(t[1]), ev(t[2])
        if t[0] == "-ln":
            return a - mp.log(b)
        ea = mp.mpf(0) if a == -mp.inf else mp.exp(a)
        if b == mp.inf:
            return -mp.inf
        if b == 0:
            return mp.inf
        return ea - mp.log(mp.conj(b) if t[0] == "eml*" else b)
    return ev(t)


# ---------- branch-and-bound ----------

def certify(t, domain, depth=10, budget=600, disc_x=None):
    """Splits the domain until every box passes. disc_x = (center, radius): boxes of x entirely inside
    this open disc are set aside (counted in "excluded"); they come under a separate proof.
    Returns {certified, boxes, absorbed, excluded, max_depth, seconds, reason}."""
    names = sorted(domain)
    if disc_x and domain["x"].in_disc(*disc_x):
        raise ValueError("the whole domain of x is inside the excluded disc: there is nothing to certify")
    count = {"absorbed": 0}
    work, ok, excluded, dmax, t0 = [(tuple(domain[n] for n in names), 0)], 0, 0, 0, time.time()
    while work:
        boxes, d = work.pop()
        env = dict(zip(names, boxes))
        if disc_x and env["x"].in_disc(*disc_x):
            excluded += 1
            continue
        try:
            r = _ev(t, {k: v.iv() for k, v in env.items()}, "root", count)
            if isinstance(r, str):
                raise Failure(f"infinite value at the root ({r})")
            ok += 1
            dmax = max(dmax, d)
        except Failure as e:
            if d >= depth or time.time() - t0 > budget:
                return {"certified": False, "boxes": ok, "absorbed": count["absorbed"], "excluded": excluded, "max_depth": d,
                        "seconds": round(time.time() - t0, 1), "reason": f"{e}; box {env}"}
            pieces = [b.split() for b in boxes]
            combos = [()]
            for m in pieces:
                combos = [c + (p,) for c in combos for p in m]
            work += [(c, d + 1) for c in combos]
    if ok == 0:                                      # nothing was verified: this is not a success
        return {"certified": False, "boxes": 0, "absorbed": 0, "excluded": excluded, "max_depth": 0,
                "seconds": round(time.time() - t0, 1), "reason": "no box accepted (is everything inside the excluded disc?)"}
    return {"certified": True, "boxes": ok, "absorbed": count["absorbed"], "excluded": excluded, "max_depth": dmax,
            "seconds": round(time.time() - t0, 1), "reason": ""}


def E_add_lemma(a, b):
    """E_add(a, b) = eml(E_ln(a), c) rewritten a − ln c by the absorption lemma (a = 0 included)."""
    return minus_ln(a, E_exp(E_neg(b)))


_ALL = []


def _report(name, r):
    _ALL.append(bool(r["certified"]))
    print(f"{name}: {'CERTIFIED' if r['certified'] else 'FAILED'} — {r['boxes']} boxes, {r['absorbed']} absorbed crossings, "
          f"{r['excluded']} in the disc, depth {r['max_depth']}, {r['seconds']} s" + (f" | {r['reason']}" if r["reason"] else ""))


if __name__ == "__main__":
    # The decimal bounds are given as text: rounded outward, the computed boxes CONTAIN
    # the stated boxes; the excluded disc is CONTAINED in the stated disc.
    print(f"Sizes: E_add {size(E_add(X, Y))}, E_mul {size(E_mul(X, Y))} symbols")
    D = Box("0.1", 3, -1, 1)
    _report("1) x + y, x,y ∈ [0.1,3]+i[−1,1]", certify(E_add(X, Y), {"x": D, "y": D}))
    _report("2) x × y (original tree), same domain, outside |x−1| < 0.05", certify(E_mul(X, Y), {"x": D, "y": D}, depth=10, budget=3000, disc_x=(1, "0.05")))
    _report("3) x × y (lemma form), same domain, disc included", certify(E_mul_lemma(X, Y), {"x": D, "y": D}))
    R = [Box(-3, "-0.1", -3, 3), Box("0.1", 3, -3, 3), Box(-3, 3, "0.1", 3), Box(-3, 3, -3, "-0.1")]
    for name, t in (("4) x + y", E_add(X, Y)), ("5) x × y (lemma form)", E_mul_lemma(X, Y))):
        tot = ab = dm = 0
        ok = True
        for bx in R:
            for by in R:
                r = certify(t, {"x": bx, "y": by}, depth=9)
                tot, ab, dm, ok = tot + r["boxes"], ab + r["absorbed"], max(dm, r["max_depth"]), ok and r["certified"]
        _ALL.append(ok)
        print(f"{name} on [−3,3]+i[−3,3] with |Re|,|Im| < 0.1 removed: {'CERTIFIED' if ok else 'FAILED'} — {tot} boxes, {ab} absorbed crossings, depth {dm}")
    P = Box(-3, 3, -3, 3)
    _report("6) x + y (lemma form) on the full square [−3,3]+i[−3,3], zero included", certify(E_add_lemma(X, Y), {"x": P, "y": P}, depth=9))
    print(f"Nonzero finite constants enclosed during these certificates: {FINITE_CONSTANTS[0]}")
    import sys
    if not all(_ALL) or len(_ALL) != 6:
        print("AT LEAST ONE EXPECTED CERTIFICATE FAILED")
        sys.exit(1)
    print("The six expected certificates passed.")
