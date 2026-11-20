"""EML side: a single operation, eml(x, y) = exp(x) − ln(y), and the number 1 are enough to build addition,
multiplication, the logarithm, the sine, i, π... (Odrzywołek 2026).

A tree: the leaf '1', a variable ('x', 'y'), or a node eml(left, right).
ln is the principal logarithm (imaginary part in ]−π, π]). Extended reals: ln 0 = −∞, exp(−∞) = 0,
exp(+∞) = +∞, ln(+∞) = +∞; ∞ − ∞ and ln(−∞) are undefined.

WARNING, the wall: each recipe is correct only in a zone (see DOMAINS and the function `measure` in
toolbox.py). Outside that zone the tree returns a value that is wrong by a whole number of turns 2πi, without warning.
"""
from functools import lru_cache

import mpmath as mp
import sympy as sp

from .tree import node

mp.mp.dps = 50
ONE = "1"
PLUS_INF, MINUS_INF = "+inf", "-inf"


class Undefined(ArithmeticError):
    pass


def eml(a, b):
    return node("eml", a, b)


def eml_star(a, b):
    """eml★(x, y) = exp(x) − ln(conj y): the same brick, with a reflection on its second ingredient."""
    return node("eml*", a, b)


# ── Recipes (the lengths in brackets are measured by size()) ─────────────────
def exp_(z):                      # [len(z) + 2]
    return eml(z, ONE)


def ln_(z):                       # [len(z) + 6]
    return eml(ONE, eml(eml(ONE, z), ONE))


ZERO = ln_(ONE)                   # e − ln(e^e) = 0
MINUS_INFINITY = ln_(ZERO)        # ln 0 = −∞


def neg(z):
    return eml(MINUS_INFINITY, exp_(z))     # exp(−∞) − ln(e^z) = −z


def sub(a, b):
    return eml(ln_(a), exp_(b))             # e^(ln a) − ln(e^b) = a − b


def add(a, b):
    return sub(a, neg(b))


def inv(z):
    return exp_(neg(ln_(z)))


def mul(a, b):
    return exp_(add(ln_(a), ln_(b)))


def div(a, b):
    return mul(a, inv(b))


def power(a, b):
    return exp_(mul(b, ln_(a)))


MINUS_ONE = neg(ONE)
TWO = add(ONE, ONE)
HALF = inv(TWO)
MINUS_I_PI = sub(eml(ZERO, MINUS_ONE), ONE)  # 1 − ln(−1) − 1 = −iπ: the only intended passage ON the wall
MINUS_I = exp_(mul(MINUS_I_PI, HALF))        # e^(−iπ/2)
I = neg(MINUS_I)
PI = mul(MINUS_I_PI, I)


def sqrt_(z):
    return exp_(mul(HALF, ln_(z)))


def sine(x):
    ix = mul(I, x)
    return div(sub(exp_(ix), exp_(neg(ix))), mul(TWO, I))


def cosine(x):
    ix = mul(I, x)
    return div(add(exp_(ix), exp_(neg(ix))), TWO)


# ── The mirror (eml★): what eml alone cannot do ──────────────────────────────
def mirror(z):
    """conj(z) in 9 symbols: 1 1 z eml★ 1 eml eml 1 eml. Correct on the whole plane (0 included, extended reals);
    9 is the minimum (report § 16)."""
    return eml(eml(ONE, eml(eml_star(ONE, z), ONE)), ONE)


def modulus_squared(z):
    return mul(z, mirror(z))


def modulus(z):
    return sqrt_(modulus_squared(z))


def real_part(z):
    return div(add(z, mirror(z)), TWO)


# ── "Wall-free" recipes: longer, but correct where the short ones go wrong ───────────────────────────────────
# Principle: multiplication is correct everywhere (away from zero), and adding 1 is always safe.
def ln_exact(z):
    """ln z, correct on the negative reals too (ln(−2) = ln 2 + iπ): we look in the mirror, take the ln,
    and look in the mirror again. Remaining hole: z = 0 returns −∞ (extended reals)."""
    return mirror(ln_(mirror(z)))


def negate(z):
    """−z by multiplication by −1: no condition on Im z (z nonzero)."""
    return mul(MINUS_ONE, z)


def add_safe(x, y):
    """x + y = y × (x/y + 1): no wall. Remaining hole: y = 0 (division by zero)."""
    return mul(y, add(div(x, y), ONE))


def sub_safe(x, y):
    return add_safe(x, negate(y))


def sqrt_safe(z):
    """Principal square root, correct on the negative reals too (√−4 = 2i)."""
    return exp_(mul(HALF, ln_exact(z)))


def real_part_safe(z):
    """Re z = ½ ln |e^z|²: no wall."""
    return mul(HALF, ln_(modulus_squared(exp_(z))))


def imag_part(z):
    """Im z = Re(−i z)."""
    return real_part_safe(mul(MINUS_I, z))


def add_total(x, y):
    """x + y by the ORTHOGONAL AXIS: we split y into real part and imaginary part (two real numbers), and we
    only ever add reals — and adding a real is always safe (the rule of addition only looks
    at the imaginary part of what is added).   x + y = i·(−i·x + Im y) + Re y.   No wall, no hole at zero."""
    return add(mul(I, add(mul(MINUS_I, x), imag_part(y))), real_part_safe(y))


def sub_total(x, y):
    """x − y = i·(−i·x − Im y) − Re y: subtracting a real is safe too."""
    return sub(mul(I, sub(mul(MINUS_I, x), imag_part(y))), real_part_safe(y))


def sine_safe(x):
    ix = mul(I, x)
    return div(sub_safe(exp_(ix), exp_(neg(ix))), mul(TWO, I))


def cosine_safe(x):
    ix = mul(I, x)
    return div(add_safe(exp_(ix), exp_(neg(ix))), TWO)


# ── Measuring ─────────────────────────────────────────────────────────────────
@lru_cache(maxsize=None)
def size(t):
    """Length in reverse Polish notation (leaves + nodes), without sharing."""
    return 1 if isinstance(t, str) else size(t.left) + size(t.right) + 1


@lru_cache(maxsize=None)
def variables(t):
    if isinstance(t, str):
        return frozenset() if t == ONE else frozenset([t])
    return variables(t.left) | variables(t.right)


@lru_cache(maxsize=None)
def constant(t):
    """EXACT value (SymPy) of a variable-free subtree, or PLUS_INF / MINUS_INF."""
    if t == ONE:
        return sp.Integer(1)
    a, b = constant(t.left), constant(t.right)
    if b == MINUS_INF:
        raise Undefined("ln(−∞)")
    if b == PLUS_INF:
        if a == PLUS_INF:
            raise Undefined("∞ − ∞")
        return MINUS_INF
    if b == 0 or a == PLUS_INF:
        return PLUS_INF
    ea = sp.Integer(0) if a == MINUS_INF else sp.exp(a)
    return sp.simplify(ea - sp.log(sp.conjugate(b) if t.op == "eml*" else b))


def exact(t, values):
    """EXACT value (SymPy) of the tree when the variables are exact numbers (rationals, i, π...).
    Serves as the judge at the points where intervals can say nothing: on the cut, or at zero."""
    memo = {}

    def ev(u):
        if isinstance(u, str):
            return sp.Integer(1) if u == ONE else sp.sympify(values[u])
        if id(u) in memo:
            return memo[id(u)]
        if not variables(u):
            v = constant(u)
        else:
            a, b = ev(u.left), ev(u.right)
            if b == MINUS_INF:
                raise Undefined("ln(−∞)")
            if b == PLUS_INF:
                if a == PLUS_INF:
                    raise Undefined("∞ − ∞")
                v = MINUS_INF
            elif b == 0 or a == PLUS_INF:
                v = PLUS_INF
            else:
                ea = sp.Integer(0) if a == MINUS_INF else sp.exp(a)
                v = sp.simplify(ea - sp.log(sp.conjugate(b) if u.op == "eml*" else b))
        memo[id(u)] = v
        return v

    return ev(t)


@lru_cache(maxsize=None)
def constants_finite(t):
    """Do all the nodes of this constant subtree have a finite value?"""
    if isinstance(t, str):
        return True
    return constant(t) not in (PLUS_INF, MINUS_INF) and constants_finite(t.left) and constants_finite(t.right)


@lru_cache(maxsize=None)
def gates_used(t):
    """Set of the node types of the tree ('eml', 'eml*')."""
    return frozenset() if isinstance(t, str) else frozenset([t.op]) | gates_used(t.left) | gates_used(t.right)


def _to_mp(v):
    if v == PLUS_INF:
        return mp.inf
    if v == MINUS_INF:
        return -mp.inf
    z = sp.N(v, mp.mp.dps + 10)
    re, im = z.as_real_imag()
    return mp.mpf(str(re)) if im == 0 else mp.mpc(mp.mpf(str(re)), mp.mpf(str(im)))


def _eml_num(a, b, star=False):
    if star and not mp.isinf(b):
        b = mp.conj(b)
    if b == 0:
        lb = -mp.inf
    elif mp.isinf(b):
        if mp.im(b) != 0 or b < 0:
            raise Undefined("ln(−∞)")
        lb = mp.inf
    else:
        lb = mp.log(b)
    ea = mp.mpf(0) if (mp.isinf(a) and mp.re(a) < 0) else mp.exp(a)
    if mp.isinf(ea) and mp.isinf(lb):
        raise Undefined("∞ − ∞")
    return ea - lb


def evaluate(t, env, strict=False):
    """Numerical value of the tree; the variable-free pieces are computed exactly first.
    strict=True: refuses (Undefined) as soon as a node, even a constant one, equals ±∞ — no passage through ln 0."""
    memo = {}

    def ev(u):
        if u in memo:
            return memo[u]
        if not variables(u):
            v = _to_mp(constant(u))
        elif isinstance(u, str):
            v = mp.mpmathify(env[u])
        else:
            v = _eml_num(ev(u.left), ev(u.right), u.op == "eml*")
        if strict and mp.isinf(v):
            raise Undefined("a node equals infinity (strict evaluation)")
        if strict and not isinstance(u, str) and not variables(u) and not constants_finite(u):
            raise Undefined("a constant passes through infinity (strict evaluation)")
        memo[u] = v
        return v

    return ev(t)


# Where each recipe is correct. Short recipes: rules certified by intervals (report § 1–15) or derived
# from them. Safe recipes: measured on grids of hard points (test_toolbox.py), not yet certified.
# In a "hole", the tool raises the Undefined error: it refuses instead of going wrong.
DOMAINS = {
    "addition": ("add", "correct if |Im y| < π", "add_safe", "correct on all measured points; holes: x = 0 or y = 0 (refusal)"),
    "subtraction": ("sub", "correct if |Im y| < π", "sub_safe", "correct on all measured points; holes: x = 0 or y = 0 (refusal)"),
    "total addition": ("add", "correct if |Im y| < π", "add_total", "correct on all measured points, zero included"),
    "multiplication": ("mul", "correct for all nonzero x, y", None, None),
    "logarithm": ("ln_", "wrong on the negative reals (−iπ instead of +iπ)", "ln_exact", "correct on all measured points (z nonzero)"),
    "square root": ("sqrt_", "wrong on the negative reals", "sqrt_safe", "correct on all measured points (z nonzero)"),
    "sine": ("sine", "correct on the real axis; complex: only if |Im x| < ln π ≈ 1.14",
             "sine_safe", "correct on the whole measured complex grid"),
    "cosine": ("cosine", "same zone as the short sine", "cosine_safe", "correct on the whole measured complex grid"),
    "real part": ("real_part", "correct if |Im z| < π", "real_part_safe", "correct on all measured points"),
    "mirror": ("mirror", "correct everywhere, 9 symbols (proved exactly, minimum)", None, None),
}


def recipe_table():
    """(name, short symbols, short zone, safe symbols or None, safe zone) for each recipe."""
    here, rows = globals(), []
    for name, (short, zone, safe, safe_zone) in DOMAINS.items():
        two = short in ("add", "sub", "mul")  # two-input recipes
        args = ("x", "y") if two else ("x",)
        rows.append((name, size(here[short](*args)), zone, size(here[safe](*args)) if safe else None, safe_zone))
    return rows
