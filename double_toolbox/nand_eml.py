"""The eml brick itself, built in NAND gates: eml(x, y) = exp(x) − Log(y), x and y complex.

Numbers: fixed point, two's complement, `integer_bits` bits before the point and `fraction` after, for each of the
four inputs (Re x, Im x, Re y, Im y) and of the two outputs.

Three workshops, all in additions, shifts and comparisons (hence in NAND):
  1. exp(a) for a real: we remove packets of ln 2 (each one = a one-bit shift), then some ln(1 + 2^-i)
     (each one = "add your own shifted self to yourself");
  2. the rotation: the CORDIC needle starts from the length exp(Re x) and turns by the angle Im x → (real part,
     imaginary part) of exp(x);
  3. Log(y): the same needle run backwards brings y back onto the axis (we read the angle = Arg y and the length = |y|),
     then ln|y| by the method of workshop 1 in reverse.

Accepted limits: no saturation (if exp(x) exceeds what the bits can write, the result overflows);
y = 0 is signalled by an "undefined" wire; no −∞ (the extended reals do not exist in fixed point).
"""
from functools import lru_cache

import mpmath as mp

from . import nand as N
from .nand_sine import _select, _shift, _fixed, _negate_if, _add_or_sub, _subtract_if_possible

LN2 = mp.log(2)


def _extend(v, n):
    return (v + [v[-1]] * n)[:n]


def _shift_left(v, k):
    return ([0] * k + v)[:len(v)]


def _or_all(wires):
    r = 0
    for b in wires:
        r = N.or_(r, b)
    return r


def _gain(steps):
    return mp.fprod(mp.sqrt(1 + mp.mpf(4) ** -i) for i in range(steps))


# ── Workshop 1: exp of a real ────────────────────────────────────────────────
def exp_real(a, f, integer_bits, offset=0):
    """exp(a + offset) > 0, on `integer_bits` + f bits. Below 2^-(f+2), returns 0."""
    n, bias_bits = len(a), f + 2
    bias = mp.mpf(offset) + bias_bits * LN2
    maxi = mp.mpf(2) ** (n - 1 - f) + abs(bias) + 1
    w = f + int(mp.ceil(mp.log(maxi, 2))) + 2
    ab = _add_or_sub(_extend(a, w), _fixed(bias, f, w), 0)
    too_small = ab[-1]
    packets = int(maxi / LN2) + 1
    m = {}
    for j in reversed(range(packets.bit_length())):
        ab, m[j] = _subtract_if_possible(ab, 2 ** j * LN2, f)
    r, e = ab[:f + 1], _fixed(1, f, f + 3)
    for i in range(1, f + 1):
        r, taken = _subtract_if_possible(r, mp.log(1 + mp.mpf(2) ** -i), f)
        e = _select(taken, _add_or_sub(e, _shift(e, i), 0), e)
    v = e + [0] * packets
    for j, taken in m.items():
        v = _select(taken, _shift_left(v, 2 ** j), v)
    keep = N.not_(too_small)
    return [N.and_(b, keep) for b in (v[bias_bits:] + [0] * (integer_bits + f))[:integer_bits + f]]


# ── Workshop 2: rotation ─────────────────────────────────────────────────────
def _fold(angle, f):
    """|angle| brought into [0, π/2]; returns (folded angle, sign to give to the sine, sign to give to the cosine)."""
    n, sign = len(angle), angle[-1]
    a, maxi = _negate_if(angle, sign), mp.mpf(2) ** (len(angle) - 1 - f)
    j = 0
    while 2 * mp.pi * 2 ** (j + 1) <= maxi:
        j += 1
    for k in range(j, -1, -1):
        if 2 * mp.pi * 2 ** k <= maxi:
            a, _ = _subtract_if_possible(a, 2 * mp.pi * 2 ** k, f)
    half = 0
    if mp.pi <= maxi:
        a, half = _subtract_if_possible(a, mp.pi, f)
    quarter = N.not_(_add_or_sub(a + [0], _fixed(mp.pi / 2, f, n + 1), 1)[-1])
    a = _select(quarter, _add_or_sub(_fixed(mp.pi, f, n), a, 1), a)
    return a, N.xor(sign, half), N.xor(half, quarter)


def exp_complex(xr, xi, f, integer_bits):
    """(Re, Im) of exp(xr + i·xi)."""
    steps = f + 1
    length = exp_real(xr, f, integer_bits + 1, offset=-mp.log(_gain(steps)))     # the gain of the needle is paid in advance
    a, s_sin, s_cos = _fold(xi, f)
    w = len(length) + 1
    x, y, z = _extend(length + [0], w), [0] * w, (a + [0] * (f + 3))[:f + 3]
    for i in range(steps):
        negative = z[-1]
        positive = N.not_(negative)
        x, y = _add_or_sub(x, _shift(y, i), positive), _add_or_sub(y, _shift(x, i), negative)
        z = _add_or_sub(z, _fixed(mp.atan(mp.mpf(2) ** -i), f, f + 3), positive)
    return _negate_if(x, s_cos), _negate_if(y, s_sin)


# ── Workshop 3: Log of a complex number ──────────────────────────────────────
def log_complex(u, v, f):
    """(ln|y|, Arg y, wire "y = 0") for y = u + i·v. Arg in ]−π, π]."""
    n = len(u)
    is_zero = N.not_(_or_all(u + v))
    su, sv = u[-1], v[-1]
    u, v = _negate_if(u, su), _negate_if(v, su)                 # half turn if Re y < 0: we bring it back to the right
    av, e = _negate_if(v, v[-1]), {}
    for j in reversed(range((n - 2).bit_length())):             # we enlarge y until it fills the bits (precision)
        k = 2 ** j
        if k > n - 2:
            continue
        free = N.not_(_or_all([N.or_(u[t], av[t]) for t in range(n - 1 - k, n - 1)]))
        u, v, av = (_select(free, _shift_left(q, k), q) for q in (u, v, av))
        e[j] = free
    w, steps = n + 2, f + 2
    x, y, z = _extend(u, w), _extend(v, w), [0] * (f + 3)
    for i in range(steps):
        negative = y[-1]
        positive = N.not_(negative)
        x, y = _add_or_sub(x, _shift(y, i), negative), _add_or_sub(y, _shift(x, i), positive)
        z = _add_or_sub(z, _fixed(mp.atan(mp.mpf(2) ** -i), f, f + 3), negative)
    half_turn = _select(sv, _fixed(-mp.pi, f, f + 3), _fixed(mp.pi, f, f + 3))
    angle = _add_or_sub(z, _select(su, half_turn, [0] * (f + 3)), 0)
    wl, halvings = f + 9, []
    for _ in range(2):                                          # mantissa m = x / 2^(n−2) brought into [1, 2[
        big = _or_all(x[n - 1:])
        x = _select(big, _shift(x, 1), x)
        halvings.append(big)
    ln_m = _fixed(LN2, f, wl)                                   # we push m up to 2: ln m = ln 2 − Σ ln(1 + 2^-i)
    for i in range(1, f + 2):
        t = _add_or_sub(x, _shift(x, i), 0)                     # m × (1 + 2^-i)
        fits = N.not_(_or_all(t[n - 1:]))                       # do we stay below 2?
        x = _select(fits, t, x)
        ln_m = _select(fits, _add_or_sub(ln_m, _fixed(mp.log(1 + mp.mpf(2) ** -i), f, wl), 1), ln_m)
    base = _fixed((n - 2 - f) * LN2 - mp.log(_gain(steps)), f, wl)
    for j, taken in e.items():
        base = _select(taken, _add_or_sub(base, _fixed(2 ** j * LN2, f, wl), 1), base)
    for big in halvings:
        base = _select(big, _add_or_sub(base, _fixed(LN2, f, wl), 0), base)
    return _add_or_sub(base, ln_m, 0), angle, is_zero


# ── The brick ────────────────────────────────────────────────────────────────
def brick(xr, xi, yr, yi, fraction, guard=4):
    """NAND circuit of eml(x, y). Returns (bits of Re, bits of Im, "undefined" wire when y = 0)."""
    n0, f = len(xr), fraction + guard
    integer_bits = n0 - fraction
    xr, xi, yr, yi = ([0] * guard + q for q in (xr, xi, yr, yi))
    c, s = exp_complex(xr, xi, f, integer_bits)
    modulus, angle, is_zero = log_complex(yr, yi, f)
    w = max(len(c), len(modulus), len(angle)) + 1
    half = _fixed(mp.mpf(2) ** -(fraction + 1), f, w) if guard else [0] * w
    outputs = []
    for a, b in ((c, modulus), (s, angle)):
        d = _add_or_sub(_add_or_sub(_extend(a, w), _extend(b, w), 1), half, 0)
        outputs.append(_extend(d[guard:], n0))
    return outputs[0], outputs[1], is_zero


@lru_cache(maxsize=None)
def _circuit(integer_bits, fraction, guard):
    n = integer_bits + fraction
    re, im, is_zero = brick(*(N.bits(q, n) for q in "abcd"), fraction, guard)
    return re, im, is_zero, N.gate_count(re + im + [is_zero])


def compute(x, y, integer_bits=8, fraction=16, guard=4):
    """(complex value returned by the circuit, undefined wire, number of gates) for two complex numbers x, y."""
    n = integer_bits + fraction
    re, im, is_zero, count = _circuit(integer_bits, fraction, guard)
    env = {}
    for q, val in zip("abcd", (mp.re(x), mp.im(x), mp.re(y), mp.im(y))):
        k = int(mp.nint(mp.mpf(val) * 2 ** fraction))
        if not -(1 << (n - 1)) <= k < 1 << (n - 1):
            raise ValueError(f"value {mp.nstr(val, 6)} outside what {n} bits can write")
        env.update(N.write_bits(q, k, n))

    def signed(bits):
        v = N.read_bits(bits, env)
        return mp.mpf(v - (1 << n) if v >= 1 << (n - 1) else v) / 2 ** fraction

    return mp.mpc(signed(re), signed(im)), N.evaluate(is_zero, env), count


def evaluate_with_nand(tree, env, integer_bits=10, fraction=16, guard=4):
    """Evaluates an EML tree using ONLY the NAND brick for each node that depends on a variable.
    Returns (value, number of bricks used, gates per brick). The constant pieces are precomputed
    (exactly) then rounded to the bits; −∞ is replaced by the smallest writable number, whose exponential
    falls exactly to 0 in the circuit. +∞ is not writable: the tree is then refused."""
    from . import eml as E
    floor, memo, count = -mp.mpf(2) ** (integer_bits - 1), {}, [0, 0]

    def ev(u):
        if isinstance(u, str) and u != E.ONE:
            return mp.mpmathify(env[u])
        if id(u) in memo:
            return memo[id(u)]
        if not E.variables(u):
            c = E.constant(u)
            if c == E.PLUS_INF:
                raise E.Undefined("+∞ is not writable in fixed point")
            v = floor if c == E.MINUS_INF else E._to_mp(c)
        else:
            d = ev(u.right)
            v, is_zero, count[1] = compute(ev(u.left), mp.conj(d) if u.op == "eml*" else d, integer_bits, fraction, guard)
            if is_zero:
                raise E.Undefined("the brick signals Log(0)")
            count[0] += 1
        if not isinstance(u, str):
            memo[id(u)] = v
        return v

    return ev(tree), count[0], count[1]
