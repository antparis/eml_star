"""The sine in NAND gates, by the method of pocket calculators (CORDIC, Volder 1959).

Idea: to find sin(angle), we turn a needle by smaller and smaller steps
(45°, 26.6°, 14°, 7.1°...), to the left or to the right, until it points at the requested angle. Each step
only needs additions, subtractions and bit shifts: nothing but NAND.

Numbers: fixed point in two's complement, `width` bits of which `fraction` after the point.
Zone: the bare method (`sine`) converges only for |angle| ≤ about 1.74 rad. `sine_full` first folds
the angle (full turns, half turn, quarter turn) and becomes correct for every representable angle. What remains is
the precision, limited by the number of bits (see `precision_curve`).
"""
import mpmath as mp

from . import nand as N


def _fixed(value, fraction, width):
    """Bits (least significant first) of the rounded integer value × 2^fraction, in two's complement."""
    n = int(mp.nint(mp.mpf(value) * 2 ** fraction))
    return [(n >> k) & 1 for k in range(width)]


def _shift(v, i):
    """Arithmetic right shift by i bits: pure wiring, zero gates."""
    return [v[k + i] if k + i < len(v) else v[-1] for k in range(len(v))]


def _add_or_sub(a, b, subtract):
    """a + b if subtract = 0, a − b if subtract = 1 (we invert b and add 1 through the carry-in)."""
    r, out = subtract, []
    for x, y in zip(a, b):
        s, r = N.full_adder(x, N.xor(y, subtract), r)
        out.append(s)
    return out


def _cordic(z, fraction, steps):
    """Core of the method: bits of sin(z) for |z| ≤ 1.74 (z in two's complement, `fraction` bits after the point)."""
    width = len(z)
    gain = mp.fprod(1 / mp.sqrt(1 + mp.mpf(4) ** -i) for i in range(steps))
    x, y, z = _fixed(gain, fraction, width), [0] * width, list(z)
    for i in range(steps):
        negative = z[-1]                     # sign bit: 1 if the needle has gone past the angle
        positive = N.not_(negative)
        x, y = _add_or_sub(x, _shift(y, i), positive), _add_or_sub(y, _shift(x, i), negative)
        z = _add_or_sub(z, _fixed(mp.atan(mp.mpf(2) ** -i), fraction, width), positive)
    return y


def sine(angle, fraction=13, steps=14):
    """Bare NAND circuit (without folding): correct only for |angle| ≤ 1.74."""
    return _cordic(angle, fraction, steps)


def _select(s, if_one, if_zero):
    """Switch: for each bit, if_one when s = 1, if_zero when s = 0."""
    ns = N.not_(s)
    return [N.nand(N.nand(a, s), N.nand(b, ns)) for a, b in zip(if_one, if_zero)]


def _negate_if(v, s):
    """−v if s = 1, v otherwise (two's complement)."""
    r, out = s, []
    for b in v:
        total, r = N.half_adder(N.xor(b, s), r)
        out.append(total)
    return out


def _subtract_if_possible(a, threshold, fraction):
    """a unsigned. Returns (a − threshold if a ≥ threshold, otherwise a; the bit "a ≥ threshold")."""
    n = len(a)
    d = _add_or_sub(a + [0], _fixed(threshold, fraction, n + 1), 1)
    enough = N.not_(d[-1])
    return _select(enough, d[:-1], a), enough


def sine_full(angle, fraction=13, guard=4, steps=None):
    """NAND circuit with folding of the angle, as in a real calculator: correct for EVERY representable
    angle. Steps: 1) sign set aside (sin(−a) = −sin a); 2) we remove the full turns of 2π;
    3) we remove a half turn (sin(a − π) = −sin a); 4) we fold onto the quarter turn (sin(π − a) = sin a);
    5) CORDIC; 6) we put the sign back. `guard` = extra bits of precision kept during the computation."""
    width, f = len(angle), fraction + guard
    if width < fraction + 2:
        raise ValueError("at least 2 bits are needed before the point")
    sign = angle[-1]
    a = [0] * guard + _negate_if(angle, sign)               # |angle|, unsigned, f bits after the point
    maxi, k = mp.mpf(2) ** (width - 1 - fraction), 0
    while 2 * mp.pi * 2 ** (k + 1) <= maxi:
        k += 1
    for j in range(k, -1, -1):
        if 2 * mp.pi * 2 ** j <= maxi:
            a, _ = _subtract_if_possible(a, 2 * mp.pi * 2 ** j, f)
    if mp.pi <= maxi:
        a, half_turn = _subtract_if_possible(a, mp.pi, f)
        sign = N.xor(sign, half_turn)
    n = len(a)
    d = _add_or_sub(a + [0], _fixed(mp.pi / 2, f, n + 1), 1)
    a = _select(N.not_(d[-1]), _add_or_sub(_fixed(mp.pi, f, n), a, 1), a)
    w = f + 3
    y = _cordic((a + [0] * w)[:w], f, steps or f + 1)
    if guard:                                              # round to nearest before dropping the guard bits
        y = _add_or_sub(y, _fixed(mp.mpf(2) ** -(fraction + 1), f, w), 0)[guard:]
    y = _negate_if(y, sign)                                # the sign last: sin(−a) = −sin(a) exactly
    return (y + [y[-1]] * width)[:width]


def compute(value, width=20, fraction=13, steps=None, full=True, guard=4):
    """(sine found by the circuit, number of gates, symbols without sharing) for an angle in radians."""
    wires = N.bits("t", width)
    circuit = sine_full(wires, fraction, guard, steps) if full else sine(wires, fraction, steps or fraction + 1)
    n = int(mp.nint(mp.mpf(value) * 2 ** fraction))
    if not -(1 << (width - 1)) <= n < 1 << (width - 1):
        raise ValueError(f"angle outside what {width} bits can write")
    raw = N.read_bits(circuit, N.write_bits("t", n, width))
    if raw >= 1 << (width - 1):
        raw -= 1 << width
    return mp.mpf(raw) / 2 ** fraction, N.gate_count(circuit), sum(N.tree_size(c) for c in circuit)


def precision_curve(settings=((12, 9), (16, 13), (24, 21), (32, 29)), points=41):
    """What precision costs: for each (bits in total, bits after the point), the number of gates and the
    worst error measured on `points` angles spread over the whole representable range."""
    rows = []
    for width, fraction in settings:
        bound = mp.mpf(2) ** (width - 1 - fraction) * mp.mpf("0.999")
        worst, count = 0, 0
        for k in range(points):
            x = -bound + 2 * bound * k / (points - 1)
            v, count, _ = compute(x, width, fraction)
            x_written = mp.nint(x * 2 ** fraction) / 2 ** fraction          # the angle actually given to the circuit
            worst = max(worst, abs(v - mp.sin(x_written)))
        rows.append((width, fraction, count, worst))
    return rows
