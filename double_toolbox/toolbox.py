"""The double toolbox: the same questions asked at the same time to the two worlds.

  NAND  : a single gate, bits 0/1, n-bit integers. Exact, with no forbidden zone.
  EML   : a single operation, eml(x, y) = exp(x) − ln(y), and the number 1. Complex numbers.
          Correct in a zone; outside the zone, wrong by a turn 2πi without warning.

`compare(a, b)` does the addition and the multiplication on both sides at once and returns the sizes.
`measure(recipe, reference, points)` tells, point by point, where an EML recipe is correct.
"""
from dataclasses import dataclass

import mpmath as mp

from . import eml as E
from . import nand as N


@dataclass(frozen=True)
class Result:
    world: str
    operation: str
    value: object
    expected: object
    correct: bool
    bricks: int           # NAND gates (with sharing) / distinct eml nodes
    symbols: int          # length of the writing without sharing


def _distinct_nodes(t):
    seen = set()

    def v(u):
        if isinstance(u, str) or u in seen:
            return
        seen.add(u)
        v(u.left)
        v(u.right)

    v(t)
    return len(seen)


def nand_operation(op, a, b, n=8):
    A, B = N.bits("a", n), N.bits("b", n)
    outputs = N.addition(A, B) if op == "addition" else N.multiplication(A, B)
    inputs = {**N.write_bits("a", a, n), **N.write_bits("b", b, n)}
    value = N.read_bits(outputs, inputs)
    expected = a + b if op == "addition" else a * b
    return Result("NAND", op, value, expected, value == expected, N.gate_count(outputs),
                  sum(N.tree_size(s) for s in outputs if s is not None))


def eml_operation(op, a, b, tolerance=1e-30):
    t = E.add("x", "y") if op == "addition" else E.mul("x", "y")
    value = E.evaluate(t, {"x": a, "y": b})
    expected = mp.mpmathify(a) + b if op == "addition" else mp.mpmathify(a) * b
    correct = abs(value - expected) <= tolerance * max(1, abs(expected))
    return Result("EML", op, value, expected, correct, _distinct_nodes(t), E.size(t))


def compare(a, b, n=8):
    """Addition and multiplication of a and b, at the same time in the two worlds."""
    return [f(op, a, b) for op in ("addition", "multiplication")
            for f in (lambda o, x, y: nand_operation(o, x, y, n), eml_operation)]


def measure(tree, reference, points, tolerance=1e-30):
    """For each point: (point, value of the tree, expected value, correct?). An undefined computation counts as wrong."""
    out = []
    for p in points:
        env = p if isinstance(p, dict) else {"x": p}
        expected = reference(**{k: mp.mpmathify(v) for k, v in env.items()})
        try:
            v = E.evaluate(tree, env)
            ok = abs(v - expected) <= tolerance * max(1, abs(expected))
        except (E.Undefined, ZeroDivisionError, ValueError):
            v, ok = None, False
        out.append((p, v, expected, bool(ok)))
    return out


def compare_sine(x, width=20, fraction=13):
    """The same sine in the two worlds: NAND (CORDIC with folding, fixed point) and EML (wall-free tree)."""
    from . import nand_sine as S
    v, n_gates, symbols = S.compute(x, width, fraction)
    x_written = mp.nint(mp.mpf(x) * 2 ** fraction) / 2 ** fraction      # NAND only sees the angle rounded to its bits
    t = E.sine_safe("x")
    w = E.evaluate(t, {"x": x})
    return [Result("NAND", "sine", v, mp.sin(x_written), abs(v - mp.sin(x_written)) < 2 ** -(fraction - 1), n_gates, symbols),
            Result("EML", "sine", w, mp.sin(x), abs(w - mp.sin(x)) < 1e-30, _distinct_nodes(t), E.size(t))]
