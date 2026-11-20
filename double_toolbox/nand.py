"""NAND side: a single logic gate, NAND(a, b) = NOT (a AND b), is enough to build all computation on bits.

A circuit is a tree: a leaf is the name of an input ('a0', 'b3', ...), a node is ('nand', left, right).
Identical sub-circuits are shared (like a wire plugged in at two places): `gate_count` counts the
real gates of the circuit, `tree_size` counts what the same computation would cost written without sharing.
"""
import sys
from functools import lru_cache

from .tree import node

sys.setrecursionlimit(max(sys.getrecursionlimit(), 50000))      # deep circuits (chained carries)


def nand(a, b):
    """A wire can also be a constant 0 or 1 (Python integer): the gate is then simplified at
    construction (NAND(0, ·) = 1, NAND(1, b) = NOT b), like a wire soldered to ground or to the power supply."""
    if isinstance(a, int) or isinstance(b, int):
        if a == 0 or b == 0:
            return 1
        if isinstance(a, int) and isinstance(b, int):
            return 0
        other = b if isinstance(a, int) else a
        return node("nand", other, other)
    return node("nand", a, b)


# ── The basic gates, all in NAND ──────────────────────────────────────────────
def not_(a):
    return nand(a, a)


def and_(a, b):
    return not_(nand(a, b))


def or_(a, b):
    return nand(not_(a), not_(b))


def xor(a, b):
    if isinstance(a, int) or isinstance(b, int):          # a constant wire: at most one NOT gate remains
        if isinstance(a, int) and isinstance(b, int):
            return a ^ b
        c, wire = (a, b) if isinstance(a, int) else (b, a)
        return not_(wire) if c else wire
    m = nand(a, b)
    return nand(nand(a, m), nand(b, m))


def half_adder(a, b):
    """(sum, carry) of two bits: 5 gates."""
    m = nand(a, b)
    return nand(nand(a, m), nand(b, m)), not_(m)


def full_adder(a, b, r):
    """(sum, carry) of three bits: 9 gates, the classic layout."""
    m1 = nand(a, b)
    s1 = nand(nand(a, m1), nand(b, m1))
    m2 = nand(s1, r)
    total = nand(nand(s1, m2), nand(r, m2))
    return total, nand(m1, m2)


# ── Integers written in n bits (bit 0 = least significant) ────────────────────
def bits(name, n):
    return [f"{name}{k}" for k in range(n)]


def addition(a, b):
    """Addition of two n-bit numbers; returns n + 1 bits."""
    total, r = [], None
    for k, (x, y) in enumerate(zip(a, b)):
        s, r = half_adder(x, y) if r is None else full_adder(x, y, r)
        total.append(s)
    return total + [r]


def multiplication(a, b):
    """Multiplication by shifts and additions (as at school, in binary); returns len(a) + len(b) bits."""
    n, m = len(a), len(b)
    acc = [and_(a[i], b[0]) for i in range(n)]
    product = [acc[0]]
    acc = acc[1:]
    for j in range(1, m):
        row = [and_(a[i], b[j]) for i in range(n)]
        upper = acc + [None]
        new, r = [], None
        for i in range(n):
            x, y = row[i], upper[i]
            if y is None:
                s, r = (x, None) if r is None else half_adder(x, r)
            elif r is None:
                s, r = half_adder(x, y)
            else:
                s, r = full_adder(x, y, r)
            new.append(s)
        product.append(new[0])
        acc = new[1:] + [r]
    return product + [x for x in acc]


# ── Evaluating and measuring ──────────────────────────────────────────────────
def evaluate(circuit, inputs, _memo=None):
    """Value (0 or 1) of a circuit for inputs {name: bit}. A missing bit (None) counts as 0."""
    memo = {} if _memo is None else _memo
    if circuit is None:
        return 0
    if isinstance(circuit, int):
        return circuit
    if isinstance(circuit, str):
        return inputs[circuit]
    key = id(circuit)
    if key not in memo:
        memo[key] = 1 - (evaluate(circuit.left, inputs, memo) & evaluate(circuit.right, inputs, memo))
    return memo[key]


def write_bits(name, value, n):
    return {f"{name}{k}": (value >> k) & 1 for k in range(n)}


def read_bits(outputs, inputs):
    memo = {}
    return sum(evaluate(s, inputs, memo) << k for k, s in enumerate(outputs))


def gate_count(outputs):
    """Number of distinct NAND gates in the circuit (with wire sharing)."""
    seen = set()

    def visit(c):
        if c is None or isinstance(c, (str, int)) or c in seen:
            return
        seen.add(c)
        visit(c.left)
        visit(c.right)

    for s in outputs:
        visit(s)
    return len(seen)


@lru_cache(maxsize=None)
def tree_size(c):
    """Number of symbols if the computation is written without sharing (leaves + gates)."""
    if c is None or isinstance(c, (str, int)):
        return 1
    return tree_size(c.left) + tree_size(c.right) + 1
