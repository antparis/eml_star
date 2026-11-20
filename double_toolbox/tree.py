"""Shared nodes ("hash-consing"): two identical subtrees are the same object in memory.

Without this, comparing or memoizing large nested trees costs a time that blows up (each comparison
walks down the whole tree again). Here equality is identity, and a dictionary guarantees that the same construction
always returns the same object.
"""


class Node:
    __slots__ = ("op", "left", "right")

    def __init__(self, op, left, right):
        self.op, self.left, self.right = op, left, right

    def __repr__(self):
        return f"{self.op}(…)"


_TABLE = {}


def _key(x):
    return ("leaf", x) if isinstance(x, str) else ("node", id(x))


def node(op, left, right):
    key = (op, _key(left), _key(right))
    n = _TABLE.get(key)
    if n is None:
        n = _TABLE[key] = Node(op, left, right)
    return n
