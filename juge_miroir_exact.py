"""Juge exact (SymPy, sans tolérance numérique) du miroir de 9 symboles et de son témoin sans étoile.
Script proposé par ChatGPT le 2026-11-20, enregistré tel quel (seul ajout : cette docstring)."""
import sympy as sp
STAR = "1 1 z eml★ 1 eml eml 1 eml".split()
CONTROL = "1 1 z eml 1 eml eml 1 eml".split()
OPS = {"eml", "eml★"}
def require(condition, message):
    if not condition:
        raise SystemExit(f"FAIL: {message}")
def parse(tokens):
    stack = []
    for token in tokens:
        if token in {"1", "z"}:
            stack.append(token)
        elif token in OPS and len(stack) >= 2:
            right = stack.pop()
            left = stack.pop()
            stack.append((token, left, right))
        else:
            raise ValueError(f"Invalid token or stack underflow: {token}")
    require(len(stack) == 1, "invalid final stack")
    return stack[0]
def evaluate(tree, z, *, extended=False, trace=None):
    if tree == "1":
        return sp.Integer(1)
    if tree == "z":
        return z
    op, left, right = tree
    a = evaluate(left, z, extended=extended, trace=trace)
    b = evaluate(right, z, extended=extended, trace=trace)
    arg = sp.conjugate(b) if op == "eml★" else b
    if arg == 0:
        require(extended, "log(0) in finite evaluation")
        log_value = -sp.oo
    else:
        require(arg != -sp.oo, "log(-infinity) is undefined")
        log_value = sp.log(arg)
    value = sp.exp(a) - log_value
    require(not value.has(sp.nan, sp.zoo), "undefined intermediate value")
    if trace is not None:
        trace.append(value)
    return value
def check_equal(actual, expected, label):
    residual = sp.simplify(actual - expected)
    require(residual == 0, f"{label}: residual = {residual}")
require(len(STAR) == len(CONTROL) == 9, "RPN length")
require(sum(t in OPS for t in STAR) == 4, "operator count")
require(STAR.count("eml★") == 1, "star count")
require(["eml" if t == "eml★" else t for t in STAR] == CONTROL,
        "control differs by more than the star")
star, control = parse(STAR), parse(CONTROL)
# complex=True, zero=False allows nonreal values.
# SymPy's nonzero=True would also assume the symbol is real.
z = sp.Symbol("z", complex=True, zero=False)
require(z.is_real is None, "generic input was accidentally restricted to reals")
require(z.is_zero is False, "generic input must exclude zero")
for tree, expected, label in (
    (star, sp.conjugate(z), "starred"),
    (control, z, "control"),
):
    check_equal(evaluate(tree, z), expected, f"{label}, generic nonzero z")
# Distinguishes the two outputs even if a symbolic assumption is changed later.
check_equal(evaluate(star, sp.I), -sp.I, "starred at i")
check_equal(evaluate(control, sp.I), sp.I, "control at i")
# On the principal cut, the chosen pointwise value is Log(-1) = i*pi.
for name, tree in (("starred", star), ("control", control)):
    check_equal(evaluate(tree, -sp.Integer(1)), -1, f"{name} on the cut")
    trace = []
    value = evaluate(tree, sp.Integer(0), extended=True, trace=trace)
    require(trace == [sp.oo, sp.oo, -sp.oo, sp.Integer(0)],
            f"{name} at zero: unexpected trace {trace}")
    check_equal(value, 0, f"{name} at zero")
print("SymPy", sp.__version__)
print("PASS: same tree shape, 9 RPN symbols, one star replacement")
print("PASS: starred = conjugate(z), control = z for finite nonzero complex z")
print("PASS: at i the outputs are -i and +i; both return -1 at -1")
print("PASS: at zero both traces are [+oo, +oo, -oo, 0] under extended evaluation")
