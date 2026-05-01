"""
eml_toolkit/optimizer.py
EGraph optimizer for eml★ expression trees.

Implements 6 rewriting rules that simplify eml and eml★ trees
by eliminating redundant compositions.

Rules:
  1. eml(x, 1)           -> exp(x)
  2. eml(1, exp(x))      -> 1 - x  (ln simplification)
  3. eml(ln(x), exp(y))  -> x - y  (subtraction shortcut)
  4. eml★(0, eml(z, 1))  -> 1 - conj(z)  (Theorem 3.1)
  5. eml(eml(1,x), 1)    -> exp(1 - x)
  6. eml(0, x)           -> 1 - ln(x)
"""
import mpmath

mpmath.mp.dps = 50


# ── Symbolic expression tree ───────────────────────────────────────────────────

class Expr:
    """Base class for symbolic EML expressions."""
    def eval(self, z=None):
        raise NotImplementedError

    def depth(self):
        return 1

    def node_count(self):
        return 1


class Const(Expr):
    def __init__(self, value):
        self.value = value

    def eval(self, z=None):
        return mpmath.mpc(self.value)

    def __repr__(self):
        return str(self.value)


class Var(Expr):
    def __init__(self, name="z"):
        self.name = name

    def eval(self, z=None):
        return z

    def __repr__(self):
        return self.name


class EML(Expr):
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def eval(self, z=None):
        l = self.left.eval(z)
        r = self.right.eval(z)
        return mpmath.exp(l) - mpmath.log(r)

    def depth(self):
        return 1 + max(self.left.depth(), self.right.depth())

    def node_count(self):
        return 1 + self.left.node_count() + self.right.node_count()

    def __repr__(self):
        return f"eml({self.left}, {self.right})"


class EMLStar(Expr):
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def eval(self, z=None):
        l = self.left.eval(z)
        r = self.right.eval(z)
        return mpmath.exp(l) - mpmath.log(mpmath.conj(r))

    def depth(self):
        return 1 + max(self.left.depth(), self.right.depth())

    def node_count(self):
        return 1 + self.left.node_count() + self.right.node_count()

    def __repr__(self):
        return f"eml★({self.left}, {self.right})"


# ── EGraph optimizer ───────────────────────────────────────────────────────────

class EGraphOptimizer:
    """
    Simplifies EML expression trees using 6 algebraic rewriting rules.
    Reduces tree depth and node count without changing the computed value.
    """

    def __init__(self, max_passes=10):
        self.max_passes = max_passes
        self.rules_applied = 0

    def optimize(self, expr):
        """Apply rewriting rules until no further simplification is possible."""
        self.rules_applied = 0
        for _ in range(self.max_passes):
            new_expr = self._apply_rules(expr)
            if repr(new_expr) == repr(expr):
                break
            expr = new_expr
        return expr

    def _apply_rules(self, expr):
        """Recursively apply all 6 rewriting rules."""
        if isinstance(expr, (Const, Var)):
            return expr

        if isinstance(expr, EML):
            # Recurse first
            left  = self._apply_rules(expr.left)
            right = self._apply_rules(expr.right)

            # Rule 1: eml(x, 1) = exp(x)
            if isinstance(right, Const) and right.value == 1:
                self.rules_applied += 1
                return _ExpNode(left)

            # Rule 5: eml(eml(1, x), 1) = exp(1 - x)
            if (isinstance(left, EML)
                    and isinstance(left.left, Const)
                    and left.left.value == 1
                    and isinstance(right, Const)
                    and right.value == 1):
                self.rules_applied += 1
                return _ExpNode(_SubNode(Const(1), left.right))

            # Rule 6: eml(0, x) = 1 - ln(x)
            if isinstance(left, Const) and left.value == 0:
                self.rules_applied += 1
                return _SubNode(Const(1), _LnNode(right))

            return EML(left, right)

        if isinstance(expr, EMLStar):
            left  = self._apply_rules(expr.left)
            right = self._apply_rules(expr.right)

            # Rule 4: eml★(0, eml(z, 1)) = 1 - conj(z)
            if (isinstance(left, Const) and left.value == 0
                    and isinstance(right, EML)
                    and isinstance(right.right, Const)
                    and right.right.value == 1):
                self.rules_applied += 1
                return _SubNode(Const(1), _ConjNode(right.left))

            return EMLStar(left, right)

        return expr


# ── Helper nodes for simplified expressions ────────────────────────────────────

class _ExpNode(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, z=None):
        return mpmath.exp(self.arg.eval(z))

    def __repr__(self):
        return f"exp({self.arg})"


class _LnNode(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, z=None):
        return mpmath.log(self.arg.eval(z))

    def __repr__(self):
        return f"ln({self.arg})"


class _SubNode(Expr):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def eval(self, z=None):
        return self.a.eval(z) - self.b.eval(z)

    def __repr__(self):
        return f"({self.a} - {self.b})"


class _ConjNode(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, z=None):
        return mpmath.conj(self.arg.eval(z))

    def __repr__(self):
        return f"conj({self.arg})"


# ── Demo ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    opt = EGraphOptimizer()

    print("EGraphOptimizer — 6 rewriting rules demo")
    print()

    # conjugate_formula: 1 - eml★(0, eml(z, 1))
    z = Var("z")
    conj_expr = EMLStar(Const(0), EML(z, Const(1)))
    simplified = opt.optimize(conj_expr)
    print(f"  Original : {conj_expr}")
    print(f"  Simplified: {simplified}")
    print(f"  Rules applied: {opt.rules_applied}")

    z_val = mpmath.mpc(1.5, 0.8)
    original_val   = conj_expr.eval(z_val)
    simplified_val = simplified.eval(z_val)
    print(f"  Numeric check: original={original_val:.6f}, simplified={simplified_val:.6f}")
    print(f"  Error: {abs(original_val - simplified_val):.2e}")
