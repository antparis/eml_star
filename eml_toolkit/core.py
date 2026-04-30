"""
eml_toolkit.core
----------------
Opérateurs de base eml et eml★.

eml(x, y)  = exp(x) - ln(y)       [Odrzywołek 2026]
eml★(x, y) = exp(x) - ln(conj(y)) [Monnerot 2026]

Théorème 3.1 (vérifié à 50 décimales) :
    z̄ = 1 - eml★(0, eml(z, 1))   pour Im(z) ∈ [-π, π)
"""

import mpmath
import sympy as sp
from typing import Any, List, Optional

# Précision par défaut
mpmath.mp.dps = 50


def eml(x: Any, y: Any) -> Any:
    """eml(x, y) = exp(x) - ln(y)  [Odrzywołek 2026]"""
    if isinstance(x, (int, float, mpmath.mpf, mpmath.mpc)) and \
       isinstance(y, (int, float, mpmath.mpf, mpmath.mpc)):
        return mpmath.exp(x) - mpmath.log(y)
    if isinstance(x, sp.Expr) or isinstance(y, sp.Expr):
        return sp.exp(x) - sp.log(y)
    raise TypeError(f"eml: types non supportés ({type(x)}, {type(y)})")


def eml_star(x: Any, y: Any) -> Any:
    """eml★(x, y) = exp(x) - ln(conj(y))  [Monnerot 2026]"""
    if isinstance(x, (int, float, mpmath.mpf, mpmath.mpc)) and \
       isinstance(y, (int, float, mpmath.mpf, mpmath.mpc)):
        return mpmath.exp(x) - mpmath.log(mpmath.conj(y))
    if isinstance(x, sp.Expr) or isinstance(y, sp.Expr):
        return sp.exp(x) - sp.log(sp.conjugate(y))
    raise TypeError(f"eml★: types non supportés ({type(x)}, {type(y)})")


def conjugate_formula(z: Any) -> Any:
    """
    z̄ = 1 - eml★(0, eml(z, 1))
    Formule de conjugaison à profondeur 2 (Théorème 3.1).
    Valide pour Im(z) ∈ [-π, π).
    """
    return 1 - eml_star(0, eml(z, 1))


class EMLExpression:
    """
    Représentation symbolique d'une expression {eml, eml★, 1}.
    Permet l'optimisation DAG + CSE et l'évaluation haute précision.
    """

    def __init__(self, expr: Any, variables: Optional[List[str]] = None):
        self.expr = expr
        self.variables = variables or []
        self._optimized = False

    @classmethod
    def from_sympy(cls, sympy_expr, variables=None):
        if isinstance(sympy_expr, str):
            sympy_expr = sp.sympify(sympy_expr)
        vars_in_expr = [str(v) for v in sympy_expr.free_symbols]
        return cls(sympy_expr, variables or vars_in_expr)

    def to_sympy(self) -> sp.Expr:
        return self.expr

    @property
    def node_count(self) -> int:
        code = str(self.expr)
        return sum(code.count(op) for op in ['+', '*', '/', '**']) + 2

    def optimize(self, method: str = "egraph") -> "EMLExpression":
        from .optimizer import optimize_eml
        optimized_expr = optimize_eml(self.expr, method=method)
        new = EMLExpression(optimized_expr, self.variables)
        new._optimized = True
        return new

    def evaluate(self, **kwargs) -> Any:
        local_dict = {
            k: mpmath.mpc(v) if isinstance(v, complex) else mpmath.mpf(v)
            for k, v in kwargs.items()
        }
        return self.expr.evalf(subs=local_dict) if hasattr(self.expr, 'evalf') else self.expr

    def __repr__(self):
        status = "optimized" if self._optimized else "raw"
        return f"EMLExpression({status}, nodes≈{self.node_count}, vars={self.variables})"

    def __str__(self):
        return str(self.expr)
