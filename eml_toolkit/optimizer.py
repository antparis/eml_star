"""
eml_toolkit.optimizer
---------------------
Moteur E-Graph pour expressions {eml, eml★, 1}.

6 règles de réécriture :
  1. rational_one_plus_delta    — réécriture fraction → 1 + δ
  2. horner_after_rewrite       — forme de Horner sur dénominateur
  3. factor_after_delta         — factorisation après réécriture
  4. cancel_common_terms        — annulation termes communs
  5. trig_identity_rule         — identités trigonométriques (sin²+cos²=1, etc.)
  6. exp_log_rule               — simplification exp/log

Auteur : Anthony Monnerot (2026), d'après collaboration avec Grok/Claude.
"""

import mpmath
import sympy as sp
from sympy import symbols, simplify, Add, Poly

s = symbols('s')

_EPS = mpmath.mpf('1e-40')


class EGraphOptimizer:
    """
    Optimiseur E-Graph pour expressions symboliques EML/eml★.

    Applique 6 règles de réécriture en boucle de saturation
    jusqu'à convergence (coût minimal atteint).
    """

    def __init__(self, eps: mpmath.mpf = _EPS):
        self.eps = eps
        self.rules = [
            self._rational_one_plus_delta,
            self._horner_after_rewrite,
            self._factor_after_delta,
            self._cancel_common_terms,
            self._trig_identity_rule,
            self._exp_log_rule,
        ]

    # ── Règles ───────────────────────────────────────────────────────────────

    def _rational_one_plus_delta(self, expr):
        """Règle 1 : fraction P/Q → 1 + δ·s/Q quand deg P = deg Q."""
        try:
            N, D = expr.as_numer_denom()
            if N.is_polynomial(s) and D.is_polynomial(s):
                pN = Poly(N, s)
                pD = Poly(D, s)
                if pN.degree() == pD.degree() == 2:
                    if abs(float(pN.LC() - pD.LC())) < float(self.eps):
                        delta_b = pN.coeff_monomial(s) - pD.coeff_monomial(s)
                        diff = N - D - delta_b * s
                        if abs(float(diff)) < float(self.eps) and delta_b != 0:
                            return 1 + (delta_b * s) / D
        except Exception:
            pass
        return expr

    def _horner_after_rewrite(self, expr):
        """Règle 2 : forme de Horner sur le dénominateur d'une fraction 1 + num/den."""
        if expr.is_Add and len(expr.args) == 2:
            one, frac = expr.args
            if one == 1 and frac.is_Mul:
                try:
                    num, den = frac.as_numer_denom()
                    if den.is_polynomial(s):
                        return 1 + num / sp.horner(den)
                except Exception:
                    pass
        return expr

    def _factor_after_delta(self, expr):
        """Règle 3 : factorisation du quotient après réécriture 1+δ."""
        if expr.is_Add and len(expr.args) == 2:
            one, frac = expr.args
            if one == 1 and frac.is_Mul:
                try:
                    num, den = frac.as_numer_denom()
                    return 1 + sp.factor(num / den)
                except Exception:
                    pass
        return expr

    def _cancel_common_terms(self, expr):
        """Règle 4 : annulation des termes communs entre numérateur et dénominateur."""
        try:
            N, D = expr.as_numer_denom()
            if N.is_Add and D.is_Add:
                common = set(N.args) & set(D.args)
                if common:
                    new_N = Add(*[t for t in N.args if t not in common])
                    new_D = Add(*[t for t in D.args if t not in common])
                    if new_D != 0:
                        return new_N / new_D
        except Exception:
            pass
        return expr

    def _trig_identity_rule(self, expr):
        """Règle 5 : identités trigonométriques (sin²+cos²=1, etc.)."""
        try:
            if expr.has(sp.sin, sp.cos):
                return sp.trigsimp(expr)
        except Exception:
            pass
        return expr

    def _exp_log_rule(self, expr):
        """Règle 6 : simplifications exp/log (exp(log(x))=x, etc.)."""
        try:
            if expr.has(sp.exp, sp.log):
                return sp.simplify(expr)
        except Exception:
            pass
        return expr

    # ── Boucle de saturation ─────────────────────────────────────────────────

    def optimize(self, expr, max_iterations: int = 30):
        """
        Applique les règles en boucle jusqu'à convergence.
        Retourne l'expression de coût minimal trouvée.
        """
        current = simplify(expr)
        best = current
        best_cost = self._cost(current)
        seen = set()

        for _ in range(max_iterations):
            improved = False
            for rule in self.rules:
                try:
                    candidate = rule(current)
                    candidate = simplify(candidate)
                except Exception:
                    continue
                key = str(candidate)
                if key not in seen:
                    seen.add(key)
                    c = self._cost(candidate)
                    if c < best_cost:
                        best = candidate
                        best_cost = c
                        current = candidate
                        improved = True
            if not improved:
                break

        return best

    def _cost(self, expr) -> int:
        """Coût = nombre d'opérations atomiques dans l'expression."""
        try:
            return sum(
                len(expr.atoms(op))
                for op in [sp.Add, sp.Mul, sp.Pow,
                           sp.sin, sp.cos, sp.exp, sp.log]
            )
        except Exception:
            return 9999


def optimize_eml(expr, method: str = "egraph"):
    """
    Interface simplifiée.
    method='egraph'  → EGraphOptimizer (défaut)
    method='basic'   → factor + horner + CSE sympy
    """
    if method == "egraph":
        return EGraphOptimizer().optimize(expr)

    # Méthode de base sans E-graph
    try:
        expr = sp.factor(expr)
    except Exception:
        pass
    expr = sp.horner(expr)
    replacements, reduced = sp.cse(expr, optimizations='basic')
    for old, new in replacements:
        reduced[0] = reduced[0].subs(old, new)
    return reduced[0]
