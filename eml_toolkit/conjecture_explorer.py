"""
eml_toolkit.conjecture_explorer
--------------------------------
Explorateur heuristique de conjectures pour {eml, eml★, 1}.

Note : ceci est une exploration aléatoire parallèle (random search),
PAS un vrai Monte Carlo Tree Search (pas d'arbre UCB, pas de
backpropagation). Le nom a été corrigé pour honnêteté mathématique.

Utile pour :
- Chercher des identités inconnues dans le système eml/eml★
- Valider numériquement des conjectures à haute précision
- Explorer l'espace des expressions de profondeur ≤ max_depth

Auteur : Anthony Monnerot (2026)
"""

import multiprocessing
import random
import mpmath
import sympy as sp
from sympy import symbols, simplify, exp, log, conjugate

from .optimizer import EGraphOptimizer

s = symbols('s')

_EPS = mpmath.mpf('1e-40')


class ConjectureExplorer:
    """
    Explorateur heuristique parallèle pour conjectures eml/eml★.

    Génère aléatoirement des expressions candidates, les optimise
    via EGraphOptimizer, et vérifie si elles approchent une cible.

    Parameters
    ----------
    max_depth    : profondeur maximale des arbres eml générés
    num_rollouts : nombre de candidats à évaluer
    num_workers  : nombre de processus parallèles
    dps          : précision mpmath (décimales)
    """

    def __init__(
        self,
        max_depth: int = 5,
        num_rollouts: int = 100,
        num_workers: int = 4,
        dps: int = 50,
    ):
        self.optimizer = EGraphOptimizer()
        self.max_depth = max_depth
        self.num_rollouts = num_rollouts
        self.num_workers = num_workers
        self.dps = dps
        self.results = []

    def _random_eml_expression(self, depth: int = 0):
        """Génère une expression eml/eml★ aléatoire jusqu'à max_depth."""
        if depth >= self.max_depth:
            return s ** 2 + symbols('b') * s + symbols('c')

        op = random.choice(['eml', 'eml_star', 'conj', 'poly'])
        if op == 'eml':
            return exp(s) - log(s ** 2 + 1)
        elif op == 'eml_star':
            return exp(s) - log(conjugate(s ** 2 + 1))
        elif op == 'conj':
            # z̄ = 1 − eml★(0, eml(z, 1))
            return 1 - (exp(0) - log(conjugate(exp(s) - log(1))))
        else:
            return s ** 2 + symbols('b') * s + symbols('c')

    def _evaluate_candidate(self, args):
        """Évalue un candidat contre la cible (appelé en sous-processus)."""
        target, dps = args
        mpmath.mp.dps = dps
        candidate = self._random_eml_expression()
        try:
            optimized = self.optimizer.optimize(candidate)
            diff = simplify(candidate - optimized)
            if diff == 0:
                err = mpmath.fabs(mpmath.re(optimized - target))
                if err < _EPS:
                    return (str(candidate), str(optimized), float(err))
        except Exception:
            pass
        return None

    def explore(self, target):
        """
        Lance l'exploration parallèle.

        Parameters
        ----------
        target : expression sympy cible à approcher

        Returns
        -------
        list of (candidate_str, optimized_str, error_float)
            Les 10 meilleurs résultats triés par erreur croissante.
        """
        mpmath.mp.dps = self.dps
        self.results = []

        args_list = [(target, self.dps)] * self.num_rollouts

        with multiprocessing.Pool(self.num_workers) as pool:
            raw = pool.map(self._evaluate_candidate, args_list)

        self.results = sorted(
            [r for r in raw if r is not None],
            key=lambda x: x[2]
        )[:10]

        return self.results

    def report(self):
        """Affiche les résultats de l'exploration."""
        if not self.results:
            print("Aucun résultat trouvé.")
            return
        print(f"{'Rang':>4}  {'Erreur':>12}  Expression optimisée")
        print("-" * 60)
        for i, (cand, opt, err) in enumerate(self.results, 1):
            print(f"{i:>4}  {err:>12.3e}  {opt}")
