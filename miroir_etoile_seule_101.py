"""Miroir strict en eml★ seule, 101 symboles — script de vérification autonome (SymPy uniquement).

Découverte : ChatGPT, 20 novembre 2026 (voir FINDINGS_20261119_branch_safety_certificate.md, § 17).
Ce fichier rebâtit l'arbre sans dépendre d'aucun autre fichier du dépôt, l'écrit en notation postfixée, puis le juge.

    S(a, b) = eml★(a, b) = exp(a) − Log(conj b),   log principal, Arg dans ]−π, π].
    Feuilles : 1 et z. Aucune autre porte, aucune constante fournie.

Évaluation STRICTE : chaque nœud doit avoir une valeur finie ; Log(0) est refusé (erreur Indefini).

Ce que le script établit : la grammaire, la longueur, l'identité J(z) = conj(z) en calcul exact sur un symbole
complexe non nul et en des points exacts, le refus en zéro, et eml★(x, J(y)) = eml(x, y).
Ce qu'il n'établit PAS : qu'aucun nœud intermédiaire ne s'annule pour TOUT z non nul. Cela vient de l'argument de
domaine du § 17 (la partie imaginaire de S(b, z) est dans [−3π, −π[), à faire relire. Ni la minimalité de 101.
Code de sortie 0 seulement si tous les contrôles passent.   Usage : python3 miroir_etoile_seule_101.py
"""
import sys

import sympy as sp


class Indefini(ArithmeticError):
    pass


def S(a, b):
    return ("S", a, b)


def X(t):                      # exp(t)
    return S(t, "1")


def Q(t):                      # e^e / conj(t)
    return X(S("1", t))


def L(t):                      # un logarithme de t : exp(L(t)) = t
    return S("1", Q(t))


def d(t, r):                   # t − r, r réel
    return S(L(t), X(r))


E_ = X("1")                                    # e
O = S("1", X(E_))                              # 0
MOINS_UN = d(S("1", E_), E_)                   # −1
K = d(S(O, MOINS_UN), "1")                     # −iπ
H = S(L(K), MOINS_UN)                          # −2πi : un tour complet
B = L(H)                                       # exp(B) = −2πi
A = L(S("1", MOINS_UN))                        # exp(exp(A)) = −e^e


def P(t):                      # −t pour t non nul
    return X(S(A, Q(t)))


def J(t):                      # conj(t) pour t non nul
    return X(P(S(B, t)))


def postfixe(t):
    return [t] if isinstance(t, str) else postfixe(t[1]) + postfixe(t[2]) + ["eml★"]


def evaluer(t, valeurs, etoile=True, _memo=None):
    """Valeur exacte (SymPy), stricte. etoile=False remplace chaque eml★ par eml (contrôle)."""
    memo = {} if _memo is None else _memo
    if isinstance(t, str):
        return sp.Integer(1) if t == "1" else valeurs[t]
    if t not in memo:
        a, b = evaluer(t[1], valeurs, etoile, memo), evaluer(t[2], valeurs, etoile, memo)
        if b == 0:
            raise Indefini("Log(0)")
        v = sp.simplify(sp.exp(a) - sp.log(sp.conjugate(b) if etoile else b))
        if not v.is_finite and v.is_finite is not None:
            raise Indefini("valeur infinie")
        memo[t] = v
    return memo[t]


def main():
    ok = []

    def controle(nom, condition):
        ok.append(bool(condition))
        print(f"[{'OK' if condition else 'ÉCHEC'}] {nom}")

    arbre = J("z")
    mots = postfixe(arbre)
    print("Arbre en notation postfixée :\n" + " ".join(mots) + "\n")
    controle(f"grammaire : seulement 1, z et eml★ (trouvé : {sorted(set(mots))})", set(mots) == {"1", "z", "eml★"})
    controle(f"longueur dépliée = {len(mots)} symboles (attendu 101)", len(mots) == 101)

    for nom, c, attendu in (("O", O, 0), ("−1", MOINS_UN, -1), ("k", K, -sp.I * sp.pi), ("h", H, -2 * sp.I * sp.pi)):
        controle(f"constante {nom} = {attendu}", sp.simplify(evaluer(c, {}) - attendu) == 0)
    controle("constante b : exp(b) = −2πi", sp.simplify(sp.exp(evaluer(B, {})) + 2 * sp.I * sp.pi) == 0)
    controle("constante a : exp(exp(a)) = −e^e", sp.simplify(sp.exp(sp.exp(evaluer(A, {}))) + sp.exp(sp.E)) == 0)

    z = sp.Symbol("z", complex=True, zero=False)
    controle("J(z) = conj(z) sur un symbole complexe non nul", sp.simplify(evaluer(arbre, {"z": z}) - sp.conjugate(z)) == 0)
    points = (sp.Integer(1), sp.Integer(-1), sp.I, sp.exp(sp.E), sp.Rational(-1, 3), sp.Integer(-50),
              sp.Rational(3, 2) + 2 * sp.I, -5 * sp.I, sp.Rational(1, 7) - 40 * sp.I)
    for q in points:
        controle(f"J({q}) = conj({q})", sp.simplify(evaluer(arbre, {"z": q}) - sp.conjugate(q)) == 0)

    try:
        evaluer(arbre, {"z": sp.Integer(0)})
        controle("refus strict en z = 0", False)
    except Indefini:
        controle("refus strict en z = 0", True)

    x, y = sp.Symbol("x", complex=True), sp.Symbol("y", complex=True, zero=False)
    porte = S("x", J("y"))
    controle(f"eml★(x, J(y)) = eml(x, y) sur des symboles ({len(postfixe(porte))} symboles)",
             sp.simplify(evaluer(porte, {"x": x, "y": y}) - (sp.exp(x) - sp.log(y))) == 0)
    for q in (sp.Integer(-2), sp.Integer(-50), sp.Rational(1, 2) + sp.I):
        controle(f"eml★(1/3, J({q})) = eml(1/3, {q})",
                 sp.simplify(evaluer(porte, {"x": sp.Rational(1, 3), "y": q}) - (sp.exp(sp.Rational(1, 3)) - sp.log(q))) == 0)

    try:                                         # contrôle négatif : sans l'étoile, l'arbre n'est PAS un miroir
        sans = evaluer(arbre, {"z": sp.I}, etoile=False)
        controle("contrôle : le même arbre sans étoile ne rend pas conj(i)", sp.simplify(sans + sp.I) != 0)
    except Indefini:
        controle("contrôle : le même arbre sans étoile ne rend pas conj(i) (il est indéfini)", True)

    print(f"\nBILAN : {sum(ok)}/{len(ok)} contrôles passés. Python {sys.version.split()[0]}, SymPy {sp.__version__}.")
    print("Non établi par ce script : absence de zéro intermédiaire pour TOUT z non nul (argument de domaine, § 17) ; minimalité.")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
