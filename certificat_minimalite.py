"""Certificat à arrondi garanti de la minimalité du miroir (2026-11-20, demandé en relecture par ChatGPT).

Point témoin RATIONNEL exact z0 = (37 + 61 i)/100. Pour chacun des 714 arbres {eml, eml★} à feuilles {1, z} de longueur
RPN ≤ 7, on encadre T(z0) − conj(z0) en arithmétique d'intervalles (mpmath.iv, arrondi vers l'extérieur) et l'arbre
est ÉLIMINÉ dès que l'intervalle de la partie réelle ou celui de la partie imaginaire exclut 0. Un arbre égal à conj
partout le serait en z0 : 714 éliminations = aucun arbre de longueur ≤ 7 ne réalise le miroir.

Convention (celle du manuscrit) : log principal, Arg ∈ (−π, π] ; évaluation étendue nœud par nœud : ln 0 = −∞,
exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞ ; ∞ − ∞ et ln(−∞) indéfinis.
- Sous-arbres CONSTANTS (sans z) : décidés exactement avec SymPy (égalité à 0, infini, réel négatif). Un réel négatif
  EXACT a pour log principal ln|b| + iπ : encadré ainsi, sans dépendre d'une coupure flottante.
- Sous-arbres dépendant de z : si l'argument d'un log est un rectangle qui contient 0 ou qui chevauche l'axe réel
  négatif, l'arbre est déclaré INDÉCIS (rien n'est deviné).
- Un arbre indéfini ou infini en z0 ne peut pas valoir conj(z0) (fini) : il est éliminé, avec ce motif."""
import sympy as sp
from mpmath import iv

iv.dps = 60
Z0_EXACT = sp.Rational(37, 100) + sp.I * sp.Rational(61, 100)
Z0 = iv.mpc(iv.mpf("0.37"), iv.mpf("0.61"))
CONJ = iv.mpc(iv.mpf("0.37"), -iv.mpf("0.61"))
PI_I = iv.mpc(0, iv.pi)


class Indecis(Exception):
    pass


class Indefini(Exception):
    pass


def exact_const(op, a, b):
    """Valeur exacte (SymPy, avec sp.oo / -sp.oo) d'un nœud dont les deux enfants sont constants."""
    if op == "eml★" and b not in (sp.oo, -sp.oo):
        b = sp.conjugate(b)
    if b == -sp.oo:
        raise Indefini("ln(−∞)")
    lb = -sp.oo if sp.simplify(b) == 0 else (sp.oo if b == sp.oo else sp.log(b))
    ea = 0 if a == -sp.oo else (sp.oo if a == sp.oo else sp.exp(a))
    if ea == sp.oo and lb == sp.oo:
        raise Indefini("∞ − ∞")
    return sp.oo if ea == sp.oo else (-sp.oo if lb == sp.oo else (sp.oo if lb == -sp.oo else sp.simplify(ea - lb)))


def _est_reel_negatif(c):
    """Décision exacte « c est un réel strictement négatif ». Si SymPy ne sait pas trancher, l'arbre est INDÉCIS
    (on ne devine pas, et on ne plante pas)."""
    im_nul = sp.im(c) == 0
    if im_nul is not True and im_nul is not False:
        raise Indecis(f"partie imaginaire d'une constante non décidée : {c}")
    if not im_nul:
        return False
    negatif = sp.re(c).is_negative
    if negatif is None:
        raise Indecis(f"signe d'une constante non décidé : {c}")
    return negatif


def ev(t):
    """Renvoie (encadrement, valeur exacte ou None). L'encadrement est un iv.mpc, ou '+inf' / '-inf'."""
    if t == "1":
        return iv.mpc(1, 0), sp.Integer(1)
    if t == "z":
        return Z0, None
    op, ga, gb = t
    A, sa = ev(ga)
    B, sb = ev(gb)
    if sa is not None and sb is not None:                      # nœud constant : exact d'abord
        s = exact_const(op, sa, sb)
        if s in (sp.oo, -sp.oo):
            return ("+inf" if s == sp.oo else "-inf"), s
    else:
        s = None
    # logarithme du second argument
    if isinstance(B, str):
        if B == "-inf":
            raise Indefini("ln(−∞)")
        LB = "+inf"
    else:
        Bc = iv.mpc(B.real, -B.imag) if op == "eml★" else B
        if sb is not None:                                      # argument constant : décision exacte
            bc = sp.conjugate(sb) if op == "eml★" else sb
            if sp.simplify(bc) == 0:
                LB = "-inf"
            elif _est_reel_negatif(bc):                           # réel négatif EXACT : Log = ln|b| + iπ
                LB = iv.log(iv.mpc(-Bc.real, 0)) + PI_I
            else:
                if Bc.real.a <= 0 and Bc.imag.a <= 0 <= Bc.imag.b:
                    raise Indecis(f"constante non décidée près de la coupure : {bc}")
                LB = iv.log(Bc)
        else:
            contient0 = Bc.real.a <= 0 <= Bc.real.b and Bc.imag.a <= 0 <= Bc.imag.b
            if contient0:
                raise Indecis("log d'un rectangle qui contient 0")
            if not (Bc.real.a > 0 or Bc.imag.a > 0 or Bc.imag.b < 0):
                raise Indecis("log d'un rectangle à cheval sur l'axe réel négatif")
            LB = iv.log(Bc)
    # exponentielle du premier argument
    if isinstance(A, str):
        EA = iv.mpc(0, 0) if A == "-inf" else "+inf"
    else:
        EA = iv.exp(A)
    if EA == "+inf":
        if LB == "+inf":
            raise Indefini("∞ − ∞")
        return "+inf", s
    if LB == "+inf":
        return "-inf", s
    if LB == "-inf":
        return "+inf", s
    return EA - LB, s


def arbres(n):
    if n == 1:
        yield "1"; yield "z"; return
    for k in range(1, n):
        for g in arbres(k):
            for d in arbres(n - k):
                for op in ("eml", "eml★"):
                    yield (op, g, d)


def separe_de_zero(d):
    """Décision d'élimination, UTILISÉE PAR juger ET PAR LE CONTRÔLE 3 : renvoie "réelle" ou "imaginaire" si la partie
    correspondante de l'intervalle d exclut strictement 0 (la partie réelle est examinée d'abord), sinon None."""
    if d.real.a > 0 or d.real.b < 0:
        return "réelle"
    if d.imag.a > 0 or d.imag.b < 0:
        return "imaginaire"
    return None


MOTIF_CONTIENT_ZERO = "INDÉCIS (les deux intervalles contiennent 0)"


def juger(t):
    try:
        v, _ = ev(t)
    except Indefini as e:
        return "éliminé (indéfini en z0 : " + str(e) + ")"
    except Indecis as e:
        return "INDÉCIS (" + str(e) + ")"
    if isinstance(v, str):
        return "éliminé (valeur infinie en z0)"
    partie = separe_de_zero(v - CONJ)
    if partie is not None:
        return f"éliminé (partie {partie} de T(z0) − conj(z0) séparée de 0)"
    return MOTIF_CONTIENT_ZERO


def controles():
    """Contrôles du code, demandés en relecture : couverture du générateur, encadrement du point exact, règle d'élimination."""
    import sympy as sp
    tous = [t for n in range(1, 5) for t in arbres(n)]
    distincts = len(set(tous))
    from math import comb
    catalan = lambda m: comb(2 * m, m) // (m + 1)
    attendu = sum(catalan(n - 1) * 2 ** n * 2 ** (n - 1) for n in range(1, 5))     # formes × feuilles × opérateurs
    print(f"Contrôle 1, générateur : {len(tous)} arbres produits, {distincts} distincts, formule de dénombrement {attendu}")
    assert len(tous) == distincts == attendu == 714
    def exact_rationnel(raw):
        """Valeur binaire EXACTE d'une borne d'intervalle, lue directement dans sa représentation interne
        (signe, mantisse, exposant), en rationnel SymPy. (Ne pas passer par mpmath.mpf : il arrondit à la précision
        courante, 53 bits par défaut, et ferait croire à tort que l'intervalle est réduit à un point.)"""
        signe, mantisse, exposant, _ = raw
        return (-1) ** signe * sp.Integer(mantisse) * sp.Integer(2) ** exposant
    for v, exact in ((Z0.real, sp.Rational(37, 100)), (Z0.imag, sp.Rational(61, 100))):
        bas, haut = exact_rationnel(v._mpi_[0]), exact_rationnel(v._mpi_[1])
        assert bas < exact < haut, "z0 hors de son intervalle"
    print("Contrôle 2, encadrement : l'intervalle de z0 contient strictement 37/100 et 61/100 (comparaison exacte en rationnels)")
    zero_reel = iv.mpc(iv.mpf([-1e-70, 1e-70]), iv.mpf([1, 2]))       # partie réelle contient 0, partie imaginaire non
    zero_imag = iv.mpc(iv.mpf([1, 2]), iv.mpf([-1e-70, 1e-70]))       # l'inverse
    zero_deux = iv.mpc(iv.mpf([-1e-70, 1e-70]), iv.mpf([-1e-70, 1e-70]))
    bord = iv.mpc(iv.mpf([0, 1]), iv.mpf([0, 1]))                      # 0 sur le bord des deux intervalles : non séparé
    assert separe_de_zero(zero_reel) == "imaginaire" and separe_de_zero(zero_imag) == "réelle"
    assert separe_de_zero(zero_deux) is None and separe_de_zero(bord) is None
    print("Contrôle 3, règle (la fonction separe_de_zero, celle qu'utilise juger) : un écart dont les deux intervalles contiennent 0, même au bord, n'est PAS éliminé ; un seul intervalle strictement séparé de 0 suffit")
    return tous


if __name__ == "__main__":
    import collections, sys
    controles()
    bilan, par_longueur, total = collections.Counter(), collections.Counter(), 0
    indecis, par_imag = [], []
    for n in range(1, 5):
        for t in arbres(n):
            total += 1
            par_longueur[2 * n - 1] += 1
            r = juger(t)
            bilan[r] += 1
            if "partie imaginaire" in r:
                par_imag.append(t)
            if r.startswith("INDÉCIS"):
                indecis.append((t, r))
    print(f"z0 = (37 + 61 i)/100 exact ; arbres énumérés : {total} ; par longueur RPN : {dict(sorted(par_longueur.items()))}")
    for k, v in sorted(bilan.items()):
        print(f"   {v:4d} × {k}")
    print(f"INDÉCIS : {len(indecis)}")
    print(f"Arbre(s) éliminé(s) par la partie imaginaire : {par_imag} (prédit en relecture : la feuille z seule, car z0 − conj(z0) = 61i/50)")
    for t, r in indecis[:10]:
        print("   ", t, r)
    F = ("eml", ("eml", "1", ("eml", ("eml★", "1", "z"), "1")), "1")
    verdict_F = juger(F)
    print("Contrôle positif, l'arbre de 9 symboles F :", verdict_F, "(exigé : « les deux intervalles contiennent 0 », puisqu'il vaut exactement conj(z0))")
    n_reel = sum(v for k, v in bilan.items() if "partie réelle" in k)
    attendus = {
        "714 arbres énumérés": total == 714,
        "0 indécis": not indecis,
        "713 éliminés par la partie réelle": n_reel == 713,
        "seule la feuille z éliminée par la partie imaginaire": par_imag == ["z"],
        "contrôle positif F non éliminé, motif exact": verdict_F == MOTIF_CONTIENT_ZERO,
    }
    for nom, vrai in attendus.items():
        print(f"   [{'OK' if vrai else 'ÉCHEC'}] {nom}")
    ok = all(attendus.values())
    print("CONCLUSION :", "714 arbres éliminés, 0 indécis : aucun arbre de longueur RPN ≤ 7 ne réalise le miroir ; 9 est minimal." if ok else "NON CONCLUANT")
    sys.exit(0 if ok else 1)
