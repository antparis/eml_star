"""Certificat de branch-safety par arithmétique d'intervalles (mpmath.iv) pour les arbres eml / eml★.

Ce que prouve un certificat réussi, sur tout le domaine donné (produit de rectangles complexes pour x et y, privé
éventuellement d'un disque autour d'un point de x) : en chaque point, l'arbre se calcule et la fonction calculée est
CONTINUE. Pour cela, chaque logarithme de l'arbre est contrôlé :

- son argument évite la coupure (−∞, 0] : rien à dire ;
- son argument traverse l'axe réel négatif (loin de 0) ET son nœud entre dans le parent par la PORTE X (le premier
  argument, celui qu'on exponentie) : le saut de 2πi est effacé par l'exponentielle du parent, car
  exp(c − log y) = e^c / y quelle que soit la branche du log. On calcule alors avec une branche tournée,
  ln(−y) + iπ, dont la coupure évite le rectangle (écart 0 ou 2πi avec le log principal) ;
- tout autre cas (traversée en porte y, traversée à la racine, 0 dans le rectangle) : échec, on subdivise ; échec
  définitif si la profondeur maximale est atteinte.

Les sous-arbres CONSTANTS (aucune variable) sont calculés à part, exactement, sur les réels étendus [−∞, +∞] :
ln 0 = −∞, exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞ ; ∞ − ∞ et ln(−∞) sont refusés. L'égalité à 0 est décidée par
SymPy (ex. E_zero = e − ln(e^e) = 0). Les constantes finies sont ensuite encadrées par un intervalle de largeur 1e−55
autour de leur valeur à 70 chiffres (evalf de SymPy, erreur contrôlée).

Arbres : copies exactes de branch_safety_final.py (E_exp, E_ln, E_zero, E_neg, E_sub, E_add, E_mul), sous forme d'arbres.
Origine du travail : docs de reprise CutLattice du 2026-11-19 (découverte de la traversée de ln(ln x) dans E_mul)."""
import time

import mpmath as mp
import sympy as sp
from mpmath import iv

iv.dps = 50
mp.mp.dps = 50

PLUS_INF, MOINS_INF = "+inf", "-inf"


class Indefini(ValueError):
    pass


class Echec(Exception):
    """Un logarithme n'a pas pu être contrôlé sur ce rectangle."""


# ---------- arbres ----------

ONE, X, Y = ("1",), ("x",), ("y",)


def eml(a, b):
    return ("eml", a, b)


def emls(a, b):
    return ("eml*", a, b)


def E_exp(z):    return eml(z, ONE)
def E_ln(z):     return eml(ONE, eml(eml(ONE, z), ONE))
def E_zero():    return eml(ONE, eml(eml(ONE, ONE), ONE))
def E_neg(z):    return eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))
def E_sub(a, b): return eml(E_ln(a), E_exp(b))
def E_add(a, b): return E_sub(a, E_neg(b))
def E_mul(a, b): return E_exp(E_add(E_ln(a), E_ln(b)))


def moins_ln(a, c):
    """Nœud « a − ln c ». Sert UNIQUEMENT à remplacer eml(E_ln(a), c), en vertu du lemme (prouvé à la main, voir
    FINDINGS_20261119_branch_safety_certificate.md) : exp(E_ln(w)) = w pour tout w ≠ 0, quelle que soit la branche, et
    = 0 = w en w = 0 avec les conventions des réels étendus."""
    return ("-ln", a, c)


def E_mul_lemme(a, b):
    """E_mul(a, b) où le seul passage par ln(ln a) est remplacé grâce au lemme : même valeur, sans ce passage."""
    return E_exp(moins_ln(E_ln(a), E_exp(E_neg(E_ln(b)))))


def variables(t):
    if t[0] in ("x", "y"):
        return {t[0]}
    if t[0] == "1":
        return set()
    return variables(t[1]) | variables(t[2])


def taille(t):
    """Nombre de symboles en notation polonaise inverse (feuilles + opérateurs)."""
    return 1 if len(t) == 1 else 1 + taille(t[1]) + taille(t[2])


def texte(t):
    if len(t) == 1:
        return t[0]
    return f"{t[0]}({texte(t[1])}, {texte(t[2])})"


# ---------- constantes exactes sur les réels étendus ----------

_cache = {}


def constante(t):
    """Valeur exacte d'un sous-arbre sans variable : expression SymPy finie, PLUS_INF ou MOINS_INF."""
    if t in _cache:
        return _cache[t]
    if t == ONE:
        v = sp.Integer(1)
    else:
        a, b = constante(t[1]), constante(t[2])
        if b == MOINS_INF:
            raise Indefini(f"ln(−∞) dans {texte(t)}")
        if isinstance(b, str):                      # b = +∞ : ln b = +∞
            if a == PLUS_INF:
                raise Indefini(f"∞ − ∞ dans {texte(t)}")
            v = MOINS_INF
        elif b == 0:                                # ln 0 = −∞
            v = PLUS_INF
        elif a == PLUS_INF:
            v = PLUS_INF
        else:
            ea = sp.Integer(0) if a == MOINS_INF else sp.exp(a)
            arg = sp.conjugate(b) if t[0] == "eml*" else b
            v = sp.simplify(ea - sp.log(arg))
    _cache[t] = v
    return v


CONSTANTES_FINIES = [0]      # compteur : combien de constantes finies non nulles ont été encadrées (aucune pour + et ×)


def _intervalle_constant(t):
    """Encadrement GARANTI d'un sous-arbre constant. ±∞ et 0 viennent du calcul exact (SymPy). Une constante finie non
    nulle est recalculée nœud par nœud en arithmétique d'intervalles (arrondi dirigé) ; si l'un de ses logarithmes
    n'est pas strictement hors de la coupure, on refuse (Echec) plutôt que de deviner."""
    v = constante(t)
    if isinstance(v, str):
        return v
    if v == 0:
        return iv.mpc(0, 0)
    if t == ONE:
        return iv.mpc(1, 0)
    CONSTANTES_FINIES[0] += 1
    a, b = _intervalle_constant(t[1]), _intervalle_constant(t[2])
    if isinstance(a, str) or isinstance(b, str):
        raise Echec(f"constante finie passant par l'infini, non traitée : {texte(t)[:80]}")
    arg = iv.mpc(b.real, -b.imag) if t[0] == "eml*" else b
    if _contient_zero(arg) or not _hors_coupure(arg):
        raise Echec(f"constante dont un logarithme touche la coupure, non traitée : {texte(t)[:80]}")
    return iv.exp(a) - iv.log(arg)


# ---------- rectangles ----------

def _borne(v, vers_le_bas):
    """Une borne décimale donnée en TEXTE ("0.1") est arrondie vers l'extérieur du rectangle, pour que le rectangle
    calculé contienne le rectangle annoncé. Un nombre Python (float, int) est pris tel quel : sa valeur binaire exacte."""
    if isinstance(v, str):
        return mp.mpf(iv.mpf(v)._mpi_[0 if vers_le_bas else 1])
    return mp.mpf(v)


class Boite:
    def __init__(self, re_min, re_max, im_min, im_max):
        self.r = (_borne(re_min, True), _borne(re_max, False))
        self.i = (_borne(im_min, True), _borne(im_max, False))

    def iv(self):
        return iv.mpc(iv.mpf([self.r[0], self.r[1]]), iv.mpf([self.i[0], self.i[1]]))

    def point(self, re, im):
        p = Boite(re, re, im, im)
        return p

    def couper(self):
        mr, mi = (self.r[0] + self.r[1]) / 2, (self.i[0] + self.i[1]) / 2
        return [Boite(self.r[0], mr, self.i[0], mi), Boite(mr, self.r[1], self.i[0], mi),
                Boite(self.r[0], mr, mi, self.i[1]), Boite(mr, self.r[1], mi, self.i[1])]

    def dans_disque(self, centre, rayon):
        """Le rectangle entier est-il, de façon GARANTIE, dans le disque ouvert ? (un disque est convexe : il suffit des
        quatre coins.) Comparaison des carrés en arithmétique d'intervalles : majorant de la distance² < minorant de
        rayon². Un rayon décimal se donne en texte ("0.05")."""
        c = mp.mpc(centre)
        r = iv.mpf(rayon) if isinstance(rayon, str) else iv.mpf(mp.mpf(rayon))
        r2 = r * r
        for a in self.r:
            for b in self.i:
                d2 = (iv.mpf(a) - iv.mpf(c.real)) ** 2 + (iv.mpf(b) - iv.mpf(c.imag)) ** 2
                if not d2.b < r2.a:
                    return False
        return True

    def __repr__(self):
        f = lambda v: mp.nstr(v, 6)
        return f"[{f(self.r[0])}, {f(self.r[1])}] + i[{f(self.i[0])}, {f(self.i[1])}]"


# ---------- évaluation par intervalles ----------

def _contient_zero(z):
    return z.real.a <= 0 <= z.real.b and z.imag.a <= 0 <= z.imag.b


def _hors_coupure(z):
    return z.real.a > 0 or z.imag.a > 0 or z.imag.b < 0


def log_tourne_point(y):
    """Branche tournée du log, coupure sur l'axe réel positif : ln(−y) + iπ."""
    return mp.log(-y) + 1j * mp.pi


def _ev(t, env, porte, compte):
    if t == ONE:
        return iv.mpc(1, 0)
    if t in (X, Y):
        return env[t[0]]
    if not variables(t):
        return _intervalle_constant(t)
    if t[0] == "-ln":                              # a − ln c : un saut de 2πik dans a devient un saut de la valeur du nœud,
        a = _ev(t[1], env, porte, compte)         # permis exactement quand le nœud lui-même est en porte x
        b = _ev(t[2], env, "y", compte)
        if isinstance(a, str) or isinstance(b, str):
            raise Echec(f"valeur infinie dans {texte(t)[:80]}")
        if _contient_zero(b):
            raise Echec(f"log d'un rectangle qui contient zéro (nœud {texte(t)[:80]})")
        if _hors_coupure(b):
            return a - iv.log(b)
        if porte == "x":
            compte["absorbees"] += 1
            return a - (iv.log(iv.mpc(-b.real, -b.imag)) + iv.mpc(0, iv.pi))
        raise Echec(f"traversée de la coupure dans a − ln c (nœud {texte(t)[:80]})")
    a = _ev(t[1], env, "x", compte)
    b = _ev(t[2], env, "y", compte)
    if isinstance(b, str):
        if b == MOINS_INF:
            raise Echec(f"ln(−∞) dans {texte(t)[:80]}")
        if a == PLUS_INF:
            raise Echec(f"∞ − ∞ dans {texte(t)[:80]}")
        return MOINS_INF
    if a == PLUS_INF:
        return PLUS_INF
    ea = iv.mpc(0, 0) if a == MOINS_INF else iv.exp(a)
    arg = iv.mpc(b.real, -b.imag) if t[0] == "eml*" else b
    if _contient_zero(arg):
        raise Echec(f"log d'un rectangle qui contient zéro (nœud {texte(t)[:80]})")
    if _hors_coupure(arg):
        return ea - iv.log(arg)
    if porte == "x":                                 # traversée absorbée par l'exp du parent
        compte["absorbees"] += 1
        return ea - (iv.log(iv.mpc(-arg.real, -arg.imag)) + iv.mpc(0, iv.pi))
    where = {"racine": "à la racine", "direct": "dans une valeur utilisée telle quelle"}.get(porte, "en porte y")
    raise Echec(f"traversée de la coupure {where} (nœud {texte(t)[:80]})")


def enclos(t, boites):
    """Encadrement de la valeur de l'arbre sur les rectangles donnés (lève Echec si un log n'est pas contrôlé)."""
    r = _ev(t, {k: v.iv() for k, v in boites.items()}, "racine", {"absorbees": 0})
    if isinstance(r, str):
        raise Echec(f"valeur infinie à la racine ({r})")
    return r


def valeur_point(t, x, y=None):
    """Valeur ordinaire (mpmath, log principal) en un point : pour vérifier que l'arbre calcule bien sa fonction."""
    def ev(t):
        if t == ONE:
            return mp.mpc(1)
        if t == X:
            return mp.mpc(x)
        if t == Y:
            return mp.mpc(y)
        if not variables(t):
            v = constante(t)
            return {PLUS_INF: mp.inf, MOINS_INF: -mp.inf}.get(v, v) if isinstance(v, str) else mp.mpc(complex(sp.N(v, 60)))
        a, b = ev(t[1]), ev(t[2])
        if t[0] == "-ln":
            return a - mp.log(b)
        ea = mp.mpf(0) if a == -mp.inf else mp.exp(a)
        if b == mp.inf:
            return -mp.inf
        if b == 0:
            return mp.inf
        return ea - mp.log(mp.conj(b) if t[0] == "eml*" else b)
    return ev(t)


# ---------- branch-and-bound ----------

def certifier(t, domaine, profondeur=10, budget=600, disque_x=None):
    """Découpe le domaine jusqu'à ce que chaque rectangle passe. disque_x = (centre, rayon) : rectangles de x entièrement
    dans ce disque ouvert laissés de côté (comptés dans « exclues ») ; ils relèvent d'une preuve à part.
    Renvoie {certifie, boites, absorbees, exclues, profondeur_max, secondes, raison}."""
    noms = sorted(domaine)
    if disque_x and domaine["x"].dans_disque(*disque_x):
        raise ValueError("tout le domaine de x est dans le disque exclu : il n'y a rien à certifier")
    compte = {"absorbees": 0}
    travail, ok, exclues, pmax, t0 = [(tuple(domaine[n] for n in noms), 0)], 0, 0, 0, time.time()
    while travail:
        boites, d = travail.pop()
        env = dict(zip(noms, boites))
        if disque_x and env["x"].dans_disque(*disque_x):
            exclues += 1
            continue
        try:
            r = _ev(t, {k: v.iv() for k, v in env.items()}, "racine", compte)
            if isinstance(r, str):
                raise Echec(f"valeur infinie à la racine ({r})")
            ok += 1
            pmax = max(pmax, d)
        except Echec as e:
            if d >= profondeur or time.time() - t0 > budget:
                return {"certifie": False, "boites": ok, "absorbees": compte["absorbees"], "exclues": exclues, "profondeur_max": d,
                        "secondes": round(time.time() - t0, 1), "raison": f"{e} ; rectangle {env}"}
            morceaux = [b.couper() for b in boites]
            combos = [()]
            for m in morceaux:
                combos = [c + (p,) for c in combos for p in m]
            travail += [(c, d + 1) for c in combos]
    if ok == 0:                                      # rien n'a été vérifié : ce n'est pas un succès
        return {"certifie": False, "boites": 0, "absorbees": 0, "exclues": exclues, "profondeur_max": 0,
                "secondes": round(time.time() - t0, 1), "raison": "aucun rectangle accepté (tout est dans le disque exclu ?)"}
    return {"certifie": True, "boites": ok, "absorbees": compte["absorbees"], "exclues": exclues, "profondeur_max": pmax,
            "secondes": round(time.time() - t0, 1), "raison": ""}


def E_add_lemme(a, b):
    """E_add(a, b) = eml(E_ln(a), c) récrit a − ln c par le lemme d'absorption (a = 0 compris)."""
    return moins_ln(a, E_exp(E_neg(b)))


_TOUS = []


def _dire(nom, r):
    _TOUS.append(bool(r["certifie"]))
    print(f"{nom} : {'CERTIFIÉ' if r['certifie'] else 'ÉCHEC'} — {r['boites']} rectangles, {r['absorbees']} traversées absorbées, "
          f"{r['exclues']} dans le disque, profondeur {r['profondeur_max']}, {r['secondes']} s" + (f" | {r['raison']}" if r["raison"] else ""))


if __name__ == "__main__":
    # Les bornes décimales sont données en texte : arrondies vers l'extérieur, les rectangles calculés CONTIENNENT
    # les rectangles annoncés ; le disque exclu est CONTENU dans le disque annoncé.
    print(f"Tailles : E_add {taille(E_add(X, Y))}, E_mul {taille(E_mul(X, Y))} symboles")
    D = Boite("0.1", 3, -1, 1)
    _dire("1) x + y, x,y ∈ [0,1;3]+i[−1,1]", certifier(E_add(X, Y), {"x": D, "y": D}))
    _dire("2) x × y (arbre d'origine), même domaine, hors |x−1| < 0,05", certifier(E_mul(X, Y), {"x": D, "y": D}, profondeur=10, budget=3000, disque_x=(1, "0.05")))
    _dire("3) x × y (forme lemme), même domaine, disque compris", certifier(E_mul_lemme(X, Y), {"x": D, "y": D}))
    R = [Boite(-3, "-0.1", -3, 3), Boite("0.1", 3, -3, 3), Boite(-3, 3, "0.1", 3), Boite(-3, 3, -3, "-0.1")]
    for nom, t in (("4) x + y", E_add(X, Y)), ("5) x × y (forme lemme)", E_mul_lemme(X, Y))):
        tot = ab = pm = 0
        ok = True
        for bx in R:
            for by in R:
                r = certifier(t, {"x": bx, "y": by}, profondeur=9)
                tot, ab, pm, ok = tot + r["boites"], ab + r["absorbees"], max(pm, r["profondeur_max"]), ok and r["certifie"]
        _TOUS.append(ok)
        print(f"{nom} sur [−3,3]+i[−3,3] privé de |Re|,|Im| < 0,1 : {'CERTIFIÉ' if ok else 'ÉCHEC'} — {tot} rectangles, {ab} traversées absorbées, profondeur {pm}")
    P = Boite(-3, 3, -3, 3)
    _dire("6) x + y (forme lemme) sur le carré plein [−3,3]+i[−3,3], zéro compris", certifier(E_add_lemme(X, Y), {"x": P, "y": P}, profondeur=9))
    print(f"Constantes finies non nulles encadrées pendant ces certificats : {CONSTANTES_FINIES[0]}")
    import sys
    if not all(_TOUS) or len(_TOUS) != 6:
        print("AU MOINS UN CERTIFICAT ATTENDU A ÉCHOUÉ")
        sys.exit(1)
    print("Les six certificats attendus sont passés.")
