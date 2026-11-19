"""Quels logarithmes de l'arbre × peut-on déplacer sans changer le résultat ? (objection de relecture du 2026-11-19 :
« pour tout choix de branche à chaque logarithme » est faux ; contre-exemple : le log 1 du dernier nœud.)
On numérote les logarithmes de l'arbre développé dans l'ordre d'évaluation, on classe chacun selon que son argument
dépend des variables ou non, et on décale UN SEUL logarithme de 2πi à la fois."""
import mpmath as mp

mp.mp.dps = 60
TOUR = 2j * mp.pi
X, Y, ONE = ("x",), ("y",), ("1",)
eml = lambda a, b: ("eml", a, b)
E_exp = lambda z: eml(z, ONE)
E_ln = lambda z: eml(ONE, eml(eml(ONE, z), ONE))
E_zero = lambda: eml(ONE, eml(eml(ONE, ONE), ONE))
E_neg = lambda z: eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))
E_sub = lambda a, b: eml(E_ln(a), E_exp(b))
E_add = lambda a, b: E_sub(a, E_neg(b))
E_mul = lambda a, b: E_exp(E_add(E_ln(a), E_ln(b)))


def a_variable(t):
    return t in (X, Y) or (len(t) == 3 and (a_variable(t[1]) or a_variable(t[2])))


def evalue(t, x, y, decale=None, tours=None):
    """Évalue l'arbre (log principal, réels étendus). `decale` = numéro du logarithme qui reçoit + 2πi.
    Renvoie (valeur, liste des logs : (numéro, argument dépend des variables ?, est-ce le log de la feuille 1 ?))."""
    logs = []

    def ev(t):
        if t == ONE:
            return mp.mpf(1)
        if t == X:
            return x
        if t == Y:
            return y
        a, b = ev(t[1]), ev(t[2])
        n = len(logs)
        logs.append((n, a_variable(t[2]), t[2] == ONE))
        if b == 0:
            lb = -mp.inf
        elif b == mp.inf:
            lb = mp.inf
        else:
            lb = mp.log(b)
        if decale == n and mp.isfinite(lb):
            lb = lb + TOUR
        ea = mp.mpf(0) if a == -mp.inf else (mp.inf if a == mp.inf else mp.exp(a))
        if abs(mp.re(ea)) < mp.mpf(10) ** -55 and mp.im(ea) == 0 and not a_variable(t[1]) and t[1] != ONE:
            pass
        return ea - lb

    # E_zero vaut 0 exactement : on force la valeur exacte de ce sous-arbre constant (en flottant : ~1e-60)
    zero = E_zero()

    def ev2(t):
        if t == zero:
            n0 = len(logs)
            for k in range(3):                                   # ses trois logs, numérotés pour être décalables
                logs.append((n0 + k, False, k != 2))
            # E_zero = eml(1, eml(eml(1, 1), 1)). Valeur EXACTE selon le logarithme décalé (calcul à la main) :
            if decale == n0:          # ln 1 du nœud le plus profond : e − 2πi, puis exp(e − 2πi) − 0 = exp(e), puis e − ln(exp e) = 0
                return mp.mpf(0)
            if decale == n0 + 1:      # ln 1 du nœud du milieu : exp(e) − 2πi, puis e − ln(exp(e) − 2πi) : fini, non nul
                return mp.e - mp.log(mp.exp(mp.e) - TOUR)
            if decale == n0 + 2:      # ln(exp e) de la racine de E_zero : e − (e + 2πi) = −2πi
                return mp.mpc(0) - TOUR
            return mp.mpf(0)
        if t in (ONE, X, Y):
            return ev(t)
        a, b = ev2(t[1]), ev2(t[2])
        if a is None or b is None:
            return None
        n = len(logs)
        logs.append((n, a_variable(t[2]), t[2] == ONE))
        lb = -mp.inf if b == 0 else (mp.inf if b == mp.inf else mp.log(b))
        if decale == n and mp.isfinite(lb):
            lb = lb + TOUR
        if tours and mp.isfinite(lb):
            lb = lb + TOUR * tours.get(n, 0)               # plusieurs logarithmes décalés EN MÊME TEMPS
        ea = mp.mpf(0) if a == -mp.inf else (mp.inf if a == mp.inf else mp.exp(a))
        positions.append((n, porte.get(id(t), "racine")))
        return ea - lb

    positions, porte = [], {}

    def marque(t):
        if len(t) == 3:
            porte[id(t[1])], porte[id(t[2])] = "x (exponentié par le parent)", "y (son log est repris par le parent)"
            marque(t[1]); marque(t[2])
    marque(t)
    v = ev2(t)
    return (v, logs, dict(positions)) if tours is not None else (v, logs)


def effet_d_un_decalage(n, x=mp.mpc(0.7, 0.3), y=mp.mpc(-1.3, 0.4)):
    """(résultat − x·y) quand le seul logarithme n° n de l'arbre × reçoit + 2πi."""
    v, _ = evalue(E_mul(X, Y), x, y, decale=n)
    return v - x * y


if __name__ == "__main__":
    x, y = mp.mpc(0.7, 0.3), mp.mpc(-1.3, 0.4)
    t = E_mul(X, Y)
    base, logs = evalue(t, x, y)
    print(f"Arbre × : {len(logs)} logarithmes. Sans décalage : écart avec x·y = {mp.nstr(abs(base - x * y), 3)}")
    bilan = {}
    for n, variable, feuille1 in logs:
        v, _ = evalue(t, x, y, decale=n)
        if not mp.isfinite(abs(v)):
            effet = "résultat infini"
        else:
            d = (v - x * y)
            effet = "aucun effet" if abs(d) < mp.mpf(10) ** -40 else (f"résultat décalé de {mp.nstr(d / TOUR, 4)} tour(s)" if abs(d / TOUR - mp.nint((d / TOUR).real)) < 1e-30 else "résultat faux")
        sorte = "argument variable" if variable else ("log de la feuille 1" if feuille1 else "argument constant")
        bilan.setdefault((sorte, effet), []).append(n)
    for (sorte, effet), ns in sorted(bilan.items()):
        print(f"   {sorte:22s} : {effet:32s} — logs n° {ns}")
    v, _ = evalue(t, x, y, decale=len(logs) - 1)
    print("Dernier nœud seul (log 1 de la racine décalé de +2πi) : résultat − x·y =", mp.nstr(v - x * y, 12), " (attendu −2πi =", mp.nstr(-TOUR, 12), ")")

    # ---- 2) où vont les 8 logarithmes libres ? et que se passe-t-il quand on les déplace TOUS ensemble ?
    import random
    _, logs, pos = evalue(t, x, y, tours={})
    libres = [n for n, variable, _ in logs if variable]
    print("Position, dans son parent, du nœud qui porte chaque logarithme à argument variable :")
    for n in libres:
        print(f"   log n° {n:2d} : porte {pos[n]}")
    rng = random.Random(1119)
    pire, n_pts = mp.mpf(0), 0
    for _ in range(3000):
        xx = mp.mpc(rng.uniform(-3, 3), rng.uniform(-6, 6)); yy = mp.mpc(rng.uniform(-3, 3), rng.uniform(-6, 6))
        k = {n: rng.randint(-3, 3) for n in libres}
        v, _, _ = evalue(t, xx, yy, tours=k)
        pire = max(pire, abs(v - xx * yy) / abs(xx * yy)); n_pts += 1
    print(f"Décalages SIMULTANÉS au hasard (−3…3 tours) sur les {len(libres)} logs libres, {n_pts} points : pire écart relatif avec x·y = {mp.nstr(pire, 3)}")
    for kk in (0, 1, -2):
        v, _, _ = evalue(t, mp.mpf(1), mp.mpc(-1.3, 0.4), tours={n: kk for n in libres})
        print(f"   en x = 1 exactement, tous les logs libres décalés de {kk:+d} tour(s) : écart avec x·y = {mp.nstr(abs(v - mp.mpc(-1.3, 0.4)), 3)}")
