"""Minimalité du miroir de 9 symboles (2026-11-20). On énumère TOUS les arbres {eml, eml★} à feuilles {1, z} de longueur
RPN 1, 3, 5, 7 (1 à 4 feuilles) et on les évalue en un seul point z0 non réel, à 60 chiffres, avec l'évaluation étendue.
Argument : si un arbre valait conj(z) pour tout z non nul, il vaudrait conj(z0) en z0. Il suffit donc qu'aucun ne
s'en approche à z0 (marge très supérieure aux erreurs d'arrondi) pour qu'aucun arbre de longueur ≤ 7 ne soit le miroir."""
import itertools
import mpmath as mp

mp.mp.dps = 60
INF = mp.inf
Z0 = mp.mpc("0.37", "0.61")


def ev(t):
    if t == "1":
        return mp.mpf(1)
    if t == "z":
        return Z0
    op, a, b = t
    a, b = ev(a), ev(b)
    if a is None or b is None:
        return None
    if op == "eml★" and b not in (INF, -INF):
        b = mp.conj(b)
    if b == -INF:
        return None                                   # ln(−∞) indéfini
    lb = -INF if b == 0 else (INF if b == INF else mp.log(b))
    ea = mp.mpf(0) if a == -INF else (INF if a == INF else mp.exp(a))
    if ea == INF and lb == INF:
        return None                                   # ∞ − ∞
    return ea - lb


def arbres(n):
    if n == 1:
        yield "1"; yield "z"; return
    for k in range(1, n):
        for g in arbres(k):
            for d in arbres(n - k):
                for op in ("eml", "eml★"):
                    yield (op, g, d)


cible = mp.conj(Z0)
total, indefinis, plus_proche = 0, 0, None
for n in range(1, 5):
    for t in arbres(n):
        total += 1
        v = ev(t)
        if v is None or not mp.isfinite(mp.re(v)) or not mp.isfinite(mp.im(v)):
            indefinis += 1
            continue
        d = abs(v - cible)
        if plus_proche is None or d < plus_proche[0]:
            plus_proche = (d, t)
print(f"Arbres de longueur RPN ≤ 7 énumérés : {total} (dont {indefinis} infinis ou indéfinis en z0 = {Z0})")
print(f"Plus petit écart à conj(z0) : {mp.nstr(plus_proche[0], 6)}, atteint par {plus_proche[1]}")
F = ("eml", ("eml", "1", ("eml", ("eml★", "1", "z"), "1")), "1")
print(f"Contrôle, l'arbre de 9 symboles en z0 : écart {mp.nstr(abs(ev(F) - cible), 3)}")
print("CONCLUSION :", "aucun arbre de longueur ≤ 7 ne réalise conj ; 9 est minimal" if plus_proche[0] > mp.mpf("1e-10") else "À EXAMINER")
