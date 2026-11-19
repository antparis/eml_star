"""Deux questions (2026-11-19) : le résultat dépend-il de l'endroit où l'on place le mur du logarithme ? et que se
passe-t-il loin du domaine certifié (grandes parties imaginaires, voisinage de 0) ?
Mutant de branche : chaque logarithme dont l'argument DÉPEND DES VARIABLES reçoit un nombre entier de tours tiré au
hasard (2πi·k). Les logarithmes de constantes (ln 1 = 0, la tête −∞ de E_neg) ne sont PAS perturbés : les perturber change
le résultat (voir branches_autorisees.py : ln 1 = 2πi à la racine donne exactement −2πi)."""
import random
import mpmath as mp

mp.mp.dps = 80
DEUX_PI_I = 2j * mp.pi


def arbres(log):
    eml = lambda x, y: mp.exp(x) - log(y)
    ONE = mp.mpc(1)
    E_exp = lambda z: eml(z, ONE)
    E_ln = lambda z: eml(ONE, eml(eml(ONE, z), ONE))
    E_neg = lambda z: -log(mp.exp(z))                 # eml(−∞, eml(z, 1)) : la constante vaut −∞, exp(−∞) = 0
    E_sub = lambda a, b: eml(E_ln(a), E_exp(b))
    E_add = lambda a, b: E_sub(a, E_neg(b))
    E_mul = lambda a, b: E_exp(E_add(E_ln(a), E_ln(b)))
    return E_add, E_mul


def log_mutant(rng):
    def log(z):
        if z == 1:
            return mp.mpc(0)                          # ln 1 = 0 : la feuille constante n'est pas mutée
        return mp.log(z) + DEUX_PI_I * rng.randint(-3, 3)
    return log


def tire(rng, im_max):
    r = 10 ** rng.uniform(-3, 2)
    return mp.mpc(r * rng.choice([1, -1]) * rng.random(), rng.uniform(-im_max, im_max))


def enroule(y):
    """y ramené dans la bande |Im| < π par des tours entiers."""
    return y - DEUX_PI_I * mp.nint(y.imag / (2 * mp.pi))


if __name__ == "__main__":
    rng = random.Random(20261119)
    # 1) mur déplacé au hasard à chaque logarithme
    add_m, mul_m = arbres(log_mutant(rng))
    pire_mul, pire_add_tours, tours_vus, n = mp.mpf(0), mp.mpf(0), set(), 0
    for _ in range(4000):
        x, y = tire(rng, 12), tire(rng, 12)
        pire_mul = max(pire_mul, abs(mul_m(x, y) - x * y) / abs(x * y))
        k = (add_m(x, y) - (x + y)) / DEUX_PI_I
        pire_add_tours = max(pire_add_tours, abs(k - mp.nint(k.real)))
        tours_vus.add(int(mp.nint(k.real)))
        n += 1
    print(f"1) mutant de branche, {n} points, |Im| jusqu'à 12 :")
    print(f"   × : pire écart relatif avec x·y = {mp.nstr(pire_mul, 3)}")
    print(f"   + : (résultat − (x+y)) / 2πi est toujours un entier (écart {mp.nstr(pire_add_tours, 3)}) ; entiers vus : {min(tours_vus)} … {max(tours_vus)}")
    # 2) logarithme principal, loin du domaine certifié
    add_p, mul_p = arbres(mp.log)
    pm, pa, hors = mp.mpf(0), mp.mpf(0), 0
    for _ in range(4000):
        x, y = tire(rng, 40), tire(rng, 40)
        pm = max(pm, abs(mul_p(x, y) - x * y) / abs(x * y))
        pa = max(pa, abs(add_p(x, y) - (x + enroule(y))))
        hors += abs(y.imag) > mp.pi
    print(f"2) log principal, 4000 points, |Im| jusqu'à 40 ({hors} avec |Im y| > π) :")
    print(f"   × : pire écart relatif avec x·y = {mp.nstr(pm, 3)}")
    print(f"   + : pire écart avec x + enroule(y) = {mp.nstr(pa, 3)}")
    print("   exemple : E_add(1, 5i) =", mp.nstr(add_p(mp.mpc(1), mp.mpc(0, 5)), 12), " ; 1 + (5 − 2π)i =", mp.nstr(1 + (5 - 2 * mp.pi) * 1j, 12))
    # 3) vers zéro
    print("3) vers zéro (log principal) :")
    for e in (5, 20, 60):
        x = mp.mpc(10) ** -e * mp.mpc(0.6, 0.8); y = mp.mpc(1.3, 0.4)
        print(f"   |x| = 1e-{e} : |E_mul(x, y)| = {mp.nstr(abs(mul_p(x, y)), 4)}, écart relatif avec x·y = {mp.nstr(abs(mul_p(x, y) - x * y) / abs(x * y), 3)}"
              f" ; |E_mul(y, x)| = {mp.nstr(abs(mul_p(y, x)), 4)}")
