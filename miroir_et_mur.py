"""Ta barre face au mur (2026-11-19) : la formule du miroir z̄ = 1 − eml★(0, eml(z, 1)) (Théorème 3.1), et ce qu'on en
construit (Re, |z|²), avec le log principal puis avec un mutant de branche (chaque log reçoit 2πi·k au hasard)."""
import random
import mpmath as mp

mp.mp.dps = 60
TOUR = 2j * mp.pi


def enroule(w):
    """w ramené dans la bande −π < Im ≤ π par des tours entiers (convention du log principal)."""
    k = mp.ceil((w.imag - mp.pi) / (2 * mp.pi))
    return w - TOUR * k


def outils(log):
    eml = lambda x, y: mp.exp(x) - log(y)
    emls = lambda x, y: mp.exp(x) - log(mp.conj(y))
    conj = lambda z: 1 - emls(mp.mpc(0), eml(z, mp.mpc(1)))
    return eml, emls, conj


def mutant(rng):
    return lambda z: mp.mpc(0) if z == 1 else mp.log(z) + TOUR * rng.randint(-3, 3)


rng = random.Random(1119)
pts = [mp.mpc(rng.uniform(-5, 5), rng.uniform(-15, 15)) for _ in range(3000)]
dans = [z for z in pts if -mp.pi <= z.imag < mp.pi]
eml, emls, conj = outils(mp.log)

print("1) Formule du miroir, log principal :")
print("   dans la bande −π ≤ Im z < π :", len(dans), "points, pire |formule − z̄| =", mp.nstr(max(abs(conj(z) - mp.conj(z)) for z in dans), 3))
print("   partout (Im jusqu'à 15) : pire |formule − enroule(z̄)| =", mp.nstr(max(abs(conj(z) - enroule(mp.conj(z))) for z in pts), 3))
z = mp.mpc(1, 5)
print("   exemple z = 1 + 5i : formule =", mp.nstr(conj(z), 10), "; z̄ =", mp.nstr(mp.conj(z), 10), "; enroule(z̄) =", mp.nstr(enroule(mp.conj(z)), 10))
print("   bords (valeur exacte attendue : 0.5 + πi dans les DEUX cas ; l'arbre brut en flottant est instable sur le bord) :")
print("      Im z = −π → arbre brut", mp.nstr(conj(mp.mpc(0.5, -mp.pi)), 8), "; Im z = +π → arbre brut", mp.nstr(conj(mp.mpc(0.5, mp.pi)), 8))

print("2) Ce qui ne dépend JAMAIS du mur :")
print("   pire |exp(formule) − exp(z̄)| / |exp(z̄)| =", mp.nstr(max(abs(mp.exp(conj(z)) - mp.exp(mp.conj(z))) / abs(mp.exp(mp.conj(z))) for z in pts), 3))
print("   pire |Re(formule) − Re(z)| =", mp.nstr(max(abs(conj(z).real - z.real) for z in pts), 3))

print("3) Mutant de branche (le mur bouge à chaque log) :")
_, _, conj_m = outils(mutant(rng))
ks = []
for z in pts:
    k = (conj_m(z) - mp.conj(z)) / TOUR
    ks.append(k)
print("   (formule − z̄) / 2πi toujours entier : pire écart", mp.nstr(max(abs(k - mp.nint(k.real)) for k in ks), 3),
      "; entiers vus", int(min(mp.nint(k.real) for k in ks)), "…", int(max(mp.nint(k.real) for k in ks)))

print("4) Ce qu'on construit avec :")
re_f = lambda z: (z + conj(z)) / 2
mod_f = lambda z: z * conj(z)
print("   Re z = (z + formule)/2 : pire écart sur la partie réelle =", mp.nstr(max(abs(re_f(z).real - z.real) for z in pts), 3),
      "; partie imaginaire parasite max =", mp.nstr(max(abs(re_f(z).imag) for z in pts), 4), "(0 dans la bande :", mp.nstr(max(abs(re_f(z).imag) for z in dans), 3), ")")
print("   |z|² = z · formule : pire écart dans la bande =", mp.nstr(max(abs(mod_f(z) - abs(z) ** 2) for z in dans), 3),
      "; hors bande, exemple z = 1 + 5i :", mp.nstr(mod_f(mp.mpc(1, 5)), 8), "au lieu de 26")

print("5) ancien eml_toolkit/core.py, modulus_squared(z) = eml(z, eml★(z, 1)) = exp(z) − Log(exp(z)) :")
for z in (mp.mpc(3, 4), mp.mpc(1, 0), mp.mpc(0.5, -0.2)):
    print("   z =", mp.nstr(z, 4), ": modulus_squared =", mp.nstr(eml(z, emls(z, mp.mpc(1))), 10), "| |z|² =", mp.nstr(abs(z) ** 2, 10), "| exp(z) − Log(exp(z)) =", mp.nstr(mp.exp(z) - mp.log(mp.exp(z)), 10))
