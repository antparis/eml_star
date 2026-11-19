"""Script du § 8 du rapport (axe orthogonal), REMIS EN FICHIER le 2026-11-19 au soir : la trace d'origine
(CONSOLE_20261119_axe_orthogonal.txt) avait été produite par ces mêmes lignes tapées directement dans le terminal,
sans fichier. Trois essais : 1) la partie réelle efface-t-elle le saut 2πi ? 2) l'addition (forme lemme) est-elle
certifiée sur le carré plein, zéro compris ? 3) témoin par points, au log principal ordinaire, collés aux coupures."""
import random
import mpmath as mp
import certify_branch_safety as cb

X, Y, B = cb.X, cb.Y, cb.Boite
# 1) porte « partie réelle » : le saut 2πi est purement imaginaire
t = cb.E_ln(cb.E_ln(X))
a, b = cb.valeur_point(t, mp.mpc(0.7, 1e-12)), cb.valeur_point(t, mp.mpc(0.7, -1e-12))
print("1) ln(ln x) autour de x=0,7 : saut total", mp.nstr(a - b, 8), "| saut de la partie réelle", mp.nstr(a.real - b.real, 3),
      "| saut de u + conj(u)", mp.nstr((a + mp.conj(a)) - (b + mp.conj(b)), 3), "| saut de la partie imaginaire", mp.nstr(a.imag - b.imag, 8))
# 2) addition avec le lemme : plus de trou autour de 0 ?
add_l = cb.moins_ln(X, cb.E_exp(cb.E_neg(Y)))
for p in [(mp.mpc(0), mp.mpc(1.2, 0.4)), (mp.mpc(-0.5, 0.2), mp.mpc(0)), (mp.mpc(0), mp.mpc(0))]:
    print("   E_add forme lemme en", p, "=", mp.nstr(cb.valeur_point(add_l, *p), 8))
r = cb.certifier(add_l, {"x": B(-3, 3, -3, 3), "y": B(-3, 3, -3, 3)}, profondeur=9, budget=300)
print("2) x + y (forme lemme) sur le carré PLEIN [-3,3]+i[-3,3], zéro compris :", "CERTIFIÉ" if r["certifie"] else "ÉCHEC", r["boites"], "rectangles", r["raison"][:150])
# 3) témoin indépendant de la logique d'intervalles
random.seed(20261119)
pire, n = 0, 0
for k in range(6000):
    if k % 3 == 0:
        x = mp.mpc(random.uniform(0.1, 1.0), random.choice([1, -1]) * 10 ** random.uniform(-30, -3))      # collé au segment (0,1)
    elif k % 3 == 1:
        x = mp.mpc(random.uniform(-3, -0.1), random.choice([1, -1]) * 10 ** random.uniform(-30, -3))    # collé à l'axe négatif
    else:
        x = mp.mpc(random.uniform(-3, 3), random.uniform(0.1, 3) * random.choice([1, -1]))
    y = mp.mpc(random.uniform(0.1, 3) * random.choice([1, -1]), random.uniform(-1, 1))
    for t, ref in ((cb.E_mul(X, Y), x * y), (cb.E_mul_lemme(X, Y), x * y), (cb.E_add(X, Y), x + y)):
        pire = max(pire, abs(cb.valeur_point(t, x, y) - ref))
        n += 1
print("3) témoin :", n, "valeurs comparées à x·y et x+y, pire écart =", mp.nstr(pire, 3))
