import numpy as np, mpmath as mp, time
mp.mp.dps = 40
# ---- 1) la proposition de Quadruplay : eml3(x, y, z) = exp(x) - ln(y) + Re(z) Im(z)
def eml3(x, y, z):
    ex = mp.mpf(0) if x == -mp.inf else mp.exp(x)
    return ex - mp.log(y) + mp.re(z) * mp.im(z)
pire = 0
import random; random.seed(1)
for _ in range(2000):
    z = mp.mpc(random.uniform(-50, 50), random.uniform(-50, 50))
    re = eml3(-mp.inf, 1, z + 1j) - eml3(-mp.inf, 1, z)
    im = eml3(-mp.inf, 1, z + 1) - eml3(-mp.inf, 1, z)
    pire = max(pire, abs(re - z.real), abs(im - z.imag))
x, y = mp.mpc(0.3, 2), mp.mpc(-1.2, 0.7)
print("1) Quadruplay : pire écart sur Re et Im, 2000 points jusqu'à |Im| = 50 :", mp.nstr(pire, 3),
      "| eml3(x, y, 1) − eml(x, y) =", mp.nstr(abs(eml3(x, y, 1) - (mp.exp(x) - mp.log(y))), 3))

# ---- 2) eml★ SEUL (sans eml), feuilles 1 et z : atteint-on conj(z) ? recherche exhaustive par nombre de feuilles
PTS = np.array([0.37 + 0.61j, -0.83 + 1.27j, 1.9 - 0.45j])
def op(a, b):                      # eml★(a, b) = exp(a) − ln(conj b), sur des tableaux (formules × points)
    with np.errstate(all="ignore"):
        return np.exp(a) - np.log(np.conj(b))
def cle(v):
    return np.round(v[:, 0].real, 9) + 1j * np.round(v[:, 0].imag, 9), np.round(v[:, 1].real, 9) + 1j * np.round(v[:, 1].imag, 9)
niveaux = {1: (np.array([np.ones(3, complex), PTS]), ["1", "z"])}
vus = {}
for v, t in zip(*niveaux[1]):
    vus[(complex(np.round(v[0], 9)), complex(np.round(v[1], 9)))] = t
cibles = {"conj(z)": np.conj(PTS), "Re(z)": PTS.real + 0j, "Im(z)": PTS.imag + 0j, "-z": -PTS, "1 - z": 1 - PTS, "i": np.full(3, 1j), "0": np.zeros(3, complex), "ln z": np.log(PTS), "z - 1": PTS - 1, "1/z": 1 / PTS, "|z|^2": (abs(PTS) ** 2) + 0j}
trouve = {}
t0 = time.time()
for n in range(2, 11):
    vals, txts = [], []
    for i in range(1, n):
        A, TA = niveaux[i]; B, TB = niveaux[n - i]
        if len(A) * len(B) > 6e7:
            print(f"   n = {n} : coupe {i}+{n-i} trop grande ({len(A)}×{len(B)}), arrêt de la recherche"); A = None; break
        R = op(A[:, None, :], B[None, :, :]).reshape(-1, 3)
        ok = ~np.any(np.isnan(R.real) | np.isnan(R.imag), axis=1)
        idx = np.nonzero(ok)[0]
        for j in idx:
            v = R[j]
            k = (complex(np.round(v[0], 9)), complex(np.round(v[1], 9)))
            if k in vus:
                continue
            ta, tb = TA[j // len(B)], TB[j % len(B)]
            txt = f"s({ta},{tb})"
            vus[k] = txt
            vals.append(v); txts.append(txt)
            for nom, c in cibles.items():
                if nom not in trouve and (np.all(np.isfinite(v)) and np.max(np.abs(v - c)) < 1e-9):
                    trouve[nom] = (n, txt)
    if A is None:
        break
    niveaux[n] = (np.array(vals), txts)
    print(f"   {n} feuilles (longueur RPN {2*n-1}) : {len(vals)} fonctions nouvelles, {time.time()-t0:.0f} s ; trouvées jusqu'ici : {sorted(trouve)}")
    if time.time() - t0 > 500:
        break
print("2) eml★ seul, feuilles {1, z}, s(a,b) = eml★(a,b) :")
for nom in cibles:
    print(f"   {nom:8s} :", f"TROUVÉ à {trouve[nom][0]} feuilles : {trouve[nom][1]}" if nom in trouve else "pas trouvé dans la recherche")
