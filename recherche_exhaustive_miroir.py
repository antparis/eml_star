"""Recherche exhaustive du plus court arbre sur {eml, eml★} avec feuilles {1, z} pour conj(z), Re, Im, |z|²…
e(a,b) = eml(a,b) ; s(a,b) = eml★(a,b). Valeurs en 3 points de la bande ; infinis admis (réels étendus de numpy)."""
import numpy as np, time
PTS = np.array([0.37 + 0.61j, -0.83 + 1.27j, 1.9 - 0.45j])
def ops(a, b):
    with np.errstate(all="ignore"):
        ea = np.exp(a)
        return ea - np.log(b), ea - np.log(np.conj(b))
niveaux = {1: (np.array([np.ones(3, complex), PTS]), ["1", "z"])}
vus = set()
for v in niveaux[1][0]:
    vus.add((complex(np.round(v[0], 9)), complex(np.round(v[1], 9))))
cibles = {"conj(z)": np.conj(PTS), "Re(z)": PTS.real + 0j, "Im(z)": PTS.imag + 0j, "|z|^2": abs(PTS) ** 2 + 0j, "-z": -PTS, "1 - z": 1 - PTS,
          "1 - conj(z)": 1 - np.conj(PTS), "-conj(z)": -np.conj(PTS), "i": np.full(3, 1j)}
trouve, t0 = {}, time.time()
for n in range(2, 10):
    vals, txts, stop = [], [], False
    for i in range(1, n):
        A, TA = niveaux[i]; B, TB = niveaux[n - i]
        if len(A) * len(B) > 4e7:
            stop = True; break
        for nom_op, R in zip("es", ops(A[:, None, :], B[None, :, :])):
            R = R.reshape(-1, 3)
            ok = ~np.any(np.isnan(R.real) | np.isnan(R.imag), axis=1)
            k0 = np.round(R[:, 0], 9); k1 = np.round(R[:, 1], 9)
            for j in np.nonzero(ok)[0]:
                k = (complex(k0[j]), complex(k1[j]))
                if k in vus:
                    continue
                vus.add(k)
                txt = f"{nom_op}({TA[j // len(B)]},{TB[j % len(B)]})"
                v = R[j]; vals.append(v); txts.append(txt)
                if np.all(np.isfinite(v)):
                    for nom, c in cibles.items():
                        if nom not in trouve and np.max(np.abs(v - c)) < 1e-9:
                            trouve[nom] = (n, txt)
    if stop:
        print(f"   {n} feuilles : trop grand, arrêt"); break
    niveaux[n] = (np.array(vals), txts)
    print(f"   {n} feuilles (longueur RPN {2*n-1}) : {len(vals)} fonctions nouvelles, {time.time()-t0:.0f} s ; trouvées : {sorted(trouve)}", flush=True)
    if time.time() - t0 > 420:
        break
for nom in cibles:
    print(f"   {nom:12s} :", f"TROUVÉ, {trouve[nom][0]} feuilles, longueur RPN {2*trouve[nom][0]-1} : {trouve[nom][1]}" if nom in trouve else "pas trouvé")
