"""Exhaustive search for the shortest tree over {eml, eml★} with leaves {1, z} for conj(z), Re, Im, |z|²…
e(a,b) = eml(a,b); s(a,b) = eml★(a,b). Values at 3 points of the strip; infinities allowed (numpy's extended reals)."""
import numpy as np, time
PTS = np.array([0.37 + 0.61j, -0.83 + 1.27j, 1.9 - 0.45j])
def ops(a, b):
    with np.errstate(all="ignore"):
        ea = np.exp(a)
        return ea - np.log(b), ea - np.log(np.conj(b))
levels = {1: (np.array([np.ones(3, complex), PTS]), ["1", "z"])}
seen = set()
for v in levels[1][0]:
    seen.add((complex(np.round(v[0], 9)), complex(np.round(v[1], 9))))
targets = {"conj(z)": np.conj(PTS), "Re(z)": PTS.real + 0j, "Im(z)": PTS.imag + 0j, "|z|^2": abs(PTS) ** 2 + 0j, "-z": -PTS, "1 - z": 1 - PTS,
           "1 - conj(z)": 1 - np.conj(PTS), "-conj(z)": -np.conj(PTS), "i": np.full(3, 1j)}
found, t0 = {}, time.time()
for n in range(2, 10):
    vals, txts, stop = [], [], False
    for i in range(1, n):
        A, TA = levels[i]; B, TB = levels[n - i]
        if len(A) * len(B) > 4e7:
            stop = True; break
        for op_name, R in zip("es", ops(A[:, None, :], B[None, :, :])):
            R = R.reshape(-1, 3)
            ok = ~np.any(np.isnan(R.real) | np.isnan(R.imag), axis=1)
            k0 = np.round(R[:, 0], 9); k1 = np.round(R[:, 1], 9)
            for j in np.nonzero(ok)[0]:
                k = (complex(k0[j]), complex(k1[j]))
                if k in seen:
                    continue
                seen.add(k)
                txt = f"{op_name}({TA[j // len(B)]},{TB[j % len(B)]})"
                v = R[j]; vals.append(v); txts.append(txt)
                if np.all(np.isfinite(v)):
                    for name, c in targets.items():
                        if name not in found and np.max(np.abs(v - c)) < 1e-9:
                            found[name] = (n, txt)
    if stop:
        print(f"   {n} leaves: too large, stopping"); break
    levels[n] = (np.array(vals), txts)
    print(f"   {n} leaves (RPN length {2*n-1}): {len(vals)} new functions, {time.time()-t0:.0f} s; found: {sorted(found)}", flush=True)
    if time.time() - t0 > 420:
        break
for name in targets:
    print(f"   {name:12s} :", f"FOUND, {found[name][0]} leaves, RPN length {2*found[name][0]-1}: {found[name][1]}" if name in found else "not found")
