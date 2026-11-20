import numpy as np, mpmath as mp, time
mp.mp.dps = 40
# ---- 1) Quadruplay's proposal: eml3(x, y, z) = exp(x) - ln(y) + Re(z) Im(z)
def eml3(x, y, z):
    ex = mp.mpf(0) if x == -mp.inf else mp.exp(x)
    return ex - mp.log(y) + mp.re(z) * mp.im(z)
worst = 0
import random; random.seed(1)
for _ in range(2000):
    z = mp.mpc(random.uniform(-50, 50), random.uniform(-50, 50))
    re = eml3(-mp.inf, 1, z + 1j) - eml3(-mp.inf, 1, z)
    im = eml3(-mp.inf, 1, z + 1) - eml3(-mp.inf, 1, z)
    worst = max(worst, abs(re - z.real), abs(im - z.imag))
x, y = mp.mpc(0.3, 2), mp.mpc(-1.2, 0.7)
print("1) Quadruplay: worst deviation on Re and Im, 2000 points up to |Im| = 50:", mp.nstr(worst, 3),
      "| eml3(x, y, 1) − eml(x, y) =", mp.nstr(abs(eml3(x, y, 1) - (mp.exp(x) - mp.log(y))), 3))

# ---- 2) eml★ ALONE (without eml), leaves 1 and z: do we reach conj(z)? exhaustive search by number of leaves
PTS = np.array([0.37 + 0.61j, -0.83 + 1.27j, 1.9 - 0.45j])
def op(a, b):                      # eml★(a, b) = exp(a) − ln(conj b), on arrays (formulas × points)
    with np.errstate(all="ignore"):
        return np.exp(a) - np.log(np.conj(b))
def key(v):
    return np.round(v[:, 0].real, 9) + 1j * np.round(v[:, 0].imag, 9), np.round(v[:, 1].real, 9) + 1j * np.round(v[:, 1].imag, 9)
levels = {1: (np.array([np.ones(3, complex), PTS]), ["1", "z"])}
seen = {}
for v, t in zip(*levels[1]):
    seen[(complex(np.round(v[0], 9)), complex(np.round(v[1], 9)))] = t
targets = {"conj(z)": np.conj(PTS), "Re(z)": PTS.real + 0j, "Im(z)": PTS.imag + 0j, "-z": -PTS, "1 - z": 1 - PTS, "i": np.full(3, 1j), "0": np.zeros(3, complex), "ln z": np.log(PTS), "z - 1": PTS - 1, "1/z": 1 / PTS, "|z|^2": (abs(PTS) ** 2) + 0j}
found = {}
t0 = time.time()
for n in range(2, 11):
    vals, txts = [], []
    for i in range(1, n):
        A, TA = levels[i]; B, TB = levels[n - i]
        if len(A) * len(B) > 6e7:
            print(f"   n = {n}: split {i}+{n-i} too large ({len(A)}×{len(B)}), search stopped"); A = None; break
        R = op(A[:, None, :], B[None, :, :]).reshape(-1, 3)
        ok = ~np.any(np.isnan(R.real) | np.isnan(R.imag), axis=1)
        idx = np.nonzero(ok)[0]
        for j in idx:
            v = R[j]
            k = (complex(np.round(v[0], 9)), complex(np.round(v[1], 9)))
            if k in seen:
                continue
            ta, tb = TA[j // len(B)], TB[j % len(B)]
            txt = f"s({ta},{tb})"
            seen[k] = txt
            vals.append(v); txts.append(txt)
            for name, c in targets.items():
                if name not in found and (np.all(np.isfinite(v)) and np.max(np.abs(v - c)) < 1e-9):
                    found[name] = (n, txt)
    if A is None:
        break
    levels[n] = (np.array(vals), txts)
    print(f"   {n} leaves (RPN length {2*n-1}): {len(vals)} new functions, {time.time()-t0:.0f} s; found so far: {sorted(found)}")
    if time.time() - t0 > 500:
        break
print("2) eml★ alone, leaves {1, z}, s(a,b) = eml★(a,b):")
for name in targets:
    print(f"   {name:8s} :", f"FOUND at {found[name][0]} leaves: {found[name][1]}" if name in found else "not found in the search")
