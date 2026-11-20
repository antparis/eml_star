"""Minimality of the 9-symbol mirror (2026-11-20). We enumerate ALL trees {eml, eml★} with leaves {1, z} of RPN
length 1, 3, 5, 7 (1 to 4 leaves) and evaluate them at a single non-real point z0, to 60 digits, with extended evaluation.
Argument: if a tree were equal to conj(z) for every nonzero z, it would be equal to conj(z0) at z0. It therefore suffices
that none comes close to it at z0 (margin far above the rounding errors) for no tree of length ≤ 7 to be the mirror."""
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
        return None                                   # ln(−∞) undefined
    lb = -INF if b == 0 else (INF if b == INF else mp.log(b))
    ea = mp.mpf(0) if a == -INF else (INF if a == INF else mp.exp(a))
    if ea == INF and lb == INF:
        return None                                   # ∞ − ∞
    return ea - lb


def trees(n):
    if n == 1:
        yield "1"; yield "z"; return
    for k in range(1, n):
        for left in trees(k):
            for right in trees(n - k):
                for op in ("eml", "eml★"):
                    yield (op, left, right)


target = mp.conj(Z0)
total, undefined_count, closest = 0, 0, None
for n in range(1, 5):
    for t in trees(n):
        total += 1
        v = ev(t)
        if v is None or not mp.isfinite(mp.re(v)) or not mp.isfinite(mp.im(v)):
            undefined_count += 1
            continue
        d = abs(v - target)
        if closest is None or d < closest[0]:
            closest = (d, t)
print(f"Trees of RPN length ≤ 7 enumerated: {total} (of which {undefined_count} infinite or undefined at z0 = {Z0})")
print(f"Smallest deviation from conj(z0): {mp.nstr(closest[0], 6)}, attained by {closest[1]}")
F = ("eml", ("eml", "1", ("eml", ("eml★", "1", "z"), "1")), "1")
print(f"Check, the 9-symbol tree at z0: deviation {mp.nstr(abs(ev(F) - target), 3)}")
print("CONCLUSION:", "no tree of length ≤ 7 realizes conj; 9 is minimal" if closest[0] > mp.mpf("1e-10") else "TO BE EXAMINED")
