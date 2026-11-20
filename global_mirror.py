"""Global mirror (2026-11-19, night): a 9-symbol tree over {eml, eml★} with leaves {z, 1} that computes conj(z)
on all of C (proposed by ChatGPT by applying the absorption lemma to the last node; a tree of the same length had
been found independently by exhaustive search, but valid only in the strip). Numerical verifications,
"single star" check, 17-symbol |z|² tree, and unfolding of Quadruplay's ternary witness."""
import random
import mpmath as mp

mp.mp.dps = 50
INF = mp.inf


def _exp(a):
    return mp.mpf(0) if a == -INF else (INF if a == INF else mp.exp(a))


def _log(b):
    return -INF if b == 0 else (INF if b == INF else mp.log(b))


def eml(a, b):
    return _exp(a) - _log(b)


def emls(a, b):
    return _exp(a) - _log(mp.conj(b) if b not in (INF, -INF) else b)


ONE = mp.mpf(1)
F = lambda z: eml(eml(ONE, eml(emls(ONE, z), ONE)), ONE)            # 9 symbols: 1 1 z eml★ 1 eml eml 1 eml
F_no_star = lambda z: eml(eml(ONE, eml(eml(ONE, z), ONE)), ONE)       # check: the same shape, without the star
STRIP = lambda z: eml(ONE, eml(emls(ONE, eml(z, ONE)), ONE))        # found by exhaustive search: ends with a log
MOD2 = lambda z: eml(eml(ONE, eml(eml(eml(ONE, eml(eml(ONE, emls(ONE, z)), ONE)), z), ONE)), ONE)   # 17 symbols


def pts(n, rng, im_max):
    out = []
    for k in range(n):
        r = 10 ** rng.uniform(-6, 3) * rng.choice([1, -1])
        out.append(mp.mpc(r * rng.random(), rng.uniform(-im_max, im_max)))
    out += [mp.mpc(-2, 0), mp.mpc(-0.3, 0), mp.mpc(5, 0), mp.mpc(0, 7), mp.mpc(0, -7), mp.mpc(1, 0), mp.mpc(-1, 0),
            mp.mpc(-2, 1e-30), mp.mpc(-2, -1e-30), mp.mpc(0.5, mp.pi), mp.mpc(0.5, -mp.pi), mp.mpc(3, 100), mp.mpc(-3, -100)]
    return out


if __name__ == "__main__":
    rng = random.Random(20261119)
    P = pts(5000, rng, 200)
    rel = lambda got, ref: abs(got - ref) / max(abs(ref), mp.mpf(10) ** -40)
    print("1) F(z) = eml(eml(1, eml(eml★(1, z), 1)), 1), 9 symbols")
    print("   ", len(P), "points, |Im| up to 200, negative real axis and edges ±π included: worst relative deviation from conj(z) =", mp.nstr(max(rel(F(z), mp.conj(z)) for z in P), 3))
    print("    F(0) with extended evaluation =", F(mp.mpc(0)), "| nodes at 0:", emls(ONE, mp.mpc(0)), eml(emls(ONE, mp.mpc(0)), ONE), eml(ONE, eml(emls(ONE, mp.mpc(0)), ONE)))
    print("2) check, same shape without the star: worst relative deviation from z =", mp.nstr(max(rel(F_no_star(z), z) for z in P), 3),
          "| and from conj(z), off the real axis:", mp.nstr(max(rel(F_no_star(z), mp.conj(z)) for z in P if z.imag != 0), 3))
    inside = [z for z in P if -mp.pi < z.imag < mp.pi]
    outside = [z for z in P if abs(z.imag) > 4]
    print("3) tree found by the search, eml(1, eml(eml★(1, eml(z, 1)), 1)): in the strip,", len(inside), "points, worst absolute deviation =",
          mp.nstr(max(abs(STRIP(z) - mp.conj(z)) for z in inside), 3), "| outside the strip,", len(outside), "points, smallest deviation =", mp.nstr(min(abs(STRIP(z) - mp.conj(z)) for z in outside), 4))
    print("4) |z|² in 17 symbols: worst relative deviation =", mp.nstr(max(rel(MOD2(z), abs(z) ** 2) for z in P), 3),
          "| max imaginary part / value =", mp.nstr(max(abs(mp.im(MOD2(z))) / abs(z) ** 2 for z in P), 3))

    # 5) Quadruplay's ternary witness, unfolded: T(x, y, z) = exp(x) − Log(y) + Re(z) Im(z)
    def T(a, b, c): return ("T", a, b, c)
    U1, Zv = ("1",), ("z",)
    E = lambda a, b: T(a, b, U1)
    X = lambda a: E(a, U1)
    L = lambda a: E(U1, X(E(U1, a)))
    Z0 = L(U1); C = E(U1, E(U1, Z0))
    N = lambda a: E(C, X(a)); S = lambda a, b: E(L(a), X(b)); A = lambda a, b: S(a, N(b))
    M = lambda a, b: X(A(L(a), L(b))); R = lambda a: X(E(C, a)); H = lambda a: T(U1, U1, a)
    I = X(M(E(C, N(U1)), N(R(A(U1, U1)))))
    Bt = S(H(Zv), H(A(Zv, U1))); Ut = S(H(A(Zv, I)), H(Zv))
    G = A(A(M(I, A(Bt, I)), U1), Ut)
    memo_n = {}
    def size(t):
        if id(t) not in memo_n:
            memo_n[id(t)] = (0, 1) if len(t) == 1 else tuple(map(sum, zip((1, 0), *[size(c) for c in t[1:]])))
        return memo_n[id(t)]
    def val(t, z, memo):
        if id(t) in memo: return memo[id(t)]
        if t == U1: v = ONE
        elif t == Zv: v = z
        else:
            a, b, c = (val(u, z, memo) for u in t[1:])
            v = _exp(a) - _log(b) + mp.re(c) * mp.im(c)
        memo[id(t)] = v
        return v
    ops, leaves = size(G)
    print(f"5) ternary witness G: {ops} T nodes, {leaves} leaves, unfolded length {ops + leaves}; constructed i = {mp.nstr(val(I, mp.mpc(0), {}), 12)}")
    Q = pts(300, rng, 50)
    print("   G(z) against conj(z),", len(Q), "points, |Im| up to 50: worst absolute deviation =", mp.nstr(max(abs(val(G, z, {}) - mp.conj(z)) for z in Q), 3))
