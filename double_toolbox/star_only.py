"""eml★ ALONE is enough: independent reconstruction of the macros proposed by ChatGPT on 2026-11-20
(PROOFS.md, archive ~/EML_replay_20261120_Vw6P7N). Nothing is copied from its script: the trees are rebuilt here
with the eml★ node of this toolbox, and judged by our own evaluators.

S(a, b) = eml★(a, b) = exp(a) − Log(conj b). Leaves: 1 and the variable. No eml node.
The idea that unlocks everything: build h = −2πi, then T(t) = h − Log(1/t). T(t) is "a logarithm of t"
(exp(T(t)) = t) shifted by a full turn, hence never zero: its imaginary part stays in [−3π, −π[.
"""
from .eml import ONE, eml_star as S


def X(t):
    return S(t, ONE)                     # exp(t)


E_ = X(ONE)                              # e
C = X(E_)                                # e^e
O = S(ONE, C)                            # 0


def Q(t):
    return X(S(ONE, t))                  # C / conj(t)


def L(t):
    return S(ONE, Q(t))                  # a logarithm of t: exp(L(t)) = t


def d(t, r):
    return S(L(t), X(r))                 # t − r for r real


N1 = d(S(ONE, E_), E_)                   # −1
K = d(S(O, N1), ONE)                     # −iπ
H = S(L(K), N1)                          # −2πi
B_ = L(H)                                # a logarithm of h


def I(t):
    return X(S(B_, t))                   # 1 / conj(t)


def T(t):
    return S(B_, I(t))                   # h − Log(1/t): never zero


def U(t):
    return L(T(t))


def B(u, v):
    return X(S(U(u), v))                 # u / conj(v)


def M(u, v):
    return B(u, I(v))                    # u × v


def J(t):
    """The mirror in eml★ alone: conj(t) for every nonzero t, all the nodes staying finite."""
    return X(M(N1, S(B_, t)))


def eml_from_star(x, y):
    """The eml gate rebuilt with eml★ alone: eml(x, y) = eml★(x, J(y)) for y nonzero."""
    return S(x, J(y))


TWO = d(ONE, N1)
FOUR = d(d(TWO, N1), N1)
MINUS_I = X(M(B(ONE, FOUR), H))
I_ = M(N1, MINUS_I)


def real_part_strict(z):
    """Re z in eml★ alone, all the nodes finite on the WHOLE plane (even if Re z = 0, even at z = 0)."""
    w = X(z)
    a = M(B(ONE, TWO), T(M(w, J(w))))    # Re z − iπ
    h1 = S(L(a), I_)                     # + iπ/2
    return S(L(h1), I_)                  # + iπ/2


def nand_bit(x, y):
    """The NAND gate on bits 0/1, in eml★ alone: t = 4 − x − y, then t(7 − t)/2 − 5."""
    def R(a, r):
        return X(S(L(L(a)), r))          # a / r for a > 1, r > 0
    integers = {1: ONE, 2: TWO}
    for k in range(3, 8):
        integers[k] = d(integers[k - 1], N1)
    t = d(d(integers[4], x), y)
    return d(R(t, R(integers[2], d(integers[7], t))), integers[5])


# ── Shortened mirror (ChatGPT, same day): a simple change of sign replaces the multiplier ─────────────────────
A_ = L(S(ONE, N1))                       # exp(A_) = e − iπ, hence exp(exp(A_)) = −e^e


def P(t):
    return X(S(A_, Q(t)))                # −t for every nonzero t


def J_short(t):
    """The mirror in eml★ alone, short version: conj(t) for every nonzero t, all the nodes finite."""
    return X(P(S(B_, t)))
