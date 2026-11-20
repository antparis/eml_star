"""The conjugation bar facing the wall (2026-11-19): the mirror formula z̄ = 1 − eml★(0, eml(z, 1)) (Theorem 3.1), and what is
built from it (Re, |z|²), with the principal log and then with a branch mutant (each log receives 2πi·k at random)."""
import random
import mpmath as mp

mp.mp.dps = 60
TURN = 2j * mp.pi


def wrap(w):
    """w brought back into the strip −π < Im ≤ π by whole turns (principal-log convention)."""
    k = mp.ceil((w.imag - mp.pi) / (2 * mp.pi))
    return w - TURN * k


def tools(log):
    eml = lambda x, y: mp.exp(x) - log(y)
    emls = lambda x, y: mp.exp(x) - log(mp.conj(y))
    conj = lambda z: 1 - emls(mp.mpc(0), eml(z, mp.mpc(1)))
    return eml, emls, conj


def mutant(rng):
    return lambda z: mp.mpc(0) if z == 1 else mp.log(z) + TURN * rng.randint(-3, 3)


rng = random.Random(1119)
pts = [mp.mpc(rng.uniform(-5, 5), rng.uniform(-15, 15)) for _ in range(3000)]
inside = [z for z in pts if -mp.pi <= z.imag < mp.pi]
eml, emls, conj = tools(mp.log)

print("1) Mirror formula, principal log:")
print("   in the strip −π ≤ Im z < π:", len(inside), "points, worst |formula − z̄| =", mp.nstr(max(abs(conj(z) - mp.conj(z)) for z in inside), 3))
print("   everywhere (Im up to 15): worst |formula − wrap(z̄)| =", mp.nstr(max(abs(conj(z) - wrap(mp.conj(z))) for z in pts), 3))
z = mp.mpc(1, 5)
print("   example z = 1 + 5i: formula =", mp.nstr(conj(z), 10), "; z̄ =", mp.nstr(mp.conj(z), 10), "; wrap(z̄) =", mp.nstr(wrap(mp.conj(z)), 10))
print("   edges (exact value expected: 0.5 + πi in BOTH cases; the raw tree in floating point is unstable on the edge):")
print("      Im z = −π → raw tree", mp.nstr(conj(mp.mpc(0.5, -mp.pi)), 8), "; Im z = +π → raw tree", mp.nstr(conj(mp.mpc(0.5, mp.pi)), 8))

print("2) What NEVER depends on the wall:")
print("   worst |exp(formula) − exp(z̄)| / |exp(z̄)| =", mp.nstr(max(abs(mp.exp(conj(z)) - mp.exp(mp.conj(z))) / abs(mp.exp(mp.conj(z))) for z in pts), 3))
print("   worst |Re(formula) − Re(z)| =", mp.nstr(max(abs(conj(z).real - z.real) for z in pts), 3))

print("3) Branch mutant (the wall moves at every log):")
_, _, conj_m = tools(mutant(rng))
ks = []
for z in pts:
    k = (conj_m(z) - mp.conj(z)) / TURN
    ks.append(k)
print("   (formula − z̄) / 2πi always an integer: worst deviation", mp.nstr(max(abs(k - mp.nint(k.real)) for k in ks), 3),
      "; integers seen", int(min(mp.nint(k.real) for k in ks)), "…", int(max(mp.nint(k.real) for k in ks)))

print("4) What is built from it:")
re_f = lambda z: (z + conj(z)) / 2
mod_f = lambda z: z * conj(z)
print("   Re z = (z + formula)/2: worst deviation on the real part =", mp.nstr(max(abs(re_f(z).real - z.real) for z in pts), 3),
      "; max spurious imaginary part =", mp.nstr(max(abs(re_f(z).imag) for z in pts), 4), "(0 in the strip:", mp.nstr(max(abs(re_f(z).imag) for z in inside), 3), ")")
print("   |z|² = z · formula: worst deviation in the strip =", mp.nstr(max(abs(mod_f(z) - abs(z) ** 2) for z in inside), 3),
      "; outside the strip, example z = 1 + 5i:", mp.nstr(mod_f(mp.mpc(1, 5)), 8), "instead of 26")

print("5) old eml_toolkit/core.py, modulus_squared(z) = eml(z, eml★(z, 1)) = exp(z) − Log(exp(z)):")
for z in (mp.mpc(3, 4), mp.mpc(1, 0), mp.mpc(0.5, -0.2)):
    print("   z =", mp.nstr(z, 4), ": modulus_squared =", mp.nstr(eml(z, emls(z, mp.mpc(1))), 10), "| |z|² =", mp.nstr(abs(z) ** 2, 10), "| exp(z) − Log(exp(z)) =", mp.nstr(mp.exp(z) - mp.log(mp.exp(z)), 10))
