"""Two questions (2026-11-19): does the result depend on where the wall of the logarithm is placed? and what
happens far from the certified domain (large imaginary parts, neighbourhood of 0)?
Branch mutant: every logarithm whose argument DEPENDS ON THE VARIABLES receives a whole number of turns drawn at
random (2πi·k). The logarithms of constants (ln 1 = 0, the −∞ head of E_neg) are NOT perturbed: perturbing them changes
the result (see allowed_branches.py: ln 1 = 2πi at the root gives exactly −2πi)."""
import random
import mpmath as mp

mp.mp.dps = 80
TWO_PI_I = 2j * mp.pi


def trees(log):
    eml = lambda x, y: mp.exp(x) - log(y)
    ONE = mp.mpc(1)
    E_exp = lambda z: eml(z, ONE)
    E_ln = lambda z: eml(ONE, eml(eml(ONE, z), ONE))
    E_neg = lambda z: -log(mp.exp(z))                 # eml(−∞, eml(z, 1)): the constant equals −∞, exp(−∞) = 0
    E_sub = lambda a, b: eml(E_ln(a), E_exp(b))
    E_add = lambda a, b: E_sub(a, E_neg(b))
    E_mul = lambda a, b: E_exp(E_add(E_ln(a), E_ln(b)))
    return E_add, E_mul


def mutant_log(rng):
    def log(z):
        if z == 1:
            return mp.mpc(0)                          # ln 1 = 0: the constant leaf is not mutated
        return mp.log(z) + TWO_PI_I * rng.randint(-3, 3)
    return log


def draw(rng, im_max):
    r = 10 ** rng.uniform(-3, 2)
    return mp.mpc(r * rng.choice([1, -1]) * rng.random(), rng.uniform(-im_max, im_max))


def wrap(y):
    """y brought back into the strip |Im| < π by whole turns."""
    return y - TWO_PI_I * mp.nint(y.imag / (2 * mp.pi))


if __name__ == "__main__":
    rng = random.Random(20261119)
    # 1) wall moved at random at every logarithm
    add_m, mul_m = trees(mutant_log(rng))
    worst_mul, worst_add_turns, turns_seen, n = mp.mpf(0), mp.mpf(0), set(), 0
    for _ in range(4000):
        x, y = draw(rng, 12), draw(rng, 12)
        worst_mul = max(worst_mul, abs(mul_m(x, y) - x * y) / abs(x * y))
        k = (add_m(x, y) - (x + y)) / TWO_PI_I
        worst_add_turns = max(worst_add_turns, abs(k - mp.nint(k.real)))
        turns_seen.add(int(mp.nint(k.real)))
        n += 1
    print(f"1) branch mutant, {n} points, |Im| up to 12:")
    print(f"   ×: worst relative difference from x·y = {mp.nstr(worst_mul, 3)}")
    print(f"   +: (result − (x+y)) / 2πi is always an integer (deviation {mp.nstr(worst_add_turns, 3)}); integers seen: {min(turns_seen)} … {max(turns_seen)}")
    # 2) principal logarithm, far from the certified domain
    add_p, mul_p = trees(mp.log)
    wm, wa, outside = mp.mpf(0), mp.mpf(0), 0
    for _ in range(4000):
        x, y = draw(rng, 40), draw(rng, 40)
        wm = max(wm, abs(mul_p(x, y) - x * y) / abs(x * y))
        wa = max(wa, abs(add_p(x, y) - (x + wrap(y))))
        outside += abs(y.imag) > mp.pi
    print(f"2) principal log, 4000 points, |Im| up to 40 ({outside} with |Im y| > π):")
    print(f"   ×: worst relative difference from x·y = {mp.nstr(wm, 3)}")
    print(f"   +: worst difference from x + wrap(y) = {mp.nstr(wa, 3)}")
    print(f"   example: E_add(1, 5i) = {mp.nstr(add_p(mp.mpc(1), mp.mpc(0, 5)), 12)}; 1 + (5 − 2π)i = {mp.nstr(1 + (5 - 2 * mp.pi) * 1j, 12)}")
    # 3) towards zero
    print("3) towards zero (principal log):")
    for e in (5, 20, 60):
        x = mp.mpc(10) ** -e * mp.mpc(0.6, 0.8); y = mp.mpc(1.3, 0.4)
        print(f"   |x| = 1e-{e}: |E_mul(x, y)| = {mp.nstr(abs(mul_p(x, y)), 4)}, relative difference from x·y = {mp.nstr(abs(mul_p(x, y) - x * y) / abs(x * y), 3)}"
              f"; |E_mul(y, x)| = {mp.nstr(abs(mul_p(y, x)), 4)}")
