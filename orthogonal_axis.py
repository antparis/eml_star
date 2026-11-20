"""Script of § 8 of the report (orthogonal axis), PUT BACK INTO A FILE on the evening of 2026-11-19: the original trace
(CONSOLE_20261119_orthogonal_axis.txt) had been produced by these same lines typed directly into the terminal,
without a file. Three tests: 1) does the real part erase the 2πi jump? 2) is the addition (lemma form)
certified on the full square, zero included? 3) witness by points, with the ordinary principal log, hugging the branch cuts."""
import random
import mpmath as mp
import certify_branch_safety as cb

X, Y, B = cb.X, cb.Y, cb.Box
# 1) "real part" gate: the 2πi jump is purely imaginary
t = cb.E_ln(cb.E_ln(X))
a, b = cb.point_value(t, mp.mpc(0.7, 1e-12)), cb.point_value(t, mp.mpc(0.7, -1e-12))
print("1) ln(ln x) around x=0.7: total jump", mp.nstr(a - b, 8), "| jump of the real part", mp.nstr(a.real - b.real, 3),
      "| jump of u + conj(u)", mp.nstr((a + mp.conj(a)) - (b + mp.conj(b)), 3), "| jump of the imaginary part", mp.nstr(a.imag - b.imag, 8))
# 2) addition with the lemma: no more hole around 0?
add_l = cb.minus_ln(X, cb.E_exp(cb.E_neg(Y)))
for p in [(mp.mpc(0), mp.mpc(1.2, 0.4)), (mp.mpc(-0.5, 0.2), mp.mpc(0)), (mp.mpc(0), mp.mpc(0))]:
    print("   E_add lemma form at", p, "=", mp.nstr(cb.point_value(add_l, *p), 8))
r = cb.certify(add_l, {"x": B(-3, 3, -3, 3), "y": B(-3, 3, -3, 3)}, depth=9, budget=300)
print("2) x + y (lemma form) on the FULL square [-3,3]+i[-3,3], zero included:", "CERTIFIED" if r["certified"] else "FAILED", r["boxes"], "boxes", r["reason"][:150])
# 3) witness independent of the interval logic
random.seed(20261119)
worst, n = 0, 0
for k in range(6000):
    if k % 3 == 0:
        x = mp.mpc(random.uniform(0.1, 1.0), random.choice([1, -1]) * 10 ** random.uniform(-30, -3))      # hugging the segment (0,1)
    elif k % 3 == 1:
        x = mp.mpc(random.uniform(-3, -0.1), random.choice([1, -1]) * 10 ** random.uniform(-30, -3))    # hugging the negative axis
    else:
        x = mp.mpc(random.uniform(-3, 3), random.uniform(0.1, 3) * random.choice([1, -1]))
    y = mp.mpc(random.uniform(0.1, 3) * random.choice([1, -1]), random.uniform(-1, 1))
    for t, ref in ((cb.E_mul(X, Y), x * y), (cb.E_mul_lemma(X, Y), x * y), (cb.E_add(X, Y), x + y)):
        worst = max(worst, abs(cb.point_value(t, x, y) - ref))
        n += 1
print("3) witness:", n, "values compared with x·y and x+y, worst difference =", mp.nstr(worst, 3))
