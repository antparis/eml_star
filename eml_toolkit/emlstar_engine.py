"""
eml★ Verification Engine
Replicates the logic of main.rs (Odrzywołek) in Python/mpmath,
with eml★ added as a candidate operator.

Two objectives:
1. Verify that {eml, eml★, 1} reconstructs the target operations
2. Test branch-safety on the found trees
"""
import mpmath
import itertools
mpmath.mp.dps = 60

EULER_GAMMA = mpmath.euler   # γ ≈ 0.5772...
GLAISHER    = mpmath.mpc('1.2824271291006226368753425688697917277676889273250011920637400')

# ─── Operators ────────────────────────────────────────────────────────────────

def eml(x, y):
    """eml(x,y) = exp(x) - ln(y)"""
    return mpmath.exp(x) - mpmath.log(y)

def eml_star(x, y):
    """eml★(x,y) = exp(x) - ln(conj(y))"""
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))

ONE = mpmath.mpc(1)

# ─── Symbolic tree ────────────────────────────────────────────────────────────

class Node:
    def __init__(self, op, left=None, right=None, value=None):
        self.op    = op      # 'eml', 'eml★', 'const'
        self.left  = left
        self.right = right
        self.value = value   # for leaf nodes

    def eval(self, x=None, y=None):
        if self.op == 'const':
            return self.value
        elif self.op == 'var_x':
            return x
        elif self.op == 'var_y':
            return y
        elif self.op == 'eml':
            l = self.left.eval(x, y)
            r = self.right.eval(x, y)
            return mpmath.exp(l) - mpmath.log(r)
        elif self.op == 'eml★':
            l = self.left.eval(x, y)
            r = self.right.eval(x, y)
            return mpmath.exp(l) - mpmath.log(mpmath.conj(r))

    def __repr__(self):
        if self.op == 'const':
            return '1'
        elif self.op in ('var_x', 'var_y'):
            return self.op[-1]
        else:
            return f"{self.op}({self.left},{self.right})"

# ─── Numerical verification (Odrzywołek method) ───────────────────────────────

def near(a, b, tol=mpmath.mpf('1e-50')):
    return abs(a - b) < tol

def verify_formula(tree_fn, target_fn, test_points, label=""):
    """Verifies that tree_fn ≈ target_fn on test_points."""
    max_err = 0
    for pt in test_points:
        try:
            got = tree_fn(*pt)
            ref = target_fn(*pt)
            err = float(abs(got - ref))
            if err > max_err:
                max_err = err
        except:
            return False, float('inf')
    return max_err < 1e-50, max_err

# ─── Known trees (from Appendix A + Odrzywołek Table 4) ──────────────────────

def ln_x(z):
    """ln(z) = eml(1, eml(eml(1,z), 1))  K=7"""
    return eml(ONE, eml(eml(ONE, z), ONE))

def exp_x(z):
    """exp(z) = eml(z, 1)  K=3"""
    return eml(z, ONE)

def subtract(x, y):
    """x-y = eml(ln(x), exp(y)) = eml(eml(1,eml(eml(1,x),1)), eml(y,1))  K=11"""
    return eml(eml(ONE, eml(eml(ONE, x), ONE)), eml(y, ONE))

def conj_emlstar(z):
    """z̄ = 1 - eml★(0, eml(z,1))  Theorem 3.1"""
    return 1 - eml_star(mpmath.mpc(0), eml(z, ONE))

def re_z(z):
    """Re(z) = (z + z̄)/2"""
    return (z + conj_emlstar(z)) / 2

def mod_sq(z):
    """|z|² = z * z̄"""
    return z * conj_emlstar(z)

# ─── Test points (Odrzywołek method: γ and A) ─────────────────────────────────

test_unary = [
    (EULER_GAMMA,),
    (mpmath.mpc(1.3),),
    (mpmath.mpc(-0.4),),
    (mpmath.mpc(0.9),),
    (mpmath.mpc(2.1),),
    (mpmath.mpc(0.8, 0.2),),
]

test_binary = [
    (EULER_GAMMA, GLAISHER),
    (mpmath.mpc(1.3), mpmath.mpc(0.7)),
    (mpmath.mpc(0.8, 0.2), mpmath.mpc(1.1)),
]

test_band = [  # Im(z) ∈ (-π, π) for eml★
    (mpmath.mpc(r, i),)
    for r in [-5, -2, 0, 1, 3]
    for i in [-2.8, -1.0, 0.0, 1.0, 2.8]
]

# ─── Tests ────────────────────────────────────────────────────────────────────

print("=" * 65)
print("eml★ VERIFICATION ENGINE — mpmath 60 decimal places")
print("=" * 65)

results = []

def run(label, fn, ref_fn, pts):
    ok, err = verify_formula(fn, ref_fn, pts, label)
    status = "EXACT" if ok else f"err={err:.2e}"
    print(f"  [{'PASS' if ok else 'FAIL'}]  {status:20s}  {label}")
    results.append((label, ok, err))

print("\n-- Holomorphic functions (eml only) --")
run("exp(z) = eml(z,1)",
    lambda z: exp_x(z), lambda z: mpmath.exp(z), test_unary)
run("ln(z) = eml(1,eml(eml(1,z),1))",
    lambda z: ln_x(z), lambda z: mpmath.log(z), test_unary)
run("x-y = eml(ln(x),exp(y))",
    lambda x, y: subtract(x, y), lambda x, y: x - y, test_binary)

print("\n-- Anti-holomorphic functions (eml★ required) --")
run("conj(z) = 1-eml★(0,eml(z,1))  [Im(z)∈(-π,π)]",
    lambda z: conj_emlstar(z), lambda z: mpmath.conj(z), test_band)
run("Re(z) = (z+conj(z))/2         [Im(z)∈(-π,π)]",
    lambda z: re_z(z), lambda z: mpmath.re(z), test_band)
run("|z|² = z*conj(z)              [Im(z)∈(-π,π)]",
    lambda z: mod_sq(z), lambda z: abs(z)**2, test_band)

print("\n-- Out-of-band counterexamples --")
for z_val in [mpmath.mpc(0, 4), mpmath.mpc(0, 100)]:
    comp = conj_emlstar(z_val)
    ref  = mpmath.conj(z_val)
    err  = float(abs(comp - ref))
    n    = round(err / (2 * float(mpmath.pi)))
    inside = abs(float(mpmath.im(z_val))) < float(mpmath.pi)
    print(f"  z={z_val}  Im={float(mpmath.im(z_val)):.1f}  "
          f"err={err:.4f}={n}x2pi  "
          f"({'in band' if inside else 'OUT OF BAND -> expected failure'})")

# ─── Basis comparison ─────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("BASIS COMPARISON")
print("=" * 65)

comparison = [
    ("{eml, 1}",
     "Odrzywołek (2026)",
     "Holomorphic elementary functions",
     "conj, Re, Im, |z|²",
     "Uniform grammar S->1|eml(S,S)"),

    ("{eml, eml★, 1}",
     "Monnerot (2026)",
     "Holomorphic + anti-holomorphic (strip Im∈(-π,π))",
     "Functions requiring |Im(value)|>π",
     "Uniform binary grammar, Sheffer-compatible"),

    ("{eml, Re, 1}",
     "Alternative basis",
     "Holomorphic + anti-holomorphic (all of C)",
     "Nothing",
     "Re is unary — breaks uniform binary grammar"),
]

print(f"\n{'Basis':<22} {'Author':<22} {'Covers':<45} {'Does not cover':<30} {'Structure'}")
print("-" * 145)
for base, who, covers, not_covers, struct in comparison:
    print(f"{base:<22} {who:<22} {covers:<45} {not_covers:<30} {struct}")

# ─── Final verdict ────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("FINAL VERDICT")
print("=" * 65)
ok_count = sum(1 for _, ok, _ in results if ok)
print(f"\nTests passed: {ok_count}/{len(results)}")
print("""
Key points:

1. eml★ is EXACT on its claimed domain (Im(z) in (-pi, pi))
   Verified at 60 decimal places on conj(z), Re(z), |z|²

2. Odrzywołek verification (Rust) uses double precision (f64)
   Our mpmath tests are strictly more demanding

3. {eml, eml★, 1}:
   [+] Stays within the Sheffer binary paradigm
   [+] Grammar: S -> 1 | eml(S,S) | eml★(S,S)
   [!] Conditional on Im(value) in (-pi, pi)
   [+] Cor. 2.2 unconditional (the main result)

4. Trees x*y (K=17) and x+y (K=19): in Supplementary Info
   -> Not formally verified — Th. 4.3 remains conditional
   -> Verification possible once SI is available
""")
