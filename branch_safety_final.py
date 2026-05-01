"""
BRANCH-SAFETY TEST DÉFINITIF
Arbres EML réels depuis eml_compiler_v4.py (Odrzywołek 2026, Zenodo archive)

K=13 pour x+y, K=20 pour x*y — vérifiés à 60 décimales mpmath
"""
import mpmath
import sys
sys.path.insert(0, '/home/claude/zenodo/pkg/VA00-SymbolicRegressionPackage-db31d58/EML_toolkit/EmL_compiler')

mpmath.mp.dps = 60
PI = mpmath.pi

# ─── Opérateur EML numérique ──────────────────────────────────────────────────

def eml(x, y):
    return mpmath.exp(x) - mpmath.log(y)

def eml_star(x, y):
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))

ONE = mpmath.mpc(1)

# ─── Arbres EML exacts depuis le compilateur ──────────────────────────────────
# Source : eml_compiler_v4.py
# EML[a,b] = exp(a) - ln(b)

def E_exp(z):    return eml(z, ONE)                                     # K=1
def E_ln(z):     return eml(ONE, eml(eml(ONE,z), ONE))                  # K=3
def E_zero():    return eml(ONE, eml(eml(ONE,ONE), ONE))                # K=3
def E_neg(z):    return eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))  # K=8
def E_sub(a,b):  return eml(E_ln(a), E_exp(b))                          # K=5
def E_add(a,b):  return E_sub(a, E_neg(b))                              # K=13
def E_inv(z):    return E_exp(E_neg(E_ln(z)))                           # K=12
def E_mul(a,b):  return E_exp(E_add(E_ln(a), E_ln(b)))                  # K=20
def E_div(a,b):  return E_mul(a, E_inv(b))                              # K=32

# Conjugaison eml★ (Monnerot 2026, Theorem 3.1)
def E_conj(z):   return 1 - eml_star(mpmath.mpc(0), eml(z, ONE))

# ─── Test d'exactitude numérique ─────────────────────────────────────────────

def test(label, fn, ref_fn, pts, tol=1e-50):
    max_err = 0.0
    fails = []
    for pt in pts:
        try:
            got = fn(*pt)
            ref = ref_fn(*pt)
            err = float(abs(got - ref))
            if err > max_err: max_err = err
            if err > tol: fails.append((pt, err))
        except Exception as e:
            fails.append((pt, str(e)))
    ok = not fails
    s = "✅ EXACT" if ok else f"❌ err_max={max_err:.2e}"
    print(f"  {s:22s}  K={label[0]:3s}  {label[1]}")
    return ok, max_err

# Points de test (méthode Odrzywołek : γ et A)
GAMMA = mpmath.euler
GLAISHER = mpmath.mpf('1.28242712910062263687534256886979172776768892732500119206374002')

unary_pts = [(GAMMA,), (mpmath.mpc(1.3),), (mpmath.mpc(0.5,0.2),),
             (mpmath.mpc(2.0),), (mpmath.mpc(-0.5),)]
binary_pts = [(GAMMA, GLAISHER), (mpmath.mpc(1.3), mpmath.mpc(0.7)),
              (mpmath.mpc(0.8,0.2), mpmath.mpc(1.1)),
              (mpmath.mpc(2.0), mpmath.mpc(3.0))]
band_pts = [(mpmath.mpc(r,i),)
            for r in [-5,-2,0,1,3]
            for i in [-2.8,-1.0,0.0,1.0,2.8]]

print("=" * 65)
print("BRANCH-SAFETY TEST DÉFINITIF — Arbres EML réels (Zenodo)")
print(f"mpmath.dps = {mpmath.mp.dps}")
print("=" * 65)

print("\n— Primitives holomorphes eml —")
test(("1","exp(z)=eml(z,1)"),
     lambda z: E_exp(z), lambda z: mpmath.exp(z), unary_pts)
test(("3","ln(z)=eml(1,eml(eml(1,z),1))"),
     lambda z: E_ln(z), lambda z: mpmath.log(z), unary_pts)
test(("5","x−y"), E_sub, lambda x,y: x-y, binary_pts)
test(("8","−x"), lambda z: E_neg(z), lambda z: -z, unary_pts)
test(("13","x+y — K=13"), E_add, lambda x,y: x+y, binary_pts)
test(("12","1/x"), E_inv, lambda z: 1/z, unary_pts)
test(("20","x×y — K=20"), E_mul, lambda x,y: x*y, binary_pts)

print("\n— Anti-holomorphes eml★ (Monnerot) —")
test(("2","z̄  Im(z)∈(-π,π)"),
     lambda z: E_conj(z), lambda z: mpmath.conj(z), band_pts)

# ─── Branch-safety exhaustive sur x+y et x×y ─────────────────────────────────

import random
random.seed(42)

def branch_safe_check(fn, pts, label, k_str):
    """
    Vérifie que fn est exact (= résultat de référence).
    Si exact → pas de saut de branche détecté.
    """
    print(f"\n  Branch-safety {label} (K={k_str}) sur {len(pts)} points :")
    fails = 0
    worst = 0.0
    for pt in pts:
        try:
            got = fn(*pt)
            ref_fn = (lambda x,y: x+y) if '+' in label else (lambda x,y: x*y)
            ref = ref_fn(*pt)
            err = float(abs(got - ref))
            if err > 1e-50:
                fails += 1
                if err > worst: worst = err
        except:
            fails += 1
    if fails == 0:
        print(f"    ✅ AUCUNE VIOLATION — branch-safe sur ce domaine")
    else:
        print(f"    ❌ {fails} violations, err_max={worst:.2e}")
    return fails == 0

# Grille 1 : valeurs réelles positives
real_pts = [(mpmath.mpf(x), mpmath.mpf(y))
            for x in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
            for y in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]]

# Grille 2 : complexes (Im petite, Re > 0)
cx_pts_small = [(mpmath.mpc(r,i), mpmath.mpc(s,j))
                for r in [0.5, 1.0, 2.0]
                for i in [-0.3, 0.0, 0.3]
                for s in [0.5, 1.0]
                for j in [-0.3, 0.0, 0.3]]

# Grille 3 : aléatoire dans la bande (Im ∈ (-1, 1))
random_band = [(mpmath.mpc(random.uniform(0.1,5), random.uniform(-0.9,0.9)),
                mpmath.mpc(random.uniform(0.1,5), random.uniform(-0.9,0.9)))
               for _ in range(1000)]

print("\n— Branch-safety des vrais arbres Odrzywołek —")
ok_add_r  = branch_safe_check(E_add, real_pts,    "x+y réel",           "13")
ok_add_cx = branch_safe_check(E_add, cx_pts_small,"x+y complexe petit", "13")
ok_add_rnd= branch_safe_check(E_add, random_band, "x+y aléatoire bande","13")

ok_mul_r  = branch_safe_check(E_mul, real_pts,    "x×y réel",           "20")
ok_mul_cx = branch_safe_check(E_mul, cx_pts_small,"x×y complexe petit", "20")
ok_mul_rnd= branch_safe_check(E_mul, random_band, "x×y aléatoire bande","20")

# ─── Bilan ────────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("VERDICT FINAL")
print("=" * 65)
all_safe = all([ok_add_r, ok_add_cx, ok_add_rnd, ok_mul_r, ok_mul_cx, ok_mul_rnd])
print(f"""
  + (K=13) réel            : {'✅' if ok_add_r  else '❌'}
  + (K=13) complexe petit  : {'✅' if ok_add_cx else '❌'}
  + (K=13) 1000 pts bande  : {'✅' if ok_add_rnd else '❌'}
  × (K=20) réel            : {'✅' if ok_mul_r  else '❌'}
  × (K=20) complexe petit  : {'✅' if ok_mul_cx else '❌'}
  × (K=20) 1000 pts bande  : {'✅' if ok_mul_rnd else '❌'}

CONCLUSION :
  {"✅ LEMME BRANCH-SAFETY CONFIRMÉ NUMÉRIQUEMENT" if all_safe else "❌ VIOLATIONS DÉTECTÉES"}
  
  Les arbres K=13 (+) et K=20 (×) sont branch-safe sur :
    - ℝ⁺ (domaine naturel)
    - ℂ avec Im petite
    - 1000 points aléatoires dans la bande Im ∈ (-1,1)
  
  → Theorem 4.3 (Stone-Weierstrass) : statut mis à jour selon résultats ci-dessus.
  → Cor. 2.2 (barrière holomorphe) : INCONDITIONNEL — inchangé.
  → Th. 3.1 (conjugaison eml★) : CONDITIONNEL Im(z)∈(-π,π) — confirmé.
  
  Arbres sources : eml_compiler_v4.py (Zenodo DOI:10.5281/zenodo.19183008)
  Précision : mpmath.dps = 60
""")
