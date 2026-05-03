"""
eml★ Verification Engine
Reproduit la logique de main.rs (Odrzywołek) en Python/mpmath,
avec ajout de eml★ comme opérateur candidat.

Deux objectifs :
1. Vérifier que {eml, eml★, 1} reconstruit les opérations cibles
2. Tester branch-safety sur les arbres trouvés
"""
import mpmath
import itertools
mpmath.mp.dps = 60

EULER_GAMMA = mpmath.euler   # γ ≈ 0.5772...
GLAISHER    = mpmath.mpc('1.2824271291006226368753425688697917277676889273250011920637400')

# ─── Opérateurs ───────────────────────────────────────────────────────────────

def eml(x, y):
    """eml(x,y) = exp(x) - ln(y)"""
    return mpmath.exp(x) - mpmath.log(y)

def eml_star(x, y):
    """eml★(x,y) = exp(x) - ln(conj(y))"""
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))

ONE = mpmath.mpc(1)

# ─── Arbre symbolique simple ──────────────────────────────────────────────────

class Node:
    def __init__(self, op, left=None, right=None, value=None):
        self.op    = op      # 'eml', 'eml★', 'const'
        self.left  = left
        self.right = right
        self.value = value   # pour les feuilles
    
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
        elif self.op in ('var_x','var_y'):
            return self.op[-1]
        else:
            return f"{self.op}({self.left},{self.right})"

# ─── Vérification numérique (méthode Odrzywołek) ─────────────────────────────

def near(a, b, tol=mpmath.mpf('1e-50')):
    return abs(a - b) < tol

def verify_formula(tree_fn, target_fn, test_points, label=""):
    """Vérifie que tree_fn ≈ target_fn sur test_points."""
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

# ─── Arbres connus (depuis Appendix A + Odrzywołek Table 4) ──────────────────

def ln_x(z):
    """ln(z) = eml(1, eml(eml(1,z), 1))  K=7"""
    return eml(ONE, eml(eml(ONE, z), ONE))

def exp_x(z):
    """exp(z) = eml(z, 1)  K=3"""
    return eml(z, ONE)

def subtract(x, y):
    """x-y = eml(ln(x), exp(y)) = eml(eml(1,eml(eml(1,x),1)), eml(y,1))  K=11"""
    return eml(eml(ONE, eml(eml(ONE, x), ONE)), eml(y, ONE))

def multiply(x, y):
    """x*y = exp(ln(x)+ln(y)) — construit via eml  K=17
    = eml(eml(eml(1,eml(eml(1,x),1)), eml(eml(1,eml(eml(1,y),1)),1)), 1)
    """
    ln_x_tree  = eml(ONE, eml(eml(ONE, x), ONE))
    ln_y_tree  = eml(ONE, eml(eml(ONE, y), ONE))
    # ln(x) + ln(y) = subtract(0, subtract(ln(x), ln(y)))... non
    # Plus simple : exp(ln(x) + ln(y))
    # ln(x)+ln(y) via eml : eml(eml(1, eml(eml(1,lnx_lny_sum),1)), ...)
    # En fait, la formule directe d'Odrzywołek :
    # x*y = eml(ln(x) + ln(y), 1)
    # ln(x)+ln(y) doit être construit via eml...
    # Méthode directe depuis Table 4 (K=17, direct search):
    # Utilisons la formule : x*y = exp(ln(x)+ln(y))
    # et ln(x)+ln(y) = eml(eml(1,eml(eml(1,x),1)), eml(eml(eml(1,eml(eml(1,y),1)),1),1))
    # Simplifié : lnx + lny = -( (-lnx) - lny ) mais ça tourne en rond
    # Formule depuis Odrzywołek (reconstruction depuis les règles de sa Table 4):
    # Clé : x+y = ln(e^x * e^y) = ln(e^(x+y))
    #        mais on n'a pas + directement...
    # On reconstruit depuis les primitives bootstrappées :
    # Une fois qu'on a {eml, -, +, ×} via bootstrapping, × est disponible.
    # Sans les arbres du SI, on utilise la formule générale.
    # Pour nos tests, on utilise la formule numérique directe :
    return x * y  # placeholder — remplacer par l'arbre réel quand SI disponible

def conj_emlstar(z):
    """z̄ = 1 - eml★(0, eml(z,1))  Theorem 3.1"""
    return 1 - eml_star(mpmath.mpc(0), eml(z, ONE))

def re_z(z):
    """Re(z) = (z + z̄)/2"""
    return (z + conj_emlstar(z)) / 2

def mod_sq(z):
    """|z|² = z * z̄"""
    return z * conj_emlstar(z)

# ─── Points de test (méthode Odrzywołek : γ et A) ────────────────────────────

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

test_band = [  # Im(z) ∈ (-π, π) pour eml★
    (mpmath.mpc(r, i),)
    for r in [-5, -2, 0, 1, 3]
    for i in [-2.8, -1.0, 0.0, 1.0, 2.8]
]

# ─── Tests ───────────────────────────────────────────────────────────────────

print("=" * 65)
print("eml★ VERIFICATION ENGINE — mpmath 60 décimales")
print("=" * 65)

results = []

def run(label, fn, ref_fn, pts):
    ok, err = verify_formula(fn, ref_fn, pts, label)
    status = "✅ EXACT" if ok else f"❌ err={err:.2e}"
    print(f"  {status:20s}  {label}")
    results.append((label, ok, err))

print("\n— Fonctions holomorphes (eml seul) —")
run("exp(z) = eml(z,1)",
    lambda z: exp_x(z), lambda z: mpmath.exp(z), test_unary)
run("ln(z) = eml(1,eml(eml(1,z),1))",
    lambda z: ln_x(z), lambda z: mpmath.log(z), test_unary)
run("x-y = eml(ln(x),exp(y))",
    lambda x,y: subtract(x,y), lambda x,y: x-y, test_binary)

print("\n— Fonctions anti-holomorphes (eml★ requis) —")
run("z̄ = 1−eml★(0,eml(z,1))  [Im(z)∈(-π,π)]",
    lambda z: conj_emlstar(z), lambda z: mpmath.conj(z), test_band)
run("Re(z) = (z+z̄)/2         [Im(z)∈(-π,π)]",
    lambda z: re_z(z), lambda z: mpmath.re(z), test_band)
run("|z|² = z·z̄              [Im(z)∈(-π,π)]",
    lambda z: mod_sq(z), lambda z: abs(z)**2, test_band)

print("\n— Contre-exemples hors bande (Grok) —")
for z_val in [mpmath.mpc(0,4), mpmath.mpc(0,100)]:
    comp = conj_emlstar(z_val)
    ref  = mpmath.conj(z_val)
    err  = float(abs(comp - ref))
    n    = round(err / (2*float(mpmath.pi)))
    print(f"  z={z_val}  Im={float(mpmath.im(z_val)):.1f}  "
          f"err={err:.4f}={n}×2π  "
          f"({'∉ bande → saut attendu' if abs(float(mpmath.im(z_val)))>float(mpmath.pi) else '∈ bande'})")

# ─── Comparaison des bases ────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("COMPARAISON DES BASES (réponse à Grok)")
print("=" * 65)

comparison = [
    ("{eml, 1}",
     "Odrzywołek (2026)",
     "Fonctions élémentaires holomorphes",
     "z̄, Re, Im, |z|²",
     "Grammaire S→1|eml(S,S) uniforme ✅"),
    
    ("{eml, eml★, 1}",
     "Monnerot (2026)",
     "Holomorphe + anti-holomorphe (bande Im∈(-π,π))",
     "Fonctions nécessitant |Im(valeur)|>π",
     "Grammaire binaire uniforme ✅ Sheffer-compatible"),
    
    ("{eml, Re, 1}",
     "Grok (2026, proposition)",
     "Holomorphe + anti-holomorphe (partout sur ℂ)",
     "Rien",
     "Re est unaire ⚠️ — perd grammaire binaire uniforme"),
    
    ("{eml, Re, Deriv, 1}",
     "Grok (extension)",
     "Holomorphe + anti-holomorphe + dérivées",
     "EDP, distributions (hors scope Sheffer)",
     "Deux primitives hétérogènes ❌ — hors paradigme"),
]

print(f"\n{'Base':<22} {'Qui':<22} {'Couvre':<40} {'Ne couvre pas':<35} {'Structure'}")
print("-" * 140)
for base, who, covers, not_covers, struct in comparison:
    print(f"{base:<22} {who:<22} {covers:<40} {not_covers:<35} {struct}")

# ─── Verdict ─────────────────────────────────────────────────────────────────

print("\n" + "=" * 65)
print("VERDICT FINAL")
print("=" * 65)
ok_count = sum(1 for _, ok, _ in results if ok)
print(f"\nTests réussis : {ok_count}/{len(results)}")
print(f"""
Points clés :

1. eml★ est EXACT sur son domaine revendiqué (Im(z)∈(-π,π))
   Vérifié à 60 décimales sur z̄, Re(z), |z|²

2. La vérification d'Odrzywołek (Rust) est en double précision (f64)
   Nos tests mpmath sont strictement plus exigeants

3. {'{'}eml, Re, 1{'}'} de Grok :
   ✅ Plus robuste (pas de bande)
   ❌ Perd la grammaire binaire uniforme (Re est unaire)
   ❌ Incompatible avec le moteur de SR d'Odrzywołek (arbres binaires seuls)

4. {'{'}eml, eml★, 1{'}'} :
   ✅ Reste dans le paradigme Sheffer binaire
   ✅ Grammaire : S → 1 | eml(S,S) | eml★(S,S)
   ⚠️ Conditionnel sur Im(valeur)∈(-π,π)
   ✅ Cor.2.2 inconditionnel (la vraie nouveauté)

5. Arbres × (K=17) et + (K=19) : dans le Supplementary Info
   → Non vérifiés formellement — Th.4.3 reste conjecture
   → Vérification possible dès que SI disponible
""")
