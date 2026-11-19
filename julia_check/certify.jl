# Deuxième avis sur certify_branch_safety.py : même règle (les deux portes, lemme A), mais un autre langage (Julia),
# une autre bibliothèque d'intervalles (IntervalArithmetic.jl sur MPFR, au lieu de mpmath.iv), et un code récrit :
# - exp et log complexes refaits à la main à partir des fonctions RÉELLES (exp, log, sin, cos, atan) ;
# - constantes exactes par un petit calcul symbolique maison (au lieu de SymPy).
# Lancer :  julia --project=. certify.jl
using IntervalArithmetic
setprecision(BigFloat, 200)
const I = Interval{BigFloat}
iv(a, b) = interval(BigFloat(a), BigFloat(b))
iv(a) = iv(a, a)
const PI = interval(BigFloat, π)

struct C; re::I; im::I; end
Base.:-(a::C, b::C) = C(a.re - b.re, a.im - b.im)
neg(a::C) = C(-a.re, -a.im)
conjc(a::C) = C(a.re, -a.im)
cexp(a::C) = (r = exp(a.re); C(r * cos(a.im), r * sin(a.im)))
lnmod(z::C) = log(z.re * z.re + z.im * z.im) / interval(BigFloat, 2)
log_droite(z::C) = C(lnmod(z), atan(z.im / z.re))                      # Re z > 0
log_haut(z::C)   = C(lnmod(z), PI / interval(BigFloat, 2) - atan(z.re / z.im))    # Im z > 0
log_bas(z::C)    = C(lnmod(z), -PI / interval(BigFloat, 2) - atan(z.re / z.im))   # Im z < 0
contient0(z::C) = inf(z.re) <= 0 <= sup(z.re) && inf(z.im) <= 0 <= sup(z.im)

struct Echec <: Exception; msg::String; end

# ---------- arbres ----------
abstract type T end
struct Un <: T end
struct Var <: T; nom::Symbol; end
struct Eml <: T; a::T; b::T; star::Bool; end
struct MoinsLn <: T; a::T; c::T; end          # a − ln c : seulement via le lemme A
const ONE = Un(); const X = Var(:x); const Y = Var(:y)
eml(a, b) = Eml(a, b, false)
E_exp(z) = eml(z, ONE)
E_ln(z) = eml(ONE, eml(eml(ONE, z), ONE))
E_zero() = eml(ONE, eml(eml(ONE, ONE), ONE))
E_neg(z) = eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))
E_sub(a, b) = eml(E_ln(a), E_exp(b))
E_add(a, b) = E_sub(a, E_neg(b))
E_mul(a, b) = E_exp(E_add(E_ln(a), E_ln(b)))
E_mul_lemme(a, b) = E_exp(MoinsLn(E_ln(a), E_exp(E_neg(E_ln(b)))))
E_add_lemme(a, b) = MoinsLn(a, E_exp(E_neg(b)))
avar(t::Un) = false; avar(t::Var) = true
avar(t::Eml) = avar(t.a) || avar(t.b); avar(t::MoinsLn) = avar(t.a) || avar(t.c)
taille(t::Un) = 1; taille(t::Var) = 1; taille(t::Eml) = 1 + taille(t.a) + taille(t.b)

# ---------- constantes exactes sur les réels étendus : :zero, :pinf, :minf ou un terme réel fini ("exp(1)", …) ----------
struct K; sorte::Symbol; terme::String; val::I; end
kz() = K(:zero, "0", iv(0)); kp() = K(:pinf, "", iv(0)); km() = K(:minf, "", iv(0))
function kexp(a::K)
    a.sorte == :minf && return kz()
    a.sorte == :pinf && return kp()
    a.sorte == :zero && return K(:fini, "1", iv(1))
    K(:fini, "exp(" * a.terme * ")", exp(a.val))
end
function klog(b::K)
    b.sorte == :minf && throw(Echec("ln(−∞)"))
    b.sorte == :zero && return km()
    b.sorte == :pinf && return kp()
    b.terme == "1" && return kz()
    m = match(r"^exp\((.*)\)$", b.terme)
    if m !== nothing && equilibre(m.captures[1])
        return K(:fini, String(m.captures[1]), log(b.val))        # ln(exp(s)) = s pour s réel
    end
    inf(b.val) > 0 || throw(Echec("log d'une constante non positive"))
    K(:fini, "ln(" * b.terme * ")", log(b.val))
end
equilibre(s) = (n = 0; for ch in s; ch == '(' && (n += 1); ch == ')' && (n -= 1); n < 0 && return false; end; n == 0)
function kmoins(e::K, l::K)
    e.sorte == :pinf && l.sorte == :pinf && throw(Echec("∞ − ∞"))
    l.sorte == :minf && return kp()
    e.sorte == :pinf && return kp()
    l.sorte == :pinf && return km()
    l.sorte == :zero && return e
    e.sorte == :fini && l.sorte == :fini && e.terme == l.terme && return kz()   # s − s = 0, exact
    e.sorte == :zero && return K(:fini, "(0 - " * l.terme * ")", -l.val)
    K(:fini, "(" * e.terme * " - " * l.terme * ")", e.val - l.val)
end
konst(t::Un) = K(:fini, "1", iv(1))
konst(t::Eml) = kmoins(kexp(konst(t.a)), klog(konst(t.b)))          # constantes réelles : la conjugaison ne change rien

# ---------- évaluation : valeur complexe finie, ou :pinf / :minf ----------
function ev(t::T, env, porte::Symbol, compte)
    t isa Un && return C(iv(1), iv(0))
    t isa Var && return env[t.nom]
    if !avar(t)
        k = konst(t)
        k.sorte == :pinf && return :pinf
        k.sorte == :minf && return :minf
        return C(k.val, iv(0))
    end
    if t isa MoinsLn
        a = ev(t.a, env, porte, compte); b = ev(t.c, env, :y, compte)
        (a isa Symbol || b isa Symbol) && throw(Echec("valeur infinie dans a − ln c"))
        return a - logc(b, porte, compte, "a − ln c")
    end
    a = ev(t.a, env, :x, compte); b = ev(t.b, env, :y, compte)
    if b isa Symbol
        b == :minf && throw(Echec("ln(−∞)"))
        a == :pinf && throw(Echec("∞ − ∞"))
        return :minf
    end
    a == :pinf && return :pinf
    ea = a == :minf ? C(iv(0), iv(0)) : cexp(a)
    ea - logc(t.star ? conjc(b) : b, porte, compte, "eml")
end
function logc(z::C, porte, compte, ou)
    contient0(z) && throw(Echec("log d'un rectangle qui contient zéro ($ou)"))
    inf(z.re) > 0 && return log_droite(z)
    inf(z.im) > 0 && return log_haut(z)
    sup(z.im) < 0 && return log_bas(z)
    if porte == :x                                  # à cheval sur l'axe négatif, absorbé par l'exp du parent
        compte[] += 1
        r = log_droite(neg(z))
        return C(r.re, r.im + PI)
    end
    throw(Echec("traversée de la coupure, porte $porte ($ou)"))
end

# ---------- découpage ----------
struct Boite; r0::BigFloat; r1::BigFloat; i0::BigFloat; i1::BigFloat; end
# Une borne décimale donnée en TEXTE ("0.1") est arrondie vers l'extérieur du rectangle ; un nombre est pris tel quel.
bas(v) = v isa AbstractString ? parse(BigFloat, v, RoundDown) : BigFloat(v)
haut(v) = v isa AbstractString ? parse(BigFloat, v, RoundUp) : BigFloat(v)
Base.parse(::Type{BigFloat}, s::AbstractString, r::RoundingMode) = setrounding(() -> parse(BigFloat, s), BigFloat, r)
B(a, b, c, d) = Boite(bas(a), haut(b), bas(c), haut(d))
civ(b::Boite) = C(interval(b.r0, b.r1), interval(b.i0, b.i1))
function couper(b::Boite)
    mr = (b.r0 + b.r1) / 2; mi = (b.i0 + b.i1) / 2
    [Boite(b.r0, mr, b.i0, mi), Boite(mr, b.r1, b.i0, mi), Boite(b.r0, mr, mi, b.i1), Boite(mr, b.r1, mi, b.i1)]
end
# Garanti : majorant de la distance² < minorant de rayon² (arithmétique d'intervalles) ; rayon décimal donné en texte.
function dans_disque(b::Boite, c, r)
    ri = r isa AbstractString ? interval(parse(BigFloat, r, RoundDown), parse(BigFloat, r, RoundUp)) : interval(BigFloat(r))
    r2 = ri * ri
    all(sup((interval(x) - interval(real(c)))^2 + (interval(y) - interval(imag(c)))^2) < inf(r2) for x in (b.r0, b.r1), y in (b.i0, b.i1))
end

function certifier(t, bx::Boite, by::Boite; prof = 10, disque = nothing)
    travail = [(bx, by, 0)]; ok = 0; exclues = 0; pmax = 0; compte = Ref(0); t0 = time()
    while !isempty(travail)
        x, y, d = pop!(travail)
        if disque !== nothing && dans_disque(x, disque[1], disque[2]); exclues += 1; continue; end
        try
            r = ev(t, Dict(:x => civ(x), :y => civ(y)), :racine, compte)
            r isa Symbol && throw(Echec("valeur infinie à la racine"))
            ok += 1; pmax = max(pmax, d)
        catch e
            e isa Echec || rethrow()
            d >= prof && return (certifie = false, boites = ok, absorbees = compte[], exclues = exclues, prof = d, s = round(time() - t0, digits = 1), raison = e.msg)
            for a in couper(x), b in couper(y); push!(travail, (a, b, d + 1)); end
        end
    end
    ok == 0 && return (certifie = false, boites = 0, absorbees = 0, exclues = exclues, prof = 0, s = round(time() - t0, digits = 1), raison = "aucun rectangle accepté (tout est dans le disque exclu ?)")
    (certifie = true, boites = ok, absorbees = compte[], exclues = exclues, prof = pmax, s = round(time() - t0, digits = 1), raison = "")
end
const TOUS = Bool[]
dire(nom, r) = (push!(TOUS, r.certifie); println(nom, " : ", r.certifie ? "CERTIFIÉ" : "ÉCHEC", " — ", r.boites, " rectangles, ", r.absorbees, " traversées absorbées, ",
                       r.exclues, " dans le disque, profondeur ", r.prof, ", ", r.s, " s", r.raison == "" ? "" : " | " * r.raison))

# ---------- contrôles de bon sens : le programme doit REFUSER ce qui saute ----------
G = B(-2, -1, -0.1, 0.1)
@assert konst(E_zero()).sorte == :zero
@assert konst(eml(ONE, E_zero())).sorte == :pinf && konst(eml(ONE, eml(eml(ONE, E_zero()), ONE))).sorte == :minf
@assert certifier(eml(eml(ONE, X), ONE), G, G; prof = 3).certifie            # porte x : absorbé
@assert !certifier(eml(ONE, eml(ONE, X)), G, G; prof = 3).certifie           # porte y : refusé
@assert !certifier(eml(ONE, X), G, G; prof = 3).certifie                     # racine : refusé
@assert !certifier(eml(eml(ONE, X), ONE), B(-0.1, 0.1, -0.1, 0.1), G; prof = 3).certifie   # log de zéro : refusé
@assert !certifier(E_mul(X, Y), B(0.9, 1.1, -0.1, 0.1), B(1.1, 1.5, -0.1, 0.1); prof = 4).certifie
@assert !certifier(E_exp(X), B(0.99, 1.01, -0.01, 0.01), G; disque = (complex(BigFloat(1)), "0.5")).certifie   # tout dans le disque : pas un succès
println("Contrôles de refus : passés. Tailles : E_add ", taille(E_add(X, Y)), ", E_mul ", taille(E_mul(X, Y)))

# ---------- les mêmes certificats que la version Python ----------
D = B("0.1", 3, -1, 1)
dire("1) x + y, x,y ∈ [0,1;3]+i[−1,1]", certifier(E_add(X, Y), D, D))
dire("2) x × y (arbre d'origine), même domaine, hors |x−1| < 0,05", certifier(E_mul(X, Y), D, D; disque = (complex(BigFloat(1)), "0.05")))
dire("3) x × y (forme lemme), même domaine, disque compris", certifier(E_mul_lemme(X, Y), D, D))
H = 3
R = [B(-3, "-0.1", -H, H), B("0.1", 3, -H, H), B(-3, 3, "0.1", H), B(-3, 3, -H, "-0.1")]
for (nom, t) in (("4) x + y", E_add(X, Y)), ("5) x × y (forme lemme)", E_mul_lemme(X, Y)))
    tot = 0; ab = 0; ok = true; pm = 0
    for bx in R, by in R
        r = certifier(t, bx, by; prof = 9); tot += r.boites; ab += r.absorbees; ok &= r.certifie; pm = max(pm, r.prof)
    end
    push!(TOUS, ok)
    println(nom, " sur [−3,3]+i[−3,3] privé de |Re|,|Im| < 0,1 : ", ok ? "CERTIFIÉ" : "ÉCHEC", " — ", tot, " rectangles, ", ab, " traversées absorbées, profondeur ", pm)
end
P = B(-3, 3, -3, 3)
dire("6) x + y (forme lemme) sur le carré plein [−3,3]+i[−3,3], zéro compris", certifier(E_add_lemme(X, Y), P, P; prof = 9))
if !(all(TOUS) && length(TOUS) == 6)
    println("AU MOINS UN CERTIFICAT ATTENDU A ÉCHOUÉ"); exit(1)
end
println("Les six certificats attendus sont passés.")
