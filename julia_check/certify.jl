# Second opinion on certify_branch_safety.py: same rule (the two gates, lemma A), but another language (Julia),
# another interval library (IntervalArithmetic.jl on MPFR, instead of mpmath.iv), and rewritten code:
# - complex exp and log redone by hand from the REAL functions (exp, log, sin, cos, atan);
# - exact constants through a small home-made symbolic computation (instead of SymPy).
# Run:  julia --project=. certify.jl
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
log_right(z::C) = C(lnmod(z), atan(z.im / z.re))                      # Re z > 0
log_upper(z::C) = C(lnmod(z), PI / interval(BigFloat, 2) - atan(z.re / z.im))    # Im z > 0
log_lower(z::C) = C(lnmod(z), -PI / interval(BigFloat, 2) - atan(z.re / z.im))   # Im z < 0
contains0(z::C) = inf(z.re) <= 0 <= sup(z.re) && inf(z.im) <= 0 <= sup(z.im)

struct Failure <: Exception; msg::String; end

# ---------- trees ----------
abstract type T end
struct One <: T end
struct Var <: T; name::Symbol; end
struct Eml <: T; a::T; b::T; star::Bool; end
struct MinusLn <: T; a::T; c::T; end          # a − ln c: only through lemma A
const ONE = One(); const X = Var(:x); const Y = Var(:y)
eml(a, b) = Eml(a, b, false)
E_exp(z) = eml(z, ONE)
E_ln(z) = eml(ONE, eml(eml(ONE, z), ONE))
E_zero() = eml(ONE, eml(eml(ONE, ONE), ONE))
E_neg(z) = eml(eml(ONE, eml(eml(ONE, E_zero()), ONE)), eml(z, ONE))
E_sub(a, b) = eml(E_ln(a), E_exp(b))
E_add(a, b) = E_sub(a, E_neg(b))
E_mul(a, b) = E_exp(E_add(E_ln(a), E_ln(b)))
E_mul_lemma(a, b) = E_exp(MinusLn(E_ln(a), E_exp(E_neg(E_ln(b)))))
E_add_lemma(a, b) = MinusLn(a, E_exp(E_neg(b)))
hasvar(t::One) = false; hasvar(t::Var) = true
hasvar(t::Eml) = hasvar(t.a) || hasvar(t.b); hasvar(t::MinusLn) = hasvar(t.a) || hasvar(t.c)
treesize(t::One) = 1; treesize(t::Var) = 1; treesize(t::Eml) = 1 + treesize(t.a) + treesize(t.b)

# ---------- exact constants on the extended reals: :zero, :pinf, :minf or a finite real term ("exp(1)", …) ----------
struct K; kind::Symbol; term::String; val::I; end
kz() = K(:zero, "0", iv(0)); kp() = K(:pinf, "", iv(0)); km() = K(:minf, "", iv(0))
function kexp(a::K)
    a.kind == :minf && return kz()
    a.kind == :pinf && return kp()
    a.kind == :zero && return K(:finite, "1", iv(1))
    K(:finite, "exp(" * a.term * ")", exp(a.val))
end
function klog(b::K)
    b.kind == :minf && throw(Failure("ln(−∞)"))
    b.kind == :zero && return km()
    b.kind == :pinf && return kp()
    b.term == "1" && return kz()
    m = match(r"^exp\((.*)\)$", b.term)
    if m !== nothing && balanced(m.captures[1])
        return K(:finite, String(m.captures[1]), log(b.val))        # ln(exp(s)) = s for real s
    end
    inf(b.val) > 0 || throw(Failure("log of a non-positive constant"))
    K(:finite, "ln(" * b.term * ")", log(b.val))
end
balanced(s) = (n = 0; for ch in s; ch == '(' && (n += 1); ch == ')' && (n -= 1); n < 0 && return false; end; n == 0)
function kminus(e::K, l::K)
    e.kind == :pinf && l.kind == :pinf && throw(Failure("∞ − ∞"))
    l.kind == :minf && return kp()
    e.kind == :pinf && return kp()
    l.kind == :pinf && return km()
    l.kind == :zero && return e
    e.kind == :finite && l.kind == :finite && e.term == l.term && return kz()   # s − s = 0, exact
    e.kind == :zero && return K(:finite, "(0 - " * l.term * ")", -l.val)
    K(:finite, "(" * e.term * " - " * l.term * ")", e.val - l.val)
end
kconst(t::One) = K(:finite, "1", iv(1))
kconst(t::Eml) = kminus(kexp(kconst(t.a)), klog(kconst(t.b)))          # real constants: conjugation changes nothing

# ---------- evaluation: finite complex value, or :pinf / :minf ----------
function ev(t::T, env, gate::Symbol, counter)
    t isa One && return C(iv(1), iv(0))
    t isa Var && return env[t.name]
    if !hasvar(t)
        k = kconst(t)
        k.kind == :pinf && return :pinf
        k.kind == :minf && return :minf
        return C(k.val, iv(0))
    end
    if t isa MinusLn
        a = ev(t.a, env, gate, counter); b = ev(t.c, env, :y, counter)
        (a isa Symbol || b isa Symbol) && throw(Failure("infinite value in a − ln c"))
        return a - logc(b, gate, counter, "a − ln c")
    end
    a = ev(t.a, env, :x, counter); b = ev(t.b, env, :y, counter)
    if b isa Symbol
        b == :minf && throw(Failure("ln(−∞)"))
        a == :pinf && throw(Failure("∞ − ∞"))
        return :minf
    end
    a == :pinf && return :pinf
    ea = a == :minf ? C(iv(0), iv(0)) : cexp(a)
    ea - logc(t.star ? conjc(b) : b, gate, counter, "eml")
end
function logc(z::C, gate, counter, site)
    contains0(z) && throw(Failure("log of a box that contains zero ($site)"))
    inf(z.re) > 0 && return log_right(z)
    inf(z.im) > 0 && return log_upper(z)
    sup(z.im) < 0 && return log_lower(z)
    if gate == :x                                   # straddling the negative axis, absorbed by the parent's exp
        counter[] += 1
        r = log_right(neg(z))
        return C(r.re, r.im + PI)
    end
    throw(Failure("branch-cut crossing, gate $gate ($site)"))
end

# ---------- subdivision ----------
struct Box; r0::BigFloat; r1::BigFloat; i0::BigFloat; i1::BigFloat; end
# A decimal bound given as TEXT ("0.1") is rounded outwards from the box; a number is taken as is.
lower(v) = v isa AbstractString ? parse(BigFloat, v, RoundDown) : BigFloat(v)
upper(v) = v isa AbstractString ? parse(BigFloat, v, RoundUp) : BigFloat(v)
Base.parse(::Type{BigFloat}, s::AbstractString, r::RoundingMode) = setrounding(() -> parse(BigFloat, s), BigFloat, r)
B(a, b, c, d) = Box(lower(a), upper(b), lower(c), upper(d))
civ(b::Box) = C(interval(b.r0, b.r1), interval(b.i0, b.i1))
function subdivide(b::Box)
    mr = (b.r0 + b.r1) / 2; mi = (b.i0 + b.i1) / 2
    [Box(b.r0, mr, b.i0, mi), Box(mr, b.r1, b.i0, mi), Box(b.r0, mr, mi, b.i1), Box(mr, b.r1, mi, b.i1)]
end
# Guaranteed: upper bound of distance² < lower bound of radius² (interval arithmetic); decimal radius given as text.
function in_disc(b::Box, c, r)
    ri = r isa AbstractString ? interval(parse(BigFloat, r, RoundDown), parse(BigFloat, r, RoundUp)) : interval(BigFloat(r))
    r2 = ri * ri
    all(sup((interval(x) - interval(real(c)))^2 + (interval(y) - interval(imag(c)))^2) < inf(r2) for x in (b.r0, b.r1), y in (b.i0, b.i1))
end

function certify(t, bx::Box, by::Box; depth = 10, disc = nothing)
    work = [(bx, by, 0)]; ok = 0; excluded = 0; dmax = 0; counter = Ref(0); t0 = time()
    while !isempty(work)
        x, y, d = pop!(work)
        if disc !== nothing && in_disc(x, disc[1], disc[2]); excluded += 1; continue; end
        try
            r = ev(t, Dict(:x => civ(x), :y => civ(y)), :root, counter)
            r isa Symbol && throw(Failure("infinite value at the root"))
            ok += 1; dmax = max(dmax, d)
        catch e
            e isa Failure || rethrow()
            d >= depth && return (certified = false, boxes = ok, absorbed = counter[], excluded = excluded, depth = d, s = round(time() - t0, digits = 1), reason = e.msg)
            for a in subdivide(x), b in subdivide(y); push!(work, (a, b, d + 1)); end
        end
    end
    ok == 0 && return (certified = false, boxes = 0, absorbed = 0, excluded = excluded, depth = 0, s = round(time() - t0, digits = 1), reason = "no box accepted (is everything inside the excluded disc?)")
    (certified = true, boxes = ok, absorbed = counter[], excluded = excluded, depth = dmax, s = round(time() - t0, digits = 1), reason = "")
end
const ALL = Bool[]
report(name, r) = (push!(ALL, r.certified); println(name, ": ", r.certified ? "CERTIFIED" : "FAILED", " — ", r.boxes, " boxes, ", r.absorbed, " absorbed crossings, ",
                       r.excluded, " in the disc, depth ", r.depth, ", ", r.s, " s", r.reason == "" ? "" : " | " * r.reason))

# ---------- sanity checks: the program must REFUSE whatever jumps ----------
G = B(-2, -1, -0.1, 0.1)
@assert kconst(E_zero()).kind == :zero
@assert kconst(eml(ONE, E_zero())).kind == :pinf && kconst(eml(ONE, eml(eml(ONE, E_zero()), ONE))).kind == :minf
@assert certify(eml(eml(ONE, X), ONE), G, G; depth = 3).certified            # x gate: absorbed
@assert !certify(eml(ONE, eml(ONE, X)), G, G; depth = 3).certified           # y gate: refused
@assert !certify(eml(ONE, X), G, G; depth = 3).certified                     # root: refused
@assert !certify(eml(eml(ONE, X), ONE), B(-0.1, 0.1, -0.1, 0.1), G; depth = 3).certified   # log of zero: refused
@assert !certify(E_mul(X, Y), B(0.9, 1.1, -0.1, 0.1), B(1.1, 1.5, -0.1, 0.1); depth = 4).certified
@assert !certify(E_exp(X), B(0.99, 1.01, -0.01, 0.01), G; disc = (complex(BigFloat(1)), "0.5")).certified   # everything in the disc: not a success
println("Refusal checks: passed. Sizes: E_add ", treesize(E_add(X, Y)), ", E_mul ", treesize(E_mul(X, Y)))

# ---------- the same certificates as the Python version ----------
D = B("0.1", 3, -1, 1)
report("1) x + y, x,y ∈ [0.1,3]+i[−1,1]", certify(E_add(X, Y), D, D))
report("2) x × y (original tree), same domain, outside |x−1| < 0.05", certify(E_mul(X, Y), D, D; disc = (complex(BigFloat(1)), "0.05")))
report("3) x × y (lemma form), same domain, disc included", certify(E_mul_lemma(X, Y), D, D))
H = 3
R = [B(-3, "-0.1", -H, H), B("0.1", 3, -H, H), B(-3, 3, "0.1", H), B(-3, 3, -H, "-0.1")]
for (name, t) in (("4) x + y", E_add(X, Y)), ("5) x × y (lemma form)", E_mul_lemma(X, Y)))
    tot = 0; ab = 0; ok = true; dm = 0
    for bx in R, by in R
        r = certify(t, bx, by; depth = 9); tot += r.boxes; ab += r.absorbed; ok &= r.certified; dm = max(dm, r.depth)
    end
    push!(ALL, ok)
    println(name, " on [−3,3]+i[−3,3] with |Re|,|Im| < 0.1 removed: ", ok ? "CERTIFIED" : "FAILED", " — ", tot, " boxes, ", ab, " absorbed crossings, depth ", dm)
end
P = B(-3, 3, -3, 3)
report("6) x + y (lemma form) on the full square [−3,3]+i[−3,3], zero included", certify(E_add_lemma(X, Y), P, P; depth = 9))
if !(all(ALL) && length(ALL) == 6)
    println("AT LEAST ONE EXPECTED CERTIFICATE FAILED"); exit(1)
end
println("The six expected certificates passed.")
