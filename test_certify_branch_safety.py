"""Tests of the branch-safety certificate (certify_branch_safety.py).
Run: python3 -m pytest test_certify_branch_safety.py -q   (or the ~/cutlattice environment: .venv/bin/python -m pytest)"""
import mpmath as mp
import pytest
import certify_branch_safety as cb

X, Y, ONE = cb.X, cb.Y, cb.ONE
eml = cb.eml


def B(a, b, c, d):
    return cb.Box(a, b, c, d)


# ---------- constants (no variable): exact values on the extended reals ----------

def test_zero_is_exact():
    assert cb.constant(cb.E_zero()) == 0


def test_the_E_neg_chain_goes_through_plus_and_minus_infinity():
    a = eml(ONE, cb.E_zero())                    # e − ln 0
    assert cb.constant(a) == cb.PLUS_INF
    assert cb.constant(eml(a, ONE)) == cb.PLUS_INF      # exp(+∞) − ln 1
    assert cb.constant(eml(ONE, eml(a, ONE))) == cb.MINUS_INF  # e − ln(+∞)


def test_constant_refuses_the_undefined():
    a = eml(ONE, cb.E_zero())                    # +∞
    with pytest.raises(cb.Undefined):
        cb.constant(eml(a, a))                   # exp(+∞) − ln(+∞) = ∞ − ∞


# ---------- the trees do compute what they say ----------

@pytest.mark.parametrize("a,b", [(mp.mpc(1.3, 0.1), mp.mpc(0.7, 0.2)), (mp.mpc(0.4, -0.5), mp.mpc(2.5, 0.3)), (mp.mpc(0.7, 1e-12), mp.mpc(1.3))])
def test_trees_correct_at_a_point(a, b):
    assert abs(cb.point_value(cb.E_add(X, Y), a, b) - (a + b)) < mp.mpf(10) ** -30
    assert abs(cb.point_value(cb.E_mul(X, Y), a, b) - a * b) < mp.mpf(10) ** -30


# ---------- the rule of the two gates ----------

LEFT = B(-2, -1, -0.1, 0.1)                    # straddling the negative real axis, away from 0


def test_crossing_at_the_x_gate_is_absorbed():
    t = eml(eml(ONE, X), ONE)                  # the node eml(1, x) enters through the x gate: exp(e − ln x) = e^e / x
    r = cb.certify(t, {"x": LEFT}, depth=4)
    assert r["certified"] and r["absorbed"] >= 1
    v = cb.enclose(t, {"x": LEFT.point(-1.5, 0.05)})
    expected = mp.exp(mp.e) / mp.mpc(-1.5, 0.05)
    assert v.real.a <= expected.real <= v.real.b and v.imag.a <= expected.imag <= v.imag.b


def test_crossing_at_the_y_gate_is_refused():
    t = eml(ONE, eml(ONE, X))                  # the node eml(1, x) enters through the y gate: its log is taken, the jump remains
    r = cb.certify(t, {"x": LEFT}, depth=4)
    assert not r["certified"] and "y gate" in r["reason"]


def test_crossing_at_the_root_is_refused():
    r = cb.certify(eml(ONE, X), {"x": LEFT}, depth=4)
    assert not r["certified"] and "root" in r["reason"]


def test_log_of_variable_zero_is_refused():
    r = cb.certify(eml(eml(ONE, X), ONE), {"x": B(-0.1, 0.1, -0.1, 0.1)}, depth=3)
    assert not r["certified"] and "zero" in r["reason"]


def test_rotation_keeps_the_half_turn():
    """ln(−y) + iπ: difference 0 or 2πi from the principal log; without the + iπ the sign would be wrong."""
    for y in (mp.mpc(-1.5, 0.05), mp.mpc(-1.5, -0.05)):
        k = (cb.rotated_log_point(y) - mp.log(y)) / (2j * mp.pi)
        assert abs(k - mp.nint(k.real)) < 1e-40 and abs(mp.exp(cb.rotated_log_point(y)) - y) < 1e-40


# ---------- the real trees ----------

def test_addition_certified():
    r = cb.certify(cb.E_add(X, Y), {"x": B(0.5, 2, -0.3, 0.3), "y": B(0.5, 2, -0.3, 0.3)}, depth=6)
    assert r["certified"], r


def test_multiplication_fails_near_x_equal_1_without_disc():
    r = cb.certify(cb.E_mul(X, Y), {"x": B(0.9, 1.1, -0.1, 0.1), "y": B(1.1, 1.5, -0.1, 0.1)}, depth=5, budget=60)
    assert not r["certified"]


def test_multiplication_certified_outside_the_disc_around_1():
    r = cb.certify(cb.E_mul(X, Y), {"x": B(0.5, 1.5, -0.3, 0.3), "y": B(1.1, 1.5, -0.1, 0.1)}, depth=8, budget=240,
                   disc_x=(1, 0.1))
    assert r["certified"], r
    assert r["excluded"] >= 1 and r["absorbed"] >= 1        # we did cross (real x < 1) and did exclude the disc


def test_the_disc_does_not_hide_everything():
    with pytest.raises(ValueError):
        cb.certify(cb.E_mul(X, Y), {"x": B(0.95, 1.05, -0.05, 0.05), "y": B(1.1, 1.5, -0.1, 0.1)}, disc_x=(1, 0.2))


# ---------- the disc around x = 1: lemma exp(E_ln(w)) = w, then the program verifies all the rest ----------

@pytest.mark.parametrize("w", [mp.mpc(-0.3, 1e-12), mp.mpc(-0.3, -1e-12), mp.mpc(2, 0.5), mp.mpc(-4, 0), mp.mpc(1e-8, -1e-9)])
def test_lemma_exp_E_ln_is_the_identity(w):
    assert abs(mp.exp(cb.point_value(cb.E_ln(X), w)) - w) < mp.mpf(10) ** -35


@pytest.mark.parametrize("a,b", [(mp.mpc(1 + 1e-6, 1e-7), mp.mpc(1.3, 0.1)), (mp.mpc(0.97, -0.02), mp.mpc(1.2, 0)), (mp.mpc(1.03, 0.04), mp.mpc(1.4, -0.05))])
def test_the_lemma_version_computes_the_same_product(a, b):
    assert abs(cb.point_value(cb.E_mul_lemma(X, Y), a, b) - a * b) < mp.mpf(10) ** -30
    assert abs(cb.point_value(cb.E_mul_lemma(X, Y), a, b) - cb.point_value(cb.E_mul(X, Y), a, b)) < mp.mpf(10) ** -30


def test_the_disc_is_certified_with_the_lemma():
    r = cb.certify(cb.E_mul_lemma(X, Y), {"x": B(0.9, 1.1, -0.1, 0.1), "y": B(1.1, 1.5, -0.1, 0.1)}, depth=6, budget=120)
    assert r["certified"], r


def test_minus_ln_a_jumping_value_is_allowed_only_at_the_x_gate():
    t = cb.minus_ln(cb.eml(ONE, X), ONE)        # a = e − ln x jumps by 2πi for negative real x
    r = cb.certify(t, {"x": LEFT}, depth=3)
    assert not r["certified"] and "root" in r["reason"]                           # at the root: refused
    r = cb.certify(cb.E_exp(t), {"x": LEFT}, depth=3)
    assert r["certified"] and r["absorbed"] >= 1                                  # under an exponential: erased
    r = cb.certify(cb.eml(ONE, t), {"x": LEFT}, depth=3)
    assert not r["certified"]                                                     # under a log: refused


def test_product_certified_for_negative_x_with_the_lemma():
    r = cb.certify(cb.E_mul_lemma(X, Y), {"x": B(-3, -0.1, -1, 1), "y": B(0.1, 3, -1, 1)}, depth=7, budget=200)
    assert r["certified"] and r["absorbed"] >= 1, r


# ---------- rigour of the bounds and of the constants (review of 2026-11-19) ----------

def test_a_decimal_bound_given_as_text_is_rounded_outward():
    b = B("0.1", "0.3", "-0.1", "0.7")
    tenth, three, seven = mp.mpf(1) / 10, mp.mpf(3) / 10, mp.mpf(7) / 10
    with mp.workdps(80):
        assert b.r[0] <= mp.mpf(1) / 10 <= mp.mpf(3) / 10 <= b.r[1] and b.i[0] <= -mp.mpf(1) / 10 and mp.mpf(7) / 10 <= b.i[1]
    assert B(0.5, 2, -1, 1).r == (mp.mpf(0.5), mp.mpf(2))            # a Python number is taken as is


def test_the_excluded_disc_is_guaranteed_inside_the_stated_disc():
    assert B(0.97, 1.03, -0.03, 0.03).in_disc(1, "0.05")
    assert not B(0.96, 1.04, -0.03, 0.03).in_disc(1, "0.05")         # corner at distance exactly 0.05: not STRICTLY inside
    assert not B(0.97, 1.03, -0.045, 0.045).in_disc(1, "0.05")


def test_the_interval_of_a_constant_contains_its_value():
    t = cb.eml(cb.eml(ONE, ONE), ONE)                                 # e^e ≈ 15.15: finite, nonzero
    z = cb._constant_interval(t)
    with mp.workdps(90):
        v = mp.exp(mp.e)
        assert mp.mpf(z.real._mpi_[0]) < v < mp.mpf(z.real._mpi_[1])
    assert mp.mpf(z.real._mpi_[1]) - mp.mpf(z.real._mpi_[0]) > 0       # nonzero width


def test_the_plus_and_times_trees_use_no_nonzero_finite_constant():
    cb.FINITE_CONSTANTS[0] = 0
    d = {"x": B(1.2, 1.4, 0.1, 0.2), "y": B(0.6, 0.8, 0.1, 0.2)}
    for t in (cb.E_add(X, Y), cb.E_mul(X, Y), cb.E_mul_lemma(X, Y), cb.E_add_lemma(X, Y)):
        assert cb.certify(t, d, depth=4)["certified"]
    assert cb.FINITE_CONSTANTS[0] == 0


def test_the_whole_domain_inside_the_disc_is_refused():
    with pytest.raises(ValueError):
        cb.certify(cb.E_exp(X), {"x": B(0.99, 1.01, -0.01, 0.01)}, disc_x=(1, "0.5"))
