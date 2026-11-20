import itertools

import mpmath as mp
import pytest

from double_toolbox import toolbox as B
from double_toolbox import eml as E
from double_toolbox import nand as N


# ── NAND ─────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("gate,table", [
    (N.not_("a"), {(0,): 1, (1,): 0}),
    (N.and_("a", "b"), {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 1}),
    (N.or_("a", "b"), {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 1}),
    (N.xor("a", "b"), {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 0}),
])
def test_truth_tables(gate, table):
    for inputs, output in table.items():
        env = dict(zip("ab", inputs))
        env.setdefault("b", 0)
        assert N.evaluate(gate, env) == output


def test_classic_gate_counts():
    assert N.gate_count([N.xor("a", "b")]) == 4
    assert N.gate_count(list(N.full_adder("a", "b", "r"))) == 9


def test_addition_and_multiplication_4_bits_exhaustive():
    A, B4 = N.bits("a", 4), N.bits("b", 4)
    s, p = N.addition(A, B4), N.multiplication(A, B4)
    for a, b in itertools.product(range(16), repeat=2):
        env = {**N.write_bits("a", a, 4), **N.write_bits("b", b, 4)}
        assert N.read_bits(s, env) == a + b
        assert N.read_bits(p, env) == a * b


# ── EML ──────────────────────────────────────────────────────────────────────
def test_known_lengths():
    assert E.size(E.add("x", "y")) == 27
    assert E.size(E.mul("x", "y")) == 41


def test_exact_constants():
    import sympy as sp
    assert E.constant(E.ZERO) == 0
    assert E.constant(E.MINUS_INFINITY) == E.MINUS_INF
    assert E.constant(E.I) == sp.I
    assert E.constant(E.PI) == sp.pi


def test_sine_cosine_on_the_reals():
    xs = [mp.mpf(k) * mp.mpf("0.37") for k in range(-60, 61)]
    assert all(ok for *_, ok in B.measure(E.sine("x"), mp.sin, xs))
    assert all(ok for *_, ok in B.measure(E.cosine("x"), mp.cos, xs))


def test_sine_near_the_real_axis():
    zs = [mp.mpc(a / 4, s * 1.1) for a in range(-40, 41) for s in (1, -1)]
    assert all(ok for *_, ok in B.measure(E.sine("x"), mp.sin, zs))


def test_the_wall_really_exists():
    # ln of a negative real: the tree returns the other side of the cut (−iπ instead of +iπ)
    v = E.evaluate(E.ln_("x"), {"x": -2})
    assert abs(v - (mp.log(2) - 1j * mp.pi)) < 1e-40
    # addition when |Im y| > π: wrong by a turn 2πi
    v = E.evaluate(E.add("x", "y"), {"x": 1, "y": mp.mpc(0, 4)})
    assert abs(v - (1 + 4j - 2j * mp.pi)) < 1e-40


def test_multiplication_without_condition():
    pts = [{"x": mp.mpc(a, b), "y": mp.mpc(c, d)} for a, b, c, d in itertools.product((-7, 2), (-9, 0.5), (-3, 5), (-8, 1))]
    assert all(ok for *_, ok in B.measure(E.mul("x", "y"), lambda x, y: x * y, pts))


def test_compare_the_two_worlds():
    for r in B.compare(13, 11):
        assert r.correct, r


# ── eml★ and the mirror ──────────────────────────────────────────────────────
def test_mirror_nine_symbols_everywhere():
    t = E.mirror("x")
    assert E.size(t) == 9
    zs = [mp.mpc(a, b) for a in (-50, -2, -0.3, 0.7, 9) for b in (-40, -mp.pi, -1, 0, 0.5, mp.pi, 7, 100)]
    assert all(ok for *_, ok in B.measure(t, mp.conj, zs))
    assert E.evaluate(t, {"x": 0}) == 0


def test_modulus_and_real_part():
    zs = [mp.mpc(a, b) for a in (-5, 0.5, 3) for b in (-20, -2, 0, 1, 30)]
    assert all(ok for *_, ok in B.measure(E.modulus_squared("x"), lambda x: abs(x) ** 2, zs))
    assert all(ok for *_, ok in B.measure(E.modulus("x"), lambda x: abs(x), zs))
    inside = [mp.mpc(a, b) for a in (-5, 0.5, 3) for b in (-3, -1, 0.2, 3)]
    assert all(ok for *_, ok in B.measure(E.real_part("x"), lambda x: mp.re(x), inside))
    outside = B.measure(E.real_part("x"), lambda x: mp.re(x), [mp.mpc(2, 5)])
    assert not outside[0][3]                     # outside the zone: wrong, as announced


# ── the sine on the NAND side ────────────────────────────────────────────────
def test_bare_nand_sine_has_a_wall():
    from double_toolbox import nand_sine as S
    assert abs(S.compute(1, 16, 13, full=False)[0] - mp.sin(1)) < 1e-3
    assert abs(S.compute(3, 16, 13, full=False)[0] - mp.sin(3)) > 0.5


def test_full_nand_sine_over_the_whole_range():
    from double_toolbox import nand_sine as S
    worst = 0
    for k in range(-63, 64):
        x = mp.mpf(k) + mp.mpf("0.37")
        if abs(x) < 63.9:
            worst = max(worst, abs(S.compute(x)[0] - mp.sin(x)))
    assert worst < 2 ** -12                      # better than 2 units of the last bit


def test_nand_sine_exactly_odd():
    from double_toolbox import nand_sine as S
    for x in (0.3, 1, 2.5, 17):
        assert S.compute(x)[0] == -S.compute(-x)[0]


def test_precision_curve():
    from double_toolbox import nand_sine as S
    rows = S.precision_curve(((12, 9), (16, 13)), points=21)
    assert all(worst * 2 ** fraction < 2 for _, fraction, _, worst in rows)
    assert rows[0][2] < rows[1][2]               # more bits = more gates


def test_compare_sine():
    for x in (1, 3, -50):
        for r in B.compare_sine(x):
            assert r.correct, r


# ── wall-free EML recipes ────────────────────────────────────────────────────
HARD = [mp.mpc(a, b) for a in (-50, -7, -2, -1, -0.3, 0.4, 1, 3, 60) for b in (-100, -7, -mp.pi, -1, 0, 0.5, mp.pi, 4, 90)]


@pytest.mark.parametrize("recipe,reference", [
    (E.ln_exact, mp.log), (E.sqrt_safe, mp.sqrt), (E.real_part_safe, mp.re),
    (E.imag_part, mp.im), (E.negate, lambda x: -x),
])
def test_safe_recipes_one_variable(recipe, reference):
    assert all(ok for *_, ok in B.measure(recipe("x"), lambda x: reference(x), HARD))


def test_safe_addition_and_subtraction():
    pairs = [{"x": a, "y": b} for a in HARD[::4] for b in HARD[::3]]
    pairs += [{"x": mp.mpf(3), "y": mp.mpf(-3)}, {"x": mp.mpc(1, 5), "y": mp.mpc(-1, -5)}]
    assert all(ok for *_, ok in B.measure(E.add_safe("x", "y"), lambda x, y: x + y, pairs))
    assert all(ok for *_, ok in B.measure(E.sub_safe("x", "y"), lambda x, y: x - y, pairs))


def test_the_holes_refuse_instead_of_lying():
    for env in ({"x": 0, "y": mp.mpc(2, 9)}, {"x": 3, "y": 0}):
        with pytest.raises(E.Undefined):
            E.evaluate(E.add_safe("x", "y"), env)


def test_safe_sine_cosine_in_the_complex_plane():
    zs = [mp.mpc(a * 0.7, b * 0.9) for a in range(-12, 13) for b in range(-8, 9)] + [mp.pi, -3 * mp.pi, mp.mpf(0)]
    assert all(ok for *_, ok in B.measure(E.sine_safe("x"), lambda x: mp.sin(x), zs))
    assert all(ok for *_, ok in B.measure(E.cosine_safe("x"), lambda x: mp.cos(x), zs + [mp.pi / 2]))


def test_recipe_table():
    t = {l[0]: l for l in E.recipe_table()}
    assert t["addition"][1] == 27 and t["multiplication"][1] == 41 and t["mirror"][1] == 9


# ── the eml brick in NAND gates ──────────────────────────────────────────────
def test_nand_brick_against_the_true_value():
    import random
    from double_toolbox import nand_eml as M
    random.seed(20261120)
    rounded = lambda z: mp.mpc(mp.nint(mp.re(z) * 65536), mp.nint(mp.im(z) * 65536)) / 65536
    cases = [(1, 1), (0, 1), (2, -3), (mp.mpc(-2, -1), mp.mpc(-0.5, 0.25)), (mp.mpc(-6, 40), mp.mpc(0.001, 0))]
    cases += [(mp.mpc(random.uniform(-6, 4.5), random.uniform(-100, 100)),
               mp.mpc(random.choice([1, -1]) * 10 ** random.uniform(-3, 2), random.choice([1, -1]) * 10 ** random.uniform(-3, 2)))
              for _ in range(25)]
    for x, y in cases:
        v, is_zero, _ = M.compute(x, y)
        true = mp.exp(rounded(x)) - mp.log(rounded(y))
        assert not is_zero
        assert abs(v - true) < 2 ** -15 * max(1, abs(mp.exp(mp.re(x)))), (x, y)


def test_nand_brick_signals_log_of_zero():
    from double_toolbox import nand_eml as M
    assert M.compute(1, 0)[1] == 1


def test_eml_trees_computed_by_nand_bricks():
    from double_toolbox import nand_eml as M
    for tree, env, true, tol in [(E.ln_("x"), {"x": 2}, mp.log(2), 1e-4), (E.add("x", "y"), {"x": 5, "y": 3}, 8, 1e-3),
                                 (E.mul("x", "y"), {"x": 13, "y": 11}, 143, 0.05),
                                 (E.mirror("x"), {"x": mp.mpc(3, 4)}, mp.mpc(3, -4), 1e-3)]:
        v, bricks, _ = M.evaluate_with_nand(tree, env)
        assert abs(v - true) < tol and bricks >= 1


# ── eml★ alone (construction proposed by ChatGPT on 2026-11-20, rebuilt here) ─
def test_mirror_in_star_only_297_symbols_strict():
    import sympy as sp
    from double_toolbox import star_only as P
    j = P.J("x")
    assert E.size(j) == 297 and E.gates_used(j) == {"eml*"}
    for z in HARD:
        if z != 0:
            assert abs(E.evaluate(j, {"x": z}, strict=True) - mp.conj(z)) < 1e-40
    for q in (sp.Rational(-1, 3), sp.Integer(-2), sp.Integer(1), -5 * sp.I):
        assert sp.simplify(E.exact(j, {"x": q}) - sp.conjugate(q)) == 0
    with pytest.raises(E.Undefined):                     # at zero, no strict mirror
        E.evaluate(j, {"x": 0}, strict=True)


def test_eml_rebuilt_from_star_only_exactly():
    import sympy as sp
    from double_toolbox import star_only as P
    t = P.eml_from_star("x", "y")
    assert E.gates_used(t) == {"eml*"}
    for x, y in ((sp.Rational(1, 2), sp.Integer(-2)), (sp.I, sp.Integer(-7)), (sp.Integer(0), sp.Rational(3, 2) + sp.I)):
        assert sp.simplify(E.exact(t, {"x": x, "y": y}) - (sp.exp(x) - sp.log(y))) == 0


def test_the_full_turn_is_necessary():
    from double_toolbox import star_only as P
    assert abs(E.evaluate(P.M("x", "y"), {"x": 1, "y": 3}, strict=True) - 3) < 1e-40
    without_turn = P.X(P.S(P.L(P.L("x")), P.I("y")))
    with pytest.raises(E.Undefined):                     # without the shift of 2πi, the input 1 passes through ln 0
        E.evaluate(without_turn, {"x": 1, "y": 3}, strict=True)


def test_strict_real_part_and_nand_in_star_only():
    from double_toolbox import star_only as P
    r = P.real_part_strict("x")
    for z in (mp.mpf(0), mp.mpc(0, 3), mp.mpc(0, 2 * mp.pi), mp.mpc(-2, 40), mp.mpc(1.5, -0.2)):
        assert abs(E.evaluate(r, {"x": z}, strict=True) - mp.re(z)) < 1e-40
    for x in (0, 1):
        for y in (0, 1):
            assert abs(E.evaluate(P.nand_bit("x", "y"), {"x": x, "y": y}, strict=True) - (1 - x * y)) < 1e-40


def test_short_mirror_in_star_only_101_symbols_strict():
    import sympy as sp
    from double_toolbox import star_only as P
    j = P.J_short("x")
    assert E.size(j) == 101 and E.gates_used(j) == {"eml*"}
    for z in HARD:
        if z != 0:
            assert abs(E.evaluate(j, {"x": z}, strict=True) - mp.conj(z)) < 1e-40
    for q in (sp.Rational(-1, 3), sp.Integer(-2), sp.Integer(1), sp.Rational(3, 2) + 2 * sp.I, -5 * sp.I):
        assert sp.simplify(E.exact(j, {"x": q}) - sp.conjugate(q)) == 0
    with pytest.raises(E.Undefined):
        E.evaluate(j, {"x": 0}, strict=True)
    for z in (mp.mpc(-3, 2), mp.mpf(-4), mp.mpf(1)):     # the "change of sign" piece alone
        assert abs(E.evaluate(P.P("x"), {"x": z}, strict=True) + z) < 1e-40
