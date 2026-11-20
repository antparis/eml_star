"""Double toolbox NAND / EML.

  python -m double_toolbox                 full demonstration
  python -m double_toolbox 200 55          addition and multiplication of two integers (0 to 255) in the two worlds
  python -m double_toolbox sine 2.5        the sine of an angle (−64 to 64) in the two worlds
  python -m double_toolbox recipes         the table of EML recipes: short / wall-free
  python -m double_toolbox curve           what precision costs on the NAND side
  python -m double_toolbox brick           the eml brick built in NAND gates, and EML trees computed with it
  python -m double_toolbox.certificates    the interval certificates of the wall-free recipes
"""
import sys

import mpmath as mp

from . import toolbox as B
from . import eml as E
from . import nand_sine as S


def _pretty(v, n=15):
    return mp.nstr(mp.chop(v, 1e-40), n)


def show_computation(a, b):
    print(f"\n=== The two worlds at the same time: {a} and {b} ===\n")
    print(f"{'world':6} {'operation':15} {'result':>10}  {'correct':7} {'bricks':>8} {'symbols without sharing':>23}")
    for r in B.compare(a, b):
        v = r.value if r.world == "NAND" else _pretty(r.value, 12)
        print(f"{r.world:6} {r.operation:15} {str(v):>10}  {'yes' if r.correct else 'NO':7} {r.bricks:>8} {r.symbols:>23,}")


def show_sine(x):
    print(f"\n=== The same sine on both sides: sin({x}) ===\n")
    for r in B.compare_sine(x):
        s = f"{r.symbols:.1e}" if r.symbols > 10 ** 15 else f"{r.symbols:,}"
        print(f"{r.world:6} {_pretty(r.value):>20}  {'correct' if r.correct else 'WRONG':7}  error {mp.nstr(abs(r.value - r.expected), 2):>8}"
              f"  bricks {r.bricks:>6}  symbols without sharing {s}")
    print("  (NAND: 20 bits of which 13 after the point; its error is counted relative to the angle rounded to its bits)")


def show_recipes():
    print("\n=== EML recipes: the short one (with its wall) and the wall-free one ===\n")
    for name, n, zone, n_safe, safe_zone in E.recipe_table():
        print(f"  {name:15} short: {n:>5} symbols — {zone}")
        if n_safe:
            print(f"  {'':15} safe : {n_safe:>5} symbols — {safe_zone}")
    print("\n  In a hole (0), the tool refuses (Undefined error) instead of returning a wrong result.")


def show_curve():
    print("\n=== NAND side: what the precision of the sine costs (folded angle, whole range) ===\n")
    print(f"  {'bits':>5} {'after the point':>17} {'gates':>8} {'worst error':>12} {'in units of the last bit':>26}")
    for width, fraction, count, worst in S.precision_curve():
        print(f"  {width:>5} {fraction:>17} {count:>8} {mp.nstr(worst, 2):>12} {mp.nstr(worst * 2 ** fraction, 2):>26}")


def show_brick():
    from . import nand_eml as M
    print("\n=== The eml brick itself, built in NAND gates ===\n")
    for x, y in ((1, 1), (mp.mpc(1, 2), mp.mpc(3, 4)), (2, -3)):
        v, _, count = M.compute(x, y, integer_bits=10)
        true = mp.exp(x) - mp.log(y)
        print(f"  eml({_pretty(mp.mpmathify(x), 4)}, {_pretty(mp.mpmathify(y), 4)}) = {_pretty(v, 8):>26}   true {_pretty(true, 8)}   error {mp.nstr(abs(v - true), 2)}")
    print(f"  → one eml brick = {count:,} NAND gates in 26 bits (10 before the point, 16 after); "
          f"{M._circuit(8, 16, 4)[3]:,} in 24 bits (8 + 16)")
    print("\n=== Closing the loop: EML trees computed only with this NAND brick ===\n")
    for name, t, env, true in (("ln(2)", E.ln_("x"), {"x": 2}, mp.log(2)), ("5 + 3", E.add("x", "y"), {"x": 5, "y": 3}, 8),
                               ("13 × 11", E.mul("x", "y"), {"x": 13, "y": 11}, 143), ("√2", E.sqrt_("x"), {"x": 2}, mp.sqrt(2)),
                               ("mirror(3+4i)", E.mirror("x"), {"x": mp.mpc(3, 4)}, mp.mpc(3, -4)),
                               ("sin(1)", E.sine("x"), {"x": 1}, mp.sin(1))):
        v, bricks, per = M.evaluate_with_nand(t, env)
        print(f"  {name:13} = {_pretty(v, 8):>24}   error {mp.nstr(abs(v - true), 2):>8}   {bricks:>3} bricks × {per:,} = {bricks * per:>10,} gates")
    print("  For comparison, in direct NAND: addition 68 gates, multiplication 600, sine 12,582.")
    print("  Limit: no saturation — a tree in which an intermediate number exceeds ±512 overflows and returns a wrong result.")
    print("  (Same 26-bit brick everywhere; the constants are precomputed; −∞ = the smallest writable number.)")


def demonstration():
    show_computation(13, 11)
    for x in (1, 3, -50):
        show_sine(x)
    print("\n=== Before / after: the walls lifted ===\n")
    v, *_ = S.compute(3, 16, 13, full=False)
    print(f"  NAND bare     sin(3)   = {_pretty(v, 6):>10}   WRONG (true: {_pretty(mp.sin(3), 6)})  → with folding: {_pretty(S.compute(3)[0], 6)}")
    for name, short, safe, z, ref in [("ln(−2)", E.ln_, E.ln_exact, mp.mpf(-2), mp.log(-2)),
                                      ("√(−4)", E.sqrt_, E.sqrt_safe, mp.mpf(-4), mp.sqrt(-4)),
                                      ("sin(1+3i)", E.sine, E.sine_safe, mp.mpc(1, 3), mp.sin(mp.mpc(1, 3)))]:
        print(f"  EML short     {name:9}= {_pretty(E.evaluate(short('x'), {'x': z}), 6):>22}   WRONG"
              f"  → wall-free: {_pretty(E.evaluate(safe('x'), {'x': z}), 6)}   (true: {_pretty(ref, 6)})")
    a, b = mp.mpc(1, 0), mp.mpc(0, 4)
    print(f"  EML short     1 + 4i   = {_pretty(E.evaluate(E.add('x', 'y'), {'x': a, 'y': b}), 6):>22}   WRONG"
          f"  → wall-free: {_pretty(E.evaluate(E.add_safe('x', 'y'), {'x': a, 'y': b}), 6)}")
    print("\n=== The eml★ mirror (9 symbols) and what it makes possible ===\n")
    for z in (mp.mpc(3, 4), mp.mpc(-2, 40), mp.mpf(-5), 0):
        print(f"  mirror({_pretty(z, 6)}) = {_pretty(E.evaluate(E.mirror('x'), {'x': z}), 6)}")
    print(f"  |3 + 4i| = {_pretty(E.evaluate(E.modulus('x'), {'x': mp.mpc(3, 4)}), 12)}   ({E.size(E.modulus('x'))} symbols)")
    print("\n=== Constants built with eml and 1 only ===\n")
    for name, t in (("i", E.I), ("π", E.PI)):
        print(f"  {name} = {_pretty(E.evaluate(t, {}))}   ({E.size(t)} symbols, exact)")
    show_recipes()
    show_curve()
    show_brick()
    print()


def main(args):
    try:
        if not args:
            demonstration()
        elif args[0] == "sine":
            show_sine(mp.mpf(args[1]))
        elif args[0] == "recipes":
            show_recipes()
        elif args[0] == "curve":
            show_curve()
        elif args[0] == "brick":
            show_brick()
        elif len(args) == 2 and all(0 <= int(a) <= 255 for a in args):
            show_computation(int(args[0]), int(args[1]))
        else:
            print(__doc__)
    except (ValueError, IndexError) as e:
        print(f"I did not understand or could not compute: {e}\n{__doc__}")


if __name__ == "__main__":
    main(sys.argv[1:])
