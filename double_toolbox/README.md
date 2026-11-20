# Double toolbox NAND / EML

Two worlds, a single brick each, queried at the same time.

- **NAND**: one logic gate. Bits, integers, fixed point. (`nand.py`, `nand_sine.py`)
- **EML**: `eml(x, y) = exp(x) − ln(y)` and the number 1; plus `eml★` for the mirror. Complex numbers. (`eml.py`)
- `toolbox.py`: `compare`, `compare_sine`, `measure`.
- `star_only.py`: the constructions that use eml★ as their only gate (mirror in 297 and 101 symbols, product,
  strict real part, NAND on bits). See the report, § 17.

This directory was written with AI assistance (Claude) on 20 November 2026 and translated from a French working
copy (a local directory named `boite_double`, never published). It has not been reviewed by a human.

## Commands (from the repository root; needs mpmath, sympy, pytest)

    python3 -m double_toolbox            # full demonstration
    python3 -m double_toolbox 200 55     # + and × of two integers (0 to 255)
    python3 -m double_toolbox sine 2.5   # sine of an angle (−64 to 64)
    python3 -m double_toolbox recipes    # EML recipes: short / wall-free
    python3 -m double_toolbox curve      # price of precision on the NAND side
    python3 -m double_toolbox brick      # the eml brick in NAND gates, and EML trees computed with it
    python3 -m double_toolbox.certificates   # interval certificates, about 20 minutes
    python3 -m pytest -q double_toolbox/test_toolbox.py   # 37 tests

Recorded outputs: `CONSOLE_20261120.txt` (tests and demonstration) and `CONSOLE_certificates_20261120.txt`
(certificates), both from runs of this English version.

## Status of what is claimed

- Certified elsewhere (report FINDINGS_20261119, § 1–16): zones of `add` and `mul`, 9-symbol mirror, minimality.
- Exact (SymPy): the constants 0, −∞, −1, 2, ½, i, π.
- **Interval certificates** (`certificates.py`): 11 expected verdicts out of 11. An expected verdict can be a
  refusal: three of the eleven are controls on the SHORT recipes, which must be found wrong, must fail, or must
  be unmasked by a witness. On the stated
  domains only, they establish the absence of branch jump for `ln_exact`, `add_safe`, `real_part_safe`, `imag_part`
  and, on a small domain, `add_total`. They do not by themselves establish equality with the target function: that
  comes from the witness points and from a derivation from the rules already established, not yet written up nor
  reviewed. Where intervals cannot decide (on the branch cut, at z = 1, at zero, on the axes), a few points are
  judged in exact arithmetic.
- **Measured only** (grids of `test_toolbox.py`): `sub_safe`, `sqrt_safe`, `sine_safe`, `cosine_safe`,
  `sub_total`, the eml★-only constructions other than the mirror, and the NAND sine.
- Known holes: `add_safe` and `sub_safe` refuse (`Undefined` error) when x = 0 or y = 0. `add_total` handles
  zero only in extended evaluation.

## The eml brick in NAND gates (`nand_eml.py`)

`brick` builds eml(x, y) = exp(x) − Log(y) in NAND: real exponential by shifts and additions, CORDIC
rotation for the cosine and the sine, inverse CORDIC for the modulus and the angle of y, logarithm by shifts and
additions. 72,513 gates in 24 bits, 77,633 in 26 bits. `evaluate_with_nand` computes a whole EML tree with this
single brick. Limits: no saturation, no +∞, y = 0 signalled by a wire, errors that accumulate from brick to brick.

## The certificates (`certify.py`, `certificates.py`)

They rely on the certifier of the repository (`certify_branch_safety.py`), on the absorption lemma, and on the
exact values of the constants by SymPy.

## Limit of the comparison

An `eml` brick already contains an exponential and a logarithm: measured, it is worth about 72,000 to 78,000
NAND gates. The comparison in gates is therefore possible (command `brick`), at fixed precision: 13 × 11 through
14 bricks costs 1,086,862 gates and returns 142.99, against 600 gates and an exact result for a direct NAND
multiplier.
