# eml★ — Completing the Sheffer Basis for Continuous Functions on ℂ

**Author:** Anthony Monnerot (independent researcher, Champigny-sur-Marne, France)  
**Date:** April 30, 2026  
**Reference:** Extension of the EML operator by Andrzej Odrzywołek (arXiv:2603.21852v2)

---

## What is this?

Odrzywołek (2026) showed that a single operator

```
eml(x, y) = exp(x) − ln(y)
```

together with the constant `1`, generates all standard elementary functions
(sin, cos, sqrt, exp, ln, …) by finite composition — the "NAND gate of
continuous mathematics."

**This repository identifies a structural limitation and proposes a minimal fix.**

The operator `eml` is holomorphic by construction.  Therefore:

- complex conjugation  `z̄`
- real part  `Re(z)`
- imaginary part  `Im(z)`
- modulus  `|z|`

are **not reachable** by any finite `eml`-composition (Theorem 2.1 + Corollary 2.2).

---

## The fix: eml★

We introduce the companion operator:

```
eml★(x, y) = exp(x) − ln( ȳ )
```

where `ȳ` is the complex conjugate of `y`.  
This is the minimal modification that introduces anti-holomorphicity.

**Key formula (Theorem 3.1):**

```
z̄  =  1  −  eml★( 0,  eml(z, 1) )       for Im(z) ∈ [−π, π)
```

Complex conjugation at **depth 2**, exact.

---

## Results (all unconditional)

| Result | Depth | Status |
|--------|-------|--------|
| `z̄ = 1 − eml★(0, eml(z,1))` | 2 | Conditional theorem (Im(z) ∈ [−π, π)) |
| `Re(z) = (z + z̄)/2` | 3 | Unconditional corollary |
| `\|z\|² = z · z̄` | 3 | Unconditional corollary |
| `\|z\| = √(z · z̄)` | 4 | Unconditional corollary |
| `{eml, eml★, 1}` dense in `C(K, ℂ)` | — | Stone–Weierstrass (K ⊂ strip [−π, π)) |

**Branch limitation (Theorem 3.2):**  
The conjugation formula holds exactly for `Im(z) ∈ [−π, π)` (half-open strip).  
The upper boundary `+π` fails; the lower boundary `−π` holds: `exp(z)` lands on the branch cut of the principal
logarithm. Interior Im(z) ∈ (−π, π) verified numerically at 50 decimal digits; boundary Im(z) = −π verified symbolically (exact SymPy proof). See `verify_theorem4.py`.

---

## Repository structure

```
eml_star/
├── README.md                      ← this file
├── paper/
│   └── eml_star_final.pdf         ← full paper (4 pages)
├── eml_toolkit/
│   ├── __init__.py
│   ├── core.py                    ← eml, eml★, conjugate_formula
│   ├── optimizer.py               ← EGraphOptimizer (6 rewrite rules)
│   └── conjecture_explorer.py     ← parallel random search
├── verify_theorem4.py             ← numerical verification (mpmath, 50 digits)
└── test_suite.py                  ← full test suite (Theorems 3.1, 3.2, Cor. 3.3)
```

---

## Quick start

```bash
pip install mpmath sympy
python verify_theorem4.py    # verify Theorem 3.1 numerically
python test_suite.py         # full test suite
```

Expected output: all tests pass, error < 10⁻⁴⁰.

---

## Paper

The full mathematical paper is in `paper/eml_star_final.pdf`.  
It is an arXiv comment on Odrzywołek's paper (2603.21852v2).

Sections:
1. Introduction
2. The holomorphic barrier (Theorem 2.1, Corollary 2.2)
3. The companion operator eml★ (Theorem 3.1, 3.2, Example 3.4)
4. Topological completeness via Stone–Weierstrass (Theorem 4.3)
5. Summary

---

## Requirements

```
mpmath >= 1.3.0
sympy  >= 1.12
python >= 3.10
```

---

## License

MIT — free to use, cite, and extend.

## Citation

```
Monnerot, A. (2026). The eml★ operator: completing the Sheffer basis
for continuous functions on ℂ. arXiv comment on 2603.21852v2.
https://github.com/[your-username]/eml_star
```

---

*This work was developed independently in April 2026,
as a direct extension of Odrzywołek (arXiv:2603.21852v2).*

---


## Caveat: intermediate branch safety

The constructions of `Re(z)`, `Im(z)`, and `|z|` rely on eml arithmetic
(addition, multiplication) from Odrzywołek [1], whose trees reach depths
of ~20–30 nodes. **It is not proven in this note** that the intermediate
`ln` arguments in those trees remain within the strip `(−π, π]` when the
inputs `z` and `z̄` are in `[−π, π)`.

If any intermediate value exits the strip, the evaluation of the full tree
may suffer a 2πi branch-cut jump even though the inputs are valid.

The clean arithmetic constructions for subtraction, negation, and small
rational constants are given explicitly in `EML_toolkit/EmL_compiler/EmL_clean_math.nb`.
These versions avoid infinity issues and improve numerical stability in the
strip [−π, π). The full addition and multiplication trees remain those of [1]
(Supplementary Information + EML compiler, available at
https://github.com/VA00/SymbolicRegressionPackage). Exhaustive verification
of intermediate nodes on arbitrary compacts is left for future work.

## Known limitations

- The conjugation formula `z̄ = 1 − eml★(0, eml(z, 1))` holds only for
  `Im(z) ∈ [−π, π)` (half-open strip [−π, π)). The upper boundary `+π` fails; the lower boundary `−π` holds due to the
  branch cut of the principal logarithm.

- Stone–Weierstrass density (Theorem 4.3) is therefore restricted to compacts
  contained in this strip. Density on arbitrary compacts in ℂ is an **open problem**:
  it would require tracking winding numbers continuously, which is impossible
  with finite compositions of analytic/anti-analytic operators (monodromy theorem).

- This does not affect the core results (Theorems 2.1, 3.1, Corollary 2.2),
  which remain valid and unconditional within their stated domains.
