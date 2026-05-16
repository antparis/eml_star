# eml★ — Minimal Anti-Holomorphic Extension of the EML Sheffer Operator

**Author:** Anthony Monnerot  
**Date:** May 2, 2026  
**Extension of:** Odrzywołek (arXiv:2603.21852v2)

---

## Background

Odrzywołek (2026) showed that a single binary operator  
`eml(x, y) = exp(x) − ln(y)` together with the constant 1  
generates all classical elementary functions (continuous Sheffer basis).

**Structural limitation identified** (Theorem 2.1 / Corollary 2.2):  
eml is holomorphic by construction → z̄, Re(z), Im(z), |z|² are **unreachable**
by any finite eml tree. This result is unconditional.

---

## Contribution: eml★

```
eml★(x, y) = exp(x) − ln(ȳ)
```

A single conjugate on the second argument is sufficient to reach:

| Result | Depth | Status |
|---|---|---|
| **Cor. 2.2**: Re, Im, conj ∉ {eml, 1} | — | **UNCONDITIONAL** |
| **Th. 3.1**: z̄ = 1 − eml★(0, eml(z,1)) | relative 2 | Conditional: Im(z) ∈ [−π, π) |
| **Th. 3.2**: exact domain [−π, π) | — | Conditional; boundary −π proved via SymPy |
| **Th. 4.3**: {eml, eml★, 1} dense in C(K, ℂ) | — | Stone–Weierstrass, conditional on branch-safety lemma; numerically confirmed at 60 dps (trees RPN-length 19 for +, RPN-length 17 for ×, Zenodo DOI:10.5281/zenodo.19183008) |

---

## Branch-Safety Verification (May 2026)

Real EML trees for + (RPN-length 19) and × (RPN-length 17) extracted from Odrzywołek's official
archive (`eml_compiler_v4.py`, Zenodo DOI:10.5281/zenodo.19183008) and tested
with mpmath at 60 decimal places:

> **Note on notation:** K denotes RPN (Reverse Polish Notation) program length,
> as defined in Odrzywołek (2026) Table 4. Tree depth is approximately K/2.

- ℝ⁺: ✅ branch-safe  
- ℂ (small Im): ✅ branch-safe  
- 1000 random points in Im ∈ (−1, 1): ✅ zero violations  

Reproducible script: `branch_safety_final.py`

---

## Choice by Use Case

- **Uniform symbolic regression (Sheffer)** → **{eml, eml★, 1}**  
  Homogeneous binary grammar `S → 1 | eml(S,S) | eml★(S,S)`.  
  Compatible with Odrzywołek's Rust/Mathematica engine without modification.

- **Physical observables without strip restriction** → **{eml, Re, 1}**  
  z̄ = 2 Re(z) − z, |z|² exact everywhere on ℂ.  
  Note: Re is a unary operator — breaks the uniform binary grammar.

Both bases are valid. The choice depends on the intended application.

---

## Honest Caveats

- Im(z) ∈ [−π, π) is a condition on the **value** of z, not on its argument.
  For normalized plane waves (|ψ| = 1), it is automatically satisfied.
- Th. 4.3 is proved via Stone–Weierstrass conditionally on a branch-safety
  lemma. The lemma is verified numerically at 60 decimal places on tested
  domains (grids and 1000 random points). A formal proof via interval
  analysis on Im(·) of every intermediate node remains an open problem.
- The system is limited to elementary functions.

---

## Project Files

| File | Role |
|---|---|
| `eml_star_final.pdf` | Full paper (4 pages + Appendix A) |
| `branch_safety_final.py` | Branch-safety test on real Zenodo trees (RPN-length 19 for +, 17 for ×) |
| `emlstar_engine.py` | Complete eml★ verification engine |
| `verify_theorem4.py` | Numerical verification of Theorem 4.3 |
| `eml_toolkit/core.py` | eml, eml★, conjugate_formula operators |
| `eml_toolkit/optimizer.py` | EGraphOptimizer (6 rewriting rules) |
| `eml_toolkit/eml_quantum.py` | Phase, geometric phase, self-testing |
| `eml_toolkit/examples_constants.py` | π/3, π/4, √2 via eml★ |
| `eml_toolkit/test_suite.py` | Complete test suite |
| `docs/PROJECT_SUMMARY.md` | Full project summary and theorem statuses |
| `reference/2603_21852v2_Odrzywolek.pdf` | Source paper (Odrzywołek 2026) |

---


---

## Application: Galaxy Rotation Curves (May 2026)

eml★ was applied to galaxy rotation curves via genetic programming symbolic regression, demonstrating its first empirical application beyond pure mathematics.

### Method
- GP engine (DEAP) with eml★ operators applied to 125 SPARC galaxies + 23 LITTLE THINGS dwarf galaxies
- Comparison protocol: fits with vs without anti-holomorphic operators
- Characterization battery on 6 synthetic complex functions

### Key Results

| Result | Status |
|--------|--------|
| eml★ = non-holomorphicity detector (0/10 holomorphic, 10/10 anti-holomorphic) | **ESTABLISHED** |
| Improvement correlates with low luminosity (rho = -0.27, p = 0.004) | **ESTABLISHED** |
| Signal replicates on independent dataset (LITTLE THINGS, 43.5%) | **ESTABLISHED** |
| MOND acceleration: not a predictor (p = 0.35) | **NEGATIVE** |
| Dark matter fraction: not a predictor (p = 0.81) | **NEGATIVE** |
| Vortex topology detection: ruled out (fake vortex 5/5) | **REFUTED** |

### Interpretation
The gravitational potential in dark-matter-dominated galaxies contains a non-holomorphic component (df/dz_bar != 0) that becomes detectable when the baryonic signal is subdominant. Luminosity is the best single predictor of eml★ response.

### Resources
- **Paper draft:** [`eml_star_paper_draft3.md`](https://github.com/antparis/oxieml-star/blob/master/eml_star_paper_draft3.md)
- **Rust engine:** [OxiEML-Star](https://github.com/antparis/oxieml-star) | [Zenodo DOI:10.5281/zenodo.20152989](https://zenodo.org/records/20152989)
- **Data:** 125 SPARC + 23 LITTLE THINGS rotation curves (public)

## License

MIT

---

## Citation

```
Monnerot, A. (2026). eml★: Minimal Anti-Holomorphic Extension of the EML
Sheffer Operator. Extension of Odrzywołek (arXiv:2603.21852v2).
GitHub: https://github.com/antparis/eml_star
Zenodo: https://doi.org/10.5281/zenodo.20102448
```


## EML-WM — EML World Model (Prototype)

`jepa_real_baseline.py` and `jepa_eml_complex_v2.py` implement a minimal prototype of **EML-WM** (EML World Model) — a complex latent predictive architecture using eml★ and eml⁰ as exact mathematical primitives in the latent predictor.

EML-WM explores whether mathematically exact complex-domain operators improve latent space stability in predictive world models. This is an original architecture — distinct from existing implementations.
