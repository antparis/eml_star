# Branch-safety certificate for the + and × trees — 2026-11-19

*Translated from French on 2026-11-20. The French original is in the repository history at revision 365d88d. The scripts were translated at the same time (comments, messages, identifiers, file names) and re-run: the console files of this version come from the English scripts and carry the same numbers as the French runs, wall-clock timings aside (timings quoted in this report are those of the original runs). Replays attributed to Anthony Monnerot below were made on the French-language versions, at or before revision c074e6a.*

Files: `certify_branch_safety.py` (certificate), `test_certify_branch_safety.py` (26 tests).
Nothing was modified in the paper, the README or the existing scripts; nothing was pushed to GitHub.

## Result in one sentence

The published assumption ("no node crosses the branch cut") is **false** for the × tree, but the property that
Theorem 4.3 actually needs (the + and × trees compute **continuous** functions on K) is now
**computer-certified** on K = [−3, 3] + i[−3, 3] minus the small square |Re|, |Im| < 0.1 around 0, for x and y,
given a three-line lemma proved by hand (Lemma A below).

## 1. What was wrong

- In `E_mul`, the subtree `E_ln(E_ln(x))` takes the logarithm of ln x. For real x in (0, 1), ln x < 0:
  this log is taken **on the branch cut**. Measured: E_ln(ln(0.7 ± 10⁻¹² i)) = −1.031 ± 3.1416 i (a jump of 2πi).
  At x = 1, it takes ln 0.
- The final result remains correct (E_mul(0.7 ± 10⁻¹² i, 1.3) = 0.91 ± 1.3·10⁻¹² i): the jump is erased higher up by
  an exponential.
- The Assumption of `eml_star_branch_safety_lemma_proof.tex` ("Im of every intermediate argument in (−π, π]") can
  remain true (ln(ln x) has Im = +π just above and ≈ −π just below): it is not the right condition.
  The sentence "every composite tree … evaluates without branch-cut crossings on K" is false.
- The README test ("1000 random points, zero violations") checks that the final value is correct, not the absence
  of crossings (mpmath computes on the branch cut without warning).

## 2. The right rule: the two gates

`eml(x, y) = exp(x) − ln(y)`. A log that crosses the negative real axis makes the value of its node jump by 2πi·k.

- Node at the **x gate** (first argument of the parent): the parent computes exp(node); exp(u + 2πik) = exp(u): erased.
- Node at the **y gate** (its log is taken in turn) or at the **root**: the jump remains: rejected.
- Log of a value that may equal **0**: rejected (except constants, §3).

**Why a successful certificate proves continuity.** On each accepted rectangle, the program computes with a
branch of the log that is continuous on that rectangle (principal if the rectangle avoids (−∞, 0]; rotated, ln(−y) + iπ, if the
rectangle straddles the negative axis and is far from 0). At the nodes at the x gate, the value computed in this way differs from the
"principal" value by 2πik, which the parent's exponential erases; elsewhere it is equal. The function at the
root therefore coincides, on each closed rectangle, with a continuous function. A finite number of closed rectangles
covers K: the function is continuous on K (pasting lemma). Correctness relies on the outward rounding
of mpmath.iv (50 digits).

## 3. The constants

Subtrees without a variable are computed exactly over the **extended reals** [−∞, +∞] (not the Riemann sphere:
exp(−∞) = 0 must be distinguished from exp(+∞) = +∞): ln 0 = −∞, exp(−∞) = 0, exp(+∞) = +∞, ln(+∞) = +∞; ∞ − ∞ and ln(−∞)
rejected. E_zero = e − ln(e^e) = 0 is decided exactly by SymPy. Chain of E_neg: eml(1, 0) = +∞,
eml(+∞, 1) = +∞, eml(1, +∞) = −∞, then exp(−∞) = 0.

## 4. Lemma A (by hand) and the point x = 1

**Lemma A.** For every w ∈ ℂ \ {0}, exp(E_ln(w)) = w.
*Proof.* E_ln(w) = e − Log(exp(e − Log w)). For every u ∈ ℂ, Log(exp u) = u − 2πi·k(u) with k(u) ∈ ℤ. Hence
E_ln(w) = Log w + 2πi·k and exp(E_ln(w)) = exp(Log w) = w. For w = 0, with the conventions of §3:
eml(1, 0) = +∞, eml(+∞, 1) = +∞, eml(1, +∞) = −∞, exp(−∞) = 0 = w. ∎

**Consequence.** eml(E_ln(a), c) = a − Log c for every a. In E_mul, the only passage through ln(ln x) is
eml(E_ln(E_ln(x)), c) = E_ln(x) − Log c: this gives the tree `E_mul_lemma`, equal to E_mul at every point (also checked
at points, down to 10⁻³⁰, including x = 1 + 10⁻⁶ + 10⁻⁷ i). It is `E_mul_lemma` that the program certifies: x = 1 no longer
needs to be excluded, and the jump of ln x for negative real x is erased by the exponential at the root.

## 5. Measurements (`python3 certify_branch_safety.py` and tests)

| Tree | Domain (for x and for y) | Result |
|---|---|---|
| E_add | [0.1, 3] + i[−1, 1] (README domain) | certified |
| E_mul (original tree) | x ∈ [0.9, 1.1] + i[−0.1, 0.1] | failure (ln(ln x)), expected |
| E_mul (original tree) | x ∈ [0.5, 1.5] + i[−0.3, 0.3] minus the disc \|x − 1\| < 0.1 | certified (absorbed crossings) |
| E_mul (original tree) | [0.1, 3] + i[−1, 1] minus the disc \|x − 1\| < 0.05 | certified: 75592 rectangles, 22784 absorbed crossings, 163 s |
| E_mul_lemma | [0.1, 3] + i[−1, 1] (disc included) | certified in one rectangle |
| E_add | [−3, 3] + i[−1, 1] minus \|Re\|, \|Im\| < 0.1 | certified: 16 rectangles, 8 absorbed crossings |
| E_mul_lemma | same | certified: 16 rectangles, 24 absorbed crossings |
| E_add | [−3, 3] + i[−3, 3] minus \|Re\|, \|Im\| < 0.1 | certified: 2056 rectangles, 584 absorbed crossings, 2 s |
| E_mul_lemma | same | certified: 16 rectangles, 24 absorbed crossings |

The tests also check that the program **rejects** what is genuinely discontinuous (a crossing at the y gate, at the
root, the log of a rectangle containing 0, a jumping value used under a log).

## 6. Proposed statement to replace the Assumption (for Anthony to decide, nothing is modified)

> **Proposition (certified continuity of the arithmetic trees).** Let K = [−3, 3] + i[−3, 3] minus
> {|Re z| < 0.1, |Im z| < 0.1}. Odrzywołek's trees E_add and E_mul, evaluated with the principal logarithm and the
> extended-real conventions for their constant subtrees, define continuous functions on K × K, equal
> to x + y and x·y. Some intermediate nodes cross the branch cut (−∞, 0] (for example ln(ln x) for x ∈ (0, 1)), but
> each crossing is followed by an exponential that erases its 2πi jump (two-gate rule). Proof: Lemma A
> and the interval certificate `certify_branch_safety.py`.

Also to be corrected: the sentence "evaluates without branch-cut crossings on K", the row of the README table
("numerically confirmed …"), and the "Branch-Safety Verification" section of the README.

## 7. Limitations

- K excludes a neighbourhood of 0: ln 0 is reached there by the variable, the multiplication goes through ln x.
- |Im| ≤ 3 < π: beyond that, ln(exp(y)) ≠ y and E_neg changes value; not tried.
- Measured sizes of the trees of `branch_safety_final.py` (leaves + operators): E_add 27, E_mul 41. The README
  announces "RPN-length 19 for +, 17 for ×" and the script itself "K=13 and K=20": these numbers do not refer to the
  same trees. To be checked against `eml_compiler_v4.py` (Zenodo) before quoting any lengths.
- Raw output of the computations: `CONSOLE_20261119_branch_safety_certificate.txt`.

## 8. Orthogonal axis (same day): what all our analyses kept fixed

`orthogonal_axis.py` is the REPRODUCTION script for this paragraph (put back into a file that same evening: the original trace
came from the same lines typed without a file). It prints its results without checking any verdict: its exit code
does not signal a failure. The command that checks the six verdicts and exits with an error on failure remains
`python certify_branch_safety.py`. Outputs: `CONSOLE_20261119_orthogonal_axis.txt` (original, with the duration) and
`CONSOLE_20261119_orthogonal_axis_replay.txt` (replay; same results printed). Scope clarified after review: what
the script OBSERVES and what is PROVED are kept separate.

- **Real part and 2πi jump.**
  *Proved (one line, without the script)*: if the two one-sided limits of the value of a node differ by 2πi·k,
  k an integer, their real parts are equal, because 2πi·k is purely imaginary; likewise for u + conj(u) = 2 Re u. The
  imaginary parts differ by 2πk; conjugation changes the gap into −2πi·k without erasing it.
  *Observed by the script, at ONE single pair of points* (x = 0.7 ± 10⁻¹² i, node ln(ln x)): difference of the two values
  close to 2πi (printed `6.2831853j`), difference of the real parts printed 0, of u + ū printed 0, of the imaginary
  parts 6.2831853. The exact jump of 2πi is that of the one-sided limits; two points at distance 10⁻¹² give
  only an approximate value of it. It is an illustration, not a verification over a domain. Remark: Re is not
  an operator of the grammar {eml, eml★, 1}; this "gate" concerns the alternative basis {eml, Re, 1} or a
  real part taken outside the tree.
- **Addition down to zero: only for the FORM REWRITTEN by the lemma.** Write P = [−3, 3] + i[−3, 3] and K = P minus
  the small open square |Re|, |Im| < 0.1.
  Case 6 of the certificates: the rewritten form x − Log(exp(E_neg(y))) is certified on P × P, zero included (256 rectangles).
  Case 4: the original tree E_add is certified on K × K (the small square is removed from each of the two domains, that
  of x and that of y).
  The link between the two forms is the absorption lemma, eml(E_ln(x), c) = x − Log c, including at x = 0 with
  the extended evaluation: it rests on the proof, not on the program.
- **Sample of points** (independent of the interval logic, same mpmath library): 6000 pairs (x, y),
  three comparisons per pair (E_mul and its rewritten form against x·y, E_add against x + y), i.e. 18000 absolute differences,
  with the ordinary principal log. x is placed right against the segment (0, 1) or the negative axis (down to 10⁻³⁰), or taken in
  [−3, 3] + i([−3, −0.1] ∪ [0.1, 3]); y satisfies |Re y| ∈ [0.1, 3] and Im y ∈ [−1, 1]. Maximum printed on this
  sample: 1.09·10⁻⁴⁹. A sample proves nothing outside the points drawn.

Leads noted at that date, with what became of them:
- Multiplication at zero ("infinity gate"). The product admits a continuous extension with value zero. But
  evaluating the raw tree at x = 0 would require an additional convention or a justified rewriting: E_ln(0) = −∞,
  then the subtree E_ln(E_ln(0)) requires log(−∞), undefined under the convention adopted; widening the bounds of the
  intervals does not resolve this operation.
- Branch invariance: proved for arbitrary simultaneous shifts of the eight logarithms with variable
  argument in the × tree, the constant logarithms keeping their convention (§12, §13). Some constant logarithms
  can also be shifted individually without effect (No. 6, §15); No. 19 cannot.
- Verification by an independent tool: done in §9, in Julia.

## 9. Second opinion: the same certificate rewritten in Julia (`julia_check/certify.jl`, output `julia_check/julia_output.txt`)

Another language (Julia 1.12), another interval library (IntervalArithmetic.jl 1.0.14 on MPFR, 200 bits), complex exp and log
redone by hand from the real functions, exact constants through a small home-made symbolic computation (without
SymPy). What remains common to the two versions: the rule (two gates, Lemma A) and its author.

| Certificate | Python (mpmath.iv) | Julia (IntervalArithmetic) |
|---|---|---|
| original ×, [0.1, 3] + i[−1, 1] outside \|x − 1\| < 0.05 | 75592 rect., 22784 absorbed, 98304 excluded | **identical**: 75592, 22784, 98304 (43 s) |
| + on [−3, 3] + i[−3, 3] minus the square around 0 | 2056 rect., 584 absorbed | **identical**: 2056, 584 |
| × lemma form, same domain | certified, 16 rect. | certified, 256 rect. |
| × lemma form, [0.1, 3] + i[−1, 1] | certified, 1 rect. | certified, 16 rect. |
| + lemma form, full square, zero included | certified, 256 rect. | certified, 256 rect. |

The two rows where the number of rectangles differs: same verdict, but the complex exp redone by hand encloses slightly more
widely than that of mpmath, so Julia splits one more time. The rejection checks (y gate, root, log of zero, × near
x = 1 without a disc) also pass in Julia.

## 10. The wall and infinity: two exact formulas (hand proof + branch mutant)

Script `wall_and_infinity.py`, output `CONSOLE_20261119_wall_and_infinity.txt`.

Only three facts are used: (i) exp(log w) = w for every branch of the log, w ≠ 0; (ii) log(exp u) = u + 2πi·(integer)
for every branch; (iii) exp(u + 2πik) = exp(u).

- E_ln(w) = e − log(exp(e − log w)) = log w + 2πi·k (w ≠ 0).
- E_neg(b) = 0 − log(exp(b)) = −b + 2πi·m (the leading constant equals −∞, exp(−∞) = 0). No condition on b.
- E_sub(a, b) = exp(E_ln(a)) − log(exp(b)) = a − b + 2πi·n (Lemma A; a = 0 included, through the extended reals).
- Hence **E_add(a, b) = a + b + 2πi·N**, N an integer, and E_mul(x, y) = exp(E_add(E_ln x, E_ln y)) = exp(log x + log y + 2πi·N').

**Theorem B (multiplication).** For all x, y ∈ ℂ \ {0}, and whatever branch is chosen for each of the
logarithms of the tree, E_mul(x, y) = x·y. (At x = 1 with the principal log, E_ln(x) = 0 and the tree passes through the
extended reals: Lemma A.) Consequences: no condition on the imaginary parts; the wall is a pure convention for ×;
the function extends continuously by 0 at x = 0 and at y = 0 (the "hole" around 0 is a removable singularity:
at x = 0 exactly, the original tree would require ln(−∞), which the extended reals refuse).

**Theorem C (addition).** With the principal log, for every x ∈ ℂ and every y ∈ ℂ, E_add(x, y) = x + wrap(y), where
wrap(y) is y brought back into the strip −π < Im ≤ π by an integer number of turns 2πi (boundary Im y = ±π to be made precise).
In particular E_add(x, y) = x + y if and only if y is in the strip |Im y| < π: the condition bears on y alone,
never on x. With another branch, the strip moves: addition, for its part, does depend on the wall.

Measurements (80 digits):
- Branch mutant (each log receives 2πi·k, k drawn at random between −3 and 3), 4,000 points, |Im| up to 12:
  ×: worst relative deviation from x·y = 5.6·10⁻⁷⁹; +: (result − (x + y)) / 2πi always an integer (deviation 1.9·10⁻⁷⁹), integers from −5 to 5.
- Principal log, 4,000 points, |Im| up to 40 (3,722 outside the strip): ×: 1.6·10⁻⁸⁰; + against x + wrap(y): 4.1·10⁻⁷⁹.
  Example: E_add(1, 5i) = 1 − 1.2832 i = 1 + (5 − 2π) i.
- Towards zero: |x| = 10⁻⁶⁰: E_mul(x, y) = x·y to within 3·10⁻⁷⁹ (relative), in both positions.

What this changes for Theorem 4.3: the continuity of × on every compact subset of (ℂ \ {0})² and of + on every compact subset of
ℂ × {|Im y| < π} follows from Theorems B and C (a proof of a few lines, to be submitted for review); the interval
certificates (Python and Julia) are its machine confirmation on the domains of §5, §8 and §9. For Stone–Weierstrass,
the functions handled will have to remain in the strip |Im| < π at the moment they are added in the y position: that is the
true hypothesis, to be written in place of the old one.

## 11. The conjugation bar facing the wall (mirror formula, Theorem 3.1)

Script `mirror_and_wall.py`, output `CONSOLE_20261119_mirror_and_wall.txt` (60 digits, 3,000 points, Im up to 15).

- **Theorem 3.1 is exact** in its strip −π ≤ Im z < π (674 points, deviation 8·10⁻⁶²), boundaries included (Im z = −π does
  give 0.5 + πi). The stated domain [−π, π) is the right one, including the closed side.
- **Outside the strip, a new exact formula**: 1 − eml★(0, eml(z, 1)) = wrap(z̄), that is, z̄ brought back into the
  strip by whole turns (example z = 1 + 5i → 1 + 1.2832 i = 1 + (2π − 5) i). Proof: the formula equals
  Log(conj(exp z)) = Log(exp(z̄)).
- **What never depends on the wall**: exp(formula) = exp(z̄) everywhere, and the real part of the formula = Re z everywhere
  (deviation 8·10⁻⁶²). With the branch mutant, formula − z̄ is always an integer number of turns 2πi (−5 … 5).
- **Consequences for what is built with it**: Re z = (z + formula)/2 always has the right real part, but outside
  the strip picks up a spurious imaginary part (up to 2π here); |z|² = z · formula is correct only in the strip
  (z = 1 + 5i gives −5.42 + 6.28 i instead of 26).

**Defect in the public code: `eml_toolkit/core.py`.**
- `modulus_squared(z) = eml(z, eml★(z, 1))` equals exp(z) − z, not |z|² (z = 3 + 4i → −16.13 − 12.92 i instead of 25;
  `modulus(3 + 4i)` = 1.51 − 4.29 i instead of 5). This file is in the GitHub repository (commit a4e2202).
- `real_part` and `imag_part` of core.py call mp.re / mp.im: not eml★ trees.
- `eml_toolkit/test_suite.py` does NOT test core.py: its import fails (core.py has no `conjugate_formula`, `eml_exp`…)
  and the suite silently switches to its own fallback definitions (which are correct). Hence "5 passed" despite the defect.
- Proposed fix (not applied): modulus_squared(z) = z * E_conj(z) (or E_mul(z, E_conj(z))), with the
  condition −π ≤ Im z < π written in the docstring, and a test suite that really imports core.py (loud failure).

**Fix applied locally on 2026-11-19 (not committed, not pushed)**: `eml_toolkit/core.py` rewritten (original kept:
`core.py.bak_2026-11-19`). The old file was also wrong for E_sub (exp(a) − b instead of a − b), E_add and E_mul
(E_mul(1.3 + 0.1i, 0.7 + 0.2i) gave −1.7·10¹⁵ + 3·10¹⁶ i instead of 0.89 + 0.33i). New tests
`eml_toolkit/test_core.py` (8, direct import); `test_suite.py` now imports core.py for real (no more fallback
definitions; original kept as `.bak_2026-11-19`).

## 12. External review (ChatGPT, same day): five corrections accepted

Verified on the machine before being accepted (`allowed_branches.py`, `CONSOLE_20261119_allowed_branches.txt`).

1. **"For every branch of each logarithm" was FALSE (§10, Theorem B).** Counterexample confirmed: ln 1 = 2πi at the
   last node eml(u, 1) alone changes the result by exactly −2πi. My mutant did not see it: it exempted ln 1.
   Measurement, one logarithm shifted at a time among the 20 of the × tree: the 8 logarithms with a variable argument change
   nothing; those of the leaf 1 change the result (shift, wrong value or undefined).
   Corrected statement: principal log; the invariance holds only for the logarithms whose argument depends on x or y.
2. **x = 1**: the extended-reals convention was granted only to the constant subtrees, whereas the node in question
   contains the variable. Corrected: extended evaluation defined node by node, at each point.
3. **Boundary of addition**: E_add(x, y) = x − P(−P(y)), P = Log∘exp; the representative has Im ∈ [−π, π) (verified: deviation 0).
   E_add(x, iπ) = x − iπ; E_add(x, −iπ) = x − iπ (correct).
4. **Theorem 4.3**: stated for A_K, the unital complex subalgebra of C(K, ℂ) generated by the tree functions
   continuous on K; it contains ℂ[z, z̄]: the question of the constant i is settled for this statement. The parenthesis
   "cut the second term into N parts" is removed (it would require a tree for g/N, including at the zeros of g).
5. **Scope**: what is withdrawn is the conclusion "without branch-cut crossing", not the bound on Im (which a negative real
   does not refute); the two-gate rule is a SUFFICIENT condition (a crossing in second position is not
   absorbed by the immediate parent: that does not prove a discontinuity at the root).

The up-to-date statements are those of `eml_star_branch_safety_lemma_proof.tex` and of the README; §10 above is kept as
written, with its error, for the record.

## 13. Simultaneous shifts (ChatGPT's question) and an argument by Gemini refuted

- **Question**: testing each logarithm separately does not prove that several can be moved together.
- **Answer**: the algebraic proof covers the 8 shifts at once (induction from the leaves to the root, each value
  known "up to 2πi"); this is now written in the proof of the × proposition, with the case of zero arguments
  (a' = Log x + 2πik = 0 forces x = 1 and k = 0).
- **Gemini's argument refuted by measurement**: "each variable logarithm enters directly into an exponential".
  False for log no. 2 (the outer log of E_ln(x)): its node is in the y gate, inside E_ln(E_ln(x)). Its shift
  is not erased by the immediate parent: it is carried unchanged (absorption lemma) up to the exponential of the
  root. The other 7 are in the x gate. Consequence: the two-gate rule is indeed only sufficient.
- **Measurement** (`allowed_branches.py`): simultaneous random shifts (−3…3 turns) on the 8 logs, 3,000 points:
  worst relative deviation 4.6·10⁻⁵⁹; at x = 1 exactly with 0, +1, −2 turns everywhere: deviation ≤ 10⁻⁵⁹.

## 14. ChatGPT's second review: three adjustments accepted

1. Shifted logarithms defined at the exceptional values: L_j(w) = Log w + 2πi k_j for finite nonzero w, L_j(0) = −∞,
   L_j(+∞) = +∞. The sentence "the value of each node up to 2πi" was inaccurate (a log applied to a shifted value
   also changes its real part): the congruence bears on the outputs of the subtrees E_ln, E_neg, E_sub, E_add.
   Counting error corrected: E_ln contains ONE log of the leaf 1, not two.
2. Remark after Theorem 4.3: representing exactly every element of A_K by a single finite tree is IMPOSSIBLE
   (the trees are countably many, A_K contains all the complex constants). The uniform density of single
   trees is a distinct question, open.
3. README line aligned: "A_K is uniformly dense"; what is withdrawn is the claim "without crossing".

ChatGPT's opinion after these adjustments: in favour of the local documentation commit; the review of core.py and of the scripts remains separate.

## 15. Code review (self-review then ChatGPT): defects fixed before the second commit

None of these defects changes a verdict. Both implementations report success on the six domains tested.
Their counters agree in cases 1, 2, 4 and 6. Cases 3 and 5 reach the same verdict with different
subdivisions (case 3: 1 rectangle in Python, 16 in Julia; case 5: 16 rectangles and 24 crossings counted in Python,
256 and 400 in Julia).

1. **Mutant, log no. 6 wrongly classified "undefined"** (the script returned `None` without computing). Exact computation: shifting the
   deepest ln 1 of E_zero gives e − 2πi, then exp(e − 2πi) = exp(e), and E_zero remains 0: no effect. No. 7
   (the middle ln 1 of E_zero) gives a wrong result. Corrected classification of the 20 logarithms: 8 with a variable argument, no
   effect; constants 8, 9, 11, no effect; leaf 1: nos. 6 and 10 no effect, no. 19 shifts by −2πi, nos. 1, 4, 7, 13, 15, 17
   make the result wrong. Regressions: `test_allowed_branches.py`.
2. **Boundary of the mirror.** Exact value of the formula: Log(exp(z̄)), hence ½ + iπ for z = ½ + iπ AND for z = ½ − iπ.
   The raw tree in floating point returns ½ − iπ at the upper boundary: on the boundary, exp(z) is a hair's breadth from the negative axis and the
   log picks a side at the whim of the rounding. `core.py`: `conjugate_formula` remains the raw tree (documented as not
   reliable ON the boundary); new `conjugate_folded`, explicit folding with a guaranteed convention at the boundaries; tests
   at both boundaries and on each side.
3. **Description of the old defect**: the old `modulus_squared` equalled exp(z) − Log(exp(z)); "exp(z) − z" (§11) is
   true only in the strip (for z = 3 + 4i the two differ, the console showed it).
4. **Constants of the certificates**: the absolute margin of 1e−55 was swallowed by the rounding to 50 digits (interval of
   zero width). Replaced by a node-by-node recomputation in interval arithmetic, with refusal if a logarithm of a
   constant touches the branch cut. A counter proves that NO finite nonzero constant is involved in the six
   certificates: the zero counter means precisely that no constant took the recomputation path of the
   finite nonzero constants other than leaf 1. The published results were not affected.
5. **Decimal bounds** (0.1; 0.05): given as text and rounded outward; excluded disc tested in intervals
   (upper bound of the distance² < lower bound of the radius²).
6. **Vacuous success**: a certificate that accepts no rectangle is now a failure (Python and Julia).
7. **Reproducibility**: `python3 certify_branch_safety.py` and `julia certify.jl` each replay the six certificates of the
   README and exit with a nonzero code if one of them fails.

53 tests pass (certificate, mutant, toolkit).

## 16. A global mirror of nine symbols (night of 19 to 20 November 2026)

Origin: while comparing eml★ with a public proposal by Quadruplay (a ternary operator
T(x, y, z) = exp(x) − Log(y) + Re(z) Im(z), discussion on the VA00/SymbolicRegressionPackage repository), ChatGPT applied
the absorption lemma to the LAST node and obtained the tree below; a local exhaustive search, launched earlier,
had found another tree of the same length, correct only in the strip.

**Proposition (global mirror; distinct from Theorem 3.1).** The tree of 9 RPN symbols, leaves z and 1,

    F(z) = eml( eml( 1, eml( eml★(1, z), 1 ) ), 1 )        postfix: 1 1 z eml★ 1 eml eml 1 eml

satisfies F(z) = conj(z) for every z ∈ ℂ \ {0} (principal log), and F(0) = 0 with extended evaluation
(outputs of the four nodes: +∞, +∞, −∞, 0). Without extended evaluation, its domain is ℂ \ {0}.

*Proof.* For z ≠ 0, set w = Log(conj z) and v = exp(e − w), which is never zero. The last two nodes give
exp(e − Log v) = e^e / v = exp(w) = conj(z). Only exp(Log u) = u (u ≠ 0) is used, never Log(exp u) = u: no
strip is involved.

*"A single star" control.* The same form with eml in place of the only eml★ gives z.

Theorem 3.1 remains exact for ITS formula, 1 − eml★(0, eml(z, 1)), in its strip −π ≤ Im z < π; outside the strip
it returns conj(z) reduced modulo 2πi (§11). It is a different formula that realises conjugation everywhere.

**Statuses.**

| Result | Status | Supporting files |
|---|---|---|
| F(z) = conj(z) for finite nonzero z; star-free control = z; at i: −i and +i; at 0: trace +∞, +∞, −∞, 0 | ESTABLISHED: short proof + SymPy 1.14.0 exact judge, without tolerance (4 PASS, replayed by Anthony) | `exact_mirror_judge.py`, `CONSOLE_20261120_exact_mirror_judge.txt` |
| Measurement: 5,013 points up to \|Im\| = 200, negative axis and boundaries ±π included, relative deviation 7·10⁻⁵¹ | OBSERVED (sample) | `global_mirror.py`, `CONSOLE_20261119_global_mirror.txt` |
| No {eml, eml★} tree with leaves {1, z} of RPN length ≤ 7 realises the mirror: 9 is minimal | CERTIFIED by machine (to be reviewed): the 714 trees (2 + 8 + 64 + 640) are enclosed at the exact rational point z0 = (37 + 61i)/100 in interval arithmetic; each has one coordinate of T(z0) − conj(z0) separated from 0 (713 by the real part, 1 by the imaginary part), 0 undecided; the constants (zero, infinity, negative real) are decided exactly by SymPy. Positive control: F remains undecided at z0, as expected. The only elimination by the imaginary part is the leaf z (predicted in review: z0 − conj(z0) = 61i/50). Checks of the code included in the script: 714 trees produced = 714 distinct = count (Catalan shapes × leaves × operators); the interval of z0 strictly contains 37/100 and 61/100 (exact comparison of the binary bounds); the decision function actually used never eliminates a deviation whose two intervals both contain 0; the exit code is 0 only if the five expected results are verified (714 trees, 0 undecided, 713 by the real part, the leaf z alone by the imaginary part, F not eliminated with the exact reason); a constant whose sign SymPy cannot decide makes the tree undecided instead of interrupting the computation. Convention: principal log, extended evaluation. The bound already holds on ℂ \ {0} (z0 ≠ 0) | `minimality_certificate.py`, `CONSOLE_20261120_minimality_certificate.txt` (and `mirror_minimality.py`, earlier floating-point version) |
| Corollary: for every compact K ⊂ ℂ, the algebra A_K is uniformly dense in C(K, ℂ) (extended evaluation if 0 ∈ K; otherwise every compact subset of ℂ \ {0}) | DERIVED, to be reviewed: the leaf gives z, F gives conj(z), continuous on K; ℂ[z, z̄]\|_K ⊂ A_K and the Stone–Weierstrass argument already reviewed applies. Does not depend on minimality | — |

**Comparison with Quadruplay's ternary operator** (same leaves z and 1, same convention). Its extraction of Re and
Im by translation is exact (verified on 2,000 points). The ternary witness examined here (built by ChatGPT, fully expanded
and verified on 313 points) builds the constant i and counts 811 symbols once fully expanded (270 nodes, 541 leaves),
against 9 for F. The minimal length in the ternary grammar remains undetermined.

**Results of the exhaustive search** (`exhaustive_search_mirror.py`, `exhaustive_search_star_only.py`,
`CONSOLE_20261120_exhaustive_search.txt`; floats, comparison at 3 points, tolerance 1e−9: HEURISTIC):
- |z|² by a single tree of 17 symbols, e(e(1,e(e(e(1,e(e(1,s(1,z)),1)),z),1)),1) (e = eml, s = eml★), measured equal to
  |z|² on 5,013 points (relative deviation 3·10⁻⁵⁰);
- Re(z), Im(z) and the constant i: no tree up to length 17;
- eml★ ALONE (without eml), leaves {1, z}: neither conj(z), nor Re, Im, −z, 1/z, i up to length 19 (4.3 million
  functions at the last level). Open question: is {eml★, 1} by itself a basis?

Still open: the uniform density of INDIVIDUAL trees (§ 4.3 of the manuscript).

## 17. eml★ alone suffices: a strict mirror of 101 symbols (20 November 2026) — not reviewed by a mathematician

This section answers the question left open at the end of § 16 ("is {eml★, 1} by itself a basis?").

**Author of the discovery: ChatGPT**, on 20 November 2026 (text `PROOFS.md`, "EML orthogonal-axis study", received
in the archive `EML_orthogonal_20261120.zip`). It first gave a mirror of 297 symbols, then, after Claude's feedback
on that length, the 101-symbol version below. Anthony Monnerot replayed the judge on his
ThinkCentre; Claude rebuilt the trees independently (without running the script received), judged them with its own
evaluators, then wrote the standalone verifier `star_only_mirror_101.py`, which Anthony replayed. No mathematician has reviewed it. No prior art has been searched for.

**Statement (sentence proposed by ChatGPT).**
"A tree of 101 symbols using eml★ alone computes complex conjugation in strict evaluation on ℂ \ {0}."

*Strict evaluation*: every node of the tree has a finite value (no passage through ln 0 = −∞). *eml★ alone*: all
the gates are eml★, the leaves are 1 and z; the tree is counted fully expanded, as in § 16.

**Construction.** Write S = eml★, S(a, b) = exp(a) − Log(conj b), principal log (Arg in (−π, π]). Each row
is a tree; the right-hand column is the proved value, not additional syntax.

| Tree | Value | Symbols |
|---|---|---|
| X(t) = S(t, 1) | exp(t) | 3 |
| C = X(X(1)) | e^e | 5 |
| O = S(1, C) | 0 | 7 |
| Q(t) = X(S(1, t)) | C / conj(t), t ≠ 0 | 5 |
| L(t) = S(1, Q(t)) | e − Log(C/t): a logarithm of t, exp(L(t)) = t, t ≠ 0 | 7 |
| d(t, r) = S(L(t), X(r)) | t − r, for t ≠ 0 and r real | — |
| n = d(S(1, X(1)), X(1)) | −1 | 17 |
| k = d(S(O, n), 1) | −iπ (since S(O, n) = 1 − Log(−1) = 1 − iπ) | 35 |
| h = S(L(k), n) | −2πi: a full turn | 59 |
| b = L(h) | a logarithm of h: exp(b) = −2πi | 65 |
| a = L(S(1, n)) | exp(a) = e − iπ, hence exp(exp(a)) = −C | 25 |
| P(t) = X(S(a, Q(t))) | −C / conj(Q(t)) = −t, t ≠ 0: sign change | 33 |
| **J(z) = X(P(S(b, z)))** | **conj(z), z ≠ 0** | **101** |

**Why it is strict on all of ℂ \ {0}.** The input given to P is w = S(b, z) = −2πi − Log(conj z). Its imaginary
part is in [−3π, −π): w is never zero, even at z = 1 and on the negative real axis. The sign change
is therefore legitimate, and J(z) = exp(−w) = exp(2πi) · conj(z) = conj(z). The computations use only exp(Log u) = u; they
never replace Log(exp u) by u. It is the "full turn" (shifting an intermediate logarithm by −2πi, which
the next exponential erases) that avoids the intermediate zero.

**Consequence: the eml gate rebuilt.** For y ≠ 0, eml★(x, J(y)) = exp(x) − Log(conj(conj y)) = eml(x, y):
103 symbols. In strict evaluation, every eml gate of a mixed tree can therefore be replaced by an assembly in
eml★ alone; the two grammars compute the same strict functions.

| Result | Status | Supporting files |
|---|---|---|
| Mirror of 297 symbols in eml★ alone (first version, with a general multiplier), J(z) = conj(z) for z ≠ 0 | ESTABLISHED by exact judge: SymPy 1.14.0 on a nonzero complex symbol and special cases; Arb enclosures (python-flint 0.9.0) at the listed points; replayed by Anthony | `~/EML_replay_20261120_Vw6P7N/`: `replay.py` (sha256 6507109b…), `PROOFS.md` (928217ad…), `replay-console.txt` (ec77fe7a…) |
| Mirror of 101 symbols: eml★ only, leaves {1, z}, fully expanded length 101, exact constants, identity J(z) = conj(z) on a nonzero complex symbol and at 9 exact points (including −1, −50, e^e), strict refusal at 0, eml★(x, J(y)) = eml(x, y), star-free control | ESTABLISHED by exact judge, twice: (1) ChatGPT's judge replayed by Anthony (output "PASS … 297 -> 101", pasted into the conversation, not archived as a file); (2) standalone script written by Claude without reusing that judge, replayed by Anthony on his ThinkCentre: Python 3.12.3, SymPy 1.14.0, 24 checks out of 24, exit code 0 (personal trace `~/EML_101_replay_NSwy5F/`, console sha256 d48f62c4…, byte-for-byte identical to that of the repository). These are two distinct implementations, both verified with SymPy | `star_only_mirror_101.py`, `CONSOLE_20261120_star_only_mirror_101.txt` |
| Independent reconstruction of the two mirrors: 297 and 101 symbols, eml★ only, no infinite node; equality with conj(z) on 89 hard points (deviation < 4·10⁻⁴⁷; negative axis, boundaries ±π, \|Im\| up to 100) and in exact SymPy computation at 5 points, two of them on the negative axis; strict refusal at 0 | OBSERVED (sample + exact points) | `double_toolbox/` (in this repository; its French working copy was a local directory named `boite_double`): `star_only.py`, `test_toolbox.py` (37 tests) |
| Strict domain = all of ℂ \ {0} (Im w ∈ [−3π, −π)) | DERIVED: short proof above, reviewed by Claude; to be submitted for review | this section |
| eml(x, y) = eml★(x, J(y)) for y ≠ 0 | DERIVED (one line) + exact judge on symbols and at 3 points including y negative real (standalone script) | `star_only_mirror_101.py` |
| Strict product, strict real part on the whole plane (2,489 symbols), NAND gate on bits (535 symbols), in eml★ alone | OBSERVED on samples (415 pairs, 93 points, 4 inputs); proofs in `PROOFS.md`, to be reviewed | `double_toolbox/star_only.py` (in this repository); `PROOFS.md` (ChatGPT's text, not in the repository) |
| Density of INDIVIDUAL trees on compact sets; strict obstruction when 0 ∈ K and exp is not injective on K (example K = {0, 2πi}) | NOT VERIFIED here: written proof by ChatGPT (`PROOFS.md` § 5–6), awaiting independent review | `PROOFS.md` (not committed to the repository) |
| Is 101 minimal in eml★ alone? | OPEN. The exhaustive search of § 16 found nothing up to length 19 (heuristic, floats): the minimum is therefore probably between 21 and 101 | — |

**Scope clarifications.**
- *At zero.* The trees of this section refuse z = 0 in strict evaluation. This is not an impossibility at the single
  point 0 (the leaf z gives the right value there). ChatGPT's statement, not verified here, is that no strict tree
  can realise the mirror both at 0 and at 2πi. The mixed 9-symbol mirror of § 16 reaches 0 in EXTENDED
  evaluation.
- *The full turn.* Removing the −2πi shift IN THIS CONSTRUCTION makes the input 1 pass through ln 0 (local check,
  not committed). This does not prove that no other assembly could do without it.
- *Floats on the branch cut.* In floating point, eml★(x, J(y)) can differ from eml(x, y) by 2πi when y is a negative
  real (observed at y = −50): J(y) falls back, to within 10⁻⁴⁹, on the wrong side of the axis. The exact computation gives zero
  deviation. This is the numerical fragility already reported at the boundaries for `conjugate_formula` (`eml_toolkit/core.py`), not a defect of the identity.
- *The two mirrors do not answer the same question.* 9 symbols: the shortest in the mixed grammar, minimal,
  correct all the way to 0 (extended). 101 symbols: a single type of gate, strict, away from 0, minimality unknown.

Still open: the minimal length of the mirror in eml★ alone; the review of the density proof; a human
review of the whole.
