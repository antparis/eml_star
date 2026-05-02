"""
eml_toolkit/eml_quantum.py
Quantum applications of the eml★ operator system.

This module provides symbolic and numerical tools for quantum information
tasks where complex conjugation appears naturally:
  - Geometric (Berry) phases
  - Self-testing of entangled states up to complex conjugation
  - Phase extraction via eml★

Mathematical basis:
  phase(z)  = Im(ln(z))  = arg(z)       [principal branch, Im(z) in (-pi, pi)]
  conj(z)   = 1 - eml★(0, eml(z, 1))   [Theorem 3.1, Monnerot 2026]

References:
  [1] Monnerot (2026). eml★: Minimal Anti-Holomorphic Extension of the EML
      Sheffer Operator. arXiv comment on 2603.21852v2.
  [2] Supic, Balanzó-Juandó, Acín et al. (2026). Self-testing of multipartite
      qubit states up to complex conjugation. QIP 2026 / Quantum Journal.
  [3] Odrzywołek (2026). All elementary functions from a single operator.
      arXiv:2603.21852v2.

Author: Anthony Monnerot (2026)
"""

import mpmath

mpmath.mp.dps = 60

# ── Core eml★ operators ────────────────────────────────────────────────────────

ONE = mpmath.mpc(1)


def eml(x, y):
    """eml(x, y) = exp(x) - ln(y)  [Odrzywołek 2026]"""
    return mpmath.exp(x) - mpmath.log(y)


def eml_star(x, y):
    """eml★(x, y) = exp(x) - ln(conj(y))  [Monnerot 2026]"""
    return mpmath.exp(x) - mpmath.log(mpmath.conj(y))


def conjugate_formula(z):
    """
    Complex conjugate via eml★ (Theorem 3.1).
    conj(z) = 1 - eml★(0, eml(z, 1))
    Valid for Im(z) in (-pi, pi).
    """
    return 1 - eml_star(mpmath.mpc(0), eml(z, ONE))


# ── Quantum phase functions ────────────────────────────────────────────────────

def phase(z):
    """
    Phase (argument) of a complex number via eml★.

    phase(z) = Im(ln(z)) = arg(z)

    This is the imaginary part of the logarithm, extracted via the
    eml★ conjugation formula at relative depth 4.

    Valid for Im(z) in (-pi, pi) and z != 0.

    Parameters
    ----------
    z : complex number (mpmath.mpc)

    Returns
    -------
    mpmath.mpf : the phase angle in (-pi, pi]
    """
    z = mpmath.mpc(z)
    if abs(z) == 0:
        raise ValueError("phase(0) is undefined")
    ln_z = mpmath.log(z)
    # Im(w) = (w - conj(w)) / (2i)
    # conj(ln_z) via eml★ (valid when Im(ln_z) = arg(z) in (-pi, pi))
    conj_ln_z = conjugate_formula(ln_z)
    return mpmath.re((ln_z - conj_ln_z) / (2 * mpmath.j))


def phase_direct(z):
    """
    Phase via mpmath.arg — reference implementation for verification.
    """
    return mpmath.arg(mpmath.mpc(z))


# ── Geometric (Berry) phase ────────────────────────────────────────────────────

def geometric_phase(states):
    """
    Discrete geometric phase for a sequence of quantum states.

    For a closed path psi_0, psi_1, ..., psi_{n-1}, psi_0, the
    discrete geometric phase (Pancharatnam phase) is:

        gamma = -Im( ln( <psi_0|psi_1> * <psi_1|psi_2> * ... * <psi_{n-1}|psi_0> ) )

    This equals the Berry phase in the continuous adiabatic limit.

    Parameters
    ----------
    states : list of lists/tuples of mpmath.mpc
        Each element is a normalized quantum state vector.
        The path is automatically closed (last state connects back to first).

    Returns
    -------
    mpmath.mpf : geometric phase in (-pi, pi]

    Example
    -------
    >>> import mpmath
    >>> # Two orthogonal states on Bloch sphere equator
    >>> psi0 = [mpmath.mpc(1, 0), mpmath.mpc(0, 0)]
    >>> psi1 = [mpmath.mpc(0, 0), mpmath.mpc(1, 0)]
    >>> geometric_phase([psi0, psi1])
    """
    n = len(states)
    if n < 2:
        raise ValueError("Need at least 2 states for a geometric phase.")

    # Compute product of overlaps along the closed path
    product = mpmath.mpc(1)
    for k in range(n):
        bra = states[k]
        ket = states[(k + 1) % n]
        overlap = _inner_product(bra, ket)
        if abs(overlap) < 1e-50:
            raise ValueError(
                f"Zero overlap between states {k} and {(k+1) % n}. "
                "Geometric phase is ill-defined."
            )
        product *= overlap

    # Geometric phase = -Im(ln(product))
    ln_product = mpmath.log(product)
    return -mpmath.im(ln_product)


def geometric_phase_emlstar(states):
    """
    Geometric phase computed entirely via eml★ (no direct Im() call).

    Uses phase(z) = Im(ln(z)) extracted via eml★ conjugation formula.
    Demonstrates that eml★ can express the full Berry phase computation.

    Parameters
    ----------
    states : list of state vectors (same format as geometric_phase)

    Returns
    -------
    mpmath.mpf : geometric phase in (-pi, pi]
    """
    n = len(states)
    if n < 2:
        raise ValueError("Need at least 2 states.")

    product = mpmath.mpc(1)
    for k in range(n):
        bra = states[k]
        ket = states[(k + 1) % n]
        product *= _inner_product(bra, ket)

    # Use eml★ to extract the imaginary part
    return -phase(product)


def _inner_product(bra, ket):
    """
    Inner product <bra|ket> = sum_i conj(bra_i) * ket_i.

    The conjugation of bra components is done via eml★.
    """
    if len(bra) != len(ket):
        raise ValueError("State vectors must have the same dimension.")
    result = mpmath.mpc(0)
    for b, k in zip(bra, ket):
        b = mpmath.mpc(b)
        k = mpmath.mpc(k)
        # conj(b) via eml★ when Im(b) in (-pi, pi) — satisfied for normalized states
        if abs(float(mpmath.im(b))) < float(mpmath.pi) - 1e-10:
            conj_b = conjugate_formula(b)
        else:
            conj_b = mpmath.conj(b)  # fallback for out-of-strip values
        result += conj_b * k
    return result


# ── Self-testing up to complex conjugation ─────────────────────────────────────

def self_test_conjugate(psi, psi_ref):
    """
    Check whether two quantum states are equal up to complex conjugation.

    A state psi self-tests (up to complex conjugation) against psi_ref if:
        psi == psi_ref   OR   psi == conj(psi_ref)   (component-wise)

    This is the standard notion from self-testing theory [2]:
    many entangled states can only be certified up to this ambiguity
    because the conjugate state gives identical measurement statistics.

    Parameters
    ----------
    psi     : list of mpmath.mpc — state to test
    psi_ref : list of mpmath.mpc — reference state

    Returns
    -------
    dict with keys:
        'equal'        : bool — psi == psi_ref (within tolerance)
        'conj_equal'   : bool — psi == conj(psi_ref) (within tolerance)
        'self_tested'  : bool — True if either equality holds
        'error_direct' : float — max component error for psi == psi_ref
        'error_conj'   : float — max component error for psi == conj(psi_ref)

    Example
    -------
    >>> import mpmath
    >>> psi     = [mpmath.mpc(1, 0.5), mpmath.mpc(0, -0.5)]
    >>> psi_ref = [mpmath.mpc(1, -0.5), mpmath.mpc(0, 0.5)]  # conj of psi
    >>> result = self_test_conjugate(psi, psi_ref)
    >>> result['self_tested']
    True
    """
    if len(psi) != len(psi_ref):
        raise ValueError("States must have the same dimension.")

    tol = mpmath.mpf("1e-50")

    # Compute conj(psi_ref) via eml★ component-wise
    conj_psi_ref = []
    for c in psi_ref:
        c = mpmath.mpc(c)
        if abs(float(mpmath.im(c))) < float(mpmath.pi) - 1e-10:
            conj_psi_ref.append(conjugate_formula(c))
        else:
            conj_psi_ref.append(mpmath.conj(c))

    # Error for psi == psi_ref
    err_direct = max(abs(mpmath.mpc(p) - mpmath.mpc(r))
                     for p, r in zip(psi, psi_ref))

    # Error for psi == conj(psi_ref)
    err_conj = max(abs(mpmath.mpc(p) - c)
                   for p, c in zip(psi, conj_psi_ref))

    equal      = err_direct < tol
    conj_equal = err_conj   < tol

    return {
        'equal':        equal,
        'conj_equal':   conj_equal,
        'self_tested':  equal or conj_equal,
        'error_direct': float(err_direct),
        'error_conj':   float(err_conj),
    }


# ── Verification ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("EML_QUANTUM.PY — Verification")
    print(f"Precision: mpmath.dps = {mpmath.mp.dps}")
    print("=" * 65)

    passed = 0
    failed = 0

    def check(label, got, expected, tol=1e-50):
        global passed, failed
        err = abs(float(got) - float(expected))
        ok = err < tol
        print(f"  [{'PASS' if ok else 'FAIL'}]  {label}  (err={err:.2e})")
        if ok: passed += 1
        else:  failed += 1

    print("\n-- phase(z) = Im(ln(z)) via eml★ --")
    test_vals = [
        mpmath.mpc(1, 0),       # phase = 0
        mpmath.mpc(0, 1),       # phase = pi/2
        mpmath.mpc(-1, 0.001),  # phase ~ pi (near branch cut)
        mpmath.mpc(1, 1) / mpmath.sqrt(2),  # phase = pi/4
        mpmath.mpc(0.5, mpmath.sqrt(3)/2),  # phase = pi/3
    ]
    for z in test_vals:
        got = phase(z)
        ref = phase_direct(z)
        check(f"phase({z})", got, ref)

    print("\n-- geometric_phase vs geometric_phase_emlstar --")
    # Four states on Bloch sphere equator (expected phase = pi)
    states_equator = [
        [mpmath.mpc(1, 0) / mpmath.sqrt(2), mpmath.mpc(1, 0) / mpmath.sqrt(2)],
        [mpmath.mpc(1, 0) / mpmath.sqrt(2), mpmath.mpc(0, 1) / mpmath.sqrt(2)],
        [mpmath.mpc(1, 0) / mpmath.sqrt(2), mpmath.mpc(-1, 0) / mpmath.sqrt(2)],
        [mpmath.mpc(1, 0) / mpmath.sqrt(2), mpmath.mpc(0, -1) / mpmath.sqrt(2)],
    ]
    gp1 = geometric_phase(states_equator)
    gp2 = geometric_phase_emlstar(states_equator)
    # The discrete geometric phase for this path is -pi (solid angle = 2pi -> phase = -pi)
    check("geometric_phase (standard)", abs(gp1), float(mpmath.pi), tol=1e-10)
    check("geometric_phase via eml★", abs(gp2), float(mpmath.pi), tol=1e-10)
    err_between = abs(float(gp1) - float(gp2))
    print(f"  [{'PASS' if err_between < 1e-50 else 'FAIL'}]  "
          f"both methods agree  (diff={err_between:.2e})")
    if err_between < 1e-50: passed += 1
    else: failed += 1

    print("\n-- self_test_conjugate --")
    # psi and conj(psi) should self-test against each other
    psi     = [mpmath.mpc(1, 0.3), mpmath.mpc(0.5, -0.2)]
    psi_ref = [mpmath.mpc(1, -0.3), mpmath.mpc(0.5, 0.2)]  # conj(psi)

    result = self_test_conjugate(psi, psi_ref)
    ok = result['self_tested'] and result['conj_equal'] and not result['equal']
    print(f"  [{'PASS' if ok else 'FAIL'}]  psi self-tests against conj(psi_ref)")
    print(f"         equal={result['equal']}, conj_equal={result['conj_equal']}")
    print(f"         err_direct={result['error_direct']:.2e}, "
          f"err_conj={result['error_conj']:.2e}")
    if ok: passed += 1
    else:  failed += 1

    # Same state — should be directly equal
    result2 = self_test_conjugate(psi, psi)
    ok2 = result2['equal'] and result2['self_tested']
    print(f"  [{'PASS' if ok2 else 'FAIL'}]  psi self-tests against itself")
    if ok2: passed += 1
    else:   failed += 1

    print("\n" + "=" * 65)
    print(f"RESULT: {passed}/{passed + failed} passed")
    if failed == 0:
        print("ALL TESTS PASSED")
    else:
        print(f"{failed} TEST(S) FAILED")
    print("=" * 65)
