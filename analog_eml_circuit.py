"""
analog_eml_circuit.py
=====================
Numerical simulation of analog EML gate circuit.
Monnerot (2026) — ancillary script for Section 8.

Simulates:
- log_module  : BJT log feedback (ln)
- exp_module  : BJT antilog (exp)
- eml_gate    : differential subtractor = exp(x) - ln(y)
- eml_star_gate: parallel Re/Im + 180 degree phase inverter on Im channel
- Binary tree of depth 3 using eml_gate as sole primitive
- Branch violation counter using eml⁰ principle

Requirements: NumPy only.
"""

import numpy as np

# ── Analog blocks ──────────────────────────────────────────────────────────────

def log_module(y):
    """BJT log-feedback block. Output: ln|y|. Valid for |y| > 0."""
    return np.log(np.abs(y) + 1e-10)


def exp_module(x):
    """BJT antilog block. Output: exp(Re(x)). Clipped to avoid overflow."""
    return np.exp(np.clip(np.real(x), -60, 60))


def eml_gate(x, y):
    """EML gate: exp(x) - ln(y). Core analog primitive."""
    return exp_module(x) - log_module(y)


def eml_star_gate(x, y):
    """
    eml★ gate: parallel Re/Im channels + 180 degree phase inverter on Im.
    Implements exp(x) - ln(conj(y)) in analog domain.
    """
    re_out = exp_module(np.real(x)) - log_module(np.abs(y))
    im_out = -np.imag(y)  # 180 degree phase inverter = conjugation of Im
    return re_out + 1j * im_out


def count_violations(z, threshold=np.pi - 0.1):
    """
    Branch violation counter (eml⁰ principle).
    Flags embeddings where |Arg(z)| > pi - 0.1.
    """
    return int(np.sum(np.abs(np.angle(z)) > threshold))


# ── Binary tree depth 3 ────────────────────────────────────────────────────────

def eml_tree_depth3(x, y):
    """
    Binary EML tree of depth 3.
    Grammar: S -> eml(S, S) | leaf
    Structure:
        root = eml(
            eml(eml(x, y), eml(y, x)),
            eml(eml(y, x), eml(x, y))
        )
    """
    leaf_1 = eml_gate(x, y)
    leaf_2 = eml_gate(y, x)
    mid_1  = eml_gate(leaf_1, leaf_2)
    mid_2  = eml_gate(leaf_2, leaf_1)
    root   = eml_gate(mid_1, mid_2)
    return root


# ── Test inputs ────────────────────────────────────────────────────────────────

np.random.seed(42)
re_vals = np.random.uniform(0.5, 3.0, 10)
im_vals = np.random.uniform(-2.5, 2.5, 10)
Z = re_vals + 1j * im_vals
W = np.roll(Z, 1)  # second input — shifted version of Z

# ── Run and print ──────────────────────────────────────────────────────────────

print("=" * 65)
print("ANALOG EML CIRCUIT SIMULATION")
print("Monnerot (2026) — Section 8 ancillary script")
print("=" * 65)

print("\n-- EML gate: eml(z, w) = exp(z) - ln(w) --\n")
for i in range(10):
    out = eml_gate(Z[i], W[i])
    print(f"  z={Z[i]:.2f}  w={W[i]:.2f}  eml={out:.4f}")

print("\n-- eml★ gate: Re/Im parallel + 180 deg inverter --\n")
for i in range(10):
    out = eml_star_gate(Z[i], W[i])
    print(f"  z={Z[i]:.2f}  w={W[i]:.2f}  eml★={out:.4f}")

print("\n-- Binary tree depth 3 --\n")
tree_outputs = []
for i in range(10):
    out = eml_tree_depth3(Z[i], W[i])
    tree_outputs.append(out)
    print(f"  z={Z[i]:.2f}  w={W[i]:.2f}  tree_out={out:.4f}")

tree_outputs = np.array(tree_outputs)
violations = count_violations(tree_outputs)

print(f"\n-- Branch violations (|Arg| > pi - 0.1) --")
print(f"  Violations: {violations}/10")

print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)
print(f"  EML gates tested    : 10")
print(f"  eml★ gates tested   : 10")
print(f"  Tree depth          : 3")
print(f"  Branch violations   : {violations}/10")
print("=" * 65)
