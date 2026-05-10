import numpy as np

def eml(x, y):
    y_safe = np.where(np.abs(y) < 1e-30, 1e-30 + 0j, y)
    return np.exp(x) - np.log(y_safe)

def eml_star(x, y):
    yc = np.complex128(y.real - 1j * y.imag)
    yc = np.where(np.abs(yc) < 1e-30, 1e-30 + 0j, yc)
    return np.exp(x) - np.log(yc)

def eml_zero(x, y):
    y_safe = np.where(np.abs(y) < 1e-30, 1e-30 + 0j, y)
    return np.exp(x) - 1j * np.angle(y_safe)

def count_violations(z, threshold=np.pi - 0.1):
    return int(np.sum(np.abs(np.angle(z)) > threshold))

omega = 1.0
t = np.linspace(0, 2 * np.pi, 100)
psi = np.exp(-1j * omega * t)

conj_psi = np.exp(1.0 - eml_star(0.0, psi))
mod_sq_emlstar = np.real(psi * conj_psi)
mod_sq_numpy = np.abs(psi)**2
mse = np.mean((mod_sq_emlstar - mod_sq_numpy)**2)
violations = count_violations(psi)

print("=" * 65)
print("QUANTUM OSCILLATION SIMULATION")
print("Using eml_star and eml_zero primitives")
print("=" * 65)
print(f"omega = {omega}")
print(f"Time points: {len(t)}")
print(f"MSE (|psi|^2 methods): {mse:.2e}")
print(f"Branch violations: {violations}/100")
print()
print("Sample (first 5 points):")
for i in range(5):
    print(f"  t={t[i]:.4f}  |psi|^2_emlstar={mod_sq_emlstar[i]:.6f}  |psi|^2_numpy={mod_sq_numpy[i]:.6f}")
print("=" * 65)
