import numpy as np

def eml(x, y):
    y_safe = np.where(np.abs(y) < 1e-30, 1e-30 + 0j, y)
    return np.exp(x) - np.log(y_safe)

def eml_star(x, y):
    yc = np.conj(y)
    yc = np.where(np.abs(yc) < 1e-30, 1e-30 + 0j, yc)
    return np.exp(x) - np.log(yc)

N = 40
x_re = np.linspace(-3.0, 3.0, N)
x_im = np.linspace(-np.pi + 0.1, np.pi - 0.1, N)
X, Y = np.meshgrid(x_re, x_im)
Z = (X + 1j * Y).ravel()
TARGET = np.conj(Z)

pred_star = 1 - eml_star(0, eml(Z, 1))
mse_star = np.mean(np.abs(pred_star - TARGET)**2)

pred_eml = 1 - eml(0, eml(Z, 1))
mse_eml = np.mean(np.abs(pred_eml - TARGET)**2)

print(f"MSE with eml* : {mse_star:.6e}")
print(f"MSE with eml  : {mse_eml:.6e}")
print(f"Ratio         : {mse_eml / mse_star:.1f}x")
