import torch
import torch.nn as nn

def eml_star(x, y):
    y_safe = torch.where(torch.abs(y) < 1e-30, torch.tensor(1e-30, dtype=y.dtype, device=y.device) + 0j, y)
    conj_y = torch.complex(y_safe.real, -y_safe.imag)
    return torch.exp(x) - torch.log(conj_y)

def eml_zero(x, y):
    y_safe = torch.where(torch.abs(y) < 1e-30, torch.tensor(1e-30, dtype=y.dtype, device=y.device) + 0j, y)
    return torch.exp(x) - 1j * torch.angle(y_safe)

def count_branch_violations(z, threshold=torch.pi - 0.1):
    return int((torch.abs(torch.angle(z)) > threshold).sum().item())

class ComplexEncoder(nn.Module):
    def __init__(self, input_dim=128, latent_dim=32):
        super().__init__()
        self.fc = nn.Linear(input_dim, latent_dim * 2)
    def forward(self, x):
        out = self.fc(x)
        real, imag = out.chunk(2, dim=-1)
        return torch.complex(real, imag)

class JEPAPredictor(nn.Module):
    def __init__(self, latent_dim=32, action_dim=8):
        super().__init__()
        self.fc = nn.Linear(latent_dim * 2 + action_dim, latent_dim * 2)
    def forward(self, z_t, action):
        z_real_imag = torch.cat([z_t.real, z_t.imag, action], dim=-1)
        out = self.fc(z_real_imag)
        real_part, imag_part = out.chunk(2, dim=-1)
        z_pred = torch.complex(real_part, imag_part)
        assert z_pred.shape == z_t.shape
        z_pred = eml_star(z_pred, z_t + 1e-8)
        _ = eml_zero(torch.zeros_like(z_pred), z_pred)
        return z_pred

encoder = ComplexEncoder()
predictor = JEPAPredictor()
optimizer = torch.optim.Adam(list(encoder.parameters()) + list(predictor.parameters()), lr=1e-3)
mse_loss = nn.MSELoss()

print("Starting EML-WM complex prototype (latent_dim=32 complex = 64 real)...")
for epoch in range(5):
    x_t = torch.randn(32, 128)
    x_t1 = torch.randn(32, 128)
    action = torch.randn(32, 8)
    z_t = encoder(x_t)
    z_t1_target = encoder(x_t1)
    z_t1_pred = predictor(z_t, action)
    loss = mse_loss(z_t1_pred.real, z_t1_target.real) + mse_loss(z_t1_pred.imag, z_t1_target.imag)
    violations = count_branch_violations(z_t1_pred)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    print(f"Epoch {epoch+1}/5 | Loss: {loss.item():.6f} | Branch violations: {violations}")
print("Done.")
