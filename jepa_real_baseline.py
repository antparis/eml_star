import torch
import torch.nn as nn

class RealEncoder(nn.Module):
    def __init__(self, input_dim=128, latent_dim=64):
        super().__init__()
        self.fc = nn.Linear(input_dim, latent_dim)
    def forward(self, x):
        return self.fc(x)

class RealPredictor(nn.Module):
    def __init__(self, latent_dim=64, action_dim=8):
        super().__init__()
        self.fc = nn.Linear(latent_dim + action_dim, latent_dim)
        self.activation = nn.ReLU()
    def forward(self, z_t, action):
        return self.activation(self.fc(torch.cat([z_t, action], dim=-1)))

def train_real_baseline(batch_size=32, latent_dim=64, action_dim=8, num_epochs=5, device="cpu"):
    encoder = RealEncoder(128, latent_dim).to(device)
    predictor = RealPredictor(latent_dim, action_dim).to(device)
    optimizer = torch.optim.Adam(list(encoder.parameters()) + list(predictor.parameters()), lr=1e-3)
    mse_loss = nn.MSELoss()
    print("Starting REAL BASELINE training (synthetic data, no eml operators)...")
    for epoch in range(num_epochs):
        x_t = torch.randn(batch_size, 128, device=device)
        x_t1 = torch.randn(batch_size, 128, device=device)
        action = torch.randn(batch_size, action_dim, device=device)
        z_t = encoder(x_t)
        z_t1_target = encoder(x_t1)
        z_t1_pred = predictor(z_t, action)
        loss = mse_loss(z_t1_pred, z_t1_target)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        print(f"Epoch {epoch+1}/{num_epochs} | Loss: {loss.item():.6f}")
    print("Real baseline training finished.")

if __name__ == "__main__":
    train_real_baseline()
