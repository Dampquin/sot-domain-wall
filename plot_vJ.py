import numpy as np, matplotlib.pyplot as plt

# --- simulation data ---
d = np.loadtxt("sot_sweep.out/vJ.txt")
J, v = d[:, 0], abs(d[:, 1])

# --- 1D model (pen-and-paper theory) ---
gamma = 1.76e11      # gyromagnetic ratio
Ms, t, alpha = 1.1e6, 0.6e-9, 0.3
D, Delta, theta = 1e-3, 5.4e-9, 0.1
hbar, e = 1.055e-34, 1.602e-19

Jm = np.linspace(0, 2.1e12, 200)
H_SH = hbar * theta * Jm / (2 * e * Ms * t)   # effective field from SOT (T)
H_D  = D / (Ms * Delta)                       # effective field from DMI (T)
vD   = np.pi / 2 * gamma * Delta * H_D        # DMI speed limit
v_model = vD / np.sqrt(1 + (alpha * H_D / np.maximum(H_SH, 1e-12))**2)

plt.plot(J / 1e11, v, "o", label="mumax3 simulation")
plt.plot(Jm / 1e11, v_model, "-", label="1D model")
plt.axhline(vD, ls="--", c="gray", label=f"DMI speed limit ≈ {vD:.0f} m/s")
plt.xlabel("Current density J (10¹¹ A/m²)")
plt.ylabel("Wall speed |v| (m/s)")
plt.title("SOT-driven Néel wall in Pt/Co (D = 1 mJ/m²)")
plt.legend(); plt.grid(alpha=0.3)
plt.savefig("vJ.png", dpi=150); plt.show()