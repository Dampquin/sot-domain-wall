import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit

dx = 2.0  # cell size in nm
def wall(x, x0, w):            # textbook wall profile
    return -np.tanh((x - x0) / w)

fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for i, (lab, ax) in enumerate(zip(["D = 0", "D = 1 mJ/m²"], axes)):
    m = np.load(f"dw_relax.out/m{i:06d}.npy")
    print("array shape:", m.shape)            # expect (3, 1, 64, 512)
    m = m[:, 0, 16:48, :].mean(axis=1)       # average central rows, skip edges
    x = np.arange(m.shape[1]) * dx
    (x0, w), _ = curve_fit(wall, x, m[2], p0=[x.mean(), 5])
    c = np.argmin(abs(x - x0))
    print(f"{lab}: width Δ = {w:.2f} nm | at centre mx = {m[0,c]:+.2f}, my = {m[1,c]:+.2f}")
    ax.plot(x - x0, m[2], label="mz")
    ax.plot(x - x0, m[0], "--", label="mx")
    ax.plot(x - x0, m[1], ":", label="my")
    ax.plot(x - x0, wall(x, x0, w) , "k", lw=0.8, alpha=0.5, label="tanh fit")
    ax.set_xlim(-40, 40); ax.set_title(lab); ax.set_xlabel("x − x₀ (nm)")
axes[0].set_ylabel("m"); axes[0].legend()
plt.tight_layout(); plt.savefig("dw_profiles.png", dpi=150); plt.show()