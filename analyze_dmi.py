import numpy as np, matplotlib.pyplot as plt

d = np.loadtxt("vs_dmi.out/vDJ.txt")      # columns: D, J, v

# ---------- 1D model ----------
gamma, Ms, t, alpha = 1.76e11, 1.1e6, 0.6e-9, 0.3
Delta, theta = 5.4e-9, 0.1
hbar, e = 1.055e-34, 1.602e-19

def model(J, D):
    H_SH = hbar * theta * J / (2 * e * Ms * t)
    H_D  = D / (Ms * Delta)
    vD   = np.pi / 2 * gamma * Delta * H_D
    return vD / np.sqrt(1 + (alpha * H_D / np.maximum(H_SH, 1e-12))**2)

# ---------- tilt angle from a snapshot ----------
def tilt_angle(k):
    mz = np.load(f"vs_dmi.out/m{k:06d}.npy")[2, 0]     # shape (64 rows, 512 columns)
    xs, ys = [], []
    for row in range(8, 56):                           # skip rows near the edges
        line = mz[row]
        i = np.where(np.diff(np.sign(line)) != 0)[0][0]    # where mz changes sign = the wall
        xw = i + line[i] / (line[i] - line[i + 1])         # fine-tune between two cells
        xs.append(xw * 2.0); ys.append(row * 2.0)          # convert to nm
    slope = np.polyfit(ys, xs, 1)[0]                   # fit a straight line through the wall
    return abs(np.degrees(np.arctan(slope)))

# ---------- plots ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
Jm = np.linspace(1e10, 2.1e12, 200)
print(" D(mJ/m²)  J(1e11)  v_sim  v_model  diff%  tilt°")
for D in np.unique(d[:, 0]):
    rows = np.where(d[:, 0] == D)[0]
    J, v = d[rows, 1], abs(d[rows, 2])
    tilts = [tilt_angle(k) for k in rows]
    p = ax1.plot(J / 1e11, v, "o", label=f"simulation, D = {D*1e3:.1f}")
    ax1.plot(Jm / 1e11, model(Jm, D), "-", color=p[0].get_color(), alpha=0.6)
    ax2.plot(J / 1e11, tilts, "o-", color=p[0].get_color(), label=f"D = {D*1e3:.1f} mJ/m²")
    for Ji, vi, ti in zip(J, v, tilts):
        vm = model(Ji, D)
        print(f"  {D*1e3:.1f}      {Ji/1e11:4.0f}   {vi:6.1f}  {vm:6.1f}  {100*(vi-vm)/vm:+5.1f}  {ti:5.1f}")

ax1.set_xlabel("J (10¹¹ A/m²)"); ax1.set_ylabel("|v| (m/s)")
ax1.set_title("Wall speed: dots = simulation, lines = 1D model"); ax1.legend(); ax1.grid(alpha=0.3)
ax2.set_xlabel("J (10¹¹ A/m²)"); ax2.set_ylabel("Tilt angle (°)")
ax2.set_title("Wall tilt while moving"); ax2.legend(); ax2.grid(alpha=0.3)
plt.tight_layout(); plt.savefig("vs_dmi.png", dpi=150); plt.show()