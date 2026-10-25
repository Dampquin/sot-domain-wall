import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import brentq

OUT, DX = "final_1nm.out", 1.0                     # output folder, cell size (nm)
d = np.loadtxt(f"{OUT}/vDJ.txt")                   # columns: D, J, v

# ---------- constants ----------
gamma, Ms, t, alpha = 1.76e11, 1.1e6, 0.6e-9, 0.3
Delta, theta = 5.49e-9, 0.1                        # wall width measured at 1 nm cells
hbar, e, mu0 = 1.055e-34, 1.602e-19, 4e-7 * np.pi
H_K = mu0 * Ms * t * np.log(2) / (np.pi * Delta)   # wall's own Bloch preference

def wall_speed(J, D, HK):                          # 1D model (HK = 0 -> simple version)
    H_SH = hbar * theta * J / (2 * e * Ms * t)
    H_D  = D / (Ms * Delta)
    bal = lambda p: (np.pi/2)*H_SH*np.cos(p)/alpha - ((np.pi/2)*H_D*np.sin(p) - HK/2*np.sin(2*p))
    p = brentq(bal, 0, np.pi/2 - 1e-9)
    return gamma * Delta * (np.pi/2) * H_SH * np.cos(p) / alpha

def tilt_angle(k):                                 # tilt of the wall line, from snapshot k
    mz = np.load(f"{OUT}/m{k:06d}.npy")[2, 0]      # shape (128 rows, 1024 columns)
    xs, ys = [], []
    for row in range(16, 112):                     # skip 16 nm at each strip edge
        line = mz[row]
        i = np.where(np.diff(np.sign(line)) != 0)[0][0]
        xs.append((i + line[i] / (line[i] - line[i + 1])) * DX)
        ys.append(row * DX)
    return abs(np.degrees(np.arctan(np.polyfit(ys, xs, 1)[0])))

# ---------- table + figure ----------
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 4.8))
Jm = np.linspace(1e10, 2.1e12, 150)
rows_out = []
print(" D    J(1e11)  v_sim  simple  diff%  extended  diff%  tilt°")
for D in np.unique(d[:, 0]):
    idx = np.where(d[:, 0] == D)[0]
    J, v = d[idx, 1], abs(d[idx, 2])
    tilt = np.array([tilt_angle(k) for k in idx])
    vs = np.array([wall_speed(j, D, 0) for j in J])
    ve = np.array([wall_speed(j, D, H_K) for j in J])
    c = a1.plot(J/1e11, v, "o", label=f"D = {D*1e3:.1f} mJ/m²")[0].get_color()
    a1.plot(Jm/1e11, [wall_speed(j, D, 0) for j in Jm], "--", color=c, alpha=0.5)
    a1.plot(Jm/1e11, [wall_speed(j, D, H_K) for j in Jm], "-", color=c)
    a2.plot(J/1e11, tilt, "o-", color=c, label=f"D = {D*1e3:.1f} mJ/m²")
    a3.plot(tilt, -100*(v - ve)/ve, "o", color=c, label=f"D = {D*1e3:.1f} mJ/m²")
    for i in range(len(J)):
        print(f"{D*1e3:.1f}   {J[i]/1e11:4.0f}   {v[i]:6.1f}  {vs[i]:6.1f}  {100*(v[i]-vs[i])/vs[i]:+5.1f}"
              f"   {ve[i]:6.1f}  {100*(v[i]-ve[i])/ve[i]:+5.1f}  {tilt[i]:5.1f}")
        rows_out.append([D*1e3, J[i], v[i], vs[i], ve[i], tilt[i]])

a1.set(xlabel="J (10¹¹ A/m²)", ylabel="|v| (m/s)", title="Speed: dots sim | dashed simple | solid extended")
a2.set(xlabel="J (10¹¹ A/m²)", ylabel="Tilt angle (°)", title="Wall tilt while moving")
a3.set(xlabel="Tilt angle (°)", ylabel="Speed deficit vs extended model (%)", title="Remaining gap vs tilt")
for a in (a1, a2, a3): a.legend(); a.grid(alpha=0.3)
plt.tight_layout(); plt.savefig("final_figure.png", dpi=200)
np.savetxt("final_results.csv", rows_out, delimiter=",",
           header="D_mJm2,J_Am2,v_sim,v_simple,v_extended,tilt_deg")
plt.show()