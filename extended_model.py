import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import brentq

d = np.loadtxt("vs_dmi.out/vDJ.txt")      # columns: D, J, v

# ---------- constants ----------
gamma, Ms, t, alpha = 1.76e11, 1.1e6, 0.6e-9, 0.3
Delta, theta = 5.4e-9, 0.1
hbar, e, mu0 = 1.055e-34, 1.602e-19, 4e-7 * np.pi

# wall's own Bloch preference (small, from the wall's internal magnetic charge)
H_K = mu0 * Ms * t * np.log(2) / (np.pi * Delta)     # ~0.034 T
print(f"mu0*H_K = {H_K:.3f} T")

def wall_speed(J, D, HK):
    H_SH = hbar * theta * J / (2 * e * Ms * t)        # SOT strength (T)
    H_D  = D / (Ms * Delta)                           # DMI strength (T)
    # speed from the SOT side minus speed from the DMI + H_K side; zero at balance
    balance = lambda phi: (np.pi/2) * H_SH * np.cos(phi) / alpha \
                          - ((np.pi/2) * H_D * np.sin(phi) - HK/2 * np.sin(2*phi))
    phi = brentq(balance, 0, np.pi/2 - 1e-9)          # find the balance angle
    v = gamma * Delta * (np.pi/2) * H_SH * np.cos(phi) / alpha
    return v, np.degrees(phi)

# ---------- compare ----------
Jm = np.linspace(1e10, 2.1e12, 200)
print(" D    J(1e11)  v_sim  simple  diff%  extended  diff%  phi°")
for D in np.unique(d[:, 0]):
    rows = np.where(d[:, 0] == D)[0]
    J, v = d[rows, 1], abs(d[rows, 2])
    p = plt.plot(J/1e11, v, "o", label=f"simulation, D = {D*1e3:.1f}")
    c = p[0].get_color()
    plt.plot(Jm/1e11, [wall_speed(j, D, 0)[0]   for j in Jm], "--", color=c, alpha=0.5)
    plt.plot(Jm/1e11, [wall_speed(j, D, H_K)[0] for j in Jm], "-",  color=c)
    for Ji, vi in zip(J, v):
        vs, _   = wall_speed(Ji, D, 0)
        ve, phi = wall_speed(Ji, D, H_K)
        print(f"{D*1e3:.1f}   {Ji/1e11:4.0f}   {vi:6.1f}  {vs:6.1f}  {100*(vi-vs)/vs:+5.1f}"
              f"   {ve:6.1f}  {100*(vi-ve)/ve:+5.1f}   {phi:4.0f}")

plt.xlabel("J (10¹¹ A/m²)"); plt.ylabel("|v| (m/s)")
plt.title("Dots: simulation | Dashed: Simple Model | Solid: + Wall Bloch preference", fontsize=9)
plt.legend(); plt.grid(alpha=0.3)
plt.savefig("extended_model.png", dpi=150); plt.show()