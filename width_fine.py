import numpy as np
from scipy.optimize import curve_fit
m = np.load("dw_fine.out/m000000.npy")[:, 0, 32:96, :].mean(axis=1)
x = np.arange(m.shape[1]) * 1.0
(x0, w), _ = curve_fit(lambda x, x0, w: -np.tanh((x - x0) / w), x, m[2], p0=[x.mean(), 5])
print(f"1 nm cells: width = {w:.2f} nm")