import numpy as np, matplotlib.pyplot as plt
d = np.loadtxt("sp4.out/table.txt")
t = d[:,0]*1e9
plt.plot(t, d[:,1], label="<mx>")
plt.plot(t, d[:,2], label="<my>")
plt.plot(t, d[:,3], label="<mz>")
plt.xlabel("t (ns)"); plt.ylabel("m"); plt.legend()
plt.savefig("sp4.png", dpi=150)
plt.show()