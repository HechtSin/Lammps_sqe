import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter


plt.rcParams.update({'font.size':18})

sqe = np.loadtxt('sqe_2k0.txt',dtype=complex)
freq = np.loadtxt('freq.txt')
energy_bin = freq[1]-freq[0] # THz
print (sqe.shape)
sqe = np.abs(sqe)
sqe += np.flip(sqe,axis=1)
sqe /= 2

fig, ax = plt.subplots(figsize=(8,6))

vmin = 1e3
vmax = 1e5
im = ax.imshow(np.rot90(sqe),norm=LogNorm(vmin=vmin,vmax=vmax),aspect='auto',cmap='viridis')

ax.set_ylim(0,200)
yticks = np.linspace(0,200,6)
yl = yticks*energy_bin
ax.set_yticks(yticks)
ax.set_yticklabels(np.round(yl,2))

xticks = np.linspace(0,20,5)
xl = np.linspace(-1,1,5)
ax.set_xticks(xticks)
ax.set_xticklabels(xl)

ax.set_ylabel('Energy (meV)')
ax.set_xlabel('2 K 0 (r.l.u.)')
fig.colorbar(im)
plt.savefig('2K0.png')
plt.show()
