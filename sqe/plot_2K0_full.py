import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams.update({'font.size':18})

sqe = np.loadtxt('sqe_2K0_full.txt',dtype=complex)
print (sqe.shape)
sqe = np.abs(sqe)
freq = np.loadtxt('freq_full.txt')
energy_bin = freq[1] - freq[0]
print (energy_bin)

#sqe = np.abs(sqe.imag)
sqe += np.flip(sqe,axis=1)
sqe /= 2

fig, ax = plt.subplots(figsize=(8,6))

vmin = 5e3
vmax = 5e5
im = ax.imshow(np.rot90(sqe),norm=LogNorm(vmin=vmin,vmax=vmax),aspect='auto',cmap='viridis',interpolation='bilinear')

ax.set_ylim(0,400)
yticks = np.linspace(0,400,6)
#yl = np.linspace(0,25,6)
yl = yticks*energy_bin
ax.set_yticks(yticks)
ax.set_yticklabels(np.round(yl,1))

xticks = np.linspace(0,20,5)
xl = np.linspace(-1,1,5)
ax.set_xticks(xticks)
ax.set_xticklabels(xl)

ax.set_ylabel('Energy (THz)')
ax.set_xlabel('2 K 0 (r.l.u.)')
fig.colorbar(im)
plt.savefig('2K0_full.png')
plt.show()
