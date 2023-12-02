import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

four = np.loadtxt('four.txt')

plt.rcParams.update({'font.size':16})
fig, ax = plt.subplots(figsize=(8,6))
H_min = -6
H_max = 6
K_min = -6
K_max = 6
ncell = 10
vmin = 1e1
vmax = 1e3
im = ax.imshow(four,origin='lower',norm=LogNorm(vmin=vmin,vmax=vmax))
xpos = np.linspace(0,ncell*(H_max-H_min),(H_max-H_min)+1,dtype=int)
xl = np.linspace(H_min,H_max,(H_max-H_min)+1,dtype=int)
ax.set_xticks(xpos)
ax.set_xticklabels(xl)
ypos = np.linspace(0,ncell*(K_max-K_min),(K_max-K_min)+1,dtype=int)
yl = np.linspace(K_min,K_max,(K_max-K_min)+1,dtype=int)
ax.set_yticks(ypos)
ax.set_yticklabels(yl)
ax.set_xlabel('H (r.l.u.)')
ax.set_ylabel('K (r.l.u.)')
fig.colorbar(im,ax=ax)
plt.savefig('elastic_diffusescattering.png',dpi=480)

