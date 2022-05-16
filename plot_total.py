import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':20})
fig, ax = plt.subplots(figsize=(10,8))


#freq = np.loadtxt('freq_test.txt')
sqe = np.loadtxt('sqe_2k0.txt',dtype=complex)
print (sqe.shape)

sqe = sqe[:,0:4096].real**2
plt.imshow(np.rot90(np.log10(sqe)),vmin=6.5,vmax=8.5,aspect='auto')


## the real energy step should be 1000/(2^13*2)*4.136 = 0.2524 meV. Here, I simply use 0.25 meV
yticks = np.arange(4096,0,-10)
yl = np.arange(0,1024,2.5)
ax.set_yticks(yticks)
ax.set_yticklabels(yl)
ax.set_ylim(4096,3996)
ax.set_ylabel('Energy (meV)')

xticks = np.linspace(0,20,11)
xl = np.linspace(-1,1,11)
ax.set_xticks(xticks)
ax.set_xticklabels(np.around(xl,decimals=1))
ax.set_xlabel('[2,K,0] (r.l.u)')

plt.colorbar()
plt.savefig('2K0_sqe.png')
#plt.show()
