import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import groupby
import collections
from matplotlib import ticker
#from mayavi import mlab
from scipy import stats
from numpy.fft import ifft
from matplotlib.colors import LogNorm

#plt.rcParams.update({'font.size': 20})
#fig,ax = plt.subplots(figsize=(10,8))

## The input files of this script is lammps structure and trajectories ##
################### Read equilibrium position from supercell ###########

path = './'

total_step = 8192 ##2 power 13
total_num = 8000 ## Total atomic number
lc = 56.3095038318 ##Angstrom
ncell = 10 ## 10*10*10 supercell
unit_lc = lc/ncell ##10*10*10 supercell

initial = np.loadtxt(path+'structure.in',skiprows=19)
print (initial.shape)

positionx = initial[:,2]
positiony = initial[:,3]
positionz = initial[:,4]

### Read trajectories ###
def generate_specific_rows(filepath): ## in order to skip every 625+9 lines
    with open(filepath,'r') as f:
        for i, line in enumerate(f):
            if i%8009 > 8:
                yield line

trajectory = np.loadtxt(generate_specific_rows(path+'nvtpositions.lammpstrj'))
## Only take first total_steps
trajectory = trajectory[0:total_step*total_num,:]
print (trajectory.shape)

trajectory = np.reshape(trajectory,(total_step,total_num,3))
print (trajectory.shape)

##################### Calculate relative displacement ############
rel_trajec = np.empty((total_step,total_num,3))

print (rel_trajec.shape)
for i in range(total_step):
    for j in range(total_num):
        rel_trajec[i,j,0] = trajectory[i,j,0]-positionx[j]
        rel_trajec[i,j,1] = trajectory[i,j,1]-positiony[j]
        rel_trajec[i,j,2] = trajectory[i,j,2]-positionz[j]
        if rel_trajec[i,j,0] > 0.5*lc:
           rel_trajec[i,j,0] -= 1*lc
           trajectory[i,j,0] -= 1*lc
        if rel_trajec[i,j,1] > 0.5*lc:
           rel_trajec[i,j,1] -= 1*lc
           trajectory[i,j,1] -= 1*lc
        if rel_trajec[i,j,2] > 0.5*lc:
           rel_trajec[i,j,2] -= 1*lc
           trajectory[i,j,2] -= 1*lc
        if rel_trajec[i,j,0] < -0.5*lc:
           rel_trajec[i,j,0] += 1*lc
           trajectory[i,j,0] += 1*lc
        if rel_trajec[i,j,1] < -0.5*lc:
           rel_trajec[i,j,1] += 1*lc
           trajectory[i,j,1] += 1*lc
        if rel_trajec[i,j,2] < -0.5*lc:
           rel_trajec[i,j,2] += 1*lc
           trajectory[i,j,2] += 1*lc

## timestep
dtau = 0.002 #ps
dt = 0.002 #ps
reciprocal_mat = np.array(([2*np.pi/unit_lc,0,0],[0,2*np.pi/unit_lc,0],[0,0,2*np.pi/unit_lc])) #to reciprocal space


#### save consecutive 2000 steps ####
struct_avg = np.zeros((total_num,3))
start_step = 1000
end_step = 3000
temporary = np.zeros((total_num,3))
for timestep in range(start_step,end_step):
    temporary[:,0] += trajectory[timestep,:,0]/(end_step-start_step)
    temporary[:,1] += trajectory[timestep,:,1]/(end_step-start_step)
    temporary[:,2] += trajectory[timestep,:,2]/(end_step-start_step)
print (temporary)
struct_avg[:,0] = temporary[:,0]/lc
struct_avg[:,1] = temporary[:,1]/lc
struct_avg[:,2] = temporary[:,2]/lc

## Fourier transform to calcualte diffuse scattring ###
Qstep = 1/ncell
H_min = -6
H_max = 6
K_min = -6
K_max = 6
L = 0
Q_2D = np.mgrid[H_min:H_max+Qstep:Qstep, K_min:K_max+Qstep:Qstep].reshape(2,-1).T
Q_2D = Q_2D.reshape((ncell*(H_max-H_min)+1,ncell*(K_max-K_min)+1,2)) # the cell is 10 10 10

## coherent neutron scattering length ##
b_Na = 3.63
b_Cl = 9.577

four = np.zeros((ncell*(H_max-H_min)+1,ncell*(K_max-K_min)+1))
for i in range(ncell*(H_max-H_min)+1):
    for j in range(ncell*(K_max-K_min)+1):
        temp = 0
        Q = np.array([Q_2D[i,j][0],Q_2D[i,j][1],L])*ncell # 10 10 10 cells
        # Na
        temp = b_Na*np.sum(np.exp(2*np.pi*1j*np.inner(Q,struct_avg[0:4000])))
        # Cl
        temp += b_Cl*np.sum(np.exp(2*np.pi*1j*np.inner(Q,struct_avg[4000:8000])))
        four[i,j] = np.abs(temp)**2

np.savetxt('four.txt',four)

plt.rcParams.update({'font.size':16})
fig, ax = plt.subplots(figsize=(8,6))
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

