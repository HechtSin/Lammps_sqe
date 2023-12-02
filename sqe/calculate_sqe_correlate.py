import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import groupby
import collections
from matplotlib import ticker
#from mayavi import mlab
from scipy import stats
from numpy.fft import ifft

#plt.rcParams.update({'font.size': 20})
#fig,ax = plt.subplots(figsize=(10,8))

## The input files of this script is lammps structure and trajectories ##
################### Read equilibrium position from supercell ###########

############################
#### START of USER INPUT ###
############################
path = './'
total_step = 8192 ##2 power 13; better use the power of 2 for FFT
total_num = 8000 ## Total atomic number
lc = 56.3095038318 ##lattice constant of supercell in Angstrom
ncell = 10 ## 10*10*10 supercell
unit_lc = lc/ncell ##10*10*10 supercell
initial_struc_file = 'structure.in' ## initial structure file
trajectory_file = 'nvtpositions.lammpstrj' ## trajectory file of atomic positions
time_bin = 0.002 ## time step of MD in unit of ps
Q_steps = 21
k_list = np.linspace([2,-1,0],[2,1,0],Q_steps) ## Q points for SQE simulations
############################
##### END of USER INPUT ####
############################

##############################################################
initial = np.loadtxt(path+initial_struc_file,skiprows=19)
print (initial.shape)

positionx = initial[:,2]
positiony = initial[:,3]
positionz = initial[:,4]

### Read trajectories ###
def generate_specific_rows(filepath): ## in order to skip every 625+9 lines
    with open(filepath,'r') as f:
        for i, line in enumerate(f):
            if i%(total_num+9) > 8:
                yield line

trajectory = np.loadtxt(generate_specific_rows(path+trajectory_file))
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
dtau = time_bin #ps
dt = time_bin #ps
reciprocal_mat = np.array(([2*np.pi/unit_lc,0,0],[0,2*np.pi/unit_lc,0],[0,0,2*np.pi/unit_lc])) #to reciprocal space

## calculate at a single k point
#k_raw = (2,0.5,0)
#k = np.matmul(k_raw,reciprocal_mat)
#rho = np.zeros(total_step,dtype=complex)
#phase = np.zeros((total_step,total_num),dtype=complex)
## Calculate phase first
#for t in range(total_step):
#    for i in range(total_num):
#        phase[t,i] = np.exp(-1j*np.dot(k,trajectory[t,i,:]))
## calculate F(k,t)
#Fkt = np.zeros(total_step,dtype=complex)
#for t in range(total_step):
#    for tau in range(0,total_step-t):
#        Fkt[t] += np.sum(phase[tau,:])*np.sum(np.conjugate(phase[t+tau,:]))
#sqe = ifft(Fkt) 
#freq = np.fft.fftfreq(total_step,d=dt)
#np.savetxt('freq_single.txt',freq)
#np.savetxt('sqe_single.txt',sqe)


### loop k points
sqe_total = np.zeros((Q_steps,total_step),dtype=complex)
sqe_total_full = np.zeros((Q_steps,total_step*2-1),dtype=complex)
for index, k_raw in enumerate(k_list):
    print (k_raw)
    k = np.matmul(k_raw,reciprocal_mat)
    rho = np.zeros(total_step,dtype=complex)

    phase = np.zeros((total_step,total_num),dtype=complex)
    ## Calculate phase first
    for t in range(total_step):
        for i in range(total_num):
            phase[t,i] = np.exp(-1j*np.dot(k,trajectory[t,i,:]))

    ## calculate F(k,t)
    Fkt = np.zeros(total_step,dtype=complex)

    ## Hard code, very slow ##
    #for t in range(total_step):
    #    for tau in range(0,total_step-t):
    #        Fkt[t] += np.sum(phase[tau,:])*np.sum(np.conjugate(phase[t+tau,:]))

    ## numpy correlate ##
    phase_sum = np.zeros(total_step,dtype=complex)

    for t in range(total_step):
        phase_sum[t] = np.sum(phase[t,:])
        
    Fkt = np.correlate(phase_sum,phase_sum,mode='same')
    Fkt_full = np.correlate(phase_sum,phase_sum,mode='full')

    sqe = ifft(Fkt)
    sqe_full = ifft(Fkt_full)
    freq = np.fft.fftfreq(total_step,d=dt)

    sqe_total[index,:] = sqe
    sqe_total_full[index,:] = sqe_full

np.savetxt('freq.txt',freq) #Thz
np.savetxt('sqe_2k0.txt',sqe_total)
np.savetxt('sqe_2k0_full.txt',sqe_total_full)
