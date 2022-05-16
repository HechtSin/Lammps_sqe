#!/opt/apps/rhel8/Anaconda3-2021.05/bin/python -u

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

path = './'

total_step = 8192 ##2 power 13
total_num = 625 ## Total atomic number
lc = 20.1271510125 ##Angstrom
unit_lc = lc/5 ##5*5*5 supercell

initial = np.loadtxt(path+'structure_sc555.in',skiprows=20)
print (initial.shape,flush=True)

positionx = initial[:,2]
positiony = initial[:,3]
positionz = initial[:,4]


### Read trajectories ###
def generate_specific_rows(filepath): ## in order to skip every 625+9 lines
    with open(filepath,'r') as f:
        for i, line in enumerate(f):
            if i%634 > 8:
                yield line

trajectory = np.loadtxt(generate_specific_rows(path+'nvtpositions.lammpstrj'))
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
k_list = np.linspace([2,-1,0],[2,1,0],21)
sqe_total = np.zeros((21,total_step),dtype=complex)
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

    sqe = ifft(Fkt)
    freq = np.fft.fftfreq(total_step,d=dt)

    sqe_total[index,:] = sqe

np.savetxt('freq.txt',freq) #Thz
np.savetxt('sqe_2k0.txt',sqe_total)
