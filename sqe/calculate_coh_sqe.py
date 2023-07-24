import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import groupby
import collections
from matplotlib import ticker
from scipy import stats
from numpy.fft import ifft

## The input files of this script is lammps structure and atomic position trajectories ##


########## Start of User Inputs ##########

path = './'
trajectory_file = 'nvtpositions.lammpstrj' # atomic position trajectory file
lc = 56.3095038318 #Angstrom
unit_lc = lc/10 #10*10*10 supercell
initial_structure = 'structure.in' # initial structural file
time_step = 0.002 #time step of lammps trajectory in ps
total_step = 8192 ## Steps from trajectories used to do the calculation. Use power of two steps to speed up
total_num = 8000 ## Total atomic number
atomic_type = 2 # number of atomic types
atomic_num = np.array([4000,4000]) # number of each atomic type in the same order of structural file
#atomic_mass = np.array([22.989769,35.453]) # atomic mass of each atomic specie in the same order of structural file
coh_length = np.array([3.63,9.577]) # Coherent neutron scattering length in the same order of atomic types. Taken from https://www.ncnr.nist.gov/resources/n-lengths/
k_list = np.linspace([2,-1,0],[2,1,0],21) ## Start Q, End Q, number of Q
sqe_file = 'sqe_2K0.txt'
sqe_file_full = 'sqe_2K0_full.txt'


########## End of User Inputs ########





################### Read equilibrium position from supercell ###########
initial = np.loadtxt(path+initial_structure,skiprows=19)
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

##################### Calculate relative displacement and correct periodic boundary condition ############
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
dtau = time_step #ps
dt = time_step #ps
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
sqe_total = np.zeros((k_list.shape[0],total_step),dtype=complex)
sqe_total_full = np.zeros((k_list.shape[0],total_step*2-1),dtype=complex)
for index, k_raw in enumerate(k_list):
    print (k_raw)
    k = np.matmul(k_raw,reciprocal_mat)
    #rho = np.zeros(total_step,dtype=complex)

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
        for atype in range(atomic_type):
            if atype == 0:
                phase_sum[t] += coh_length[0]*np.sum(phase[t,0:atomic_num[0]])
            else:
                phase_sum[t] += coh_length[atype]*np.sum(phase[t,atomic_num[atype-1]:atomic_num[atype]])
        
    Fkt = np.correlate(phase_sum,phase_sum,mode='same')
    Fkt_full = np.correlate(phase_sum,phase_sum,mode='full')

    sqe = ifft(Fkt)
    sqe_full = ifft(Fkt_full)
    freq = np.fft.fftfreq(total_step,d=dt)
    freq_full = np.fft.fftfreq(total_step*2-1,d=dt)

    sqe_total[index,:] = sqe
    sqe_total_full[index,:] = sqe_full

np.savetxt('freq.txt',freq) #Thz
np.savetxt('freq_full.txt',freq_full) #Thz
np.savetxt(sqe_file,sqe_total)
np.savetxt(sqe_file_full,sqe_total_full)
