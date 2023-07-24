This is a simple script to calculate S(Q,E) from lammps trajectories using Fourier transform and correlation functions. It should work for any trajectory files as long as they contain atomic positions of each time step. The equations I used are from ref "Witte, B. B. L., et al. "Ab initio simulations of the dynamic ion structure factor of warm dense lithium." Physical Review B 95.14 (2017): 144105.", but any general reference should be okay.

1. Please use lampps\_input folder to generate lammps trajectories..

2. Then you go to sqe folder to calculate dynamical structure factor
``` 
python calculate_coh_sqe.py
```
Or submit a job by
```
sbatch post.sh
```
I did this on a cluster, which took about 60 mins.

3. Step 2 should generate `freq.txt` and `sqe_2K0_full.txt`. Then you run
```
python plot_2K0_full.py
```
to obtain the SQE figure.

Comments:
1. Gamma points also have singularities, which I don't know why yet.
2. Will check later when I have more time.
3. Once I have more time, I will also do phonon spectral energy density (SED) and phonon density of states (DOS) calculations
