This is a simple script to calculate S(Q,E) from lammps trajectories using Fourier transform and correlation functions. It should work for any trajectory files as long as they contain atomic positions of each time step. The equations I used are from ref "Witte, B. B. L., et al. "Ab initio simulations of the dynamic ion structure factor of warm dense lithium." Physical Review B 95.14 (2017): 144105.", but any general reference should be okay.

1. Please unzip trajectory.zip first, which would provide the initial structure and trajectories containing atomic positions of each timestep.

2. Then you run
``` 
python calculate_sqe_correlate.py
```

I did this on a cluster, which took about 10 mins.

3. Step 2 should generate `freq.txt` and `sqe_2k0.txt`. Then you run
```
python plot_total.py
```

to obtain the SQE figure.

Comments:
1. My math is not good, but I hope this script is correct.
2. I didn't include neutron/x-ray scattering length, so the intensities are not directly comparable with experiments. But can be included easily.
3. The supercell size is 5*5*5 so incommensurate Q points have signularity.
4. Gamma points also have singularities, which I don't know why yet.
5. Will check later when I have more time.
