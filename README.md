*mps* - 3D incompressible fluid dynamics solver by using moving particle semi-implicit(MPS) method with Fortran.  
        the code is parallelized by OpenMP.
======
moving particle semi-implicit(MPS) method: http://li.mit.edu/Stuff/CNSE/Paper/Koshizuka96Oka.pdf  
  
how to run:  
~~~~
$ vi mps.F90 # adjust the parameters  
$ make  
$ ./mps  
$ ./create_anime.sh  
~~~~
  
TODO:  
---
- create the animation by matplotlib
- implement ICCG method and replace CG method with it
- create a function to find neighbor particles in the effective radius(using hash table is better?)
- continue the verification
- parallelize the code by MPI
- still some particles are running away from the box, so modify the code not to calculate the escaped particles

performance comparison:
---
parameter setting for the performance comparison:  
~~~~
nparticles_fluid = (/8, 16,16/)  
nparticles_box   = (/32,16,16/)  
dist_particles   = 7.5d-2  
tol = 1.0d-10  
dt  = 3.0d-3  
tstep_max  = 20
freq_write = 20  
~~~~
Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz, 28 cores, with AVX2  
wtih 1 core  : 214.4s  
with 28 cores: 29.3s  
---
a calculation result(water collapse in a box) is below.  
in this case, size of the box is 2.4x1.2x1.2(m), size of the water in the initial state is 0.6x1.2x1.2(m), time is 0.0~6.0(s)(2000 steps).  
it took about 1 hour with 28 cores to calculate the process.  
![Alt text](./water_collapse.gif?raw=true "water collapse")
