*mps* - 3D incompressive fluid dynamics solver by using moving particle semi-implicit(MPS) method with Fortran.  
        the code is parallelized by OpenMP.
======
moving particle semi-implicit(MPS) method: http://li.mit.edu/Stuff/CNSE/Paper/Koshizuka96Oka.pdf  
currently intel compiler is required to build it.  
  
how to run:
    
    $ vi mps.F90 # adjust the parameters  
    $ make  
    $ ./mps
    $ ./create_anime.sh
  
TODO:  
---
- implement ICCG method and replace CG method with it
- create a function to find neighber particles in the effective radius(using hash table is better?)
- continue the verification
- parallelize the code by MPI  

performance comparison:
---


---
a calculation result(water collapse in a box) is below.  
in this case, size of the box is 2.4x1.2x1.2(m), size of the water in the initial state is 1.2x1.2x1.2(m), time is 0~1.2(s)(400 steps).  
![Alt text](./water_collapse.gif?raw=true "water collapse")
