*mps* - 3D incompressive fluid solver by using moving partigle semi-implicit method with Fortran.
======
moving particle semi-implicit(MPS) method: http://li.mit.edu/Stuff/CNSE/Paper/Koshizuka96Oka.pdf  
currently intel compiler is required to test it.  
  
how to run:
    
    $ vi mps.F90 # adjust the parameters  
    $ make  
    $ ./mps  
  
TODO:  
---
- implement ICCG method and replace CG method with it
- create a function to find neighber particles in the effective radius(using hash table is better?)
- continue the verification
- parallelization by OpenMP, MPI  

---
water collapse in a box.
![Alt text](./water_collapse.gif?raw=true "water collapse")
