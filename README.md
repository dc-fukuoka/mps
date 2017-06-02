*mps* - 3D incompressive fluid solver by using moving partigle semi-implicit method with Fortran. still in progress...
======
implementing http://li.mit.edu/Stuff/CNSE/Paper/Koshizuka96Oka.pdf  
  
TODO:  
---
- implement ICCG method and replace CG method with it
- create a function to find neighber particles in the effective radius(using hash table is better?)
- continue the verification
- parallelization by OpenMP, MPI  

![Alt text](./water_collapse.gif?raw=true "water collapse")
