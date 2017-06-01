*mps* - 3D incompressive fluid solver by using moving partigle semi-implicit method with Fortran. still in progress...
======
implementing http://li.mit.edu/Stuff/CNSE/Paper/Koshizuka96Oka.pdf  
  
TODO:  
---
- still some particles are escaping from the box, find the reason  
- do some more tests  
- parallelization by OpenMP, MPI  

current result is the following:  
still the behavior is not natural...
  
![Alt text](./water_collapse.gif?raw=true "water collapse")
