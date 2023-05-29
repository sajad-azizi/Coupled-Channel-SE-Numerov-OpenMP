# Coupled-Channel Schrödinger Equation

This code solves a Coupled-Channel Schrödinger Equation in 2 dimenstional and it is simply can be upgraded to a full 3D.


Here we used Johnson renormalized Numerov Method [1] , which is a convinient method, to solve this equation and caculate bound state accuratly.
However, we modified the algorithm for the boundary condition which makes it very accurate in caculating bound state energy. For more information regarding this method and the way we implemented it, the interested reader referes to the Appendix A of my PhD thesis.







[1] B. R. Johnson; The renormalized Numerov method applied to calculating bound states of the coupled‐channel Schroedinger equation. J. Chem. Phys. 15 November 1978; 69 (10): 4678–4688. https://doi.org/10.1063/1.436421
