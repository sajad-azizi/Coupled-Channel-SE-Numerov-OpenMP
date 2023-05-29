# Coupled-Channel Schrödinger Equation

This code solves a coupled-channel Schrödinger equation in two dimensions, and it can simply be upgraded to full 3D.


Here we used Johnson's renormalized Numerov Method [1], which is a convenient method, to solve this equation and calculate the bound state accurately.
However, we modified the algorithm for the boundary condition, which makes it very accurate in caculating bound state energy. For more information regarding this method and the way we implemented it, the interested reader is referred to Appendix A of my PhD thesis.


This code calculates bound states and continuum states for an arbitrary anisotropic potential.










[1] B. R. Johnson; The renormalized Numerov method applied to calculating bound states of the coupled‐channel Schroedinger equation. J. Chem. Phys. 15 November 1978; 69 (10): 4678–4688. https://doi.org/10.1063/1.436421
