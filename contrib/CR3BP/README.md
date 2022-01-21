# Circular Restricted Three-Body Problem (CR3BP) Package

## Currently Implemented Functionality
- natural dynamics propagation
- dynamics + State-Transition-Matrix propagation
- compute characteristic values of CR3BP

## Functionality Left To Implement
- CR3BP rotating to inertial frame convertion
- inertial to CR3BP rotating frame convertion
- convertion to and from Ephemeris frame with SPICE API
- differential correction methods (single/multiple shooting)
- compute periodic orbits with initial guess
- analytical approximate solutions (first order/third order) for orbits like Halo and Lyapanov
- analytical computation of eigenvalues / eigenvectors of Monodromy matrix and stability of orbits
- stable and unstable manifold estimation

## Functionality I Required Help With
- Integration of code with astropy constants and units without breaking JIT Compilation
- Integration of ploting with existing poliastro interactive plotter (currently in the test file I use matplotlib and mplot3d)
