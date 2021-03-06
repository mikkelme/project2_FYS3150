# FYS 3150 - Project 2
This repository contains the work of Fredrik Hoftun and Mikkel Metzsch Jensen for project 2 in Computational Physics FYS3150. 

To use the compiler just type `make` in the commandline in the `/code` directory.
The compiler will generate a file `Eigen_solver.exe` that can be used in the command line to generate files for eigenvalues, eigenvectors, number of iterations, time used and the mean error.

Usage of `Eigen_solver.exe` in Bash: `$ ./Eigen_solver.exe N V rho_max Omega^2`, where N defines the matrix dimensionality `NxN`, V is the chosen potential where `V=1` gives the one electron potential, `V=2` gives the two electron potential with Coulomb interaction and any other value gives the buckling beam, no potential, problem. 

Example: `$ ./Eigen_solver.exe 200 2 18 1`
