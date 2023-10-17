Comparative analysis of parallel computing standards (MPI, OpenMP, CUDA) for accelerating Poisson solver performance.
BC: u=0 on all boundaries.
Solving Poisson equation is the most time-consuming part of an incompressible flow solver.
The equation is solved in a three-dimensional scenario.
Parallelized using MPI_CART_CREATE in all 3 directions.
2nd order central-difference is used for discretization.
Gauss-Seidel with SOR is used for solving the linear system.
