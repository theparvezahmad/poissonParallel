#!/bin/bash
echo "====================================="
echo "  Poisson Solver: Serial | Fortran   "
echo "====================================="
echo "Clean Rebuild"
(cd ../build && make clean)
(cd ../build && make)
echo "====================================="
# echo "Starting execution using 8 procs"
(cd ../build && ./mpi.x)
# nohup ./source/fsi3l1.x &
# gnuplot script.gnu