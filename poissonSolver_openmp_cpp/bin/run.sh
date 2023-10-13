#!/bin/bash
echo "====================================="
echo "     Poisson Solver: OpenMP | C++    "
echo "====================================="
echo "Clean Rebuild"
(cd ../build && make clean)
(cd ../build && make)
echo "====================================="
export OMP_NUM_THREADS=6
# echo "Starting execution using 8 procs"
(cd ../build && time ./openmp_cpp.x)
# nohup ./source/fsi3l1.x &
# gnuplot script.gnu