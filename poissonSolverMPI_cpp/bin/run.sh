#!/bin/bash
echo "Clean Rebuild"
(cd ../build && make clean)
(cd ../build && make)
echo "====================================="
echo "Starting execution using 8 procs"
(cd ../build && mpirun -np 8 ./mpi.x)
# nohup ./source/fsi3l1.x &
# gnuplot script.gnu