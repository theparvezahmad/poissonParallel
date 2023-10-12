#!/bin/bash
echo "====================================="
echo "   Poisson Solver: Serial | Python   "
echo "====================================="
# echo "Starting execution using 8 procs"
(cd ../src && python main_serial.py)
# nohup ./source/fsi3l1.x &
# gnuplot script.gnu