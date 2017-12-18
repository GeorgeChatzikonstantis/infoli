#!/bin/bash

# Get the Directories
BASE_DIR=$(dirname $SCRIPT_DIR)
BIN_DIR="$BASE_DIR/bin"

echo $BIN_DIR 

# Print experiment parameters
echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Ranks= ${RANKS}|Working Dir=$WORK_DIR"

# prepare experiment: 
cd $WORK_DIR

MYJOB="$BIN_DIR/infoli.x -n $SIZE -p ${PROB} -t ${STIME}"
export OMP_NUM_THREADS=$THREADSNUM
echo "mpiexec -genvall -np $RANKS $MYJOB"

/opt/intel/impi/2017.3.196/bin64/mpiexec -genvall -np $RANKS $MYJOB

echo "Simulation Finished, results in \"$WORK_DIR\" folder.\n"
echo "-----------------------------------\n"
