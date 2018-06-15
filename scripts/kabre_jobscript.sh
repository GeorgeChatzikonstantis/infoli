#!/bin/bash

# Get the Directories
BASE_DIR=$(dirname $SCRIPT_DIR)
BIN_DIR="$BASE_DIR/bin"

echo $BIN_DIR 

echo ""
numactl --hardware
echo ""

# Print experiment parameters
echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Ranks= ${RANKS}|Working Dir=$WORK_DIR"

# prepare experiment: 
cd $WORK_DIR

#MYJOB="numactl --membind 1 $BIN_DIR/infoli.x -n $SIZE -p ${PROB} -t ${STIME}"
MYJOB="$BIN_DIR/infoli.x -n $SIZE -p ${PROB} -t ${STIME}"
export OMP_NUM_THREADS=$THREADSNUM
echo "mpiexec -genvall -np $RANKS $MYJOB"

/usr/bin/time -v /opt/intel/impi/2017.3.196/bin64/mpiexec -genvall -np $RANKS $MYJOB

echo "Simulation Finished, results in \"$WORK_DIR\" folder.\n"
echo "-----------------------------------\n"
