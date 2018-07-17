#!/bin/bash

source /etc/profile.d/modules.sh
module load intel/18.0.128 intel_mpi/18.0.128

cd $HCBASE
export MYHOME=`pwd`
export OPENMPI_ROOT=/lustre/scafellpike/local/apps/gcc/openmpi/2.1.1
export PATH=$MPIROOT/bin:$PATH
export LD_LIBRARY_PATH=$MPIROOT/lib:$MPIROOT/lib/openmpi:$LD_LIBRARY_PATH

# Get the Directories
BASE_DIR=$(dirname $SCRIPT_DIR)
BIN_DIR="$BASE_DIR/bin"

echo $BIN_DIR 

# Print experiment parameters
echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Ranks= ${RANKS}|Threads= $THREADSNUM|Working Dir=$WORK_DIR"

# prepare experiment: 
cd $WORK_DIR

MYJOB="$BIN_DIR/infoli.x -n $SIZE -p ${PROB} -t ${STIME}"
export OMP_NUM_THREADS=$THREADSNUM

echo ""
echo "mpirun -genvall -np $RANKS $MYJOB"
echo ""

mpirun -genvall -np $RANKS $MYJOB

echo "Simulation Finished, results in \"$WORK_DIR\" folder.\n"
echo "-----------------------------------\n"
