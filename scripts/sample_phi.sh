if [ $# -ne 4 ]
  then
    echo "Script Usage: <script_name> <Network_Size>\
 <Network_Density ([0,1])> <Simulation_Time (ms)> <#Threads_Used>"
  exit 1
fi

size=$1
density=$2
simtime=$3
threads=$4

# root directory of InfOli dir
rootdir="~/InfOliFull"

# set name of attached MIC
MICNAME="mic0"

# prepare experiment: clean up the working room
cd ../run
rm -rf input/*

# prep input: compiling executable
cd ../src
make omp_phi
mv OpenMP/infoli.x ../run/input

# prep input: copy runtime lib to input
cd ../run/
cp runtime_libs/hartree/libiomp5.so input

# preparations complete, specify how threads will be placed and conduct the experiment
cd ../run/input
export MYJOB="./infoli.x $size $density $simtime"
/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" ssh ${MICNAME} "\
export LD_LIBRARY_PATH=~/InfOliFull/run/input:$LD_LIBRARY_PATH; \
export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=57c,4t; \
export OMP_NUM_THREADS=${threads}; cd $rootdir/run/input/; ${MYJOB}"

# the experiment is complete, clean up input and move results to output
mv InferiorOlive_Output* ../output
rm -rf ../input/*
