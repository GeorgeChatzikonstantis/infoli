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
rootdir=$(cd ..; pwd)

# prepare experiment: clean up the working room
cd ../run
rm -rf input/*

# prep input: compiling executable
cd ../src
make omp_xeon
mv OpenMP/infoli.x ../run/input

# prep input: copy runtime lib to input
cd ../run/
cp runtime_libs/hartree/libiomp5.so input

# specify how threads will be placed and topology
export KMP_AFFINITY=balanced
export KMP_PLACE_THREADS=20c,1t

# preparations complete, conduct the experiment
cd ../run/input
export MYJOB="./infoli.x $size $density $simtime"
export OMP_NUM_THREADS=$threads
/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" ${MYJOB}

# the experiment is complete, clean up input and move results to output
mv InferiorOlive_Output* ../output
rm -rf ../input/*
