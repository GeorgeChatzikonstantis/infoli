#!/bin/bash 

# set some script parameters
max_reps=1
simtime=1000

echo "Fixed Experiments Time to $simtime ms."

# set name of attached MIC
MICNAME="mic0"

# root directory of InfOli dir
rootdir="/home/georgec/InfOliFull"

# prepare experiment: clean up the working room
cd $rootdir/run
rm -rf input/*

# prep input: compiling executable
cd $rootdir/src
make omp_phi
mv OpenMP/infoli.x ../run/input

# prep input: copy runtime lib to input
cd $rootdir/run
cp runtime_libs/microlab/phaethon/libiomp5.so input
cp runtime_libs/power/lib* input

# copy necessities to the mic
cd input
export tag=$RANDOM
ssh $MICNAME "mkdir experiment$tag"
scp * $MICNAME:~/experiment$tag

# preparations complete, conduct the experiment
for (( reps=1; reps<=$max_reps; reps++ ))
do

echo "Repetition $reps/$max_reps Commencing."

	for size in 1000 2000 4000 6000 8000 10000
	do

	echo "NW Size = $size neurons."

		for synapses in 0 250 500 750 1000
		do

			echo "Synapses = $synapses"

			density=$(bc -l <<< "scale=6; $synapses/$size")
			export MYJOB="./infoli.x -n $size -p $density -t $simtime"

			for threads in 200
			do
				export THREADSN=$threads

				export vtuneCommand="amplxe-cl -c general-exploration \
				-result-dir vtune_report${threads}_${size}_${synapses} \
				-target-system=mic-native:0 -knob enable-vpu-metrics=true \
				-knob collect-memory-bandwidth=true -knob enable-tlb-metrics=true --"

#				export vtuneCommand="amplxe-cl -c memory-access \
#				-result-dir vtune_report${threads}_${size} -target-system=mic-native:0 -- "

#				$vtuneCommand "export \
#				LD_LIBRARY_PATH=/home/georgec/experiment$tag:$LD_LIBRARY_PATH; \
#				export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=57c,4t; \
#				export OMP_NUM_THREADS=${THREADSN}; cd ~/experiment$tag; ${MYJOB}"

				ssh ${MICNAME} "\
				export LD_LIBRARY_PATH=/home/georgec/experiment$tag:$LD_LIBRARY_PATH; \
				export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=57c,4t; \
				export OMP_NUM_THREADS=${THREADSN}; cd ~/experiment$tag; ${MYJOB}"

			done
		done
	done
done

echo "Done!"

# the experiment is complete, clean up input and mic-trash and move results to output
ssh $MICNAME "rm -rf ~/experiment$tag"
mv vtune_report* $rootdir/run/output
#rm -rf $rootdir/run/input/*
