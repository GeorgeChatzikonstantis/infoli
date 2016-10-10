# set name of attached MIC
MICNAME="mic0"

# prepare experiment: clean up the working room
cd ~/infoli/run
rm -rf input/*

# prep input: compiling executable
cd ~/infoli/src
make omp_phi
mv OpenMP/infoli.x ~/infoli/run/input

# prep input: compile desired connectivity generator
cd ~/infoli/tools
make count
mv connectivity/3d_synapse_count/conn_generator.x ~/infoli/run/input

# prep input: copy runtime lib to input
cd ~/infoli/run
cp runtime_libs/phaethon/libiomp5.so input

# copy necessities to the mic
export tag=$RANDOM
ssh $MICNAME "mkdir experiment$tag"
cd ~/infoli/run/input
scp infoli.x libiomp5.so $MICNAME:~/experiment$tag

# preparations complete, conduct the experiment
cd ~/infoli/workspace/input

for size in 1000 2000 5000 10000
do
	for synapses in 10 20 50 100 200 500
	do

		export MYJOB="./infoli.x $size"

		# generate connectivity map and move it to the mic
		./conn_generator.x $((size/100)) 10 10 $synapses 1
		scp cellConnections.txt $MICNAME:~/experiment$tag

		for threads in 50 100 200
		do
			export THREADSN=$threads

			export vtuneCommand="amplxe-cl -c general-exploration \
			-result-dir vtune_report${threads}_${size}_${synapses} \
			-target-system=mic-native:0 -knob enable-vpu-metrics=true \
			-knob collect-memory-bandwidth=true -knob enable-tlb-metrics=true --"

#			export vtuneCommand="amplxe-cl -c memory-access \
#			-result-dir vtune_report${threads}_${size} -target-system=mic-native:0 -- "

			$vtuneCommand "export \
			LD_LIBRARY_PATH=/home/georgec/experiment$tag:$LD_LIBRARY_PATH; \
			export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=57c,4t; \
			export OMP_NUM_THREADS=${THREADSN}; cd ~/experiment$tag; ${MYJOB}"

#			/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" ssh ${MICNAME} "\
#			export LD_LIBRARY_PATH=/home/georgec/experiment$tag:$LD_LIBRARY_PATH; \
#			export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=57c,4t; \
#			export OMP_NUM_THREADS=${THREADSN}; cd ~/experiment$tag; ${MYJOB}"

		done
	done
done

# the experiment is complete, clean up input and mic-trash and move results to output
ssh $MICNAME "rm -rf ~/experiment$tag"
mv vtune_report* ~/infoli/run/output
rm -rf ~/infoli/run/input/*
