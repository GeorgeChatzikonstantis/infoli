#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/InfOliFull/run/output/georgecjob.out
#BSUB -e /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/InfOliFull/run/output/georgecjob.err
#BSUB -R "span[ptile=24]"
#BSUB -W 24:00
#BSUB -n 1
#BSUB -q phiq

# prepare Hartree's modules
source /etc/profile.d/modules.sh
module purge
module load intel_mpi/5.0.3 intel/15.2.164 intel_vtune

# prepare experiment: clean up the working room
cd ~/InfOliFull/run
rm -rf input/*

# prep input: compiling executable
cd ~/InfOliFull/src
make omp_xeon
mv OpenMP/infoli.x ~/InfOliFull/run/input

# prep input: compile desired connectivity generator
cd ~/InfOliFull/tools
make count
mv connectivity/3d_synapse_count/conn_generator.x ~/InfOliFull/run/input

# prep input: copy runtime lib to input
cd ~/InfOliFull/run/
cp runtime_libs/hartree/libiomp5.so input

# specify how threads will be placed and topology
export KMP_AFFINITY=balanced
export KMP_PLACE_THREADS=20c,1t

# preparations complete, specify how threads will be placed and conduct the experiment
cd ~/InfOliFull/run/input
for size in 500000
do

	for synapses in 500
	do

		./conn_generator.x $((size/100)) 10 10 $synapses 1
		export MYJOB="/gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input/infoli.x $size"

		for reps in 1 2 3 4
		do

			for threads in 20
			do
				export OMP_NUM_THREADS=$threads
#				export vtuneCommand="amplxe-cl -c general-exploration \
#				-result-dir vtune_report_xeon${size}_${synapses}_${reps} \
#				--target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#				$vtuneCommand ${MYJOB}
#				echo "$size $synapses $reps"
				/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" ${MYJOB}
			done
		done
	done
done

# the experiment is complete, clean up input and move results to output
mv vtune_report* /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output
rm -rf /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input/*
