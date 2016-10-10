#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output/georgecjob.out
#BSUB -e /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output/georgecjob.err
#BSUB -R "span[ptile=1]"
#BSUB -W 6:00
#BSUB -n 1
#BSUB -q phiq

# prepare Hartree's modules
source /etc/profile.d/modules.sh
module purge
module load intel_mpi/5.0.3_mic intel/15.2.164_mic intel_vtune

# get name of attached MIC
string=$HOSTNAME
MICNAME=${string/ib0/mic0}
echo "$MICNAME"

# prepare experiment: clean up the working room
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run
rm -rf input/*

# prep input: compiling executable
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/src
make omp_phi
mv OpenMP/infoli.x /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input

# prep input: compile desired connectivity generator
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/tools
make count
mv connectivity/3d_synapse_count/conn_generator.x /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input

# prep input: copy runtime lib to input
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/
cp runtime_libs/hartree/libiomp5.so input

# preparations complete, conduct the experiment
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input
for size in 1000000
do

	for synapses in 100
	do

		./conn_generator.x $((size/100)) 10 10 $synapses 1
		export MYJOB="/gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input/infoli.x $size"

		for times in 1 2
		do
			for threads in 200
			do
				export THREADSN=$threads
				export vtuneCommand="amplxe-cl -c general-exploration -result-dir vtune_report${size}_${synapses}_${times} \
				-target-system=mic-native:0 -knob enable-vpu-metrics=true \
				--target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#				$vtuneCommand "export \
#				LD_LIBRARY_PATH=/gpfs/stfc/local/systems/sja01/ljj24-sja01/testjob/intel/clck/2.2.1.003/share/mic/:$LD_LIBRARY_PATH; \
#				export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=60c,4t; export OMP_NUM_THREADS=${THREADSN}; \
#				cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/infoli_g/workspace/input; ${MYJOB}"
				/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" ssh ${MICNAME} "export \
				LD_LIBRARY_PATH=/gpfs/stfc/local/systems/sja01/ljj24-sja01/testjob/intel/clck/2.2.1.003/share/mic/:$LD_LIBRARY_PATH; \
				export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=60c,4t; export OMP_NUM_THREADS=${THREADSN}; \
				cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input; ${MYJOB}"
			done
		done
	done
done

# the experiment is complete, clean up input and move results to output
mv vtune_report* /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output
rm -rf /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input/*
