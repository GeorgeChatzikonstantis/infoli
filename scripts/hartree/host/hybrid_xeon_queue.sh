#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output/georgecjob.out
#BSUB -e /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output/georgecjob.err
#BSUB -R "span[ptile=24]"
#BSUB -W 6:00
#BSUB -n 24
#BSUB -q phiq

# prepare Hartree's modules
source /etc/profile.d/modules.sh
module purge
module load intel_mpi/5.0.3 intel/15.2.164 intel_vtune

# library paths and necessary flags
export I_MPI_FABRICS=shm:tcp
export OFFLOAD_INIT=on_start
export I_MPI_PIN_DOMAIN=omp
export KMP_AFFINITY=balanced
export KMP_PLACE_THREADS=24c,1t
export LD_LIBRARY_PATH=/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/compiler/lib/intel64:\
/gpfs/stfc/local/apps/intel/intel_mpi/5.0.2.044/intel64/lib:\
/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:\
/gpfs/stfc/local/apps/intel/intel_mpi/4.1.3.049/intel64/lib:\
/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/mpirt/lib/intel64:\
/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/ipp/../compiler/lib/intel64:\
/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/ipp/lib/intel64:\
/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/mkl/lib/intel64:\
/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/tbb/lib/intel64/gcc4.4:\
/lib:/usr/lib:/usr/X11R6/lib:/usr/local/lib:/usr/local/lib:\
/gpfs/stfc/local/apps/intel/intel_mpi/5.0.2.044/intel64/lib:$LD_LIBRARY_PATH;

# prepare experiment: clean up the working room
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/
rm -rf input/*

# prep input: compiling executable
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/src
make hybrid_xeon
mv Hybrid/infoli.x /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input

# prep input: compile desired connectivity generator
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/tools
make count
mv connectivity/3d_synapse_count/conn_generator.x /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input

# prep input: copy runtime lib to input
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/ 
cp runtime_libs/hartree/libiomp5.so input

# specify vtune analysis type
export analysis_type="general-exploration"

# preparations complete, conduct the experiment
cd /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input
for size in 10000
do

	for synapses in 10 20 50 100 200 500 1000
	do

		./conn_generator.x $((size/100)) 10 10 $synapses 1
		export MYJOB="/gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input/infoli.x $size"

		export OMP_NUM_THREADS=20
		export vtuneCommand="amplxe-cl -c $analysis_type -result-dir vtune_report1_$size --target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#		$vtuneCommand mpirun -np 1 -genvall ${MYJOB}
		/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np 1 -genvall ${MYJOB}

		export OMP_NUM_THREADS=10
		export vtuneCommand="amplxe-cl -c $analysis_type -result-dir vtune_report2_$size --target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#		$vtuneCommand mpirun -np 2 -genvall ${MYJOB}
		/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np 2 -genvall ${MYJOB}

		export OMP_NUM_THREADS=4
		export vtuneCommand="amplxe-cl -c $analysis_type -result-dir vtune_report5_$size --target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#		$vtuneCommand mpirun -np 5 -genvall ${MYJOB}
		/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np 5 -genvall ${MYJOB}

		export OMP_NUM_THREADS=2
		export vtuneCommand="amplxe-cl -c $analysis_type -result-dir vtune_report10_$size --target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#		$vtuneCommand mpirun -np 10 -genvall ${MYJOB}
		/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np 10 -genvall ${MYJOB}

		export OMP_NUM_THREADS=1
		export vtuneCommand="amplxe-cl -c $analysis_type -result-dir vtune_report20_$size --target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#		$vtuneCommand mpirun -np 20 -genvall ${MYJOB}
		/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np 20 -genvall ${MYJOB}

	done
done

# the experiment is complete, clean up input and move results to output
mv vtune_report* /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/output
rm -rf /gpfs/stfc/local/HCPhi005/ddr01/gxc30-ddr01/InfOliFull/run/input/*
