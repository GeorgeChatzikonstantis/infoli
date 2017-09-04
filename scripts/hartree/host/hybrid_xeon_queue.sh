#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/InfOliFull/run/output/georgecjob.out
#BSUB -e /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/InfOliFull/run/output/georgecjob.err
#BSUB -R "span[ptile=2]"
#BSUB -W 24:00
#BSUB -n 10
#BSUB -q phiq

# number of ranks used
ranks=10

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

# root directory of InfOli dir
rootdir=$(cd ../../../; pwd)

# prepare experiment: clean up the working room
cd $rootdir/run/
rm -rf input/*

# prep input: compiling executable
cd $rootdir/src
make hybrid_xeon
mv Hybrid_New/infoli.x $rootdir/run/input

# prep input: copy runtime lib to input
cd $rootdir/run/ 
cp runtime_libs/hartree/libiomp5.so input

# specify vtune analysis type
export analysis_type="general-exploration"

# preparations complete, conduct the experiment
cd $rootdir/run/input
ls
for size in 50000
do
	for pct in 0.1
	do
		for threads in 10
		do

			export MYJOB="$rootdir/run/input/infoli.x $size $pct 100"
			export OMP_NUM_THREADS=$threads
			export vtuneCommand="amplxe-cl -c $analysis_type -result-dir vtune_report10_$size --target-install-dir=/gpfs/stfc/local/apps/intel/intel_cs/vtune_amplifier_xe_2015 --"
#			$vtuneCommand mpirun -np 10 -genvall ${MYJOB}
			/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np $ranks -genvall ${MYJOB}
		done
	done
done

# the experiment is complete, clean up input and move results to output
mv vtune_report* $rootdir/run/output
#rm -rf $rootidr/run/input/*
