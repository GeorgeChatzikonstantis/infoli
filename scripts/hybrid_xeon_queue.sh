#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/job.out
#BSUB -e /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/job.err
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

# set directories
rootdir=$(cd ../; pwd)
bindir=$(cd ../bin; pwd)

# prep input: compiling executable
cd $rootdir/src
make hybrid_xeon

# preparations complete, conduct the experiment
cd $HOME
for size in 50000
do
	for pct in 0.1
	do
		for threads in 10
		do

			export MYJOB="$bindir/infoli.x $size $pct 10"
			export OMP_NUM_THREADS=$threads

			/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -np $ranks -genvall ${MYJOB}
		done
	done
done
