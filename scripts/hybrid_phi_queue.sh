#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/job.out
#BSUB -e /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/job.err
#BSUB -R "span[ptile=1]"
#BSUB -W 6:00
#BSUB -n 1
#BSUB -q phiq

# number of ranks used
ranks=1
ranks_per_node=1

# prepare Hartree's modules
source /etc/profile.d/modules.sh
module purge
module load intel_mpi/5.1.3_mic intel/16.2.062_mic intel_vtune

# library paths and necessary flags
export I_MPI_FABRICS=shm:tcp
export OFFLOAD_INIT=on_start
export I_MPI_PIN_DOMAIN=omp
export KMP_AFFINITY=balanced
export KMP_PLACE_THREADS=60c,4t
export I_MPI_MIC=enable
export LD_LIBRARY_PATH=/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/compiler/lib/mic:\
/gpfs/stfc/local/apps/intel/intel_cs/2015.1.133/composer_xe_2015.1.133/compiler/lib/intel64:\
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
export SINK_LD_LIBRARY_PATH=/gpfs/stfc/local/apps/intel/intel_mpi/5.1.3.181/mic/lib:$MIC_LD_LIBRARY_PATH;

# set directories
rootdir=$(cd ../; pwd)
bindir=$(cd ../bin; pwd)

# get name of attached MIC
string=$HOSTNAME
MICNAME=${string/ib0/mic0}
echo "$MICNAME"

# prep input: compiling executable
cd $rootdir/src
make hybrid_phi
micnativeloadex $bindir/infoli.x -l

# preparations complete, conduct the experiment
cd $HOME
for size in 1000
do

	for pct in 0.1
	do

		export MYJOB="$bindir/infoli.x -n $size -p $pct -t 10"
		export OMP_NUM_THREADS=$(bc -l <<< "scale=0; 200/$ranks_per_node")

		micnativeloadex $bindir/infoli.x -a "1000 0.1 10" -e "export OMP_NUM_THREADS=200"
#		/usr/bin/time -f "Total Time:\t%E\tMem Usage:\t%MkB" mpirun -envall -n $ranks -host $MICNAME $MYJOB

	done
done
