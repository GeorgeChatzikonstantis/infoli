#BSUB -J georgecjob
#BSUB -o /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/InfOliFull/run/output/georgecjob.out
#BSUB -e /gpfs/stfc/local/HCEEC005/nnm17/gxc30-nnm17/InfOliFull/run/output/georgecjob.err
#BSUB -R "span[ptile=1]"
#BSUB -W 24:00
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

# library paths and necessary flags
export KMP_AFFINITY=balanced
export KMP_PLACE_THREADS=60c,4t
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

# root directory of InfOli dir
rootdir="~/InfOliFull"

# prepare experiment: clean up the working room
cd ~/InfOliFull/run
rm -rf input/*
rm -rf output/georgecjob.* output/*.txt

# prep input: compiling executable
cd ~/InfOliFull/src
make omp_phi
mv OpenMP/infoli.x ~/InfOliFull/run/input

# prep input: copy runtime lib to input
cd ~/InfOliFull/run/
cp runtime_libs/hartree/libiomp5.so input

# preparations complete, conduct the experiment
cd ~/InfOliFull/run/input
for size in 1000
do

	for synapses in 0 250
	do

		density=$(bc -l <<< "scale=6; $synapses/$size")
		export MYJOB="$rootdir/run/input/infoli.x $size $density 100"

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
				ssh ${MICNAME} "export \
				LD_LIBRARY_PATH=/gpfs/stfc/local/systems/sja01/ljj24-sja01/testjob/intel/clck/2.2.1.003/share/mic/:$LD_LIBRARY_PATH; \
				export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=60c,4t; export OMP_NUM_THREADS=${THREADSN}; \
				cd $rootdir/run/input; ls; ${MYJOB}"
			done
		done
	done
done

# the experiment is complete, clean up input and move results to output
mv vtune_report* ~/InfOliFull/run/output
#rm -rf ~/InfOliFull/run/input/*
