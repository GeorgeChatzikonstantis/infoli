#!/bin/bash

# library paths and necessary flags
export I_MPI_FABRICS=shm:tcp
export OFFLOAD_INIT=on_start
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

MODEL=$1
for size in 100
do
	cp ./VersionVerification/cellConnections.txt cellConnections.txt
	cp ./VersionVerification/gcal_file.txt gcal_file.txt
	cp ./VersionVerification/libiomp5.so libiomp5.so
	#run the stable version
	./VersionVerification/stable_version.x $size
	mv InferiorOlive_Output.txt correct_output.txt

	#check the MPI implementation
	if [ "$MODEL" = "mpi_xeon" ]; then 
		#run the new version
		mpirun -np 1 -genvall ./MPI/infoli.x $size
		cp InferiorOlive_Output0.txt new_output.txt 
		#compare the outputs
		diff --ignore-all-space new_output.txt correct_output.txt > differences.txt
		if [ -s differences.txt ]; then
			echo "The new version is not consistent with the stable version."
		else 
			echo "The new version is consistent with the stable version."
		fi
		rm InferiorOlive_Output0.txt
	else
		#check the OpenMP implementation
		if [ "$MODEL" = "omp_xeon" ]; then 
			export KMP_AFFINITY=balanced
			export KMP_PLACE_THREADS=20c,1t
			#run the new version
                	export OMP_NUM_THREADS=10
			./OpenMP/infoli.x $size

			cp InferiorOlive_Output.txt new_output.txt  
			#compare the outputs
			diff new_output.txt correct_output.txt > differences.txt
			if [ -s differences.txt ]; then
				echo "The new version is not consistent with the stable version."
			else 
				echo "The new version is consistent with the stable version."
			fi
			rm InferiorOlive_Output.txt
		else
			#check the hybrid implementation
			export I_MPI_PIN_DOMAIN=omp
			export KMP_AFFINITY=balanced
			export KMP_PLACE_THREADS=24c,1t

			#run the new version
			export OMP_NUM_THREADS=20
			mpirun -np 1 -genvall ./Hybrid/infoli.x $size
			cp InferiorOlive_Output0.txt new_output.txt 

			#compare the outputs
			diff new_output.txt correct_output.txt > differences.txt
			if [ -s differences.txt ]; then
				echo "The new version is not consistent with the stable version."
			else 
				echo "The new version is consistent with the stable version."
			fi
			rm InferiorOlive_Output0.txt
		fi
	fi

	cp correct_output.txt ./VersionVerification/correct_output.txt
	cp new_output.txt ./VersionVerification/new_output.txt
	cp differences.txt ./VersionVerification/differences.txt
	rm correct_output.txt new_output.txt differences.txt
	rm cellConnections.txt gcal_file.txt libiomp5.so
done
