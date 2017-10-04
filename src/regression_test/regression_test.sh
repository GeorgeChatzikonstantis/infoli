#!/bin/bash

echo "Runs a regression test on the current build of infoli simulator source code."
echo "Based on pre-stored verified simulation output found in scripts/regression_test."
echo "Usage: /path/to/regression_test.sh -src <omp or hybrid> -str <standard or full>"
echo "If no argument is given it runs a default standard (no stimulus) test on openmp."

# Get the Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR=$SCRIPT_DIR/../..
SRC_DIR=$REPO_DIR/src
BIN_DIR=$REPO_DIR/bin
CUR_DIR=$(pwd)
HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"

BUILD="omp"
STRENGTH="standard"

while [[ $# -gt 1 ]]
do
	key="$1"

	case $key in
		-src|-build)
			BUILD="$2"
			shift
			;;

		-str|-version)
			STRENGTH="$2"
			shift
			;;

		*)
		# unknown option
			echo "Invalid argument: $1"
			exit 1
	esac
	shift
done

cd $SRC_DIR
if [ "$BUILD" == 'omp' ]; then
	make omp_test
else
	if [ "$BUILD" == 'hybrid' ]; then
		make hybrid_test
	fi
fi

cd $SCRIPT_DIR
make
cp checker.x $BIN_DIR
make clean

cd $BIN_DIR
cp $SCRIPT_DIR/InferiorOlive_Output_Baseline.txt .
cp $SCRIPT_DIR/gcal_file.txt .
export OMP_NUM_THREADS=1

neurons=10
density=1
if [ "$STRENGTH" == 'standard' ]; then
simtime=500
else
	if [ "$STRENGTH" == 'full' ]; then
		simtime=2000
	fi
fi

if [ "$BUILD" == 'omp' ]; then
	./infoli.x $neurons $density $simtime 
else
	if [ "$BUILD" == 'hybrid' ]; then
		mpirun -np 1 ./infoli.x $neurons $density $simtime
		mv InferiorOlive_Output0.txt InferiorOlive_Output.txt
	fi
fi

./checker.x $neurons

rm -f InferiorOlive_Output.txt InferiorOlive_Output_Baseline.txt
rm -f gcal_file.txt infoli.x checker.x
cd $CUR_DIR
