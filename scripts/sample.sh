#!/bin/bash

echo "Runs a job with the infoli simulator from current Directory."
echo "It creates a new folder or a given folder outputs the results there"
echo ""
echo "Usage: ./sample.sh -dir <dir to put results in, leave blank for current dir>"
echo "-pl <platform selection: x86 or phi> -bld <build: omp or hybrid>"
echo "-ns <network size> -pb <network density> -st <simulation time>"
echo "-np <number of MPI Ranks (if hybrid build)> -th <threads num>"
echo ""
echo "If no argument is given it runs on default values:"
echo "An OpenMP x86 build is chosen, 100 neurons with 20% density"
echo "is simulated for 10ms on 4 threads."

# Get the Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CUR_DIR=$(pwd)
BASE_DIR=$(dirname $SCRIPT_DIR)
SRC_DIR=$BASE_DIR/src
BIN_DIR=$BASE_DIR/bin
HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"

if [ "$CUR_DIR" == "$HOME" ]
then
    mkdir -p "$CUR_DIR/brainframe-results"
    CUR_DIR="$CUR_DIR/brainframe-results"
fi

DATE=`date +%T`
DNAME="testRun"_$DATE

# Thread options
MICNAME="mic0"
THREADSNUM=4
RANKSNUM=2
PLATFORM="x86"
BUILD="omp"

# Experiment variables
SIZE=100
PROB=0.2
STIME=10


while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -dir|-working_dir)
            DNAME="$2"
            shift # past argument
            ;;
	-pl|-platform)
            PLATFORM="$2"
            shift # past argument
            ;;
	-bld|-build)
            BUILD="$2"
            shift # past argument
            ;;
        -ns|-net_size)
            SIZE="$2"
            shift # past argument
            ;;
        -pb|-probability)
            PROB="$2"
            shift # past argument
            ;;
        -st|-sim_time)
            STIME="$2"
            shift # past argument
            ;;
        -th|-threads)
            THREADSNUM="$2"
            shift # past argument
            ;;
	-np|-n)
	    RANKSNUM="$2"
	    shift # past argument
	    ;;

        *)
            # unknown option
            echo "Invalid argument: $1"
            exit 1
    esac
    shift # past argument or value
done

# make for the correct target
if [[ "$PLATFORM" == "x86" ]]; then
	if [[ "$BUILD" == "omp" ]]; then
		cd $SRC_DIR
		make omp_xeon || exit 1
		cd -
	elif [[ "$BUILD" == "hybrid" ]]; then
		cd $SRC_DIR
		make hybrid_xeon || exit 1
		cd -
	else
		echo "Unable to make for target.Please review build options.\n"
		exit 1
	fi
elif [[ "$PLATFORM" == "phi" ]]; then
	if [[ "$BUILD" == "omp" ]]; then
		cd $SRC_DIR
		make omp_phi || exit 1
		cd -
	elif [[ "$BUILD" == "hybrid" ]]; then
		cd $SRC_DIR
		make hybrid_phi || exit 1
		cd -
	else
		echo "Unable to make for target.Please review build options.\n"
		exit 1
	fi
else
	echo "Unable to make for target.Please review platform options.\n"
	exit 1
fi

echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Threads=$THREADSNUM DNameP:$DNAME"
mkdir -p $CUR_DIR/$DNAME

if [[ "$BUILD" == "omp" ]]; then
	MYJOB="$BIN_DIR/infoli.x -n ${SIZE} -p ${PROB} -t ${STIME}"
elif [[ "$BUILD" == "hybrid" ]]; then
	MYJOB="mpirun -np ${RANKSNUM} $BIN_DIR/infoli.x -n ${SIZE} -p ${PROB} -t ${STIME}"
fi
echo $MYJOB

if [[ "$PLATFORM" == "phi" ]]; then
	echo ""
	echo "Support for KNC executions through this script is pending."
	echo "Please review scripts/microlab_top_level.sh for similar functionality."
	exit 1
fi

export KMP_AFFINITY=balanced;
export OMP_NUM_THREADS=$THREADSNUM;
echo ""
echo "Executing for $BUILD build on $PLATFORM ..."
$MYJOB
echo "Moving Results ..."
mv InferiorOlive_Output*.txt $CUR_DIR/${DNAME}

echo "Simulation Finished, results in \"$CUR_DIR/$DNAME/\" folder.\n"
echo "-----------------------------------\n"

echo "Everything ok!"
