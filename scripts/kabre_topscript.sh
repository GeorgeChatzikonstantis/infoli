#!/bin/bash

echo "Runs a job with the infoli simulator from current Directory."
echo "It creates a new folder or a given folder outputs the results there"
echo "Usage: ./script_name.sh -dir <dir to put results (doesn't need"
echo "to exist)> -ns <network size> -pb <network density> -st <simulation time>"
echo "-np <number of MPI ranks to spawn> -nd <number of nodes to be used>"
echo "If no argument is given it runs on default values!"

# Get the Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CUR_DIR=$(pwd)
HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"

if [ "$CUR_DIR" == "$HOME" ]
then
    mkdir -p "$CUR_DIR/brainframe-res"
    CUR_DIR="$CUR_DIR/brainframe-res"
fi

echo $SCRIPT_DIR
echo $CUR_DIR
echo $NAME

# LSF queue system options

WALL="walltime=72:00:00"
BQUE=phi-n18h72
NODES=1
RANKS=4
THREADSNUM=50

# Experiment variables
SIZE=100
PROB=0.2
STIME=10


while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -dir|-working_dir)
            NAME="$2"
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
	-np|-ranks)
	    RANKS="$2"
	    shift # past argument
	    ;;
	-nd|-nodes)
	    NODES="$2"
	    shift # past argument
	    ;;

        *)
            # unknown option
            echo "Invalid argument: $1"
            exit 1
    esac
    shift # past argument or value
done

TOPOLOGY="nodes=$NODES:ppn=4"

DATE=`date +%T`
NAME="exp"_${SIZE}_${PROB}_$DATE

SYNAPSES=$(bc -l <<< "scale=0; $SIZE*$PROB")
JOB="exp"_${SIZE}_${SYNAPSES}
B_OUT=$CUR_DIR/$NAME/job.out
B_ERR=$CUR_DIR/$NAME/job.err

echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}"
WORK_DIR=$CUR_DIR/$NAME
mkdir -p $WORK_DIR

CMD="qsub -N $JOB -o $B_OUT -e $B_ERR -q $BQUE -l $TOPOLOGY -l $WALL -v SCRIPT_DIR=$SCRIPT_DIR,WORK_DIR=$WORK_DIR,SIZE=$SIZE,PROB=$PROB,STIME=$STIME,RANKS=$RANKS,THREADSNUM=$THREADSNUM $SCRIPT_DIR/kabre_jobscript.sh"
echo $CMD
echo $JOB $B_OUT $B_ERR

eval $CMD
