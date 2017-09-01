#!/bin/bash

echo "Runs a job with the infoli simulator from current Directory."
echo "It creates a new folder or a given folder outputs the results there"
echo "Usage: ./hartree_top_level.sh -dir <dir to put results (doesn't need
to exists)> -ns <network size> -pb <network density> -st <simulation time>
-req <resource requirements>\n"
echo "If no argument is given it runs on default values!"

# Get the Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CUR_DIR=$(pwd)
DATE=`date +%T`
NAME="phiH"_$DATE

echo $SCRIPT_DIR
echo $CUR_DIR
echo $NAME

# LSF queue system options

WALL="6:00"
REQ="hname=idb1a03"
NUM_TASKS=1
BQUE=phiq

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
        -req|-request)
            REQ="$2"
            shift # past argument
            ;;

        *)
            # unknown option
            echo "Invalid argument: $1"
            exit 1
    esac
    shift # past argument or value
done

JOB=$NAME"-job"
B_OUT=$CUR_DIR/$NAME/job.out
B_ERR=$CUR_DIR/$NAME/job.err

echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Request=$REQ"
mkdir -p $CUR_DIR/$NAME



CMD="/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub -J $JOB -o $B_OUT -e $B_ERR -R \"$REQ\" -W $WALL -n $NUM_TASKS -q $BQUE \"$SCRIPT_DIR/hartree_queue_jobscript.sh -dir $CUR_DIR/$NAME -ns $SIZE -pb $PROB -st $STIME\""
echo $CMD
#echo $JOB $B_OUT $B_ERR
export LSF_ENVDIR=/opt/lsf/conf
export LSF_SERVERDIR=/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/etc

eval $CMD

#Ta parakato tha ginontai pio prin


