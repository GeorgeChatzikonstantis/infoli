#!/bin/bash

echo "Runs a job with the infoli simulator from current Directory."
echo "It creates a new folder or a given folder outputs the results there"
echo "Usage: ./microlab_top_level.sh -dir <dir to put results (doesn't need
to exists)> -ns <network size> -pb <network density> -st <simulation time> -th <threads num>\n"
echo "If no argument is given it runs on default values!"

# Get the Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CUR_DIR=$(pwd)
BASE_DIR=$(dirname $SCRIPT_DIR)
HOME="$(getent passwd $USER | awk -F ':' '{print $6}')"

if [ "$CUR_DIR" == "$HOME" ]
then
    mkdir -p "$CUR_DIR/brainframe-reses"
    CUR_DIR="$CUR_DIR/brainframe-reses"
fi

DATE=`date +%T`
DNAME="phiH"_$DATE

echo $SCRIPT_DIR
echo $CUR_DIR
echo $DNAME

# Thread options
MICNAME="mic0"
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
            DNAME="$2"
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

        *)
            # unknown option
            echo "Invalid argument: $1"
            exit 1
    esac
    shift # past argument or value
done

echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Threads=$THREADSNUM DNameP:$DNAME"
mkdir -p $CUR_DIR/$DNAME

MYJOB="../infoli.x ${SIZE} ${PROB} ${STIME}"

ssh $MICNAME "mkdir -p $DNAME"


ssh ${MICNAME} "export LD_LIBRARY_PATH=~:$LD_LIBRARY_PATH; \
    export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=57c,4t;\
    export OMP_NUM_THREADS=$THREADSNUM;cd ${DNAME};${MYJOB}; \
    scp InferiorOlive_Output.txt $USER@phaethon.microlab.ntua.gr:$CUR_DIR/${DNAME}"


echo "Simulation Finished, results in \"$CUR_DIR/$DNAME/\" folder.\n"
echo "-----------------------------------\n"
# the experiment is complete, clean up input and mic-trash and move results to output
echo "Cleaning shared memory."
ssh $MICNAME "rm -rf ~/$DNAME"
echo "Everything ok!"

