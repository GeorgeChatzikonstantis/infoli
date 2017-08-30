#!/bin/bash

# Get the Directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR=$(dirname $SCRIPT_DIR)

echo $BASE_DIR 

# get name of attached MIC
string=$HOSTNAME
MICNAME=${string/ib0/mic0}
echo "$MICNAME"

# Get experiment parameters
THREADSNUM=200
WORK_DIR=$(pwd)

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
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
        -dir|-working_dir)
            WORK_DIR="$2"
            shift # past argument
            ;;

        *)
            # unknown option
            echo "Invalid argument: $1"
            exit 1
    esac
    shift # past argument or value
done
echo "SIZE= ${SIZE}|Density= ${PROB}|Simtime= ${STIME}|Working Dir=$WORK_DIR"


# prepare experiment: 
cd $WORK_DIR
MYJOB="$BASE_DIR/bin/infoli.x $SIZE ${PROB} ${STIME}"
echo $MYJOB

for threads in $THREADSNUM 
do
    export THREADSN=$threads

				ssh ${MICNAME} "export \
				LD_LIBRARY_PATH=/gpfs/stfc/local/systems/sja01/ljj24-sja01/testjob/intel/clck/2.2.1.003/share/mic/:$LD_LIBRARY_PATH; \
				export KMP_AFFINITY=balanced; export KMP_PLACE_THREADS=60c,4t; export OMP_NUM_THREADS=${THREADSN}; \
				cd $WORK_DIR; ${MYJOB}"
done

echo "Simulation Finished, results in \"$WORK_DIR\" folder.\n"
echo "-----------------------------------\n"
