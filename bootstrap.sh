#!/bin/bash

# determine environment for setup
hostname=$HOSTNAME

# set root directory of repository
root_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $hostname == "phaethon.microlab.ntua.gr" ];
then
	cd $root_dir/src
	make omp_phi
	cd -
	scp $root_dir/bin/infoli.x mic0:~
	scp $root_dir/bin/default.conf mic0:~
	scp $root_dir/lib/lib* mic0:~
else if [[ $hostname == "crb"* ]];
	then
		source /etc/profile.d/modules.sh
		module load intel_mpi/5.0.3_mic intel/15.2.164_mic intel_vtune
		bsub -q phiq -o $root_dir/bootstrap_log.out -e $root_dir/bootstrap_log.err "cd $root_dir/src; make omp_phi"
	fi
fi
