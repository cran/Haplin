#!/bin/bash

#------------------------
############## README #############
## This script shows how to run an R script in a parallel fashion.
## It is assumed that the cluster uses SLURM queueing system and has modules
## that load specific software.
## Lines that start with #SBATCH specify input parameters for the queueing system,
## so be careful to check the requirements and limitations on your cluster and
## adjust the values accordingly.
## The parts marked with CHANGE HERE needs to be updated before each run.
#------------------------

#SBATCH --job-name=my_haplin_run
#SBATCH --output=haplin_cluster_run.out
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=8
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=100
#SBATCH --mail-user=user@domain.com
#SBATCH --mail-type=ALL

# To enable use of R and openMPI, load the modules (this part is cluster-specific,
# so check the guidelines for your cluster!)
module load R
module load openmpi

# Printout names of nodes assigned to the job
echo "nodes: $SLURM_JOB_NODELIST"

# !!! CHANGE HERE !!!
# 'NODES' should be the same as in #SBATCH --nodes above
NODES=2
# 'PPN' should be the same as in #SBATCH --ntasks-per-node above
PPN=8
# !!!

# name of the file that will contain the list of nodes
myhostfile="cur_nodes.dat"

#-------------------------
# Write the list of available nodes to the hostfile
if [ ! -e $myhostfile ]
then
    touch $myhostfile
else
    rm $myhostfile
fi

if [ `echo $SLURM_JOB_NODELIST | cut -d"[" -f 1- --output-delimiter=" " | wc -w` -gt 1 ]
then
    # first part of the nodes names:
    nodes_first=`echo $SLURM_JOB_NODELIST | cut -d"[" -f 1`
    # list of nodes:
    nodes_list=`echo $SLURM_JOB_NODELIST | cut -d"[" -f 2 | cut -d"]" -f 1`
    
	# nodelist is in format: [0,4]                                                                                                      
	if [ `echo $nodes_list | cut -d"," -f 1- --output-delimiter=" " | wc -w` -gt 1 ]
	then
		nodes_list=$(echo $nodes_list | cut -d"," -f 1- | sed 's/,/\n/g')
	fi
	# nodelist is in format: [0-3]                                                                                                          
	if [ `echo $nodes_list | cut -d"-" -f 1- --output-delimiter=" " | wc -w` -gt 1 ]
	then
		new_nodes_list=()
		for line in ${nodes_list[*]}
		do
		if [ `echo $line | cut -d"-" -f 1- --output-delimiter=" " | wc -w` -gt 1 ]
			then
			first_core=$(echo $line | cut -d"-" -f 1)
			last_core=$(echo $line | cut -d"-" -f 2)
			line=$(seq $first_core $last_core)
			fi
		new_nodes_list=($new_nodes_list $line)
		done
	fi
	echo "range of nodes: " ${new_nodes_list[*]}

    for nn in ${new_nodes_list[*]}
    do
	for core in `seq $PPN`
	do
	    echo ${nodes_first}$nn >> $myhostfile
	done
    done
else
    echo "one node"
    for nn in $SLURM_JOB_NODELIST
    do
	for core in `seq $PPN`
	do
	    echo $nn >> $myhostfile
	done
    done
fi
#-------------------------

echo "----STARTING THE JOB----"
date
echo "------------------------"

# !!! CHANGE HERE !!!
# 'haplin_cluster_run' should be changed to the name of your R script
# you can also change 'mpi_run.out' to anything you wish
mpiexec --hostfile $myhostfile -n 1 R --save < haplin_cluster_run.r >& mpi_run.out
# !!!

exit_status=$?
echo "----JOB EXITED WITH STATUS---: $exit_status"
exit $exit_status
echo "----DONE----"
