#!/bin/sh
#
# Your job name
#$ -N NBODY_OMP_MPI
#
# Use current working directory
#$ -cwd
#
# pe (Parallel environment) request. Set your number of processors here.
#$ -pe openmpi_6x1 48
#
# Run job through bash shell
#$ -S /bin/bash

# If modules are needed, source modules environment:
. /etc/profile.d/modules.sh

# Add any modules you might require:
module add shared openmpi/gcc/64/1.8.8
module add shared gcc/7.2.0

echo $PE_HOSTFILE

export OMP_NUM_THREADS=6

# Run your application
mpirun -npernode 1 nbody_mpi -f NBodyInput/input_16384.txt -o false > nbody_mpi_output_file.txt
