#!/bin/sh
#
# Your job name
#$ -N NBODY_MPI
#
# Use current working directory
#$ -cwd
#
# pe (Parallel environment) request. Set your number of processors here.
#$ -pe openmpi 2
#
# Run job through bash shell
#$ -S /bin/bash

# If modules are needed, source modules environment:
. /etc/profile.d/modules.sh

# Add any modules you might require:
module add shared openmpi/gcc/64/1.8.8
module add shared gcc/7.2.0

# The following output will show in the output file
echo "Got $NSLOTS processors."

# Run your application
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false > nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt
mpirun nbody_mpi -f NBodyInput/input_64.txt -o false >> nbody_mpi_output_file_64.txt