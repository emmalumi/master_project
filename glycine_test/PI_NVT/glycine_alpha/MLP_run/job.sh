#!/bin/bash

#SBATCH --time=96:00:00
#SBATCH --partition=Slim
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=24
#SBATCH --nodelist=cosmosrv6


rm /tmp/ipi_dftb_nvt /tmp/ipi_gap_d_nvt

i-pi input.xml > log.ipi &
sleep 30

mpirun -np 12 i-pi-py_driver -u -a gap_d_nvt -m rascal -o ../../CSD_GAP_model.json,init.extxyz > log.gap &

mpirun -np 12 /home/lumiaro/dftb+/bin/dftb+ > log.dftb &
wait


