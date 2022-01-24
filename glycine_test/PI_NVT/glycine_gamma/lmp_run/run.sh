#!/bin/bash

rm /tmp/ipi_lmp-base-gg 
rm /tmp/ipi_lmp-delta-gg


i-pi input.xml > log.ipi &
sleep 30
#mpirun -n 2 /home/lumiaro/Documents/CODE/lammps/src/lmp_mpi < lmp_pbe_ts.in > log.lmp_pbe_ts &
/home/lumiaro/Documents/CODE/lammps/src/lmp_mpi < lmp_pbe_ts.in > log.lmp_pbe_ts &
/home/lumiaro/Documents/CODE/lammps/src/lmp_mpi < lmp_pbe0_mbd.in > log.lmp_pbe0_mbd &

#mpirun -n 2 /home/lumiaro/Documents/CODE/lammps/src/lmp_mpi < lmp_pbe0_mbd.in > log.lmp_pbe0_mbd &
wait
