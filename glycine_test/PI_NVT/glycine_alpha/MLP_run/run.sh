#!/bin/bash

# remove the driver file if it already exists
rm /tmp/ipi_dftb_nvt /tmp/ipi_gap_d_nvt

# initialize the socket and set up the simulation
i-pi input.xml > log.ipi & 
sleep 30
i-pi-py_driver -u -a gap_d_nvt -m rascal -o ../../CSD_GAP_model.json,init.extxyz > log.gap &
/home/lumiaro/dftb+/bin/dftb+ > log.dftb &
wait

