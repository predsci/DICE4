#!/bin/bash

trap ctrl_c INT

function ctrl_c() {
  echo "==> Caught CTRL-C, shutting down!"
  kill $(jobs -p)
  wait
  exit 1
}

NPROC=0

MAXPROC=4
NUM_JOBS=4

for (( f=18; f<=43; f++ ))
do
   NPROC=$(($NPROC+1))
   if [ "$NPROC" -gt "$MAXPROC" ]; then
     wait
     NPROC=0
     echo "Submitted up to ${f}!"
   fi

Rscript sd-example.R year 2017 nMCMC 2e6 epi_model 2 model 3 isingle 1 Temp 2 nfit ${f} > log${f} &

done
wait
echo "All Done!"
exit 0

