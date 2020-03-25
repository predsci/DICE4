#!/bin/bash

trap ctrl_c INT

function ctrl_c() {
  echo "==> Caught CTRL-C, shutting down!"
  kill $(jobs -p)
  wait
  exit 1
}

NPROC=0

MAXPROC=5
NUM_JOBS=5

allNames=("AW" "BO" "AR" "EC" "PE")
for f in ${allNames[@]};
do
   NPROC=$(($NPROC+1))
   if [ "$NPROC" -gt "$MAXPROC" ]; then
     wait
     NPROC=0
     echo "Submitted up to ${f}!"
   fi

Rscript zika-example.R year 2017 nMCMC 2e6 nreal 1 epi_model 2 model 4 isingle 1 Temp 10 nfit 80 mod_name ${f} > log${f} &

done
wait
echo "All Done!"
exit 0

