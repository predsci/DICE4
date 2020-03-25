#!/bin/bash

trap ctrl_c INT

function ctrl_c() {
  echo "==> Caught CTRL-C, shutting down!"
  kill $(jobs -p)
  wait
  exit 1
}

NPROC=0

MAXPROC=14
NUM_JOBS=14

allNames=("AG" "CW" "DM" "GP" "JM" "PR" "BL" "KN" "TT" "VI" "EC" "CR" "MX" "NI")
for f in ${allNames[@]};
do
   NPROC=$(($NPROC+1))
   if [ "$NPROC" -gt "$MAXPROC" ]; then
     wait
     NPROC=0
     echo "Submitted up to ${f}!"
   fi

Rscript zika-example.R year 2016 nMCMC 2e6 nreal 1 epi_model 2 model 4 isingle 1 Temp 10 nfit 80 mod_name ${f} > log${f} &

done
wait
echo "All Done!"
exit 0

