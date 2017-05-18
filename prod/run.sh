#!/bin/bash


### Firt parameter is the energy in MeV
declare -i m

for m in `seq 1 100`
 do
  ./MainRecoEventMC 3 $1 $m 1 Test2  
done
