# ALanalysis

########################
SIMULATION:
########################

aesoplite.inp: Geometry file for FLUKA
source.f: Define the beam of source particle
nmgdraw.f: Define the boundary crossing. 

To Compile:
$FLUPRO/flutil/lfluka -m fluka -o ALfluka source.f nmgdraw.f
It creates the executable ALfluka

To run the first cycle:
$FLUPRO/flutil/rfluka -e ALfluka -N0 -M1 aesoplite.inp

One output file is created per cycle, e.g., aesoplite001_fort.99

########################
CONVERSION TO ROOT FORMAT
########################

write99toroot.C: Convert the output file ".99" to the ROOT format
The text file is read line by line.
Multiple text files can be converted and merged into a single ROOT file. 
More details are available in the comments of the file.

To run:
root -l -q -x write99toroot.C


