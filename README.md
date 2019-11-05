# ANALYSIS_SOFT
Analysis package for AESOPLITE: MC simulations with FLUKA, pattern recognition, track reconstructions with Kalman filter 

Pre-requisite packages for installation:

ROOT, libx11-dev, libxpm-dev 

Instructions to install and use the package:

**Download the package, in the power directory, type 
 source ./setup
 
**Install the libraries and compile source code by typing:
cd src
make

**Create the main executable file by typing:
xmkmf
make

**Executable file will be in the directory /prod
cd prod
./Main

All paths have been changed to relative ones, there should be no modifications necessary in the makefile. 

