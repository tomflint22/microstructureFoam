#!/bin/bash --login
# Runs in current dir by default
#SBATCH -p multicore  # (or --partition=) 
#SBATCH -n 8 # -c (or --cpus-per-task=) Number of cores to use for OpenMP (2--40)

# Can load modulefiles
module load openfoam/10-foss-2021a

source $FOAM_BASH

cp -r initial 0
blockMesh &> blockMesh.log
setFields &> setFields.log
TesselateFoam &> TesselateFoam.log
decomposePar &> decomposePar.log
mpirun -np 8 microstructureFoam -parallel &> microstructureFoam.log
reconstructPar &> reconstructPar.log
