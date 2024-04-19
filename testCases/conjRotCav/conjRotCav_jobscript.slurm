#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --exclusive
#SBATCH --time=12:00:00

#SBATCH --partition=standard
#SBATCH --qos=standard
####lowpriority
####standard

#SBATCH --account=ec203

# Load the openfoam module and source the bashrc file

export FOAM_INSTALL_PATH='/work/ec203/ec203/jp905/OpenFOAM-v2006-Mix'

source /work/ec203/ec203/jp905/site/modules.sh
source ${FOAM_INSTALL_PATH}/etc/bashrc

# Compose OpenFOAM work in the ususal way, except that parallel
# executables are launched via srun. For example:

srun chtMultiRegionSRFFoam -parallel
