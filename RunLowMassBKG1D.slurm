#!/bin/bash
#SBATCH --job-name=LowMassBKG1D18
#SBATCH --time=00:40:00
#SBATCH --ntasks=1
#SBATCH --mem=2000M
#SBATCH --output=LowMassBKG1D.out.%j
#SBATCH --account=122747014489

module load cctools
export PARROT_CVMFS_ALIEN_CACHE=/scratch/group/mitchcomp/CVMFS_cache
cmsSite=/scratch/group/mitchcomp/CMS/LOCAL_TAMU_HPRC
cmsMount=--mount=/cvmfs/cms.cern.ch/SITECONF/local=$cmsSite

parrot_run $cmsMount $SHELL << EOF
shopt -s expand_aliases
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
root -l -b -q LowMassBKGFit1D18.C
EOF
