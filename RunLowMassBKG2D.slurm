#!/bin/bash
#SBATCH --job-name=LowMassBKG2D18
#SBATCH --time=02:30:00
#SBATCH --ntasks=1
#SBATCH --mem=700M
#SBATCH --output=LowMassBKG2D.out.%j
#SBATCH --account=122747014089

module load cctools
export PARROT_CVMFS_ALIEN_CACHE=/scratch/group/mitchcomp/CVMFS_cache
cmsSite=/scratch/group/mitchcomp/CMS/LOCAL_TAMU_HPRC
cmsMount=--mount=/cvmfs/cms.cern.ch/SITECONF/local=$cmsSite

parrot_run $cmsMount $SHELL << EOF
shopt -s expand_aliases
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
root -l -b -q LowMassBKGPlot2D18.C
EOF
