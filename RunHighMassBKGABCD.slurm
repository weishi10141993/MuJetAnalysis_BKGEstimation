#!/bin/bash
#SBATCH --job-name=ABCD18
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10000M
#SBATCH --output=ABCD18.out.%j
#SBATCH --account=122747014089

module load cctools
export PARROT_CVMFS_ALIEN_CACHE=/scratch/group/mitchcomp/CVMFS_cache
cmsSite=/scratch/group/mitchcomp/CMS/LOCAL_TAMU_HPRC
cmsMount=--mount=/cvmfs/cms.cern.ch/SITECONF/local=$cmsSite

parrot_run $cmsMount $SHELL << EOF
shopt -s expand_aliases
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
root -l -b -q HighMassBKGABCD18.C
EOF
