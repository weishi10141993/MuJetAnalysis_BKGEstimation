# MuJetAnalysis Background Estimation
This instruction is divided into two sections. First section is dedicated to low mass regime (0.25-9GeV). The second section is for high mass regime (11-59GeV).
General setup is:
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
git clone -b test git@github.com:weishi10141993/MuJetAnalysis_bbBarEstimation
```
## User Config
Edit Config.h to change the year (default is 2018) you want to do the estimation. This configuration be picked up in later macros. Many constants used in later macros are also defined in Constants.h in the main directory.

## Low mass
A data driven method is used to estimate the SM yield.
```
root -l -b -q LowMassBKGFit1D18.C #For plot and fit 1D template and save to work space
root -l -b -q LowMassBKGPlot2D18.C  #For 2D template and draw 2D data points
```

## High mass
A MC vs data comparison at the control region. And it also gives the MC distributions at signal region.
```
root -l -b -q HighMassBKGShape18.C
```

A sanity check is done with ABCD method, which is purely data-driven. It gives estimated yield in signal region A. The stability of this method is also checked.
```
root -l -b -q HighMassBKGABCD18.C
```
