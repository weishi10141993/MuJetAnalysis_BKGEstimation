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

## Low mass
A data driven method is used to estimate the SM yield, which mostly come from QCD processes.
```
root -l -b -q LowMassBKGFit1D.C #For plot and fit 1D template and save to work space
root -l -b -q LowMassBKGPlot2D.C  #For 2D template and draw 2D data points

```

## High mass
A MC vs data comparison at the control region is used. Meanwhile, a sanity check is done with ABCD method, which is purely data-driven.
```
root -l -b -q HighMassBKGShape.C #this compares MC to DATA at control region and fit MC in SR to get a bkg shape
root -l -b -q HighMassBKGABCD.C #this is datadriven ABCD method and checks bkg yield and the stability of the method
```
