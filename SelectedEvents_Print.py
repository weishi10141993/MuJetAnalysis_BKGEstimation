import ROOT, array, os, re, math, random
from math import *

tfile = ROOT.TFile("FitNtuple_Run2012_DoubleMu_v1.root")
#ttree = tfile.Get("FitNtupleTrackerMuons/dimudimu")
ttree_dimudimu   = tfile.Get("FitNtuplePFMuons/dimudimu")
ttree_dimuorphan = tfile.Get("FitNtuplePFMuons/dimuorphan")

n_tot = 0
n_sel_m1 = 0
n_sel_m2 = 0
for t in ttree_dimuorphan:
  n_tot = n_tot + 1
  # containstrig2 > 0 && mass > 0.25 && mass < 3.55
  if t.containstrig2 > 0 and t.mass > 0.25 and t.mass < 3.55:
    n_sel_m1 = n_sel_m1 + 1
  # containstrig  > 0 && mass > 0.25 && mass < 3.55
  if t.containstrig > 0 and t.mass > 0.25 and t.mass < 3.55:
    n_sel_m2 = n_sel_m2 + 1

print "dimuorphan n_tot:", n_tot
print "           n_sel_m1:", n_sel_m1
print "           n_sel_m2:", n_sel_m2

n_tot = 0
n_sel_signal = 0
n_sel_offDiagonal = 0
n_sel_control_offDiagonal = 0
for t in ttree_dimudimu:
  n_tot = n_tot + 1
  # isoC_1mm<2. && isoF_1mm<2. && abs(massC-massF) <= (0.13 + 0.065*(massC+massF)/2.)
  if t.isoC_1mm < 2. and t.isoF_1mm < 2. and abs(t.massC - t.massF) <= (0.13 + 0.065*(t.massC + t.massF)/2.0):
    print "Signal:", t.massC, t.massF, t.isoC_1mm, t.isoF_1mm, t.run, t.lumi, t.event
    n_sel_signal = n_sel_signal + 1
  if t.isoC_1mm < 2. and t.isoF_1mm < 2. and abs(t.massC - t.massF) > (0.13 + 0.065*(t.massC + t.massF)/2.0) and t.massC>0.25 and t.massC<3.55 and t.massF>0.25 and t.massF<3.55:
    print "bbBar:", t.massC, t.massF, t.isoC_1mm, t.isoF_1mm, t.run, t.lumi, t.event
    n_sel_offDiagonal = n_sel_offDiagonal + 1
  # abs(massC-massF) > (0.13 + 0.065*(massC+massF)/2.) && massC>0.25 && massC<3.55 && massF>0.25 && massF<3.55
  if abs(t.massC-t.massF) > (0.13 + 0.065*(t.massC+t.massF)/2.) and t.massC>0.25 and t.massC<3.55 and t.massF>0.25 and t.massF<3.55:
    n_sel_control_offDiagonal = n_sel_control_offDiagonal + 1

print "dimudimu n_tot:", n_tot
print "         n_sel_signal:             ", n_sel_signal
print "         n_sel_offDiagonal:        ", n_sel_offDiagonal
print "         n_sel_control_offDiagonal:", n_sel_control_offDiagonal


