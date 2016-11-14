#!bin/bash

files_B="$(ls file_B_*)"
files_D="$(ls file_D_*)"
files_E="$(ls file_E_*)"
files_G="$(ls file_G_*)"

Num=0
inc=1
for sigle_file in $files_B
do
  echo hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/crab_Run2016B-PromptReco-v2_${Num}.root  @$sigle_file
  #hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/crab_Run2016B-PromptReco-v2_${Num}.root  @$sigle_file
  Num=$(( Num + inc ))
done
hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/crab_Run2016B-PromptReco-v2.root /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/crab_Run2016B-PromptReco-v2_*.root

#Num=0
#for sigle_file in $files_C
#do
#  echo hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016C-PromptReco-v2/crab_Run2016C-PromptReco-v2_${Num}.root  @$sigle_file
#  hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016C-PromptReco-v2/crab_Run2016C-PromptReco-v2_${Num}.root  @$sigle_file
#  Num=$(( Num + inc ))
#done
#hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016C-PromptReco-v2/crab_Run2016C-PromptReco-v2.root /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016C-PromptReco-v2/crab_Run2016C-PromptReco-v2_*.root

Num=0
for sigle_file in $files_D
do
  echo hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/crab_Run2016D-PromptReco-v2_${Num}.root  @$sigle_file
  hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/crab_Run2016D-PromptReco-v2_${Num}.root  @$sigle_file
  Num=$(( Num + inc ))
done
hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/crab_Run2016D-PromptReco-v2.root /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/crab_Run2016D-PromptReco-v2_*.root

Num=0
for sigle_file in $files_E
do
  echo hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/crab_Run2016E-PromptReco-v2_${Num}.root  @$sigle_file
  hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/crab_Run2016E-PromptReco-v2_${Num}.root  @$sigle_file
  Num=$(( Num + inc ))
done
hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/crab_Run2016E-PromptReco-v2.root /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/crab_Run2016E-PromptReco-v2_*.root

#Num=0
#for sigle_file in $files_F
#do
#  echo hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016F-PromptReco-v1/crab_Run2016F-PromptReco-v1_${Num}.root  @$sigle_file
#  #hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016F-PromptReco-v1/crab_Run2016F-PromptReco-v1_${Num}.root  @$sigle_file
#  Num=$(( Num + inc ))
#done
#hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016F-PromptReco-v1/crab_Run2016F-PromptReco-v1.root /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016F-PromptReco-v1/crab_Run2016F-PromptReco-v1_*.root

Num=0
for sigle_file in $files_G
do
  echo hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/crab_Run2016G-PromptReco-v1_${Num}.root  @$sigle_file
  hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/crab_Run2016G-PromptReco-v1_${Num}.root  @$sigle_file
  Num=$(( Num + inc ))
done
hadd -f /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/crab_Run2016G-PromptReco-v1.root /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/crab_Run2016G-PromptReco-v1_*.root

ls -ltr /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/crab_Run2016B-PromptReco-v2*.root
ls -ltr /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/crab_Run2016D-PromptReco-v2*.root
ls -ltr /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/crab_Run2016E-PromptReco-v2*.root
ls -ltr /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/crab_Run2016G-PromptReco-v1*.root
