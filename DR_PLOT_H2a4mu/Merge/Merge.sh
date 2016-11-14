rm -rf file*.txt
find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/161110_164632/ | grep root | grep ana | grep -v log | grep -v ailed > fileB.txt
find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/161110_164709/ | grep root | grep ana | grep -v log | grep -v ailed > fileD.txt
find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/161110_164728/ | grep root | grep ana | grep -v log | grep -v ailed > fileE.txt
find /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/161110_164805/ | grep root | grep ana | grep -v log | grep -v ailed > fileG.txt

rm -rf /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016B-PromptReco-v2/crab_Run2016B-PromptReco-v2_*root
rm -rf /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016D-PromptReco-v2/crab_Run2016D-PromptReco-v2_*root
rm -rf /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016E-PromptReco-v2/crab_Run2016E-PromptReco-v2_*root
rm -rf /fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016G-PromptReco-v1/crab_Run2016G-PromptReco-v1_*root

split -l 500 fileB.txt file_B_
split -l 500 fileD.txt file_D_
split -l 500 fileE.txt file_E_
split -l 500 fileG.txt file_G_

bash Systematic_Merge.sh
