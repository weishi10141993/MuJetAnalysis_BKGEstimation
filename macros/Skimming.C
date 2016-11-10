#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <string> 

//.x Skimming.C+("Input_2015D_v3_ONE.txt","BB_reduced.root")
void Skimming( TString fileList, TString outName ){

  //Doing the chain
  cout<<"Starting... "<<endl;
  TChain chain_data("cutFlowAnalyzer/Events");
  TChain chain_databb("cutFlowAnalyzer/Events_orphan");
  std::ifstream Myfile( fileList.Data() );
  std::string Line;
  if( !Myfile ) std::cout<<"ERROR opening "<<fileList<<std::endl;
  while (std::getline(Myfile, Line)){
    TString Line2(Line);
    if( Line2.Contains("root") ){
	chain_data.Add(Line2.Data());
	chain_databb.Add(Line2.Data());
    }
  }
  cout<<"Chain done!"<<endl;
  //Variables needed
  bool orph_passOffLineSel_bb, orph_passOffLineSelPt_bb, orph_passOffLineSelPt1788_bb, orph_FiredTrig_bb, orph_FiredTrig_pt_bb, orph_FiredTrig_ptColl_bb;
  int containstrig2_bb, containstrig_bb;
  float orph_dimu_mass_bb, orph_dimu_isoTk_bb;
  chain_databb.SetBranchAddress("orph_passOffLineSel",&orph_passOffLineSel_bb);
  chain_databb.SetBranchAddress("orph_passOffLineSelPt",&orph_passOffLineSelPt_bb);
  chain_databb.SetBranchAddress("orph_passOffLineSelPt1788",&orph_passOffLineSelPt1788_bb);
  chain_databb.SetBranchAddress("orph_FiredTrig",&orph_FiredTrig_bb);
  chain_databb.SetBranchAddress("orph_FiredTrig_pt",&orph_FiredTrig_pt_bb);
  chain_databb.SetBranchAddress("orph_FiredTrig_ptColl",&orph_FiredTrig_ptColl_bb);
  chain_databb.SetBranchAddress("orph_FiredTrig",&orph_FiredTrig_bb);
  chain_databb.SetBranchAddress("orph_dimu_isoTk",&orph_dimu_isoTk_bb);
  chain_databb.SetBranchAddress("containstrig2",&containstrig2_bb);
  chain_databb.SetBranchAddress("orph_dimu_mass",&orph_dimu_mass_bb);
  chain_databb.SetBranchAddress("containstrig",&containstrig_bb);
  float massC_mu, massF_mu, isoC_1mm_mu, isoF_1mm_mu;
  chain_data.SetBranchAddress("massC",&massC_mu);
  chain_data.SetBranchAddress("massF",&massF_mu);
  chain_data.SetBranchAddress("isoC_1mm",&isoC_1mm_mu);
  chain_data.SetBranchAddress("isoF_1mm",&isoF_1mm_mu);
  //My file
  TFile *f = new TFile(outName.Data(),"RECREATE");
  f->cd();
  TTree Events("Events","");
  float massC, massF, isoC_1mm, isoF_1mm;
  Events.Branch("massC",&massC,"massC/F");
  Events.Branch("massF",&massF,"massF/F");
  Events.Branch("isoC_1mm",&isoC_1mm,"isoC_1mm/F");
  Events.Branch("isoF_1mm",&isoF_1mm,"isoF_1mm/F");
  bool orph_passOffLineSel, orph_passOffLineSelPt, orph_passOffLineSelPt1788, orph_FiredTrig, orph_FiredTrig_pt, orph_FiredTrig_ptColl;
  int containstrig2, containstrig;
  float orph_dimu_mass, orph_dimu_isoTk;
  TTree Events_orphan("Events_orphan","");
  Events_orphan.Branch("orph_passOffLineSel",&orph_passOffLineSel,"orph_passOffLineSel/O");
  Events_orphan.Branch("orph_passOffLineSelPt",&orph_passOffLineSelPt,"orph_passOffLineSelPt/O");
  Events_orphan.Branch("orph_passOffLineSelPt1788",&orph_passOffLineSelPt1788,"orph_passOffLineSelPt1788/O");
  Events_orphan.Branch("orph_FiredTrig",&orph_FiredTrig,"orph_FiredTrig/O");
  Events_orphan.Branch("orph_FiredTrig_pt",&orph_FiredTrig_pt,"orph_FiredTrig_pt/O");
  Events_orphan.Branch("orph_FiredTrig_ptColl",&orph_FiredTrig_ptColl,"orph_FiredTrig_ptColl/O");
  Events_orphan.Branch("containstrig2",&containstrig2,"containstrig2/I");
  Events_orphan.Branch("containstrig",&containstrig,"containstrig/I");
  Events_orphan.Branch("orph_dimu_isoTk",&orph_dimu_isoTk,"orph_dimu_isoTk/F");
  Events_orphan.Branch("orph_dimu_mass",&orph_dimu_mass,"orph_dimu_mass/F");
  //Loops
  cout<<"Now starting loop on bb."<<endl;
  Long64_t entries = chain_databb.GetEntries(); 
  for(Long64_t i=0; i<entries; i++){ 
    if(i%10000==0) cout<<"i = "<<i<<" / "<<entries<<endl;
    chain_databb.GetEntry(i);
    if(orph_dimu_mass>-999.){
	orph_passOffLineSel = orph_passOffLineSel_bb;
	orph_passOffLineSelPt = orph_passOffLineSelPt_bb;
	orph_passOffLineSelPt1788 = orph_passOffLineSelPt1788_bb;
	orph_FiredTrig = orph_FiredTrig_bb;
	orph_FiredTrig_pt = orph_FiredTrig_pt_bb;
	orph_FiredTrig_ptColl = orph_FiredTrig_ptColl_bb;
	containstrig2 = containstrig2_bb;
	containstrig = containstrig_bb;
	orph_dimu_isoTk = orph_dimu_isoTk_bb;
	orph_dimu_mass = orph_dimu_mass_bb;
	Events_orphan.Fill();
    }
  }
  cout<<"Now starting loop on mu."<<endl;
  entries = chain_data.GetEntries(); 
  for(Long64_t i=0; i<entries; i++){ 
    if(i%10000==0) cout<<"i = "<<i<<" / "<<entries<<endl;
    chain_data.GetEntry(i);
    if(massC_mu>-999. || massF_mu>-999.){
	massC = massC_mu;
	massF = massF_mu;
	isoC_1mm = isoC_1mm_mu;
	isoF_1mm = isoF_1mm_mu;
	Events.Fill();
    }
  }
  cout<<"CLosing files"<<endl;
  f->Write();
  f->Close();
  cout<<"The end."<<endl;
}
