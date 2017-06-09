#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1D.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TFormula.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMarker.h>
#include <TChain.h>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include "TTree.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TEventList.h"
#include <iostream>
#include <sstream>  
#include <fstream>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

void DR(){ 
  TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
  //Input file
  TChain chain_data_dimuorphan("cutFlowAnalyzerPXBL2PXFL2_Data/Events_orphan");
  std::ifstream Myfile( "../Input_2016BCDEFG_60GeV.txt" );
  std::string Line;
  if( !Myfile ) std::cout<<"ERROR opening Myfile."<<std::endl;
  while (std::getline(Myfile, Line)){
    TString Line2(Line);
    if( Line2.Contains("root") ){
	TString FileName = "../" + Line2;
	chain_data_dimuorphan.Add(FileName.Data());
    }
  }

  TTreeReader myReader(&chain_data_dimuorphan);
  TTreeReaderValue<float> orph_PtMu1(myReader,      "orph_PtMu1");
  TTreeReaderValue<float> orph_EtaMu1(myReader,     "orph_EtaMu1");
  TTreeReaderValue<float> orph_PhiMu1(myReader,     "orph_PhiMu1");
  TTreeReaderValue<float> orph_PtMu0(myReader,      "orph_PtMu0");
  TTreeReaderValue<float> orph_EtaMu0(myReader,     "orph_EtaMu0");
  TTreeReaderValue<float> orph_PhiMu0(myReader,     "orph_PhiMu0");
  TTreeReaderValue<float> orph_PtOrph(myReader,     "orph_PtOrph");
  TTreeReaderValue<float> orph_EtaOrph(myReader,    "orph_EtaOrph");
  TTreeReaderValue<float> orph_PhiOrph(myReader,    "orph_PhiOrph");
  TTreeReaderValue<float> orph_dimu_mass(myReader,  "orph_dimu_mass");
  TTreeReaderValue<float> orph_dimu_isoTk(myReader, "orph_dimu_isoTk");
  TTreeReaderValue<float> orph_isoTk(myReader,      "orph_isoTk");
  TTreeReaderValue<int>   containstrig(myReader,    "containstrig");
  TTreeReaderValue<int>   containstrig2(myReader,   "containstrig2");
  TTreeReaderValue<int>   NPATJet(myReader,         "NPATJet");
  TTreeReaderArray<float> PAT_jet_pt(myReader,      "PAT_jet_pt");
  TTreeReaderArray<float> PAT_jet_eta(myReader,     "PAT_jet_eta");
  TTreeReaderArray<float> PAT_jet_phi(myReader,     "PAT_jet_phi");
  TTreeReaderArray<float> PAT_jet_en(myReader,      "PAT_jet_en");
  TTreeReaderArray<float> PAT_jet_Btag1(myReader,   "PAT_jet_Btag1");
  TTreeReaderArray<float> PAT_jet_Btag2(myReader,   "PAT_jet_Btag2");
  TTreeReaderArray<float> PAT_jet_Btag3(myReader,   "PAT_jet_Btag3");
  //Histos
  TH1F *h_DR_min                    = new TH1F("h_DR_min","", 100, 0., 6.); h_DR_min->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_noCut              = new TH1F("h_DR_min_noCut","", 100, 0., 6.); h_DR_min_noCut->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_cut                = new TH1F("h_DR_min_cut","", 100, 0., 6.); h_DR_min_cut->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_2Jets              = new TH1F("h_DR_min_2Jets","", 100, 0., 6.); h_DR_min_2Jets->GetXaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_mass               = new TH2F("h_DR_min_mass","", 20, 0., 10., 50, 0., 6.); h_DR_min_mass->GetXaxis()->SetTitle("Mass [GeV]"); h_DR_min_mass->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_DR                 = new TH2F("h_DR_min_DR","", 50, 0., 1.4, 50, 0., 6.); h_DR_min_DR->GetXaxis()->SetTitle("#Delta R(#mu1 #mu 2)"); h_DR_min_DR->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_Iso                = new TH2F("h_DR_min_Iso","", 50, 0., 2., 50, 0., 6.); h_DR_min_Iso->GetXaxis()->SetTitle("Iso."); h_DR_min_Iso->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_Iso_nocut          = new TH2F("h_DR_min_Iso_nocut","", 50, 0., 50., 50, 0., 6.); h_DR_min_Iso_nocut->GetXaxis()->SetTitle("Iso."); h_DR_min_Iso_nocut->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_IsoOr_nocut        = new TH2F("h_DR_min_IsoOr_nocut","", 50, 0., 50., 50, 0., 6.); h_DR_min_IsoOr_nocut->GetXaxis()->SetTitle("Iso."); h_DR_min_IsoOr_nocut->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DRmumu_mass               = new TH2F("h_DRmumu_mass","", 10, 0., 10., 50, 0., 6.); h_DRmumu_mass->GetXaxis()->SetTitle("Mass [GeV]"); h_DRmumu_mass->GetYaxis()->SetTitle("#Delta R(#mu 1,#mu 2)");
  TH1F *h_JetPt_lead                = new TH1F("h_JetPt_lead","", 100, 0., 300.); h_JetPt_lead->GetXaxis()->SetTitle("Pt [GeV]");
  TH1F *h_JetPt_Sublead             = new TH1F("h_JetPt_Sublead","", 100, 0., 300.); h_JetPt_Sublead->GetXaxis()->SetTitle("Pt [GeV]");
  TH1F *h_Njet_0                    = new TH1F("h_Njet_0", "", 30, -0.5, 30.5); h_Njet_0->GetXaxis()->SetTitle("#Jet");
  TH1F *h_Njet_20                   = new TH1F("h_Njet_20","", 10, -0.5, 10.5); h_Njet_20->GetXaxis()->SetTitle("#Jet (Pt>20 GeV)");
  TH1F *h_Njet_30                   = new TH1F("h_Njet_30","", 10, -0.5, 10.5); h_Njet_30->GetXaxis()->SetTitle("#Jet (Pt>30 GeV)");
  TH2F *h_DR_min_Njet_20            = new TH2F("h_DR_min_Njet_20","",10, -0.5, 10.5, 50, 0., 6.); h_DR_min_Njet_20->GetXaxis()->SetTitle("#Jet (Pt>20 GeV)"); h_DR_min_Njet_20->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_Njet_30            = new TH2F("h_DR_min_Njet_30","",10, -0.5, 10.5, 50, 0., 6.); h_DR_min_Njet_30->GetXaxis()->SetTitle("#Jet (Pt>20 GeV)"); h_DR_min_Njet_30->GetYaxis()->SetTitle("#Delta R min");
  TH1F *h_Jet_Btag2_lowptAll        = new TH1F("h_Jet_Btag2_lowptAll","",100, 0., 1.); h_Jet_Btag2_lowptAll->GetXaxis()->SetTitle("Btag 2 (low Pt)");
  TH1F *h_Jet_Btag2_pt20All         = new TH1F("h_Jet_Btag2_pt20All","",100, 0., 1.); h_Jet_Btag2_pt20All->GetXaxis()->SetTitle("Btag 2 (Pt>20 GeV)");
  TH1F *h_Jet_Btag2_pt30All         = new TH1F("h_Jet_Btag2_pt30All","",100, 0., 1.); h_Jet_Btag2_pt30All->GetXaxis()->SetTitle("Btag 2 (Pt>30 GeV)");
  TH1F *h_DRLeadJetMu               = new TH1F("h_DRLeadJetMu","", 50, 0., 0.6); h_DRLeadJetMu->GetXaxis()->SetTitle("#Delta R(LeadJet-Mu)");
  TH1F *h_DRSubLeadJetMu            = new TH1F("h_DRSubLeadJetMu","", 50, 0., 0.6); h_DRSubLeadJetMu->GetXaxis()->SetTitle("#Delta R(SubLeadJet-Mu)");
  TH1F *h_LeadJet_Btag2             = new TH1F("h_LeadJet_Btag2","", 50, 0., 1.); h_LeadJet_Btag2->GetXaxis()->SetTitle("B-Tag2 Lead Jet");
  TH1F *h_SubLeadJet_Btag2          = new TH1F("h_SubLeadJet_Btag2","", 50, 0., 1.); h_SubLeadJet_Btag2->GetXaxis()->SetTitle("B-Tag2 SubLead Jet");
  TH1F *h_LeadJet_Btag2_matchCut    = new TH1F("h_LeadJet_Btag2_matchCut","", 50, 0., 1.); h_LeadJet_Btag2_matchCut->GetXaxis()->SetTitle("B-Tag2 Lead Jet");
  TH1F *h_SubLeadJet_Btag2_matchCut = new TH1F("h_SubLeadJet_Btag2_matchCut","", 50, 0., 1.); h_SubLeadJet_Btag2_matchCut->GetXaxis()->SetTitle("B-Tag2 SubLead Jet");
  TH2F *h_DR_min_BtagLead           = new TH2F("h_DR_min_BtagLead","", 50, 0., 1., 50, 0., 6.); h_DR_min_BtagLead->GetXaxis()->SetTitle("B-Tag2 Lead Jet"); h_DR_min_BtagLead->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_BtagSubLead        = new TH2F("h_DR_min_BtagSubLead","", 50, 0., 1., 50, 0., 6.); h_DR_min_BtagSubLead->GetXaxis()->SetTitle("B-Tag2 SubLead Jet"); h_DR_min_BtagSubLead->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_BtagLead_matchCut  = new TH2F("h_DR_min_BtagLead_matchCut","", 50, 0., 1., 50, 0., 6.); h_DR_min_BtagLead_matchCut->GetXaxis()->SetTitle("B-Tag2 Lead Jet"); h_DR_min_BtagLead_matchCut->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_BtagSubLead_matchCut= new TH2F("h_DR_min_BtagSubLead_matchCut","", 50, 0., 1., 50, 0., 6.); h_DR_min_BtagSubLead_matchCut->GetXaxis()->SetTitle("B-Tag2 SubLead Jet"); h_DR_min_BtagSubLead_matchCut->GetYaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_HadLeadSubLead     = new TH1F("h_DR_min_HadLeadSubLead","", 100, 0., 6.); h_DR_min_HadLeadSubLead->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_matchCut           = new TH1F("h_DR_min_matchCut","", 100, 0., 6.); h_DR_min_matchCut->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_BtagCutOrph        = new TH1F("h_DR_min_BtagCutOrph","", 100, 0., 6.); h_DR_min_BtagCutOrph->GetXaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_BtagOrph           = new TH2F("h_DR_min_BtagOrph","", 50, 0., 1., 50, 0., 6.); h_DR_min_BtagOrph->GetXaxis()->SetTitle("B-tag Orph."); h_DR_min_BtagOrph->GetYaxis()->SetTitle("#Delta R min");
  TH1F *h_bbMass17_2MatchJet        = new TH1F("h_bbMass17_2MatchJet","", 210, 0., 60.); h_bbMass17_2MatchJet->GetXaxis()->SetTitle("Mass [GeV]");
  TH1F *h_bbMassMix_2MatchJet       = new TH1F("h_bbMassMix_2MatchJet","", 210, 0., 60.); h_bbMassMix_2MatchJet->GetXaxis()->SetTitle("Mass [GeV]");
  TH1F *h_bbMass17_2MatchJet_10     = new TH1F("h_bbMass17_2MatchJet_10","", 210, 0., 10.); h_bbMass17_2MatchJet_10->GetXaxis()->SetTitle("Mass [GeV]");
  TH1F *h_bbMassMix_2MatchJet_10    = new TH1F("h_bbMassMix_2MatchJet_10","", 210, 0., 10.); h_bbMassMix_2MatchJet_10->GetXaxis()->SetTitle("Mass [GeV]");
  TH1F *h_bbMass17_0MatchJet_10     = new TH1F("h_bbMass17_0MatchJet_10","", 210, 0., 10.); h_bbMass17_0MatchJet_10->GetXaxis()->SetTitle("Mass [GeV]");
  TH1F *h_bbMassMix_0MatchJet_10    = new TH1F("h_bbMassMix_0MatchJet_10","", 210, 0., 10.); h_bbMassMix_0MatchJet_10->GetXaxis()->SetTitle("Mass [GeV]");
  TH1F *h_bbMassMix_0MatchJet_DRcut = new TH1F("h_bbMassMix_0MatchJet_DRcut","", 210, 0., 60.); h_bbMassMix_0MatchJet_DRcut->GetXaxis()->SetTitle("Mass [GeV]");

  // Loop over all entries
  while (myReader.Next()) {
    if(fabs(*orph_EtaMu0)<2.4 && fabs(*orph_EtaMu1)<2.4 && fabs(*orph_EtaOrph)<2.4 && (*containstrig>0 || *containstrig2>0) && *orph_dimu_mass>0.){
	bool Orphan_to_Leading = true;
	if(*orph_dimu_mass<9){
	  // Di-muons - orphan
	  TLorentzVector Mu0, Mu1, MuOr, diMu;
	  Mu0.SetPtEtaPhiM(*orph_PtMu0,*orph_EtaMu0,*orph_PhiMu0,0.);
	  Mu1.SetPtEtaPhiM(*orph_PtMu1,*orph_EtaMu1,*orph_PhiMu1,0.);
	  MuOr.SetPtEtaPhiM(*orph_PtOrph,*orph_EtaOrph,*orph_PhiOrph,0.);
	  diMu = Mu0 + Mu1;
	  float DR0 = Mu0.DeltaR(MuOr);
	  float DR1 = Mu1.DeltaR(MuOr);
	  float DR_min = DR0 < DR1 ? DR0 : DR1;
	  // First Plots
	  h_DR_min_noCut->Fill(DR_min); //Bump at 3 and at 0. Why at 0?
	  h_DR_min_Iso_nocut->Fill(*orph_dimu_isoTk, DR_min); //Bump at 0 partially for non ISO dimuon
	  h_DR_min_IsoOr_nocut->Fill(*orph_isoTk, DR_min,DR_min); //Bump at 0 not due to ISO orphan
	  // Apply isolation
	  if(*orph_dimu_isoTk<2 && *orph_dimu_isoTk>=0){
	    h_DR_min->Fill(DR_min);  // Still events at 0
	    if(DR_min>Mu0.DeltaR(Mu1)) h_DR_min_cut->Fill(DR_min); // Small bump disappear. Still a tail toward 0.
	    h_DR_min_mass->Fill(*orph_dimu_mass,DR_min); // Not due to mass
	    h_DRmumu_mass->Fill(*orph_dimu_mass,Mu0.DeltaR(Mu1));  // Larger mass means larger DR between mu0 and mu1
	    h_DR_min_DR->Fill(Mu0.DeltaR(Mu1), DR_min); // Not due to dimuon angle. Means when mu0 and m1 are close, also orphan is close!
	    h_DR_min_Iso->Fill(*orph_dimu_isoTk, DR_min); //Bump at 0 not due to ISO orphan
	    // CHECK the number of jets in the event
	    int nJet_0=0, nJet_20=0, nJet_30=0;
	    for( int i=0; i<*NPATJet; i++ ){
		if( fabs(PAT_jet_eta[i])<2.4 ){
		  nJet_0++;
		  if( PAT_jet_pt[i]>20. ){
		    h_Jet_Btag2_pt20All->Fill(PAT_jet_Btag2[i]);
		    nJet_20++;
		  }
		  if( PAT_jet_pt[i]>30. ){
		    h_Jet_Btag2_pt30All->Fill(PAT_jet_Btag2[i]);
		    nJet_30++;
		  }
		  if( PAT_jet_pt[i]<30. ) h_Jet_Btag2_lowptAll->Fill(PAT_jet_Btag2[i]);
		}
	    }
	    h_Njet_0->Fill(nJet_0); h_Njet_20->Fill(nJet_20); h_Njet_30->Fill(nJet_30); // Mostly 2 jets, but also more sometimes
	    h_DR_min_Njet_20->Fill(nJet_20, DR_min); // Maybe asking 2 jets remove a bit the tail at 0
	    h_DR_min_Njet_30->Fill(nJet_30, DR_min);
	    if(nJet_20==2) h_DR_min_2Jets->Fill(DR_min); // maybe
	    //Look for the Leading and SubLeading jets
	    float minPt = 20; int id_Lead = -1, id_subLead=-1;
	    for( int i=0; i<*NPATJet; i++ ){
		float thisPt = PAT_jet_pt[i];
		if( thisPt > minPt  ){
		  minPt = thisPt;
		  id_Lead = i;
		}
	    }
	    minPt = 20;
	    for( int i=0; i<*NPATJet; i++ ){
		float thisPt = PAT_jet_pt[i];
		if( thisPt > minPt && i!=id_Lead ){
		  minPt = thisPt;
		  id_subLead = i;
		}
	    }
	    if(id_Lead!=-1)    h_JetPt_lead->Fill( PAT_jet_pt[id_Lead] );
	    else               h_JetPt_lead->Fill( -50 );
	    if(id_subLead!=-1) h_JetPt_Sublead->Fill( PAT_jet_pt[id_subLead] );
	    else               h_JetPt_Sublead->Fill( -50 );
	    //Matching with b-jets
	    TLorentzVector LeadJet, SubLeadJet;
	    LeadJet.SetPtEtaPhiM(PAT_jet_pt[id_Lead],PAT_jet_eta[id_Lead],PAT_jet_phi[id_Lead],PAT_jet_en[id_Lead]);
	    SubLeadJet.SetPtEtaPhiM(PAT_jet_pt[id_subLead],PAT_jet_eta[id_subLead],PAT_jet_phi[id_subLead],PAT_jet_en[id_subLead]);
	    float minDR_Lead=-1, minDR_subLead = -1;
	    if(LeadJet.DeltaR(MuOr) < LeadJet.DeltaR(diMu)){
		if(SubLeadJet.DeltaR(MuOr) < LeadJet.DeltaR(MuOr)){
		  minDR_subLead = SubLeadJet.DeltaR(MuOr);
		  minDR_Lead = LeadJet.DeltaR(diMu);
		  Orphan_to_Leading = false;
		}
		else{
		  minDR_Lead = LeadJet.DeltaR(MuOr);
		  minDR_subLead = SubLeadJet.DeltaR(diMu);
		  Orphan_to_Leading = true;
		}
	    }
	    else{
		if(SubLeadJet.DeltaR(diMu) < LeadJet.DeltaR(diMu)){
		  minDR_subLead = SubLeadJet.DeltaR(diMu);
		  minDR_Lead = LeadJet.DeltaR(MuOr);
		  Orphan_to_Leading = true;
		}
		else{
		  minDR_Lead = LeadJet.DeltaR(diMu);
		  minDR_subLead = SubLeadJet.DeltaR(MuOr);
		  Orphan_to_Leading = false;
		}
	    }
	    h_DRLeadJetMu->Fill(minDR_Lead);
	    h_DRSubLeadJetMu->Fill(minDR_subLead);
	    //B-TAG
	    if( id_Lead!=-1 and id_subLead!=-1 ){
		h_DR_min_HadLeadSubLead->Fill(DR_min);
		h_LeadJet_Btag2->Fill(PAT_jet_Btag2[id_Lead]);
		h_SubLeadJet_Btag2->Fill(PAT_jet_Btag2[id_subLead]);
		h_DR_min_BtagLead->Fill(PAT_jet_Btag2[id_Lead],DR_min);
		h_DR_min_BtagSubLead->Fill(PAT_jet_Btag2[id_subLead],DR_min);
		if(minDR_Lead<0.1 && minDR_subLead<0.1){
		  h_LeadJet_Btag2_matchCut->Fill(PAT_jet_Btag2[id_Lead]);
		  h_SubLeadJet_Btag2_matchCut->Fill(PAT_jet_Btag2[id_subLead]);
		  h_DR_min_BtagLead_matchCut->Fill(PAT_jet_Btag2[id_Lead],DR_min);
		  h_DR_min_BtagSubLead_matchCut->Fill(PAT_jet_Btag2[id_subLead],DR_min);
		  h_DR_min_matchCut->Fill(DR_min); //Seems working!
		  if(Orphan_to_Leading )  h_DR_min_BtagOrph->Fill(PAT_jet_Btag2[id_Lead],DR_min);
		  if(!Orphan_to_Leading ) h_DR_min_BtagOrph->Fill(PAT_jet_Btag2[id_subLead],DR_min);
		  if(Orphan_to_Leading  && PAT_jet_Btag2[id_Lead]>0.5 ) h_DR_min_BtagCutOrph->Fill(DR_min);
		  if(!Orphan_to_Leading && PAT_jet_Btag2[id_subLead]>0.5 ) h_DR_min_BtagCutOrph->Fill(DR_min);
		}
	    }
	  }// Only For isolated dimuons
	}// Control region bb
	// Now study the mass
	if(*orph_dimu_mass<60 && *orph_dimu_isoTk<2 && *orph_dimu_isoTk>=0){
	  float minPt = 30; int id_Lead = -1; int id_subLead=-1;
	  for( int i=0; i<*NPATJet; i++ ){
	    float thisPt = PAT_jet_pt[i];
	    if( thisPt > minPt  ){
		minPt = thisPt;
		id_Lead = i;
	    }
	  }
	  minPt = 30;
	  for( int i=0; i<*NPATJet; i++ ){
	    float thisPt = PAT_jet_pt[i];
	    if( thisPt > minPt && i!=id_Lead ){
		minPt = thisPt;
		id_subLead = i;
	    }
	  }
	  //Matching with b-jets
	  TLorentzVector Mu0, Mu1, MuOr, diMu;
	  Mu0.SetPtEtaPhiM(*orph_PtMu0,*orph_EtaMu0,*orph_PhiMu0,0.);
	  Mu1.SetPtEtaPhiM(*orph_PtMu1,*orph_EtaMu1,*orph_PhiMu1,0.);
	  MuOr.SetPtEtaPhiM(*orph_PtOrph,*orph_EtaOrph,*orph_PhiOrph,0.);
	  diMu = Mu0 + Mu1;
	  float DR0 = Mu0.DeltaR(MuOr);
	  float DR1 = Mu1.DeltaR(MuOr);
	  float DR_min = DR0<DR1 ? DR0 : DR1;
	  TLorentzVector LeadJet, SubLeadJet;
	  LeadJet.SetPtEtaPhiM(PAT_jet_pt[id_Lead],PAT_jet_eta[id_Lead],PAT_jet_phi[id_Lead],PAT_jet_en[id_Lead]);
	  SubLeadJet.SetPtEtaPhiM(PAT_jet_pt[id_subLead],PAT_jet_eta[id_subLead],PAT_jet_phi[id_subLead],PAT_jet_en[id_subLead]);
	  float minDR_Lead=-1, minDR_subLead = -1;
	  if(LeadJet.DeltaR(MuOr) < LeadJet.DeltaR(diMu)){
	    if(SubLeadJet.DeltaR(MuOr) < LeadJet.DeltaR(MuOr)){
		minDR_subLead = SubLeadJet.DeltaR(MuOr);
		minDR_Lead = LeadJet.DeltaR(diMu);
		Orphan_to_Leading = false;
	    }
	    else{
		minDR_Lead = LeadJet.DeltaR(MuOr);
		minDR_subLead = SubLeadJet.DeltaR(diMu);
		Orphan_to_Leading = true;
	    }
	  }
	  else{
	    if(SubLeadJet.DeltaR(diMu) < LeadJet.DeltaR(diMu)){
		minDR_subLead = SubLeadJet.DeltaR(diMu);
		minDR_Lead = LeadJet.DeltaR(MuOr);
		Orphan_to_Leading = true;
	    }
	    else{
		minDR_Lead = LeadJet.DeltaR(diMu);
		minDR_subLead = SubLeadJet.DeltaR(MuOr); 
		Orphan_to_Leading = false;
	    }
	  }
	  if(id_Lead!=-1 && id_subLead!=-1 && minDR_Lead<0.1 && minDR_subLead<0.1){ // Seems to work!
	    if(*containstrig2 > 0) h_bbMassMix_2MatchJet->Fill(*orph_dimu_mass);
	    if(*containstrig > 0)  h_bbMass17_2MatchJet->Fill(*orph_dimu_mass);
	    if(*containstrig2 > 0) h_bbMassMix_2MatchJet_10->Fill(*orph_dimu_mass);
	    if(*containstrig > 0)  h_bbMass17_2MatchJet_10->Fill(*orph_dimu_mass);
	  }
	  if(*containstrig2 > 0){
	    h_bbMassMix_0MatchJet_10->Fill(*orph_dimu_mass);
	    if(DR_min>2 and DR_min<4){ h_bbMassMix_0MatchJet_DRcut->Fill(*orph_dimu_mass); }
	  }
	  if(*containstrig > 0)  h_bbMass17_0MatchJet_10->Fill(*orph_dimu_mass);
	}//Studying mass
    }
  }
  gStyle->SetOptStat(0);
  h_DR_min_noCut->Draw();                       myc1->SaveAs("dr_studies/h_NOCUT_DR_min.pdf");                delete h_DR_min_noCut;
  h_DR_min_Iso_nocut->Draw("colz");             myc1->SaveAs("dr_studies/h_NOCUT_DRmin_Iso.pdf");             delete h_DR_min_Iso_nocut;
  h_DR_min_IsoOr_nocut->Draw("colz");           myc1->SaveAs("dr_studies/h_NOCUT_DRmin_IsoOr.pdf");           delete h_DR_min_IsoOr_nocut;
  h_DR_min->Draw();                             myc1->SaveAs("dr_studies/h_ISOCUT_DR_min.pdf");
  h_DR_min_cut->SetLineColor(2); h_DR_min_cut->Draw("same");
  h_DR_min_2Jets->SetLineColor(kGreen); h_DR_min_2Jets->Draw("same"); h_DR_min_HadLeadSubLead->SetLineColor(kViolet); h_DR_min_HadLeadSubLead->Draw("same"); h_DR_min_matchCut->SetLineColor(kOrange); h_DR_min_matchCut->Draw("same");
  h_DR_min_BtagCutOrph->SetLineColor(kBlack); h_DR_min_BtagCutOrph->Draw("same"); myc1->SaveAs("dr_studies/h_ISOCUT_DR_min_IfGreatDimuonAngle.pdf");
  h_DR_min_BtagOrph->Draw("colz");              myc1->SaveAs("dr_studies/h_JET_DRmin_BtagOrph.pdf");              delete h_DR_min_BtagOrph;
  h_DR_min_mass->Draw("colz");                  myc1->SaveAs("dr_studies/h_ISOCUT_DR_min_mass.pdf"); TProfile *p_DR_min_mass = h_DR_min_mass->ProfileX(); p_DR_min_mass->SetMinimum(0);   p_DR_min_mass->Draw();  myc1->SaveAs("dr_studies/p_DR_min_mass.pdf");
  h_DR_min_DR->Draw("colz");                    myc1->SaveAs("dr_studies/h_ISOCUT_DR_min_DR.pdf");            delete h_DR_min_DR;
  h_DR_min_Iso->Draw("colz");                   myc1->SaveAs("dr_studies/h_ISOCUT_DR_min_Iso.pdf");           delete h_DR_min_Iso;
  h_DRmumu_mass->Draw("colz");                  myc1->SaveAs("dr_studies/h_ISOCUT_DRmumu_mass.pdf");          delete h_DRmumu_mass;
  gStyle->SetOptStat(111111);
  h_Njet_0->Draw();                             myc1->SaveAs("dr_studies/h_JET_Njet_0.pdf");                  delete h_Njet_0;
  h_Njet_20->Draw();                            myc1->SaveAs("dr_studies/h_JET_Njet_20.pdf");                 delete h_Njet_20;
  h_Njet_30->Draw();                            myc1->SaveAs("dr_studies/h_JET_Njet_30.pdf");                 delete h_Njet_30;
  h_DR_min_Njet_20->Draw("colz");               myc1->SaveAs("dr_studies/h_JET_DRmin_Njet20.pdf");            delete h_DR_min_Njet_20;
  h_DR_min_Njet_30->Draw("colz");               myc1->SaveAs("dr_studies/h_JET_DRmin_Njet30.pdf");            delete h_DR_min_Njet_30;
  h_JetPt_lead->Draw();                         myc1->SaveAs("dr_studies/h_JET_JetPt_lead.pdf");              delete h_JetPt_lead;
  h_JetPt_Sublead->Draw();                      myc1->SaveAs("dr_studies/h_JET_JetPt_Sublead.pdf");           delete h_JetPt_Sublead;
  gStyle->SetOptStat(0);
  h_Jet_Btag2_lowptAll->Draw();                 myc1->SaveAs("dr_studies/h_Jet_Btag2_lowptAll.pdf");          delete h_Jet_Btag2_lowptAll;
  h_Jet_Btag2_pt20All->Draw();                  myc1->SaveAs("dr_studies/h_Jet_Btag2_pt20All.pdf");           delete h_Jet_Btag2_pt20All;
  h_Jet_Btag2_pt30All->Draw();                  myc1->SaveAs("dr_studies/h_Jet_Btag2_pt30All.pdf");           delete h_Jet_Btag2_pt30All;
  h_DRLeadJetMu->Draw();                        myc1->SaveAs("dr_studies/h_JET_Match_LeadJetMu.pdf");             delete h_DRLeadJetMu;
  h_DRSubLeadJetMu->Draw();                     myc1->SaveAs("dr_studies/h_JET_match_SubLeadJetMu.pdf");          delete h_DRSubLeadJetMu;
  h_LeadJet_Btag2->Draw();                      myc1->SaveAs("dr_studies/h_Btag2_LeadJet.pdf");               delete h_LeadJet_Btag2;
  h_SubLeadJet_Btag2->Draw();                   myc1->SaveAs("dr_studies/h_Brag2_SubLeadJet.pdf");            delete h_SubLeadJet_Btag2;
  h_LeadJet_Btag2_matchCut->Draw();             myc1->SaveAs("dr_studies/h_JET_Btag2_LeadJet_matchCut.pdf");      delete h_LeadJet_Btag2_matchCut;
  h_SubLeadJet_Btag2_matchCut->Draw();          myc1->SaveAs("dr_studies/h_JET_Btag2_SubLeadJet_matchCut.pdf");   delete h_SubLeadJet_Btag2_matchCut;
  h_DR_min_BtagLead->Draw("colz");              myc1->SaveAs("dr_studies/h_DR_min_BtagLead.pdf");             delete h_DR_min_BtagLead;
  h_DR_min_BtagSubLead->Draw("colz");           myc1->SaveAs("dr_studies/h_DR_min_BtagSubLead.pdf");          delete h_DR_min_BtagSubLead;
  h_DR_min_BtagLead_matchCut->Draw("colz");     myc1->SaveAs("dr_studies/h_JET_DR_min_BtagLead_matchCut.pdf");    delete h_DR_min_BtagLead_matchCut;
  h_DR_min_BtagSubLead_matchCut->Draw("colz");  myc1->SaveAs("dr_studies/h_JET_DR_min_BtagSubLead_matchCut.pdf"); delete h_DR_min_BtagSubLead_matchCut;
  h_bbMass17_2MatchJet->Draw();                 myc1->SaveAs("dr_studies/h_JET60_bbMass17_2MatchJet.pdf");          delete h_bbMass17_2MatchJet;
  h_bbMassMix_2MatchJet->Draw();                myc1->SaveAs("dr_studies/h_JET60_bbMassMix_2MatchJet.pdf");         delete h_bbMassMix_2MatchJet;
  h_bbMass17_2MatchJet_10->Draw();              myc1->SaveAs("dr_studies/h_JET60_bbMass17_2MatchJet_10.pdf");       
  h_bbMassMix_2MatchJet_10->Draw();             myc1->SaveAs("dr_studies/h_JET60_bbMassMix_2MatchJet_10.pdf");
  h_bbMass17_0MatchJet_10->Draw();              myc1->SaveAs("dr_studies/h_JET60_bbMass17_0MatchJet_10.pdf");
  h_bbMassMix_0MatchJet_10->Draw();             myc1->SaveAs("dr_studies/h_JET60_bbMassMix_0MatchJet_10.pdf");
  h_bbMass17_2MatchJet_10->Scale(1./h_bbMass17_2MatchJet_10->GetEntries()); h_bbMass17_2MatchJet_10->Draw();
  h_bbMass17_0MatchJet_10->Scale(3./h_bbMass17_0MatchJet_10->GetEntries()); h_bbMass17_0MatchJet_10->SetLineColor(2); h_bbMass17_0MatchJet_10->Draw("same"); myc1->SaveAs("dr_studies/h_JET60_bbMass17_WithWithout_10.pdf");
  h_bbMassMix_2MatchJet_10->Scale(1./h_bbMassMix_2MatchJet_10->GetEntries()); h_bbMassMix_2MatchJet_10->Draw();
  h_bbMassMix_0MatchJet_10->Scale(3./h_bbMassMix_0MatchJet_10->GetEntries()); h_bbMassMix_0MatchJet_10->SetLineColor(2); h_bbMassMix_0MatchJet_10->Draw("same"); myc1->SaveAs("dr_studies/h_JET60_bbMassMix_WithWithout_10.pdf");
  delete h_bbMass17_2MatchJet_10; delete h_bbMassMix_2MatchJet_10; delete h_bbMass17_0MatchJet_10; delete h_bbMassMix_0MatchJet_10;
  h_bbMassMix_0MatchJet_DRcut->Draw();          myc1->SaveAs("dr_studies/h_JET60_bbMassMix_0MatchJet_DRcut.pdf"); delete h_bbMassMix_0MatchJet_DRcut;
  delete myc1;
}
