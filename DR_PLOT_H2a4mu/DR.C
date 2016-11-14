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
  //Files and variables
  TString Nfile = "/fdata/hepx/store/user/lpernie/DoubleMuon/crab_Run2016BDEG-PromptReco.root";
  TFile* file = TFile::Open(Nfile.Data());
  if( !file ) cout<<"Waring! File: "<<file<<" not present!"<<endl;
  TTreeReader myReader("cutFlowAnalyzerPXBL2PXFL2/Events_orphan", file);
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
  TTreeReaderArray<float> PAT_jet_Btag1(myReader,   "PAT_jet_Btag1");
  TTreeReaderArray<float> PAT_jet_Btag2(myReader,   "PAT_jet_Btag2");
  TTreeReaderArray<float> PAT_jet_Btag3(myReader,   "PAT_jet_Btag3");
  //Histos
  TH1F *h_JetPt_lead         = new TH1F("h_JetPt_lead","", 100, 0., 300.); h_JetPt_lead->GetXaxis()->SetTitle("Pt [GeV]");
  TH1F *h_CloseJetsN         = new TH1F("h_CloseJetsN","", 10, -0.5, 10.5); h_CloseJetsN->GetXaxis()->SetTitle("N Close Jets");
  TH1F *h_DrMin_ifJet        = new TH1F("h_DrMin_ifJet","", 100, 0., 6.); h_DrMin_ifJet->GetXaxis()->SetTitle("DR_min");
  TH1F *h_CloseJetDR         = new TH1F("h_CloseJetDR","", 50, 0., 0.5); h_CloseJetDR->GetXaxis()->SetTitle("DR closest jet");
  TH1F *h_CloseJetBtag2      = new TH1F("h_CloseJetBtag2","", 50, 0., 1); h_CloseJetBtag2->GetXaxis()->SetTitle("Btag 2");
  TH2F *h_DR_min_Btag2Closer = new TH2F("h_DR_min_Btag2Closer","", 50, 0., 1, 50, 0., 6.); h_DR_min_Btag2Closer->GetXaxis()->SetTitle("Btag 2"); h_DR_min_Btag2Closer->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_Btag3Closer = new TH2F("h_DR_min_Btag3Closer","", 50, 0., 1, 50, 0., 6.); h_DR_min_Btag3Closer->GetXaxis()->SetTitle("Btag 3"); h_DR_min_Btag3Closer->GetYaxis()->SetTitle("#Delta R min");
  TH1F *h_Njet_0             = new TH1F("h_Njet_0", "", 30, -0.5, 30.5); h_Njet_0->GetXaxis()->SetTitle("#Jet");
  TH1F *h_Njet_20            = new TH1F("h_Njet_20","", 10, -0.5, 10.5); h_Njet_20->GetXaxis()->SetTitle("#Jet (Pt>20 GeV)");
  TH1F *h_Njet_30            = new TH1F("h_Njet_30","", 10, -0.5, 10.5); h_Njet_30->GetXaxis()->SetTitle("#Jet (Pt>30 GeV)");
  TH1F *h_HigBtag2           = new TH1F("h_HigBtag2","", 100, 0., 6.); h_HigBtag2->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_LowBtag2           = new TH1F("h_LowBtag2","", 100, 0., 6.); h_LowBtag2->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_Jet_Btag2_lowptAll = new TH1F("h_Jet_Btag2_lowptAll","",100, 0., 1.); h_Jet_Btag2_lowptAll->GetXaxis()->SetTitle("Btag 2 (low Pt)");
  TH1F *h_Jet_Btag2_pt20All  = new TH1F("h_Jet_Btag2_pt20All","",100, 0., 1.); h_Jet_Btag2_pt20All->GetXaxis()->SetTitle("Btag 2 (Pt>20 GeV)");
  TH1F *h_Jet_Btag2_pt30All  = new TH1F("h_Jet_Btag2_pt30All","",100, 0., 1.); h_Jet_Btag2_pt30All->GetXaxis()->SetTitle("Btag 2 (Pt>30 GeV)");
  TH1F *h_Jet_Btag1          = new TH1F("h_Jet_Btag1","", 50, 0., 1.); h_Jet_Btag1->GetXaxis()->SetTitle("B-tag 1");
  TH2F *h_Jet_Btag1_eta      = new TH2F("h_Jet_Btag1_eta","", 50, 0., 1., 100, -2.4, 2.4); h_Jet_Btag1_eta->GetXaxis()->SetTitle("B-tag 1"); h_Jet_Btag1_eta->GetYaxis()->SetTitle("#eta");
  TH2F *h_Jet_Btag1_pt       = new TH2F("h_Jet_Btag1_pt","", 50, 0., 1., 100, 0., 100); h_Jet_Btag1_pt->GetXaxis()->SetTitle("B-tag 1"); h_Jet_Btag1_pt->GetYaxis()->SetTitle("pt [GeV]");
  TH1F *h_Jet_Btag2          = new TH1F("h_Jet_Btag2","", 50, 0., 1.); h_Jet_Btag2->GetXaxis()->SetTitle("B-tag 2");
  TH2F *h_Jet_Btag2_eta      = new TH2F("h_Jet_Btag2_eta","", 50, 0., 1., 100, -2.4, 2.4); h_Jet_Btag2_eta->GetXaxis()->SetTitle("B-tag 2"); h_Jet_Btag2_eta->GetYaxis()->SetTitle("#eta");
  TH2F *h_Jet_Btag2_pt       = new TH2F("h_Jet_Btag2_pt","", 50, 0., 1., 100, 0., 100); h_Jet_Btag2_pt->GetXaxis()->SetTitle("B-tag 2"); h_Jet_Btag2_pt->GetYaxis()->SetTitle("pt [GeV]");
  TH1F *h_Jet_Btag3          = new TH1F("h_Jet_Btag3","", 50, 0., 1.); h_Jet_Btag3->GetXaxis()->SetTitle("B-tag 3");
  TH2F *h_Jet_Btag3_eta      = new TH2F("h_Jet_Btag3_eta","", 50, 0., 1., 100, -2.4, 2.4); h_Jet_Btag3_eta->GetXaxis()->SetTitle("B-tag 3"); h_Jet_Btag3_eta->GetYaxis()->SetTitle("#eta");
  TH2F *h_Jet_Btag3_pt       = new TH2F("h_Jet_Btag3_pt","", 50, 0., 1., 100, 0., 100); h_Jet_Btag3_pt->GetXaxis()->SetTitle("B-tag 3"); h_Jet_Btag3_pt->GetYaxis()->SetTitle("pt [GeV]");
  TH1F *h_DR_min             = new TH1F("h_DR_min","", 100, 0., 6.); h_DR_min->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_DR_min_cut         = new TH1F("h_DR_min_cut","", 100, 0., 6.); h_DR_min_cut->GetXaxis()->SetTitle("#Delta R min");
  TH1F *h_isoDR0             = new TH1F("h_isoDR0","", 100, 0., 2.); h_isoDR0->GetXaxis()->SetTitle("Iso.");
  TH2F *h_DR_min_mass        = new TH2F("h_DR_min_mass","", 12, 0., 12., 50, 0., 6.); h_DR_min_mass->GetXaxis()->SetTitle("Mass [GeV]"); h_DR_min_mass->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_DR          = new TH2F("h_DR_min_DR","", 50, 0., 2., 50, 0., 6.); h_DR_min_DR->GetXaxis()->SetTitle("#Delta R(#mu1 #mu 2)"); h_DR_min_DR->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_Iso         = new TH2F("h_DR_min_Iso","", 50, 0., 2., 50, 0., 6.); h_DR_min_Iso->GetXaxis()->SetTitle("Iso."); h_DR_min_Iso->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_Iso_nocut   = new TH2F("h_DR_min_Iso_nocut","", 200, 0., 100., 50, 0., 6.); h_DR_min_Iso_nocut->GetXaxis()->SetTitle("Iso."); h_DR_min_Iso_nocut->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DR_min_IsoOrp      = new TH2F("h_DR_min_IsoOrp","", 100, 0., 50., 50, 0., 6.); h_DR_min_IsoOrp->GetXaxis()->SetTitle("Iso. Orphan"); h_DR_min_IsoOrp->GetYaxis()->SetTitle("#Delta R min");
  TH2F *h_DRmumuor_mass      = new TH2F("h_DRmumuor_mass","", 12, 0., 12., 50, 0., 6.); h_DRmumuor_mass->GetXaxis()->SetTitle("Mass [GeV]"); h_DRmumuor_mass->GetYaxis()->SetTitle("#Delta R(#mu,di-#mu)");
  TH2F *h_DRmumu_mass        = new TH2F("h_DRmumu_mass","", 12, 0., 12., 50, 0., 6.); h_DRmumu_mass->GetXaxis()->SetTitle("Mass [GeV]"); h_DRmumu_mass->GetYaxis()->SetTitle("#Delta R(#mu 1,#mu 2)");

  // Loop over all entries
  while (myReader.Next()) {
    if(fabs(*orph_EtaMu0)<2.4 && fabs(*orph_EtaMu1)<2.4 && fabs(*orph_EtaOrph)<2.4 && (*containstrig>0 || *containstrig2>0) && *orph_dimu_mass<12 && *orph_dimu_mass>0.1){
	// Di-muons - orphan
	TLorentzVector Mu0, Mu1, MuOr, diMu;
	Mu0.SetPtEtaPhiM(*orph_PtMu0,*orph_EtaMu0,*orph_PhiMu0,0.);
	Mu1.SetPtEtaPhiM(*orph_PtMu1,*orph_EtaMu1,*orph_PhiMu1,0.);
	MuOr.SetPtEtaPhiM(*orph_PtOrph,*orph_EtaOrph,*orph_PhiOrph,0.);
	diMu = Mu0 + Mu1;
	float DR0 = Mu0.DeltaR(MuOr);
	float DR1 = Mu1.DeltaR(MuOr);
	float DR_min = DR0<DR1 ? DR0 : DR1;
	h_DR_min_Iso_nocut->Fill(*orph_dimu_isoTk,DR_min);
	if(*orph_dimu_isoTk<2 && *orph_dimu_isoTk>=0) h_DR_min_IsoOrp->Fill(*orph_isoTk,DR_min);
	if(*orph_dimu_isoTk<2 && *orph_dimu_isoTk>=0 && *orph_isoTk<2 && *orph_isoTk>=0){
	  h_DR_min->Fill(DR_min);
	  if(DR_min>Mu0.DeltaR(Mu1)) h_DR_min_cut->Fill(DR_min);
	  h_DR_min_mass->Fill(*orph_dimu_mass,DR_min);
	  h_DRmumuor_mass->Fill(*orph_dimu_mass,diMu.DeltaR(MuOr));
	  h_DRmumu_mass->Fill(*orph_dimu_mass,Mu0.DeltaR(Mu1));
	  h_DR_min_DR->Fill(Mu0.DeltaR(Mu1),DR_min);
	  h_DR_min_Iso->Fill(*orph_dimu_isoTk,DR_min);
	  if( DR_min<0.1 && Mu0.DeltaR(Mu1)<0.1 ) h_isoDR0->Fill(*orph_dimu_isoTk);
	}
	// b-jets
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
	h_Njet_0->Fill(nJet_0); h_Njet_20->Fill(nJet_20); h_Njet_30->Fill(nJet_30);
	if(*orph_dimu_isoTk<2 && *orph_dimu_isoTk>=0 && *orph_isoTk<2 && *orph_isoTk>=0){
	  float mimPt = -1.; int id_Lead = -1;
	  float DR_min2 = 0.3; int NJetCloseOrph=0, CloserJetID=-1;
	  for( int i=0; i<*NPATJet; i++ ){
	    if( fabs(PAT_jet_eta[i])<2.4 && PAT_jet_pt[i]>25. ){
		h_Jet_Btag1->Fill(PAT_jet_Btag1[i]);
		h_Jet_Btag1_pt->Fill(PAT_jet_Btag1[i],PAT_jet_pt[i]);
		h_Jet_Btag1_eta->Fill(PAT_jet_Btag1[i],PAT_jet_eta[i]);
		h_Jet_Btag2->Fill(PAT_jet_Btag2[i]);
		h_Jet_Btag2_pt->Fill(PAT_jet_Btag2[i],PAT_jet_pt[i]);
		h_Jet_Btag2_eta->Fill(PAT_jet_Btag2[i],PAT_jet_eta[i]);
		h_Jet_Btag3->Fill(PAT_jet_Btag3[i]);
		h_Jet_Btag3_pt->Fill(PAT_jet_Btag3[i],PAT_jet_pt[i]);
		h_Jet_Btag3_eta->Fill(PAT_jet_Btag3[i],PAT_jet_eta[i]);
		float thisPt = PAT_jet_pt[i];
		if( thisPt > mimPt  ){
		  mimPt = thisPt;
		  id_Lead = i;
		}
		TLorentzVector ThisJet;
		ThisJet.SetPtEtaPhiM(PAT_jet_pt[i],PAT_jet_eta[i],PAT_jet_phi[i],0.);
		float DR_jet_orph= ThisJet.DeltaR(MuOr);
		if(DR_jet_orph<DR_min2){
		  NJetCloseOrph++;
		  DR_min2=DR_jet_orph;
		  CloserJetID=i;
		}
	    }//Eta && Pt requirement
	  }//For all Jets
	  if(id_Lead>-1) h_JetPt_lead->Fill( PAT_jet_pt[id_Lead] );
	  h_CloseJetsN->Fill(NJetCloseOrph);
	  if(CloserJetID>-1){
	    h_DrMin_ifJet->Fill(DR_min);
	    h_CloseJetDR->Fill(DR_min2);
	    h_CloseJetBtag2->Fill(PAT_jet_Btag2[CloserJetID]);
	    h_DR_min_Btag2Closer->Fill(PAT_jet_Btag2[CloserJetID],DR_min);
	    h_DR_min_Btag3Closer->Fill(PAT_jet_Btag3[CloserJetID],DR_min);
	    if(PAT_jet_Btag2[CloserJetID]>0.8) h_HigBtag2->Fill(DR_min);
	    if(PAT_jet_Btag2[CloserJetID]<0.8) h_LowBtag2->Fill(DR_min);
	  }
	}// Only For isolated dimuons
    }// Control region bb
  }
  gStyle->SetOptStat(0);
  h_Njet_0->Draw();                    myc1->SaveAs("figures/h_Njet_0.pdf");              delete h_Njet_0;
  h_Njet_20->Draw();                   myc1->SaveAs("figures/h_Njet_20.pdf");             delete h_Njet_20;
  h_Njet_30->Draw();                   myc1->SaveAs("figures/h_Njet_30.pdf");             delete h_Njet_30;
  h_JetPt_lead->Draw();                myc1->SaveAs("figures/h_JetPt_lead.pdf");          delete h_JetPt_lead;
  h_Jet_Btag2_lowptAll->Draw();        myc1->SaveAs("figures/h_Jet_Btag2_lowptAll.pdf");  delete h_Jet_Btag2_lowptAll;
  h_Jet_Btag2_pt20All->Draw();         myc1->SaveAs("figures/h_Jet_Btag2_pt20All.pdf");   delete h_Jet_Btag2_pt20All;
  h_Jet_Btag2_pt30All->Draw();         myc1->SaveAs("figures/h_Jet_Btag2_pt30All.pdf");   delete h_Jet_Btag2_pt30All;
  h_Jet_Btag1->Draw();                 myc1->SaveAs("figures/h_Jet_Btag1.pdf");           delete h_Jet_Btag1;
  h_Jet_Btag1_eta->Draw("colz");       myc1->SaveAs("figures/h_Jet_Btag1_eta.pdf");       delete h_Jet_Btag1_eta;
  h_Jet_Btag1_pt->Draw("colz");        myc1->SaveAs("figures/h_Jet_Btag1_pt.pdf");        delete h_Jet_Btag1_pt;
  h_Jet_Btag2->Draw();                 myc1->SaveAs("figures/h_Jet_Btag2.pdf");           delete h_Jet_Btag2;
  h_Jet_Btag2_eta->Draw("colz");       myc1->SaveAs("figures/h_Jet_Btag2_eta.pdf");       delete h_Jet_Btag2_eta;
  h_Jet_Btag2_pt->Draw("colz");        myc1->SaveAs("figures/h_Jet_Btag2_pt.pdf");        delete h_Jet_Btag2_pt;
  h_Jet_Btag3->Draw();                 myc1->SaveAs("figures/h_Jet_Btag3.pdf");           delete h_Jet_Btag3;
  h_Jet_Btag3_eta->Draw("colz");       myc1->SaveAs("figures/h_Jet_Btag3_eta.pdf");       delete h_Jet_Btag3_eta;
  h_Jet_Btag3_pt->Draw("colz");        myc1->SaveAs("figures/h_Jet_Btag3_pt.pdf");        delete h_Jet_Btag3_pt;
  h_CloseJetsN->Draw();                myc1->SaveAs("figures/h_CloseJetsN.pdf");          delete h_CloseJetsN;
  h_DrMin_ifJet->Draw();               myc1->SaveAs("figures/h_DrMin_ifJet.pdf");         delete h_DrMin_ifJet;
  h_CloseJetDR->Draw();                myc1->SaveAs("figures/h_CloseJetDR.pdf");          delete h_CloseJetDR;
  h_CloseJetBtag2->Draw();             myc1->SaveAs("figures/h_CloseJetBtag2.pdf");       delete h_CloseJetBtag2;
  h_DR_min_Btag2Closer->Draw("colz");  myc1->SaveAs("figures/h_DR_min_Btag2Closer.pdf");  delete h_DR_min_Btag2Closer;
  h_DR_min_Btag3Closer->Draw("colz");  myc1->SaveAs("figures/h_DR_min_Btag3Closer.pdf");  delete h_DR_min_Btag3Closer;
  h_HigBtag2->Draw();                  myc1->SaveAs("figures/h_HigBtag2.pdf");            delete h_HigBtag2;
  h_LowBtag2->Draw();                  myc1->SaveAs("figures/h_LowBtag2.pdf");            delete h_LowBtag2;
  h_DR_min->Draw();                    myc1->SaveAs("figures/h_DR_min.pdf"); h_DR_min_cut->SetLineColor(2); h_DR_min_cut->Draw("same");  myc1->SaveAs("figures/h_DR_min_cut.pdf");
  h_isoDR0->Draw();                    myc1->SaveAs("figures/h_isoDR0.pdf");              delete h_isoDR0;
  h_DR_min_mass->Draw("colz");         myc1->SaveAs("figures/h_DR_min_mass.pdf");         TProfile *p_DR_min_mass = h_DR_min_mass->ProfileX(); p_DR_min_mass->SetMinimum(0);   p_DR_min_mass->Draw();  myc1->SaveAs("figures/p_DR_min_mass.pdf");
  h_DR_min_DR->Draw("colz");           myc1->SaveAs("figures/h_DR_min_DR.pdf");           delete h_DR_min_DR;
  h_DR_min_Iso->Draw("colz");          myc1->SaveAs("figures/h_DR_min_Iso.pdf");          delete h_DR_min_Iso;
  h_DR_min_Iso_nocut->Draw("colz");    myc1->SaveAs("figures/h_DR_min_Iso_nocut.pdf");    delete h_DR_min_Iso_nocut;
  h_DR_min_IsoOrp->Draw("colz");       myc1->SaveAs("figures/h_DR_min_IsoOrp.pdf");       delete h_DR_min_IsoOrp;
  h_DRmumuor_mass->Draw("colz");       myc1->SaveAs("figures/h_DRmumuor_mass.pdf");       delete h_DRmumuor_mass;
  h_DRmumu_mass->Draw("colz");         myc1->SaveAs("figures/h_DRmumu_mass.pdf");         delete h_DRmumu_mass;
  delete myc1;
}
