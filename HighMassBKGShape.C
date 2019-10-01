//***********************************************************************************************
//* For background shape modeling by comparing MC to DATA at Control Region                     *
//* Source most recent root version:                                                            *
//* . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh *
//* Run it as: root -l -b -q HighMassBKGShape.C                                                 *
//*                    Wei Shi @Sep 24, 2019, Rice U.                                           *
//***********************************************************************************************
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>
#include "stdio.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TChain.h"
#include "TBranch.h"

//****************************
//* USER modify section      *
//****************************
//The order is important here: DY 0-2J, ZZTo4L, TTJets_DiLept, ggHToZZTo4L
//Same order as ScaleFactors
TString MC_files[6] = {
  "HighMassBKGShape_DYToLL_0J.root",
  "HighMassBKGShape_DYToLL_1J.root",
  "HighMassBKGShape_DYToLL_2J.root",
  "HighMassBKGShape_ZZTo4L.root",
  "HighMassBKGShape_TTJets_DiLept.root",
  "HighMassBKGShape_ggHToZZTo4L.root"
};

//Scale MC events to Data, not including analysis SF like trigger, muon id etc
Float_t MC_ScaleFactors[6]={2.2398524E+00, 3.7937427E-01, 2.4801970E-01, 3.0278414E-03, 7.0805939E-02, 4.6355268E-04};
Color_t MC_Colors[6]={20, 30, 40, 7, 8, 9};

TString DATA_files[4] = {
  "HighMassBKGShape_2017C.root",
  "HighMassBKGShape_2017D.root",
  "HighMassBKGShape_2017E.root",
  "HighMassBKGShape_2017F.root"
};

TString store = "/fdata/hepx/store/user/wshi/SMBKGatHighMass"; //input dir
const double       m_min  = 11.0;
const double       m_max  = 59.0;
const unsigned int m_bins = 24;
//****************************
//* USER modify above ONLY   *
//****************************

void HighMassBKGShape()
{
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);
  THStack *MC_hs_CR_m1 = new THStack("MC_hs_CR_m1", "");
  THStack *MC_hs_CR_m2 = new THStack("MC_hs_CR_m2", "");
  THStack *MC_hs_SR_m1 = new THStack("MC_hs_SR_m1", "");
  THStack *MC_hs_SR_m2 = new THStack("MC_hs_SR_m2", "");
  //For plotting summed error for above stacked plots
  TH1F *MC_CR_m1 = new TH1F("MC_CR_m1", "", m_bins, m_min, m_max);
  TH1F *MC_CR_m2 = new TH1F("MC_CR_m2", "", m_bins, m_min, m_max);
  TH1F *MC_SR_m1 = new TH1F("MC_SR_m1", "", m_bins, m_min, m_max);
  TH1F *MC_SR_m2 = new TH1F("MC_SR_m2", "", m_bins, m_min, m_max);

  //SM BKG MC
  TString MC_file_name;
  for (int i = 0; i < 6; i++) {
    MC_file_name.Form( "%s/%s", store.Data(), MC_files[i].Data() );
    std::cout << "Opening file #"<< i+1 << ": " << MC_file_name.Data() << std::endl;
    file_tmp = TFile::Open(MC_file_name);
    if (!file_tmp) {
      std::cout << "ERROR: could not open file " << MC_file_name.Data() << std::endl;
      return;
    }

    TH1F *MCBKGShapeCRmassC = (TH1F*)file_tmp->Get("BKGShapeCRmassC")->Clone("MCBKGShapeCRmassC");
    MCBKGShapeCRmassC->Scale(MC_ScaleFactors[i]);
    MCBKGShapeCRmassC->SetLineColor(MC_Colors[i]);
    MCBKGShapeCRmassC->SetFillColor(MC_Colors[i]);
    MCBKGShapeCRmassC->SetMarkerColor(MC_Colors[i]);
    MC_hs_CR_m1->Add(MCBKGShapeCRmassC);
    MC_CR_m1->Add(MCBKGShapeCRmassC);

    TH1F *MCBKGShapeCRmassF = (TH1F*)file_tmp->Get("BKGShapeCRmassF")->Clone("MCBKGShapeCRmassF");
    MCBKGShapeCRmassF->Scale(MC_ScaleFactors[i]);
    MCBKGShapeCRmassF->SetLineColor(MC_Colors[i]);
    MCBKGShapeCRmassF->SetFillColor(MC_Colors[i]);
    MCBKGShapeCRmassF->SetMarkerColor(MC_Colors[i]);
    MC_hs_CR_m2->Add(MCBKGShapeCRmassF);
    MC_CR_m2->Add(MCBKGShapeCRmassF);

    TH1F *MCBKGShapeSRmassC = (TH1F*)file_tmp->Get("BKGShapeSRmassC")->Clone("MCBKGShapeSRmassC");
    MCBKGShapeSRmassC->Scale(MC_ScaleFactors[i]);
    MCBKGShapeSRmassC->SetLineColor(MC_Colors[i]);
    MCBKGShapeSRmassC->SetFillColor(MC_Colors[i]);
    MCBKGShapeSRmassC->SetMarkerColor(MC_Colors[i]);
    MC_hs_SR_m1->Add(MCBKGShapeSRmassC);
    MC_SR_m1->Add(MCBKGShapeSRmassC);

    TH1F *MCBKGShapeSRmassF = (TH1F*)file_tmp->Get("BKGShapeSRmassF")->Clone("MCBKGShapeSRmassF");
    MCBKGShapeSRmassF->Scale(MC_ScaleFactors[i]);
    MCBKGShapeSRmassF->SetLineColor(MC_Colors[i]);
    MCBKGShapeSRmassF->SetFillColor(MC_Colors[i]);
    MCBKGShapeSRmassF->SetMarkerColor(MC_Colors[i]);
    MC_hs_SR_m2->Add(MCBKGShapeSRmassF);
    MC_SR_m2->Add(MCBKGShapeSRmassF);
  }

  //DATA
  TString DATA_file_name;
  TH1F *DATA_CR_m1 = new TH1F("DATA_CR_m1", "", m_bins, m_min, m_max);
  TH1F *DATA_CR_m2 = new TH1F("DATA_CR_m2", "", m_bins, m_min, m_max);
  TH1F *DATA_SR_m1 = new TH1F("DATA_SR_m1", "", m_bins, m_min, m_max);
  TH1F *DATA_SR_m2 = new TH1F("DATA_SR_m2", "", m_bins, m_min, m_max);

  for (int j = 0; j < 4; j++) {
    DATA_file_name.Form( "%s/%s", store.Data(), DATA_files[j].Data() );
    std::cout << "Opening file #"<< j+1 << ": " << DATA_file_name.Data() << std::endl;
    file_tmp = TFile::Open(DATA_file_name);
    if (!file_tmp) {
      std::cout << "ERROR: could not open file " << DATA_file_name.Data() << std::endl;
      return;
    }

    TH1F *DATABKGShapeCRmassC = (TH1F*)file_tmp->Get("BKGShapeCRmassC")->Clone("DATABKGShapeCRmassC"); DATA_CR_m1->Add(DATABKGShapeCRmassC);
    TH1F *DATABKGShapeCRmassF = (TH1F*)file_tmp->Get("BKGShapeCRmassF")->Clone("DATABKGShapeCRmassF"); DATA_CR_m2->Add(DATABKGShapeCRmassF);
    TH1F *DATABKGShapeSRmassC = (TH1F*)file_tmp->Get("BKGShapeSRmassC")->Clone("DATABKGShapeSRmassC"); DATA_SR_m1->Add(DATABKGShapeSRmassC);
    TH1F *DATABKGShapeSRmassF = (TH1F*)file_tmp->Get("BKGShapeSRmassF")->Clone("DATABKGShapeSRmassF"); DATA_SR_m2->Add(DATABKGShapeSRmassF);
  }

  //write to output file
  TFile myPlot("HighMassBKGShape_FINAL.root","RECREATE");

  TCanvas *CR1=new TCanvas("CR1","CR m1",700,500);
  CR1->Clear();
  //MC vs DATA
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  MC_hs_CR_m1->Draw("HIST");
  MC_hs_CR_m1->SetMaximum(50);
  MC_hs_CR_m1->GetYaxis()->SetTitle("Events/2GeV");
  //Plotting error
  MC_CR_m1->SetLineColor(2);
  MC_CR_m1->SetFillColor(2);
  MC_CR_m1->SetFillStyle(3004);
  MC_CR_m1->Draw("E2 SAME");
  Double_t MC_CR_m1_error;
  Double_t MC_CR_m1_integral = MC_CR_m1->IntegralAndError(1, m_bins, MC_CR_m1_error, "");
  std::cout << "MC CR m1 integral = " << MC_CR_m1_integral << " +/- " << MC_CR_m1_error << std::endl;
  DATA_CR_m1->SetFillColor(1);
  DATA_CR_m1->SetLineColor(1);
  DATA_CR_m1->SetMarkerStyle(20);
  DATA_CR_m1->Draw("E1 X0 SAME");//Draw Error bars
  Double_t DATA_CR_m1_error;
  Double_t DATA_CR_m1_integral = DATA_CR_m1->IntegralAndError(1, m_bins, DATA_CR_m1_error, ""); // "": width
  std::cout << "DATA CR m1 integral = " << DATA_CR_m1_integral << " +/- " << DATA_CR_m1_error << std::endl;
  gPad->RedrawAxis();
  CR1->cd(); CR1->Update();
  //For plotting pull distribution
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  //fill pull histogram
  TH1F *pull_CR_m1 = new TH1F("pull_CR_m1","", m_bins, m_min, m_max);
  for(unsigned int iB=1; iB<=m_bins; iB++){
    std::cout << "Bin #" << iB << "; Bin center: " << DATA_CR_m1->GetXaxis()->GetBinCenter(iB)<<" GeV" << std::endl;//Debug
    float pull_CR_m1_iB = ( DATA_CR_m1->GetBinContent(iB) - MC_CR_m1->GetBinContent(iB) ) / DATA_CR_m1->GetBinError(iB);//pull definition
    pull_CR_m1->SetBinContent(iB, pull_CR_m1_iB );
  }
  pull_CR_m1->GetYaxis()->SetTitle("Pull");
  pull_CR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  pull_CR_m1->SetMinimum(-4);
  pull_CR_m1->SetMaximum(4);
  pull_CR_m1->SetMarkerStyle(20);
  pull_CR_m1->Draw("P");
  CR1->Write();

  TCanvas *CR2=new TCanvas("CR2","CR m2",700,500);
  CR2->cd();
  MC_hs_CR_m2->Draw("HIST");
  MC_hs_CR_m2->SetMaximum(50);
  MC_hs_CR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  MC_hs_CR_m2->GetYaxis()->SetTitle("Events/2GeV");
  MC_CR_m2->SetLineColor(2);
  MC_CR_m2->SetFillColor(2);
  MC_CR_m2->SetFillStyle(3004);
  MC_CR_m2->Draw("E2 SAME");
  Double_t MC_CR_m2_error;
  Double_t MC_CR_m2_integral = MC_CR_m2->IntegralAndError(1, m_bins, MC_CR_m2_error, "");
  std::cout << "MC CR m2 integral = " << MC_CR_m2_integral << " +/- " << MC_CR_m2_error << std::endl;
  DATA_CR_m2->SetFillColor(1);
  DATA_CR_m2->SetLineColor(1);
  DATA_CR_m2->SetMarkerStyle(20);
  DATA_CR_m2->Draw("E1 X0 SAME");
  Double_t DATA_CR_m2_error;
  Double_t DATA_CR_m2_integral = DATA_CR_m2->IntegralAndError(1, m_bins, DATA_CR_m2_error, "");
  std::cout << "DATA CR m2 integral = " << DATA_CR_m2_integral << " +/- " << DATA_CR_m2_error << std::endl;
  gPad->RedrawAxis();
  CR2->Write();

  //Don't Draw Data at SR, blinded until approval
  TCanvas *SR1=new TCanvas("SR1","SR m1",700,500);
  SR1->cd();
  MC_hs_SR_m1->Draw("HIST");
  MC_hs_SR_m1->SetMaximum(10);
  MC_hs_SR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  MC_hs_SR_m1->GetYaxis()->SetTitle("Events/2GeV");
  MC_SR_m1->SetLineColor(2);
  MC_SR_m1->SetFillColor(2);
  MC_SR_m1->SetFillStyle(3004);
  MC_SR_m1->Draw("E2 SAME");
  Double_t MC_SR_m1_error;
  Double_t MC_SR_m1_integral = MC_SR_m1->IntegralAndError(1, m_bins, MC_SR_m1_error, "");
  std::cout << "MC SR m1 integral = " << MC_SR_m1_integral << " +/- " << MC_SR_m1_error << std::endl;
  gPad->RedrawAxis();
  SR1->Write();

  TCanvas *SR2=new TCanvas("SR2","SR m2",700,500);
  SR2->cd();
  MC_hs_SR_m2->Draw("HIST");
  MC_hs_SR_m2->SetMaximum(10);
  MC_hs_SR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  MC_hs_SR_m2->GetYaxis()->SetTitle("Events/2GeV");
  MC_SR_m2->SetLineColor(2);
  MC_SR_m2->SetFillColor(2);
  MC_SR_m2->SetFillStyle(3004);
  MC_SR_m2->Draw("E2 SAME");
  Double_t MC_SR_m2_error;
  Double_t MC_SR_m2_integral = MC_SR_m2->IntegralAndError(1, m_bins, MC_SR_m2_error, "");
  std::cout << "MC SR m2 integral = " << MC_SR_m2_integral << " +/- " << MC_SR_m2_error << std::endl;
  gPad->RedrawAxis();
  SR2->Write();

  myPlot.Close();

} // End function: void
