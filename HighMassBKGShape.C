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

//****************************
//* USER modify above ONLY   *
//****************************

void HighMassBKGShape()
{
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);
  THStack *MC_hs_CR_m1 = new THStack("MC_hs_CR_m1","");
  THStack *MC_hs_CR_m2 = new THStack("MC_hs_CR_m2","");
  THStack *MC_hs_SR_m1 = new THStack("MC_hs_SR_m1","");
  THStack *MC_hs_SR_m2 = new THStack("MC_hs_SR_m2","");
  //For plotting summed error for above stacked plots
  TH1F *MC_CR_m1 = new TH1F("MC_CR_m1","",30,0.0,60.0);
  TH1F *MC_CR_m2 = new TH1F("MC_CR_m2","",30,0.0,60.0);
  TH1F *MC_SR_m1 = new TH1F("MC_SR_m1","",30,0.0,60.0);
  TH1F *MC_SR_m2 = new TH1F("MC_SR_m2","",30,0.0,60.0);

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
    MCBKGShapeCRmassC->GetXaxis()->SetRangeUser(10., 60.);
    MC_hs_CR_m1->Add(MCBKGShapeCRmassC);
    MC_CR_m1->Add(MCBKGShapeCRmassC);

    TH1F *MCBKGShapeCRmassF = (TH1F*)file_tmp->Get("BKGShapeCRmassF")->Clone("MCBKGShapeCRmassF");
    MCBKGShapeCRmassF->Scale(MC_ScaleFactors[i]);
    MCBKGShapeCRmassF->SetLineColor(MC_Colors[i]);
    MCBKGShapeCRmassF->SetFillColor(MC_Colors[i]);
    MCBKGShapeCRmassF->SetMarkerColor(MC_Colors[i]);
    MCBKGShapeCRmassF->GetXaxis()->SetRangeUser(10., 60.);
    MC_hs_CR_m2->Add(MCBKGShapeCRmassF);
    MC_CR_m2->Add(MCBKGShapeCRmassF);

    TH1F *MCBKGShapeSRmassC = (TH1F*)file_tmp->Get("BKGShapeSRmassC")->Clone("MCBKGShapeSRmassC");
    MCBKGShapeSRmassC->Scale(MC_ScaleFactors[i]);
    MCBKGShapeSRmassC->SetLineColor(MC_Colors[i]);
    MCBKGShapeSRmassC->SetFillColor(MC_Colors[i]);
    MCBKGShapeSRmassC->SetMarkerColor(MC_Colors[i]);
    MCBKGShapeSRmassC->GetXaxis()->SetRangeUser(10., 60.);
    MC_hs_SR_m1->Add(MCBKGShapeSRmassC);
    MC_SR_m1->Add(MCBKGShapeSRmassC);

    TH1F *MCBKGShapeSRmassF = (TH1F*)file_tmp->Get("BKGShapeSRmassF")->Clone("MCBKGShapeSRmassF");
    MCBKGShapeSRmassF->Scale(MC_ScaleFactors[i]);
    MCBKGShapeSRmassF->SetLineColor(MC_Colors[i]);
    MCBKGShapeSRmassF->SetFillColor(MC_Colors[i]);
    MCBKGShapeSRmassF->SetMarkerColor(MC_Colors[i]);
    MCBKGShapeSRmassF->GetXaxis()->SetRangeUser(10., 60.);
    MC_hs_SR_m2->Add(MCBKGShapeSRmassF);
    MC_SR_m2->Add(MCBKGShapeSRmassF);
  }

  //DATA
  TString DATA_file_name;
  TH1F *DATA_CR_m1 = new TH1F("DATA_CR_m1","",30,0.0,60.0);
  TH1F *DATA_CR_m2 = new TH1F("DATA_CR_m2","",30,0.0,60.0);
  TH1F *DATA_SR_m1 = new TH1F("DATA_SR_m1","",30,0.0,60.0);
  TH1F *DATA_SR_m2 = new TH1F("DATA_SR_m2","",30,0.0,60.0);

  for (int j = 0; j < 4; j++) {
    DATA_file_name.Form( "%s/%s", store.Data(), DATA_files[j].Data() );
    std::cout << "Opening file #"<< j+1 << ": " << DATA_file_name.Data() << std::endl;
    file_tmp = TFile::Open(DATA_file_name);
    if (!file_tmp) {
      std::cout << "ERROR: could not open file " << DATA_file_name.Data() << std::endl;
      return;
    }

    TH1F *DATABKGShapeCRmassC = (TH1F*)file_tmp->Get("BKGShapeCRmassC")->Clone("DATABKGShapeCRmassC"); DATABKGShapeCRmassC->GetXaxis()->SetRangeUser(10., 60.); DATA_CR_m1->Add(DATABKGShapeCRmassC);
    TH1F *DATABKGShapeCRmassF = (TH1F*)file_tmp->Get("BKGShapeCRmassF")->Clone("DATABKGShapeCRmassF"); DATABKGShapeCRmassF->GetXaxis()->SetRangeUser(10., 60.); DATA_CR_m2->Add(DATABKGShapeCRmassF);
    TH1F *DATABKGShapeSRmassC = (TH1F*)file_tmp->Get("BKGShapeSRmassC")->Clone("DATABKGShapeSRmassC"); DATABKGShapeSRmassC->GetXaxis()->SetRangeUser(10., 60.); DATA_SR_m1->Add(DATABKGShapeSRmassC);
    TH1F *DATABKGShapeSRmassF = (TH1F*)file_tmp->Get("BKGShapeSRmassF")->Clone("DATABKGShapeSRmassF"); DATABKGShapeSRmassF->GetXaxis()->SetRangeUser(10., 60.); DATA_SR_m2->Add(DATABKGShapeSRmassF);
  }

  //write to output file
  TFile myPlot("HighMassBKGShape_FINAL.root","RECREATE");

  TCanvas *CR1=new TCanvas("CR1","CR m1",700,500);
  CR1->cd();
  MC_hs_CR_m1->Draw("HIST");
  MC_hs_CR_m1->SetMaximum(50);
  MC_hs_CR_m1->GetXaxis()->SetRangeUser(10, 60);
  MC_hs_CR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  MC_hs_CR_m1->GetYaxis()->SetTitle("Events/2GeV");
  //Draw rectangle error
  MC_CR_m1->SetLineColor(2);
  MC_CR_m1->SetFillColor(2);
  MC_CR_m1->SetFillStyle(3004);
  MC_CR_m1->Draw("E2 SAME");
  Double_t MC_CR_m1_error;
  Double_t MC_CR_m1_integral = MC_CR_m1->IntegralAndError(6, 30, MC_CR_m1_error, ""); //nbins=30, 1,2...30, 6-30 integrates from 10-60 GeV
  std::cout << "MC CR m1 integral = " << MC_CR_m1_integral << " +/- " << MC_CR_m1_error << std::endl;
  DATA_CR_m1->SetFillColor(1);
  DATA_CR_m1->SetLineColor(1);
  DATA_CR_m1->SetMarkerStyle(20);
  DATA_CR_m1->Draw("E1 X0 SAME");//Draw Error bars
  DATA_CR_m1->GetXaxis()->SetRangeUser(10, 60);
  Double_t DATA_CR_m1_error;
  Double_t DATA_CR_m1_integral = DATA_CR_m1->IntegralAndError(6, 30, DATA_CR_m1_error, ""); // "" ... or ... "width"
  std::cout << "DATA CR m1 integral = " << DATA_CR_m1_integral << " +/- " << DATA_CR_m1_error << std::endl;
  gPad->RedrawAxis();
  CR1->Write();

  TCanvas *CR2=new TCanvas("CR2","CR m2",700,500);
  CR2->cd();
  MC_hs_CR_m2->Draw("HIST");
  MC_hs_CR_m2->SetMaximum(50);
  MC_hs_CR_m2->GetXaxis()->SetRangeUser(10, 60);
  MC_hs_CR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  MC_hs_CR_m2->GetYaxis()->SetTitle("Events/2GeV");
  MC_CR_m2->SetLineColor(2);
  MC_CR_m2->SetFillColor(2);
  MC_CR_m2->SetFillStyle(3004);
  MC_CR_m2->Draw("E2 SAME");
  Double_t MC_CR_m2_error;
  Double_t MC_CR_m2_integral = MC_CR_m2->IntegralAndError(6, 30, MC_CR_m2_error, ""); // "" ... or ... "width"
  std::cout << "MC CR m2 integral = " << MC_CR_m2_integral << " +/- " << MC_CR_m2_error << std::endl;
  DATA_CR_m2->SetFillColor(1);
  DATA_CR_m2->SetLineColor(1);
  DATA_CR_m2->SetMarkerStyle(20);
  DATA_CR_m2->Draw("E1 X0 SAME");
  DATA_CR_m2->GetXaxis()->SetRangeUser(10, 60);
  Double_t DATA_CR_m2_error;
  Double_t DATA_CR_m2_integral = DATA_CR_m2->IntegralAndError(6, 30, DATA_CR_m2_error, ""); // "" ... or ... "width"
  std::cout << "DATA CR m2 integral = " << DATA_CR_m2_integral << " +/- " << DATA_CR_m2_error << std::endl;
  gPad->RedrawAxis();
  CR2->Write();

  //Don't Draw Data at SR, blinded until approval
  TCanvas *SR1=new TCanvas("SR1","SR m1",700,500);
  SR1->cd();
  MC_hs_SR_m1->Draw("HIST");
  MC_hs_SR_m1->SetMaximum(10);
  MC_hs_SR_m1->GetXaxis()->SetRangeUser(10, 60);
  MC_hs_SR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  MC_hs_SR_m1->GetYaxis()->SetTitle("Events/2GeV");
  MC_SR_m1->SetLineColor(2);
  MC_SR_m1->SetFillColor(2);
  MC_SR_m1->SetFillStyle(3004);
  MC_SR_m1->Draw("E2 SAME");
  Double_t MC_SR_m1_error;
  Double_t MC_SR_m1_integral = MC_SR_m1->IntegralAndError(6, 30, MC_SR_m1_error, ""); // "" ... or ... "width"
  std::cout << "MC SR m1 integral = " << MC_SR_m1_integral << " +/- " << MC_SR_m1_error << std::endl;
  gPad->RedrawAxis();
  SR1->Write();

  TCanvas *SR2=new TCanvas("SR2","SR m2",700,500);
  SR2->cd();
  MC_hs_SR_m2->Draw("HIST");
  MC_hs_SR_m2->SetMaximum(10);
  MC_hs_SR_m2->GetXaxis()->SetRangeUser(10, 60);
  MC_hs_SR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  MC_hs_SR_m2->GetYaxis()->SetTitle("Events/2GeV");
  MC_SR_m2->SetLineColor(2);
  MC_SR_m2->SetFillColor(2);
  MC_SR_m2->SetFillStyle(3004);
  MC_SR_m2->Draw("E2 SAME");
  Double_t MC_SR_m2_error;
  Double_t MC_SR_m2_integral = MC_SR_m2->IntegralAndError(6, 30, MC_SR_m2_error, ""); // "" ... or ... "width"
  std::cout << "MC SR m2 integral = " << MC_SR_m2_integral << " +/- " << MC_SR_m2_error << std::endl;
  gPad->RedrawAxis();
  SR2->Write();

  myPlot.Close();

} // End function: void
