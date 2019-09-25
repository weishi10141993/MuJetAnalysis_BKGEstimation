//*****************************************************************************************************
//* For background shape modeling by comparing MC to DATA at Control Region                           *
//* Use: . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh  *
//*                                       Wei Shi @Sep 24, 2019, Rice U.                              *
//*****************************************************************************************************
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>
#include "stdio.h"
#include "math.h"
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
TString SM_BKG_MC_files[6] = {
  "HighMassBKGShape_DYToLL_0J.root",
  "HighMassBKGShape_DYToLL_1J.root",
  "HighMassBKGShape_DYToLL_2J.root",
  "HighMassBKGShape_ZZTo4L.root",
  "HighMassBKGShape_TTJets_DiLept.root",
  "HighMassBKGShape_ggHToZZTo4L.root"
};

//Scale MC events to Data, not including analysis SF like trigger, muon id etc
Float_t SM_BKG_MC_ScaleFactors[6]={2.2398524E+00, 3.7937427E-01, 2.4801970E-01, 3.0278414E-03, 7.0805939E-02, 4.6355268E-04};
Color_t SM_BKG_MC_Colors[6]={20, 30, 40, 7, 8, 9};

TString DATA_files[4] = {
  "HighMassBKGShape_2017C.root",
  "HighMassBKGShape_2017D.root",
  "HighMassBKGShape_2017E.root",
  "HighMassBKGShape_2017F.root"
};

TString store = "/fdata/hepx/store/user/wshi/SMBKGatHighMass/"; //main dir
TString outFile = "/fdata/hepx/store/user/wshi/SMBKGatHighMass/"; //output

//****************************
//* USER modify above ONLY   *
//****************************

void HighMassBKGShape()
{
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);
  TString DATA_file_name;
  TString SM_BKG_MC_file_name;
  THStack *SM_BKG_MC_hs_CR_m1 = new THStack("SM_BKG_MC_hs_CR_m1","");
  THStack *SM_BKG_MC_hs_CR_m2 = new THStack("SM_BKG_MC_hs_CR_m2","");
  THStack *SM_BKG_MC_hs_SR_m1 = new THStack("SM_BKG_MC_hs_SR_m1","");
  THStack *SM_BKG_MC_hs_SR_m2 = new THStack("SM_BKG_MC_hs_SR_m2","");
  TH1F *DATA_CR_m1 = new TH1F("DATA_CR_m1","",30,0.0,60.0);
  TH1F *DATA_CR_m2 = new TH1F("DATA_CR_m1","",30,0.0,60.0);
  TH1F *DATA_SR_m1 = new TH1F("DATA_SR_m1","",30,0.0,60.0);
  TH1F *DATA_SR_m2 = new TH1F("DATA_SR_m2","",30,0.0,60.0);

  //SM BKG MC
  for (int i = 0; i < 6; i++) {
    SM_BKG_MC_file_name.Form( "%s/%s", store.Data(), SM_BKG_MC_files[i].Data() );
    std::cout << "Adding file #"<< i+1 << ": " << SM_BKG_MC_file_name.Data() << std::endl;
    file_tmp = TFile::Open(SM_BKG_MC_file_name);

    BKGShapeCRmassC->Scale(SM_BKG_MC_ScaleFactors[i]);
    BKGShapeCRmassC->SetFillColor(SM_BKG_MC_Colors[i]);
    BKGShapeCRmassC->SetLineColor(SM_BKG_MC_Colors[i]);
    SM_BKG_MC_hs_CR_m1->Add(BKGShapeCRmassC);

    BKGShapeCRmassF->Scale(SM_BKG_MC_ScaleFactors[i]);
    BKGShapeCRmassF->SetFillColor(SM_BKG_MC_Colors[i]);
    BKGShapeCRmassF->SetLineColor(SM_BKG_MC_Colors[i]);
    SM_BKG_MC_hs_CR_m2->Add(BKGShapeCRmassF);

    BKGShapeSRmassC->Scale(SM_BKG_MC_ScaleFactors[i]);
    BKGShapeSRmassC->SetFillColor(SM_BKG_MC_Colors[i]);
    BKGShapeSRmassC->SetLineColor(SM_BKG_MC_Colors[i]);
    SM_BKG_MC_hs_SR_m1->Add(BKGShapeSRmassC);

    BKGShapeSRmassF->Scale(SM_BKG_MC_ScaleFactors[i]);
    BKGShapeSRmassF->SetFillColor(SM_BKG_MC_Colors[i]);
    BKGShapeSRmassF->SetLineColor(SM_BKG_MC_Colors[i]);
    SM_BKG_MC_hs_SR_m2->Add(BKGShapeSRmassF);

    file_tmp->Close();
  }

  //DATA
  for (int j = 0; j < 4; j++) {
    DATA_file_name.Form( "%s/%s", store.Data(), DATA_files[j].Data() );
    std::cout << "Adding file #"<< j+1 << ": " << DATA_file_name.Data() << std::endl;
    file_tmp = TFile::Open(DATA_file_name);

    DATA_CR_m1->Add(BKGShapeCRmassC);
    DATA_CR_m2->Add(BKGShapeCRmassF);
    DATA_SR_m1->Add(BKGShapeSRmassC);
    DATA_SR_m2->Add(BKGShapeSRmassF);

    file_tmp->Close();

  }

  //write to output file
  outFile = outFile + "HighMassBKGShape_FINAL.root";
  TFile myPlot(outFile,"RECREATE");

  TCanvas *CR1=new TCanvas("CR1","CR m1",700,500);
  CR1->cd();
  SM_BKG_MC_hs_CR_m1->Draw();
  DATA_CR_m1->Draw("ESAME");//Draw Error bars
  CR1->Write();

  TCanvas *CR2=new TCanvas("CR2","CR m2",700,500);
  CR2->cd();
  SM_BKG_MC_hs_CR_m2->Draw();
  DATA_CR_m2->Draw("ESAME");
  CR2->Write();

  TCanvas *SR1=new TCanvas("SR1","SR m1",700,500);
  SR1->cd();
  SM_BKG_MC_hs_SR_m1->Draw();
  DATA_SR_m1->Draw("ESAME");
  SR1->Write();

  TCanvas *SR2=new TCanvas("SR2","SR m2",700,500);
  SR2->cd();
  SM_BKG_MC_hs_SR_m2->Draw();
  DATA_SR_m2->Draw("ESAME");
  SR2->Write();

} // End function: void
