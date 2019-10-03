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
TString MC_files[7] = {
  "HighMassBKGShape_DYToLL_0J.root",
  "HighMassBKGShape_DYToLL_1J.root",
  "HighMassBKGShape_DYToLL_2J.root",
  "HighMassBKGShape_ZZTo4L.root",
  "HighMassBKGShape_TTJets_DiLept.root",
  "HighMassBKGShape_ggHToZZTo4L.root",
  "HighMassBKGShape_ggToZZTo4mu.root"
};

//Scale MC events to Data, not including analysis SF like trigger, muon id etc
Float_t MC_ScaleFactors[7]={2.2398524E+00, 3.7937427E-01, 2.4801970E-01, 3.0278414E-03, 7.0805939E-02, 4.6355268E-04, 6.3245328E-05};
Color_t MC_Colors[7]={20, 30, 40, 9, 8, 7, 6};
TString DATA_files[4] = {
  "HighMassBKGShape_2017C.root",
  "HighMassBKGShape_2017D.root",
  "HighMassBKGShape_2017E.root",
  "HighMassBKGShape_2017F.root"
};

TString store = "/fdata/hepx/store/user/wshi/SMBKGatHighMass"; //input dir
const double       m_min  = 11.0;
const double       m_max  = 59.0;
const unsigned int m_bins = 12;//bin size 4GeV
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
  for (int i = 0; i < 7; i++) {
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

  //***************
  //* For m1 at CR*
  //***************
  TCanvas *CR1=new TCanvas("CR1","CR m1",700,500); CR1->Clear();
  TPad *CR1pad1 = new TPad("CR1pad1", "CR1pad1", 0, 0.3, 1, 1.0);//xlow, ylow, xup, yup
  CR1pad1->SetBottomMargin(0); CR1pad1->Draw();
  TPad *CR1pad2 = new TPad("CR1pad2", "CR1pad2", 0, 0.0, 1, 0.29);
  CR1pad2->SetTopMargin(0); CR1pad2->SetBottomMargin(0.3); CR1pad2->SetGridy(); CR1pad2->Draw();
  //MC vs DATA
  CR1pad1->cd();
  //Plot stacked histogram from MC
  MC_hs_CR_m1->Draw("HIST"); MC_hs_CR_m1->SetMaximum(50); MC_hs_CR_m1->GetYaxis()->SetTitle("Events/4GeV");
  //Plot MC error
  MC_CR_m1->SetLineColor(2); MC_CR_m1->SetFillColor(2); MC_CR_m1->SetFillStyle(3004); MC_CR_m1->Draw("E2 SAME");
  Double_t MC_CR_m1_error;
  Double_t MC_CR_m1_integral = MC_CR_m1->IntegralAndError(1, m_bins, MC_CR_m1_error, "");
  std::cout << "MC CR m1 integral = " << MC_CR_m1_integral << " +/- " << MC_CR_m1_error << std::endl;
  //Overlay data
  DATA_CR_m1->SetFillColor(1); DATA_CR_m1->SetLineColor(1); DATA_CR_m1->SetMarkerStyle(20); DATA_CR_m1->Draw("E1 X0 SAME");//Draw Error bars
  Double_t DATA_CR_m1_error;
  Double_t DATA_CR_m1_integral = DATA_CR_m1->IntegralAndError(1, m_bins, DATA_CR_m1_error, ""); // "": width
  std::cout << "DATA CR m1 integral = " << DATA_CR_m1_integral << " +/- " << DATA_CR_m1_error << std::endl;
  //Build Legend
  TLegend* CR1pad1L = CR1pad1->BuildLegend();
  CR1pad1L->SetBorderSize(0); CR1pad1L->SetFillStyle(0); CR1pad1L->SetNColumns(2);
  TList *CR1pad1P = CR1pad1L->GetListOfPrimitives();
  TIter CR1pad1next(CR1pad1P);
  TObject *CR1pad1obj;
  TLegendEntry *CR1pad1li;
  int CR1pad1iEntry = 0;
  while ((CR1pad1obj = CR1pad1next())) {
    CR1pad1li = (TLegendEntry*)CR1pad1obj;
    CR1pad1iEntry++;
    if (CR1pad1iEntry==1) CR1pad1li->SetLabel("DYToLL (0J)");
    if (CR1pad1iEntry==2) CR1pad1li->SetLabel("DYToLL (1J)");
    if (CR1pad1iEntry==3) CR1pad1li->SetLabel("DYToLL (2J)");
    if (CR1pad1iEntry==4) CR1pad1li->SetLabel("qqToZZTo4L");
    if (CR1pad1iEntry==5) CR1pad1li->SetLabel("TTJetsToLL");
    if (CR1pad1iEntry==6) CR1pad1li->SetLabel("ggHToZZTo4L");
    if (CR1pad1iEntry==7) CR1pad1li->SetLabel("ggToZZTo4mu");
    if (CR1pad1iEntry==8) {CR1pad1li->SetLabel("MC Error"); CR1pad1li->SetOption("f");}
    if (CR1pad1iEntry==9) {CR1pad1li->SetLabel("Data"); CR1pad1li->SetOption("ep");}
  }
  CR1pad1->Update(); CR1pad1L->SetX1NDC(0.15); CR1pad1L->SetX2NDC(0.5); CR1pad1L->SetY1NDC(0.5); CR1pad1L->SetY2NDC(0.9); CR1pad1->Modified();
  gPad->RedrawAxis();
  CR1->cd(); CR1->Update();
  //Plot pull distribution
  CR1pad2->cd();
  //fill pull histogram
  TH1F *pull_CR_m1 = new TH1F("pull_CR_m1","", m_bins, m_min, m_max);
  for(unsigned int iB=1; iB<=m_bins; iB++){
    float pull_CR_m1_iB = ( DATA_CR_m1->GetBinContent(iB) - MC_CR_m1->GetBinContent(iB) ) / DATA_CR_m1->GetBinError(iB);//pull definition
    pull_CR_m1->SetBinContent(iB, pull_CR_m1_iB );//iB starts from #1
  }
  pull_CR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  pull_CR_m1->GetXaxis()->SetTitleSize(15);
  pull_CR_m1->GetXaxis()->SetTitleFont(43);
  pull_CR_m1->GetXaxis()->SetTitleOffset(3.0);
  pull_CR_m1->GetXaxis()->SetLabelSize(15);// labels will be 15 pixels
  pull_CR_m1->GetXaxis()->SetLabelFont(43);// Absolute font size in pixel (precision 3)
  pull_CR_m1->GetYaxis()->SetTitle("Pull");
  pull_CR_m1->GetYaxis()->CenterTitle();
  pull_CR_m1->GetYaxis()->SetTitleSize(15);
  pull_CR_m1->GetYaxis()->SetTitleFont(43);
  pull_CR_m1->GetYaxis()->SetTitleOffset(.9);
  pull_CR_m1->GetYaxis()->SetLabelSize(15);
  pull_CR_m1->GetYaxis()->SetLabelFont(43);
  pull_CR_m1->SetMinimum(-3.5);
  pull_CR_m1->SetMaximum(3.5);
  pull_CR_m1->SetStats(0);
  pull_CR_m1->SetMarkerStyle(20);
  pull_CR_m1->Draw("P");
  CR1->Write();

  //***************
  //* For m2 at CR*
  //***************
  TCanvas *CR2=new TCanvas("CR2","CR m2",700,500); CR2->Clear();
  TPad *CR2pad1 = new TPad("CR2pad1", "CR2pad1", 0, 0.3, 1, 1.0);//xlow, ylow, xup, yup
  CR2pad1->SetBottomMargin(0); CR2pad1->Draw();
  TPad *CR2pad2 = new TPad("CR2pad2", "CR2pad2", 0, 0.0, 1, 0.29);
  CR2pad2->SetTopMargin(0); CR2pad2->SetBottomMargin(0.3); CR2pad2->SetGridy(); CR2pad2->Draw();
  //MC vs DATA
  CR2pad1->cd();
  //Plot stacked histogram from MC
  MC_hs_CR_m2->Draw("HIST"); MC_hs_CR_m2->SetMaximum(50); MC_hs_CR_m2->GetYaxis()->SetTitle("Events/4GeV");
  //Plot MC error
  MC_CR_m2->SetLineColor(2); MC_CR_m2->SetFillColor(2); MC_CR_m2->SetFillStyle(3004); MC_CR_m2->Draw("E2 SAME");
  Double_t MC_CR_m2_error;
  Double_t MC_CR_m2_integral = MC_CR_m2->IntegralAndError(1, m_bins, MC_CR_m2_error, "");
  std::cout << "MC CR m2 integral = " << MC_CR_m2_integral << " +/- " << MC_CR_m2_error << std::endl;
  //Overlay data
  DATA_CR_m2->SetFillColor(1); DATA_CR_m2->SetLineColor(1); DATA_CR_m2->SetMarkerStyle(20); DATA_CR_m2->Draw("E1 X0 SAME");
  Double_t DATA_CR_m2_error;
  Double_t DATA_CR_m2_integral = DATA_CR_m2->IntegralAndError(1, m_bins, DATA_CR_m2_error, "");
  std::cout << "DATA CR m2 integral = " << DATA_CR_m2_integral << " +/- " << DATA_CR_m2_error << std::endl;
  //Build Legend
  TLegend* CR2pad1L = CR2pad1->BuildLegend();
  CR2pad1L->SetBorderSize(0); CR2pad1L->SetFillStyle(0); CR2pad1L->SetNColumns(2);
  TList *CR2pad1P = CR2pad1L->GetListOfPrimitives();
  TIter CR2pad1next(CR2pad1P);
  TObject *CR2pad1obj;
  TLegendEntry *CR2pad1li;
  int CR2pad1iEntry = 0;
  while ((CR2pad1obj = CR2pad1next())) {
    CR2pad1li = (TLegendEntry*)CR2pad1obj;
    CR2pad1iEntry++;
    if (CR2pad1iEntry==1) CR2pad1li->SetLabel("DYToLL (0J)");
    if (CR2pad1iEntry==2) CR2pad1li->SetLabel("DYToLL (1J)");
    if (CR2pad1iEntry==3) CR2pad1li->SetLabel("DYToLL (2J)");
    if (CR2pad1iEntry==4) CR2pad1li->SetLabel("qqToZZTo4L");
    if (CR2pad1iEntry==5) CR2pad1li->SetLabel("TTJetsToLL");
    if (CR2pad1iEntry==6) CR2pad1li->SetLabel("ggHToZZTo4L");
    if (CR2pad1iEntry==7) CR2pad1li->SetLabel("ggToZZTo4mu");
    if (CR2pad1iEntry==8) {CR2pad1li->SetLabel("MC Error"); CR2pad1li->SetOption("f");}
    if (CR2pad1iEntry==9) {CR2pad1li->SetLabel("Data"); CR2pad1li->SetOption("ep");}
  }
  CR2pad1->Update(); CR2pad1L->SetX1NDC(0.15); CR2pad1L->SetX2NDC(0.5); CR2pad1L->SetY1NDC(0.65); CR2pad1L->SetY2NDC(0.9); CR2pad1->Modified();
  gPad->RedrawAxis();
  CR2->cd(); CR2->Update();
  //Plot pull distribution
  CR2pad2->cd();
  //fill pull histogram
  TH1F *pull_CR_m2 = new TH1F("pull_CR_m2","", m_bins, m_min, m_max);
  for(unsigned int iB=1; iB<=m_bins; iB++){
    float pull_CR_m2_iB = ( DATA_CR_m2->GetBinContent(iB) - MC_CR_m2->GetBinContent(iB) ) / DATA_CR_m2->GetBinError(iB);//pull definition
    pull_CR_m2->SetBinContent(iB, pull_CR_m2_iB );
  }
  pull_CR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  pull_CR_m2->GetXaxis()->SetTitleSize(15);
  pull_CR_m2->GetXaxis()->SetTitleFont(43);
  pull_CR_m2->GetXaxis()->SetTitleOffset(3.0);
  pull_CR_m2->GetXaxis()->SetLabelSize(15);
  pull_CR_m2->GetXaxis()->SetLabelFont(43);//text size in unit of pixel, not the size of the pad
  pull_CR_m2->GetYaxis()->SetTitle("Pull");
  pull_CR_m2->GetYaxis()->CenterTitle();
  pull_CR_m2->GetYaxis()->SetTitleSize(15);
  pull_CR_m2->GetYaxis()->SetTitleFont(43);
  pull_CR_m2->GetYaxis()->SetTitleOffset(.9);
  pull_CR_m2->GetYaxis()->SetLabelSize(15);
  pull_CR_m2->GetYaxis()->SetLabelFont(43);
  pull_CR_m2->SetMinimum(-3.5);
  pull_CR_m2->SetMaximum(3.5);
  pull_CR_m2->SetStats(0);
  pull_CR_m2->SetMarkerStyle(20);
  pull_CR_m2->Draw("P");
  CR2->Write();

  //***************
  //* For m1 at SR*
  //***************
  //Consider to do a simple poly 0/1-fit on the shape
  //Data blinded until approval
  TCanvas *SR1=new TCanvas("SR1","SR m1",700,500);
  SR1->cd();
  MC_hs_SR_m1->Draw("HIST");
  MC_hs_SR_m1->SetMaximum(10);
  MC_hs_SR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  MC_hs_SR_m1->GetYaxis()->SetTitle("Events/4GeV");
  MC_SR_m1->SetLineColor(2);
  MC_SR_m1->SetFillColor(2);
  MC_SR_m1->SetFillStyle(3004);
  MC_SR_m1->Draw("E2 SAME");
  Double_t MC_SR_m1_error;
  Double_t MC_SR_m1_integral = MC_SR_m1->IntegralAndError(1, m_bins, MC_SR_m1_error, "");
  std::cout << "MC SR m1 integral = " << MC_SR_m1_integral << " +/- " << MC_SR_m1_error << std::endl;
  //Build Legend
  TLegend* SR1L = SR1->BuildLegend();
  SR1L->SetBorderSize(0); SR1L->SetFillStyle(0); SR1L->SetNColumns(2);
  TList *SR1P = SR1L->GetListOfPrimitives();
  TIter SR1next(SR1P);
  TObject *SR1obj;
  TLegendEntry *SR1li;
  int SR1iEntry = 0;
  while ((SR1obj = SR1next())) {
    SR1li = (TLegendEntry*)SR1obj;
    SR1iEntry++;
    if (SR1iEntry==1) SR1li->SetLabel("DYToLL (0J)");
    if (SR1iEntry==2) SR1li->SetLabel("DYToLL (1J)");
    if (SR1iEntry==3) SR1li->SetLabel("DYToLL (2J)");
    if (SR1iEntry==4) SR1li->SetLabel("qqToZZTo4L");
    if (SR1iEntry==5) SR1li->SetLabel("TTJetsToLL");
    if (SR1iEntry==6) SR1li->SetLabel("ggHToZZTo4L");
    if (SR1iEntry==7) SR1li->SetLabel("ggToZZTo4mu");
    if (SR1iEntry==8) {SR1li->SetLabel("MC Error"); SR1li->SetOption("f");}
  }
  SR1->Update(); SR1L->SetX1NDC(0.15); SR1L->SetX2NDC(0.5); SR1L->SetY1NDC(0.5); SR1L->SetY2NDC(0.9); SR1->Modified();
  gPad->RedrawAxis();
  SR1->Write();

  //***************
  //* For m2 at SR*
  //***************
  TCanvas *SR2=new TCanvas("SR2","SR m2",700,500);
  SR2->cd();
  MC_hs_SR_m2->Draw("HIST");
  MC_hs_SR_m2->SetMaximum(10);
  MC_hs_SR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  MC_hs_SR_m2->GetYaxis()->SetTitle("Events/4GeV");
  MC_SR_m2->SetLineColor(2);
  MC_SR_m2->SetFillColor(2);
  MC_SR_m2->SetFillStyle(3004);
  MC_SR_m2->Draw("E2 SAME");
  Double_t MC_SR_m2_error;
  Double_t MC_SR_m2_integral = MC_SR_m2->IntegralAndError(1, m_bins, MC_SR_m2_error, "");
  std::cout << "MC SR m2 integral = " << MC_SR_m2_integral << " +/- " << MC_SR_m2_error << std::endl;
  //Build Legend
  TLegend* SR2L = SR2->BuildLegend();
  SR2L->SetBorderSize(0); SR2L->SetFillStyle(0); SR2L->SetNColumns(2);
  TList *SR2P = SR2L->GetListOfPrimitives();
  TIter SR2next(SR2P);
  TObject *SR2obj;
  TLegendEntry *SR2li;
  int SR2iEntry = 0;
  while ((SR2obj = SR2next())) {
    SR2li = (TLegendEntry*)SR2obj;
    SR2iEntry++;
    if (SR2iEntry==1) SR2li->SetLabel("DYToLL (0J)");
    if (SR2iEntry==2) SR2li->SetLabel("DYToLL (1J)");
    if (SR2iEntry==3) SR2li->SetLabel("DYToLL (2J)");
    if (SR2iEntry==4) SR2li->SetLabel("qqToZZTo4L");
    if (SR2iEntry==5) SR2li->SetLabel("TTJetsToLL");
    if (SR2iEntry==6) SR2li->SetLabel("ggHToZZTo4L");
    if (SR2iEntry==7) SR2li->SetLabel("ggToZZTo4mu");
    if (SR2iEntry==8) {SR2li->SetLabel("MC Error"); SR2li->SetOption("f");}
  }
  SR2->Update(); SR2L->SetX1NDC(0.15); SR2L->SetX2NDC(0.5); SR2L->SetY1NDC(0.5); SR2L->SetY2NDC(0.9); SR2->Modified();
  gPad->RedrawAxis();
  SR2->Write();

  myPlot.Close();

} // End function: void
