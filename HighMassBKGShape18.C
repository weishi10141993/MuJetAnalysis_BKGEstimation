//=========================================================================
//= cmsenv                                                                =
//= Run it as: root -l -b -q HighMassBKGShape18.C                         =
//=          Wei Shi @Nov 20, 2019, Rice U.                               =
//=========================================================================
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

#include "Constants.h"//local const file
#include "Config.h"

void HighMassBKGShape18()
{
  // Configure inputs for year
  BKG_cfg::ConfigureInput(year);

  TLegend *txtHeader1 = new TLegend(0.09, 0.905, 0.89, 0.94);
  txtHeader1->SetFillColor(kWhite);
  txtHeader1->SetFillStyle(0);
  txtHeader1->SetBorderSize(0);
  txtHeader1->SetTextFont(42);
  txtHeader1->SetTextSize(0.045);
  txtHeader1->SetTextAlign(22);
  txtHeader1->SetHeader(header1);

  TLegend *txtHeader = new TLegend(0.09, 0.905, 0.89, 0.94);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader(header);

  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);
  THStack *MC_hs_CR_m1 = new THStack("MC_hs_CR_m1", "");
  THStack *MC_hs_CR_m2 = new THStack("MC_hs_CR_m2", "");
  THStack *MC_hs_CR_Lxy1 = new THStack("MC_hs_CR_Lxy1", "");
  THStack *MC_hs_CR_Lxy2 = new THStack("MC_hs_CR_Lxy2", "");
  THStack *MC_hs_SR_m1 = new THStack("MC_hs_SR_m1", "");
  THStack *MC_hs_SR_m2 = new THStack("MC_hs_SR_m2", "");

  // For plotting summed error for above stacked plots
  TH1F *MC_CR_m1 = new TH1F("MC_CR_m1", "", HM_m_bins, HM_m_min, HM_m_max);
  TH1F *MC_CR_m2 = new TH1F("MC_CR_m2", "", HM_m_bins, HM_m_min, HM_m_max);
  TH1F *MC_CR_Lxy1 = new TH1F("MC_CR_Lxy1", "", HM_Lxy_bins, HM_Lxy_min, HM_Lxy_max);
  TH1F *MC_CR_Lxy2 = new TH1F("MC_CR_Lxy2", "", HM_Lxy_bins, HM_Lxy_min, HM_Lxy_max);

  // used for plot error
  TH1F *MC_SR_m1 = new TH1F("MC_SR_m1", "", HM_m_bins, HM_m_min, HM_m_max);
  TH1F *MC_SR_m2 = new TH1F("MC_SR_m2", "", HM_m_bins, HM_m_min, HM_m_max);

  // SM BKG MC
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

    TH1F *MCBKGShapeCRLxyC = (TH1F*)file_tmp->Get("BKGShapeCRLxyC")->Clone("MCBKGShapeCRLxyC");
    MCBKGShapeCRLxyC->Scale(MC_ScaleFactors[i]);
    MCBKGShapeCRLxyC->SetLineColor(MC_Colors[i]);
    MCBKGShapeCRLxyC->SetFillColor(MC_Colors[i]);
    MCBKGShapeCRLxyC->SetMarkerColor(MC_Colors[i]);
    MC_hs_CR_Lxy1->Add(MCBKGShapeCRLxyC);
    MC_CR_Lxy1->Add(MCBKGShapeCRLxyC);

    TH1F *MCBKGShapeCRLxyF = (TH1F*)file_tmp->Get("BKGShapeCRLxyF")->Clone("MCBKGShapeCRLxyF");
    MCBKGShapeCRLxyF->Scale(MC_ScaleFactors[i]);
    MCBKGShapeCRLxyF->SetLineColor(MC_Colors[i]);
    MCBKGShapeCRLxyF->SetFillColor(MC_Colors[i]);
    MCBKGShapeCRLxyF->SetMarkerColor(MC_Colors[i]);
    MC_hs_CR_Lxy2->Add(MCBKGShapeCRLxyF);
    MC_CR_Lxy2->Add(MCBKGShapeCRLxyF);

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

  // DATA
  TString DATA_file_name;
  TH1F *DATA_CR_m1 = new TH1F("DATA_CR_m1", "", HM_m_bins, HM_m_min, HM_m_max);
  TH1F *DATA_CR_m2 = new TH1F("DATA_CR_m2", "", HM_m_bins, HM_m_min, HM_m_max);
  TH1F *DATA_CR_Lxy1 = new TH1F("DATA_CR_Lxy1", "", HM_Lxy_bins, HM_Lxy_min, HM_Lxy_max);
  TH1F *DATA_CR_Lxy2 = new TH1F("DATA_CR_Lxy2", "", HM_Lxy_bins, HM_Lxy_min, HM_Lxy_max);
  TH1F *DATA_SR_m1 = new TH1F("DATA_SR_m1", "", HM_m_bins, HM_m_min, HM_m_max);
  TH1F *DATA_SR_m2 = new TH1F("DATA_SR_m2", "", HM_m_bins, HM_m_min, HM_m_max);

  for (int j = 0; j < 1; j++) {
    DATA_file_name.Form( "%s/%s", store.Data(), DATA_files[j].Data() );
    std::cout << "Opening file #"<< j+1 << ": " << DATA_file_name.Data() << std::endl;
    file_tmp = TFile::Open(DATA_file_name);
    if (!file_tmp) {
      std::cout << "ERROR: could not open file " << DATA_file_name.Data() << std::endl;
      return;
    }

    TH1F *DATABKGShapeCRmassC = (TH1F*)file_tmp->Get("BKGShapeCRmassC")->Clone("DATABKGShapeCRmassC"); DATA_CR_m1->Add(DATABKGShapeCRmassC);
    TH1F *DATABKGShapeCRmassF = (TH1F*)file_tmp->Get("BKGShapeCRmassF")->Clone("DATABKGShapeCRmassF"); DATA_CR_m2->Add(DATABKGShapeCRmassF);
    TH1F *DATABKGShapeCRLxyC = (TH1F*)file_tmp->Get("BKGShapeCRLxyC")->Clone("DATABKGShapeCRLxyC"); DATA_CR_Lxy1->Add(DATABKGShapeCRLxyC);
    TH1F *DATABKGShapeCRLxyF = (TH1F*)file_tmp->Get("BKGShapeCRLxyF")->Clone("DATABKGShapeCRLxyF"); DATA_CR_Lxy2->Add(DATABKGShapeCRLxyF);
    TH1F *DATABKGShapeSRmassC = (TH1F*)file_tmp->Get("BKGShapeSRmassC")->Clone("DATABKGShapeSRmassC"); DATA_SR_m1->Add(DATABKGShapeSRmassC);
    TH1F *DATABKGShapeSRmassF = (TH1F*)file_tmp->Get("BKGShapeSRmassF")->Clone("DATABKGShapeSRmassF"); DATA_SR_m2->Add(DATABKGShapeSRmassF);
  }

  // write to output file
  TFile myPlot(outFileHMShape, "RECREATE");

  //===============
  //=   m1 at CR
  //===============
  TCanvas *CR1=new TCanvas("CR1", "CR m1", 700, 500); CR1->Clear();
  TPad *CR1pad1 = new TPad("CR1pad1", "CR1pad1", 0, 0.3, 1, 1.0);//xlow, ylow, xup, yup
  CR1pad1->SetBottomMargin(0); CR1pad1->Draw();
  TPad *CR1pad2 = new TPad("CR1pad2", "CR1pad2", 0, 0.0, 1, 0.29);
  CR1pad2->SetTopMargin(0); CR1pad2->SetBottomMargin(0.3); CR1pad2->SetGridy(); CR1pad2->Draw();
  // MC vs DATA
  CR1pad1->cd();
  // Plot stacked histogram from MC
  MC_hs_CR_m1->Draw("HIST"); MC_hs_CR_m1->SetMaximum(60); MC_hs_CR_m1->GetYaxis()->SetTitle("Events/3.5GeV");
  // Plot MC error
  MC_CR_m1->SetLineColor(2); MC_CR_m1->SetFillColor(2); MC_CR_m1->SetFillStyle(3004); MC_CR_m1->Draw("E2 SAME");
  Double_t MC_CR_m1_error;
  Double_t MC_CR_m1_integral = MC_CR_m1->IntegralAndError(1, HM_m_bins, MC_CR_m1_error, "");
  std::cout << "MC CR m1 integral = " << MC_CR_m1_integral << " +/- " << MC_CR_m1_error << std::endl;
  // Overlay data
  DATA_CR_m1->SetFillColor(1); DATA_CR_m1->SetLineColor(1); DATA_CR_m1->SetMarkerStyle(20); DATA_CR_m1->SetMarkerSize(0.6); DATA_CR_m1->Draw("E1 X0 SAME"); txtHeader1->Draw("SAME");//Draw Error bars
  Double_t DATA_CR_m1_error;
  Double_t DATA_CR_m1_integral = DATA_CR_m1->IntegralAndError(1, HM_m_bins, DATA_CR_m1_error, ""); // "": width
  std::cout << "DATA CR m1 integral = " << DATA_CR_m1_integral << " +/- " << DATA_CR_m1_error << std::endl;
  // Build Legend
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
  // Plot pull distribution
  CR1pad2->cd();
  // fill pull histogram
  TH1F *pull_CR_m1 = new TH1F("pull_CR_m1","", HM_m_bins, HM_m_min, HM_m_max);
  for(unsigned int iB=1; iB<=HM_m_bins; iB++){
    // pull definition: considering data and MC error
    float pull_CR_m1_iB = ( DATA_CR_m1->GetBinContent(iB) - MC_CR_m1->GetBinContent(iB) ) / sqrt( pow(DATA_CR_m1->GetBinError(iB), 2) + pow(MC_CR_m1->GetBinError(iB), 2) );
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
  pull_CR_m1->SetMarkerSize(0.6);
  pull_CR_m1->Draw("P");
  CR1->Write();
  CR1->SaveAs("HighMassShape/CRm1.pdf");

  //===============
  //=   m2 at CR  =
  //===============
  TCanvas *CR2=new TCanvas("CR2", "CR m2", 700, 500); CR2->Clear();
  TPad *CR2pad1 = new TPad("CR2pad1", "CR2pad1", 0, 0.3, 1, 1.0);//xlow, ylow, xup, yup
  CR2pad1->SetBottomMargin(0); CR2pad1->Draw();
  TPad *CR2pad2 = new TPad("CR2pad2", "CR2pad2", 0, 0.0, 1, 0.29);
  CR2pad2->SetTopMargin(0); CR2pad2->SetBottomMargin(0.3); CR2pad2->SetGridy(); CR2pad2->Draw();
  // MC vs DATA
  CR2pad1->cd();
  // Plot stacked histogram from MC
  MC_hs_CR_m2->Draw("HIST"); MC_hs_CR_m2->SetMaximum(60); MC_hs_CR_m2->GetYaxis()->SetTitle("Events/3.5GeV");
  // Plot MC error
  MC_CR_m2->SetLineColor(2); MC_CR_m2->SetFillColor(2); MC_CR_m2->SetFillStyle(3004); MC_CR_m2->Draw("E2 SAME");
  Double_t MC_CR_m2_error;
  Double_t MC_CR_m2_integral = MC_CR_m2->IntegralAndError(1, HM_m_bins, MC_CR_m2_error, "");
  std::cout << "MC CR m2 integral = " << MC_CR_m2_integral << " +/- " << MC_CR_m2_error << std::endl;
  // Overlay data
  DATA_CR_m2->SetFillColor(1); DATA_CR_m2->SetLineColor(1); DATA_CR_m2->SetMarkerStyle(20); DATA_CR_m2->SetMarkerSize(0.6); DATA_CR_m2->Draw("E1 X0 SAME"); txtHeader1->Draw("SAME");
  Double_t DATA_CR_m2_error;
  Double_t DATA_CR_m2_integral = DATA_CR_m2->IntegralAndError(1, HM_m_bins, DATA_CR_m2_error, "");
  std::cout << "DATA CR m2 integral = " << DATA_CR_m2_integral << " +/- " << DATA_CR_m2_error << std::endl;
  // Build Legend
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
  CR2pad1->Update(); CR2pad1L->SetX1NDC(0.15); CR2pad1L->SetX2NDC(0.5); CR2pad1L->SetY1NDC(0.5); CR2pad1L->SetY2NDC(0.9); CR2pad1->Modified();
  gPad->RedrawAxis();
  CR2->cd(); CR2->Update();
  // Plot pull distribution
  CR2pad2->cd();
  // fill pull histogram
  TH1F *pull_CR_m2 = new TH1F("pull_CR_m2","", HM_m_bins, HM_m_min, HM_m_max);
  for(unsigned int iB=1; iB<=HM_m_bins; iB++){
    float pull_CR_m2_iB = ( DATA_CR_m2->GetBinContent(iB) - MC_CR_m2->GetBinContent(iB) ) / sqrt( pow(DATA_CR_m2->GetBinError(iB), 2) + pow(MC_CR_m2->GetBinError(iB), 2) );
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
  pull_CR_m2->SetMarkerSize(0.6);
  pull_CR_m2->Draw("P");
  CR2->Write();
  CR2->SaveAs("HighMassShape/CRm2.pdf");

  //++++++++++++++++++++++++++++++++++++++++++++
  //+ Cross check data-MC agreement in Lxy at CR
  //++++++++++++++++++++++++++++++++++++++++++++

  //=========================
  //=  Lxy for mumu1 at CR  =
  //=========================
  TCanvas *CRLxy1=new TCanvas("CRLxy1", "CR Lxy1", 700, 500); CRLxy1->Clear();
  TPad *CRLxy1pad1 = new TPad("CRLxy1pad1", "CRLxy1pad1", 0, 0.3, 1, 1.0);//xlow, ylow, xup, yup
  CRLxy1pad1->SetBottomMargin(0); CRLxy1pad1->Draw();
  TPad *CRLxy1pad2 = new TPad("CRLxy1pad2", "CRLxy1pad2", 0, 0.0, 1, 0.29);
  CRLxy1pad2->SetTopMargin(0); CRLxy1pad2->SetBottomMargin(0.3); CRLxy1pad2->SetGridy(); CRLxy1pad2->Draw();
  // MC vs DATA
  CRLxy1pad1->cd();
  // Plot stacked histogram from MC
  MC_hs_CR_Lxy1->Draw("HIST"); MC_hs_CR_Lxy1->SetMaximum(150); MC_hs_CR_Lxy1->GetYaxis()->SetTitle("Events/cm");
  // Plot MC error
  MC_CR_Lxy1->SetLineColor(2); MC_CR_Lxy1->SetFillColor(2); MC_CR_Lxy1->SetFillStyle(3004); MC_CR_Lxy1->Draw("E2 SAME");
  Double_t MC_CR_Lxy1_error;
  Double_t MC_CR_Lxy1_integral = MC_CR_Lxy1->IntegralAndError(1, HM_Lxy_bins, MC_CR_Lxy1_error, "");
  std::cout << "MC CR Lxy1 integral = " << MC_CR_Lxy1_integral << " +/- " << MC_CR_Lxy1_error << std::endl;
  // Overlay data
  DATA_CR_Lxy1->SetFillColor(1); DATA_CR_Lxy1->SetLineColor(1); DATA_CR_Lxy1->SetMarkerStyle(20); DATA_CR_Lxy1->SetMarkerSize(0.6); DATA_CR_Lxy1->Draw("E1 X0 SAME"); txtHeader1->Draw("SAME");//Draw Error bars
  Double_t DATA_CR_Lxy1_error;
  Double_t DATA_CR_Lxy1_integral = DATA_CR_Lxy1->IntegralAndError(1, HM_Lxy_bins, DATA_CR_Lxy1_error, "");
  std::cout << "DATA CR Lxy1 integral = " << DATA_CR_Lxy1_integral << " +/- " << DATA_CR_Lxy1_error << std::endl;
  // Build Legend
  TLegend* CRLxy1pad1L = CRLxy1pad1->BuildLegend();
  CRLxy1pad1L->SetBorderSize(0); CRLxy1pad1L->SetFillStyle(0); CRLxy1pad1L->SetNColumns(2);
  TList *CRLxy1pad1P = CRLxy1pad1L->GetListOfPrimitives();
  TIter CRLxy1pad1next(CRLxy1pad1P);
  TObject *CRLxy1pad1obj;
  TLegendEntry *CRLxy1pad1li;
  int CRLxy1pad1iEntry = 0;
  while ((CRLxy1pad1obj = CRLxy1pad1next())) {
    CRLxy1pad1li = (TLegendEntry*)CRLxy1pad1obj;
    CRLxy1pad1iEntry++;
    if (CRLxy1pad1iEntry==1) CRLxy1pad1li->SetLabel("DYToLL (0J)");
    if (CRLxy1pad1iEntry==2) CRLxy1pad1li->SetLabel("DYToLL (1J)");
    if (CRLxy1pad1iEntry==3) CRLxy1pad1li->SetLabel("DYToLL (2J)");
    if (CRLxy1pad1iEntry==4) CRLxy1pad1li->SetLabel("qqToZZTo4L");
    if (CRLxy1pad1iEntry==5) CRLxy1pad1li->SetLabel("TTJetsToLL");
    if (CRLxy1pad1iEntry==6) CRLxy1pad1li->SetLabel("ggHToZZTo4L");
    if (CRLxy1pad1iEntry==7) CRLxy1pad1li->SetLabel("ggToZZTo4mu");
    if (CRLxy1pad1iEntry==8) {CRLxy1pad1li->SetLabel("MC Error"); CRLxy1pad1li->SetOption("f");}
    if (CRLxy1pad1iEntry==9) {CRLxy1pad1li->SetLabel("Data"); CRLxy1pad1li->SetOption("ep");}
  }
  CRLxy1pad1->Update(); CRLxy1pad1L->SetX1NDC(0.5); CRLxy1pad1L->SetX2NDC(0.9); CRLxy1pad1L->SetY1NDC(0.5); CRLxy1pad1L->SetY2NDC(0.9); CRLxy1pad1->Modified();
  gPad->RedrawAxis();
  CRLxy1->cd(); CRLxy1->Update();
  // Plot pull distribution
  CRLxy1pad2->cd();
  // fill pull histogram
  TH1F *pull_CR_Lxy1 = new TH1F("pull_CR_Lxy1", "", HM_Lxy_bins, HM_Lxy_min, HM_Lxy_max);
  for(unsigned int iB=1; iB<=HM_Lxy_bins; iB++){
    // pull definition: considering data and MC error
    // check if error is zero
    if ( DATA_CR_Lxy1->GetBinContent(iB) != 0 && MC_CR_Lxy1->GetBinContent(iB) != 0 ){
      float pull_CR_Lxy1_iB = ( DATA_CR_Lxy1->GetBinContent(iB) - MC_CR_Lxy1->GetBinContent(iB) ) / sqrt( pow(DATA_CR_Lxy1->GetBinError(iB), 2) + pow(MC_CR_Lxy1->GetBinError(iB), 2) );
      pull_CR_Lxy1->SetBinContent(iB, pull_CR_Lxy1_iB );//iB starts from #1
    }
  }
  pull_CR_Lxy1->GetXaxis()->SetTitle("L_{xy} of #mu#mu1 [cm]");
  pull_CR_Lxy1->GetXaxis()->SetTitleSize(15);
  pull_CR_Lxy1->GetXaxis()->SetTitleFont(43);
  pull_CR_Lxy1->GetXaxis()->SetTitleOffset(3.0);
  pull_CR_Lxy1->GetXaxis()->SetLabelSize(15);// labels will be 15 pixels
  pull_CR_Lxy1->GetXaxis()->SetLabelFont(43);// Absolute font size in pixel (precision 3)
  pull_CR_Lxy1->GetYaxis()->SetTitle("Pull");
  pull_CR_Lxy1->GetYaxis()->CenterTitle();
  pull_CR_Lxy1->GetYaxis()->SetTitleSize(15);
  pull_CR_Lxy1->GetYaxis()->SetTitleFont(43);
  pull_CR_Lxy1->GetYaxis()->SetTitleOffset(.9);
  pull_CR_Lxy1->GetYaxis()->SetLabelSize(15);
  pull_CR_Lxy1->GetYaxis()->SetLabelFont(43);
  pull_CR_Lxy1->SetMinimum(-3.5);
  pull_CR_Lxy1->SetMaximum(3.5);
  pull_CR_Lxy1->SetStats(0);
  pull_CR_Lxy1->SetMarkerStyle(20);
  pull_CR_Lxy1->SetMarkerSize(0.6);
  pull_CR_Lxy1->Draw("P");
  CRLxy1->Write();
  CRLxy1->SaveAs("HighMassShape/CRLxy1.pdf");

  //=========================
  //=  Lxy for mumu2 at CR  =
  //=========================
  TCanvas *CRLxy2=new TCanvas("CRLxy2", "CR Lxy2", 700, 500); CRLxy2->Clear();
  TPad *CRLxy2pad1 = new TPad("CRLxy2pad1", "CRLxy2pad1", 0, 0.3, 1, 1.0);//xlow, ylow, xup, yup
  CRLxy2pad1->SetBottomMargin(0); CRLxy2pad1->Draw();
  TPad *CRLxy2pad2 = new TPad("CRLxy2pad2", "CRLxy2pad2", 0, 0.0, 1, 0.29);
  CRLxy2pad2->SetTopMargin(0); CRLxy2pad2->SetBottomMargin(0.3); CRLxy2pad2->SetGridy(); CRLxy2pad2->Draw();
  // MC vs DATA
  CRLxy2pad1->cd();
  // Plot stacked histogram from MC
  MC_hs_CR_Lxy2->Draw("HIST"); MC_hs_CR_Lxy2->SetMaximum(150); MC_hs_CR_Lxy2->GetYaxis()->SetTitle("Events/cm");
  // Plot MC error
  MC_CR_Lxy2->SetLineColor(2); MC_CR_Lxy2->SetFillColor(2); MC_CR_Lxy2->SetFillStyle(3004); MC_CR_Lxy2->Draw("E2 SAME");
  Double_t MC_CR_Lxy2_error;
  Double_t MC_CR_Lxy2_integral = MC_CR_Lxy2->IntegralAndError(1, HM_Lxy_bins, MC_CR_Lxy2_error, "");
  std::cout << "MC CR Lxy2 integral = " << MC_CR_Lxy2_integral << " +/- " << MC_CR_Lxy2_error << std::endl;
  // Overlay data
  DATA_CR_Lxy2->SetFillColor(1); DATA_CR_Lxy2->SetLineColor(1); DATA_CR_Lxy2->SetMarkerStyle(20); DATA_CR_Lxy2->SetMarkerSize(0.6); DATA_CR_Lxy2->Draw("E1 X0 SAME"); txtHeader1->Draw("SAME");//Draw Error bars
  Double_t DATA_CR_Lxy2_error;
  Double_t DATA_CR_Lxy2_integral = DATA_CR_Lxy2->IntegralAndError(1, HM_Lxy_bins, DATA_CR_Lxy2_error, ""); // "": width
  std::cout << "DATA CR Lxy2 integral = " << DATA_CR_Lxy2_integral << " +/- " << DATA_CR_Lxy2_error << std::endl;
  // Build Legend
  TLegend* CRLxy2pad1L = CRLxy2pad1->BuildLegend();
  CRLxy2pad1L->SetBorderSize(0); CRLxy2pad1L->SetFillStyle(0); CRLxy2pad1L->SetNColumns(2);
  TList *CRLxy2pad1P = CRLxy2pad1L->GetListOfPrimitives();
  TIter CRLxy2pad1next(CRLxy2pad1P);
  TObject *CRLxy2pad1obj;
  TLegendEntry *CRLxy2pad1li;
  int CRLxy2pad1iEntry = 0;
  while ((CRLxy2pad1obj = CRLxy2pad1next())) {
    CRLxy2pad1li = (TLegendEntry*)CRLxy2pad1obj;
    CRLxy2pad1iEntry++;
    if (CRLxy2pad1iEntry==1) CRLxy2pad1li->SetLabel("DYToLL (0J)");
    if (CRLxy2pad1iEntry==2) CRLxy2pad1li->SetLabel("DYToLL (1J)");
    if (CRLxy2pad1iEntry==3) CRLxy2pad1li->SetLabel("DYToLL (2J)");
    if (CRLxy2pad1iEntry==4) CRLxy2pad1li->SetLabel("qqToZZTo4L");
    if (CRLxy2pad1iEntry==5) CRLxy2pad1li->SetLabel("TTJetsToLL");
    if (CRLxy2pad1iEntry==6) CRLxy2pad1li->SetLabel("ggHToZZTo4L");
    if (CRLxy2pad1iEntry==7) CRLxy2pad1li->SetLabel("ggToZZTo4mu");
    if (CRLxy2pad1iEntry==8) {CRLxy2pad1li->SetLabel("MC Error"); CRLxy2pad1li->SetOption("f");}
    if (CRLxy2pad1iEntry==9) {CRLxy2pad1li->SetLabel("Data"); CRLxy2pad1li->SetOption("ep");}
  }
  CRLxy2pad1->Update(); CRLxy2pad1L->SetX1NDC(0.5); CRLxy2pad1L->SetX2NDC(0.9); CRLxy2pad1L->SetY1NDC(0.5); CRLxy2pad1L->SetY2NDC(0.9); CRLxy2pad1->Modified();
  gPad->RedrawAxis();
  CRLxy2->cd(); CRLxy2->Update();
  //Plot pull distribution
  CRLxy2pad2->cd();
  //fill pull histogram
  TH1F *pull_CR_Lxy2 = new TH1F("pull_CR_Lxy2", "", HM_Lxy_bins, HM_Lxy_min, HM_Lxy_max);
  for(unsigned int iB=1; iB<=HM_Lxy_bins; iB++){
    //pull definition: considering data and MC error
    if ( DATA_CR_Lxy2->GetBinContent(iB) != 0 && MC_CR_Lxy2->GetBinContent(iB) != 0 ){
      float pull_CR_Lxy2_iB = ( DATA_CR_Lxy2->GetBinContent(iB) - MC_CR_Lxy2->GetBinContent(iB) ) / sqrt( pow(DATA_CR_Lxy2->GetBinError(iB), 2) + pow(MC_CR_Lxy2->GetBinError(iB), 2) );
      pull_CR_Lxy2->SetBinContent(iB, pull_CR_Lxy2_iB );//iB starts from #1
    }
  }
  pull_CR_Lxy2->GetXaxis()->SetTitle("L_{xy} of #mu#mu2 [cm]");
  pull_CR_Lxy2->GetXaxis()->SetTitleSize(15);
  pull_CR_Lxy2->GetXaxis()->SetTitleFont(43);
  pull_CR_Lxy2->GetXaxis()->SetTitleOffset(3.0);
  pull_CR_Lxy2->GetXaxis()->SetLabelSize(15);// labels will be 15 pixels
  pull_CR_Lxy2->GetXaxis()->SetLabelFont(43);// Absolute font size in pixel (precision 3)
  pull_CR_Lxy2->GetYaxis()->SetTitle("Pull");
  pull_CR_Lxy2->GetYaxis()->CenterTitle();
  pull_CR_Lxy2->GetYaxis()->SetTitleSize(15);
  pull_CR_Lxy2->GetYaxis()->SetTitleFont(43);
  pull_CR_Lxy2->GetYaxis()->SetTitleOffset(.9);
  pull_CR_Lxy2->GetYaxis()->SetLabelSize(15);
  pull_CR_Lxy2->GetYaxis()->SetLabelFont(43);
  pull_CR_Lxy2->SetMinimum(-3.5);
  pull_CR_Lxy2->SetMaximum(3.5);
  pull_CR_Lxy2->SetStats(0);
  pull_CR_Lxy2->SetMarkerStyle(20);
  pull_CR_Lxy2->SetMarkerSize(0.6);
  pull_CR_Lxy2->Draw("P");
  CRLxy2->Write();
  CRLxy2->SaveAs("HighMassShape/CRLxy2.pdf");

  //++++++++++++++++++++++++++++++++++++++++++++
  //+ End cross check in Lxy at CR
  //++++++++++++++++++++++++++++++++++++++++++++

  //=====================================
  //=  m1 at SR: data blinded version   =
  //=====================================

  TCanvas *SR1=new TCanvas("SR1", "SR m1", 700, 500);
  SR1->cd();
  MC_hs_SR_m1->Draw("HIST");
  MC_hs_SR_m1->SetMaximum(9);
  MC_hs_SR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  MC_hs_SR_m1->GetYaxis()->SetTitle("Events/3.5GeV");
  MC_SR_m1->SetLineColor(2);
  MC_SR_m1->SetFillColor(2);
  MC_SR_m1->SetFillStyle(3004);
  MC_SR_m1->Draw("E2 SAME"); txtHeader->Draw("SAME");
  Double_t MC_SR_m1_error;
  Double_t MC_SR_m1_integral = MC_SR_m1->IntegralAndError(1, HM_m_bins, MC_SR_m1_error, "");
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
  SR1->SaveAs("HighMassShape/SRm1.pdf");

  //---------------------------------------------------------------------------------------------
  //- Print bin content and bin error for all bins to be used for background pdf shape variation
  //- Bin #0 contains the underflow. Normally starts from #1
  //---------------------------------------------------------------------------------------------
  std::cout << "=== Print bin content and bin error: MC SR m1 ===" << std::endl;
  std::cout << "bin # * " << "bin content * " << "bin error" << std::endl;
  for(unsigned int iB=1; iB<=HM_m_bins; iB++){ std::cout << left << setw(8) << iB << left << setw(14) << MC_SR_m1->GetBinContent(iB) << left << setw(9) << MC_SR_m1->GetBinError(iB) << std::endl; }
  //-----------------------------------------------------------------
  //- Draw interpolation lines between bins using nominal bin content
  //-----------------------------------------------------------------
  double plotMCm1[16];
  // lower and upper bound is fixed
  plotMCm1[0] = 11.0; plotMCm1[15] = 60.0;
  // nominal, +1/-1 sigma cases: 14+2points, 2 is for the lower and upper bounds
  double plotMCm1contentNom[16];
  double plotMCm1contentUp[16];
  double plotMCm1contentDn[16];
  double plotMCm1contentBraidA[16];
  double plotMCm1contentBraidB[16];
  for (int i = 0; i < 14; i++) {
    plotMCm1[i+1] = MCBinCenterMass[i];
    plotMCm1contentNom[i+1] = MCBinContentm1[i];
    plotMCm1contentUp[i+1] = MCBinContentm1[i] + MCBinErrm1[i];
    plotMCm1contentDn[i+1] = MCBinContentm1[i] - MCBinErrm1[i];
    if ( (i+1) % 2 == 0 ) { plotMCm1contentBraidA[i+1] = MCBinContentm1[i] - MCBinErrm1[i]; }
    else { plotMCm1contentBraidA[i+1] = MCBinContentm1[i] + MCBinErrm1[i]; }
    if ( (i+1) % 2 == 0 ) { plotMCm1contentBraidB[i+1] = MCBinContentm1[i] + MCBinErrm1[i]; }
    else { plotMCm1contentBraidB[i+1] = MCBinContentm1[i] - MCBinErrm1[i]; }
  }
  plotMCm1contentNom[0]    = BKG_cfg::My_BKGShapem1(11.0);            plotMCm1contentNom[15]    = BKG_cfg::My_BKGShapem1(60.0);
  plotMCm1contentUp[0]     = BKG_cfg::My_BKGShapem1SigmaUp(11.0);     plotMCm1contentUp[15]     = BKG_cfg::My_BKGShapem1SigmaUp(60.0);
  plotMCm1contentDn[0]     = BKG_cfg::My_BKGShapem1SigmaDn(11.0);     plotMCm1contentDn[15]     = BKG_cfg::My_BKGShapem1SigmaDn(60.0);
  plotMCm1contentBraidA[0] = BKG_cfg::My_BKGShapem1SigmaBraidA(11.0); plotMCm1contentBraidA[15] = BKG_cfg::My_BKGShapem1SigmaBraidA(60.0);
  plotMCm1contentBraidB[0] = BKG_cfg::My_BKGShapem1SigmaBraidB(11.0); plotMCm1contentBraidB[15] = BKG_cfg::My_BKGShapem1SigmaBraidB(60.0);
  //plot
  TGraph* BkgPdfm1Nom = new TGraph(16, plotMCm1, plotMCm1contentNom);
  TGraph* BkgPdfm1Up  = new TGraph(16, plotMCm1, plotMCm1contentUp);
  TGraph* BkgPdfm1Dn  = new TGraph(16, plotMCm1, plotMCm1contentDn);
  TGraph* BkgPdfm1BraidA  = new TGraph(16, plotMCm1, plotMCm1contentBraidA);
  TGraph* BkgPdfm1BraidB  = new TGraph(16, plotMCm1, plotMCm1contentBraidB);
  BkgPdfm1Nom->SetLineColor(2); BkgPdfm1Nom->SetLineStyle(1); BkgPdfm1Nom->SetLineWidth(2);
  BkgPdfm1Up->SetLineColor(2); BkgPdfm1Up->SetLineStyle(2); BkgPdfm1Up->SetLineWidth(2);
  BkgPdfm1Dn->SetLineColor(2); BkgPdfm1Dn->SetLineStyle(3); BkgPdfm1Dn->SetLineWidth(2);
  BkgPdfm1BraidA->SetLineColor(4); BkgPdfm1BraidA->SetLineStyle(4); BkgPdfm1BraidA->SetLineWidth(2);
  BkgPdfm1BraidB->SetLineColor(28); BkgPdfm1BraidB->SetLineStyle(5); BkgPdfm1BraidB->SetLineWidth(2);
  BkgPdfm1Nom->Draw("L"); BkgPdfm1Up->Draw("L"); BkgPdfm1Dn->Draw("L"); BkgPdfm1BraidA->Draw("L"); BkgPdfm1BraidB->Draw("L");
  auto BkgPdfm1Legend = new TLegend(0.6, 0.5, 0.9, 0.9);
  BkgPdfm1Legend->SetHeader("Shape Interpolation: m_{#mu#mu1}", "C");
  BkgPdfm1Legend->AddEntry(BkgPdfm1Nom,    "Nominal", "l");
  BkgPdfm1Legend->AddEntry(BkgPdfm1Up,     "Nominal + 1 #sigma", "l");
  BkgPdfm1Legend->AddEntry(BkgPdfm1Dn,     "Nominal - 1 #sigma", "l");
  BkgPdfm1Legend->AddEntry(BkgPdfm1BraidA, "Braid A", "l");
  BkgPdfm1Legend->AddEntry(BkgPdfm1BraidB, "Braid B", "l");
  BkgPdfm1Legend->Draw();
  SR1->SaveAs("HighMassShape/SRm1Variations.pdf");

  //========================================
  //= For m2 at SR: data blinded version   =
  //========================================
  TCanvas *SR2=new TCanvas("SR2", "SR m2", 700, 500);
  SR2->cd();
  MC_hs_SR_m2->Draw("HIST");
  MC_hs_SR_m2->SetMaximum(9);
  MC_hs_SR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  MC_hs_SR_m2->GetYaxis()->SetTitle("Events/3.5GeV");
  MC_SR_m2->SetLineColor(2);
  MC_SR_m2->SetFillColor(2);
  MC_SR_m2->SetFillStyle(3004);
  MC_SR_m2->Draw("E2 SAME"); txtHeader->Draw("SAME");
  Double_t MC_SR_m2_error;
  Double_t MC_SR_m2_integral = MC_SR_m2->IntegralAndError(1, HM_m_bins, MC_SR_m2_error, "");
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
  SR2->SaveAs("HighMassShape/SRm2.pdf");

  // Print bin content and bin error for all bins to be used for background pdf shape variation
  std::cout << "=== Print bin content and bin error: MC SR m2 ===" << std::endl;
  std::cout << "bin # * " << "bin content * " << "bin error" << std::endl;
  for(unsigned int iB=1; iB<=HM_m_bins; iB++){ std::cout << left << setw(8) << iB << left << setw(14) << MC_SR_m2->GetBinContent(iB) << left << setw(9) << MC_SR_m2->GetBinError(iB) << std::endl; }
  //-----------------------------------------------------------------
  //- Draw interpolation lines between bins using nominal bin content
  //-----------------------------------------------------------------
  // nominal case: 14+2points, 2 is for the 11 GeV and 60 GeV bound
  double plotMCm2[16];
  plotMCm2[0] = 11.0; plotMCm2[15] = 60.0;
  double plotMCm2contentNom[16];
  double plotMCm2contentUp[16];
  double plotMCm2contentDn[16];
  double plotMCm2contentBraidA[16];
  double plotMCm2contentBraidB[16];
  for (int i = 0; i < 14; i++) {
    plotMCm2[i+1] = MCBinCenterMass[i];
    plotMCm2contentNom[i+1] = MCBinContentm2[i];
    plotMCm2contentUp[i+1] = MCBinContentm2[i] + MCBinErrm2[i];
    plotMCm2contentDn[i+1] = MCBinContentm2[i] - MCBinErrm2[i];
    if ( (i+1) % 2 == 0 ) { plotMCm2contentBraidA[i+1] = MCBinContentm2[i] - MCBinErrm2[i]; }
    else { plotMCm2contentBraidA[i+1] = MCBinContentm2[i] + MCBinErrm2[i]; }
    if ( (i+1) % 2 == 0 ) { plotMCm2contentBraidB[i+1] = MCBinContentm2[i] + MCBinErrm2[i]; }
    else { plotMCm2contentBraidB[i+1] = MCBinContentm2[i] - MCBinErrm2[i]; }
  }
  plotMCm2contentNom[0]    = BKG_cfg::My_BKGShapem2(11.0);            plotMCm2contentNom[15]    = BKG_cfg::My_BKGShapem2(60.0);
  plotMCm2contentUp[0]     = BKG_cfg::My_BKGShapem2SigmaUp(11.0);     plotMCm2contentUp[15]     = BKG_cfg::My_BKGShapem2SigmaUp(60.0);
  plotMCm2contentDn[0]     = BKG_cfg::My_BKGShapem2SigmaDn(11.0);     plotMCm2contentDn[15]     = BKG_cfg::My_BKGShapem2SigmaDn(60.0);
  plotMCm2contentBraidA[0] = BKG_cfg::My_BKGShapem2SigmaBraidA(11.0); plotMCm2contentBraidA[15] = BKG_cfg::My_BKGShapem2SigmaBraidA(60.0);
  plotMCm2contentBraidB[0] = BKG_cfg::My_BKGShapem2SigmaBraidB(11.0); plotMCm2contentBraidB[15] = BKG_cfg::My_BKGShapem2SigmaBraidB(60.0);
  //plot
  TGraph* BkgPdfm2Nom = new TGraph(16, plotMCm2, plotMCm2contentNom);
  TGraph* BkgPdfm2Up  = new TGraph(16, plotMCm2, plotMCm2contentUp);
  TGraph* BkgPdfm2Dn  = new TGraph(16, plotMCm2, plotMCm2contentDn);
  TGraph* BkgPdfm2BraidA  = new TGraph(16, plotMCm2, plotMCm2contentBraidA);
  TGraph* BkgPdfm2BraidB  = new TGraph(16, plotMCm2, plotMCm2contentBraidB);
  BkgPdfm2Nom->SetLineColor(2); BkgPdfm2Nom->SetLineStyle(1); BkgPdfm2Nom->SetLineWidth(2);
  BkgPdfm2Up->SetLineColor(2); BkgPdfm2Up->SetLineStyle(2); BkgPdfm2Up->SetLineWidth(2);
  BkgPdfm2Dn->SetLineColor(2); BkgPdfm2Dn->SetLineStyle(3); BkgPdfm2Dn->SetLineWidth(2);
  BkgPdfm2BraidA->SetLineColor(4); BkgPdfm2BraidA->SetLineStyle(4); BkgPdfm2BraidA->SetLineWidth(2);
  BkgPdfm2BraidB->SetLineColor(28); BkgPdfm2BraidB->SetLineStyle(5); BkgPdfm2BraidB->SetLineWidth(2);
  BkgPdfm2Nom->Draw("L"); BkgPdfm2Up->Draw("L"); BkgPdfm2Dn->Draw("L"); BkgPdfm2BraidA->Draw("L"); BkgPdfm2BraidB->Draw("L");
  auto BkgPdfm2Legend = new TLegend(0.6, 0.5, 0.9, 0.9);
  BkgPdfm2Legend->SetHeader("Shape Interpolation: m_{#mu#mu2}", "C");
  BkgPdfm2Legend->AddEntry(BkgPdfm2Nom,    "Nominal", "l");
  BkgPdfm2Legend->AddEntry(BkgPdfm2Up,     "Nominal + 1 #sigma", "l");
  BkgPdfm2Legend->AddEntry(BkgPdfm2Dn,     "Nominal - 1 #sigma", "l");
  BkgPdfm2Legend->AddEntry(BkgPdfm2BraidA, "Braid A", "l");
  BkgPdfm2Legend->AddEntry(BkgPdfm2BraidB, "Braid B", "l");
  BkgPdfm2Legend->Draw();
  SR2->SaveAs("HighMassShape/SRm2Variations.pdf");

  std::cout << "*******************************" << std::endl;
  std::cout << "*          Unblinding         *" << std::endl;
  std::cout << "*******************************" << std::endl;
  //==========================================================================
  //=  m1 at SR: data unblinded version (after green light from conveners)   =
  //==========================================================================
  TCanvas *SR1_Unblind = new TCanvas("SR1_Unblind", "SR unblinded m1", 700, 500); SR1_Unblind->Clear();
  TPad *SR1_Unblindpad1 = new TPad("SR1_Unblindpad1", "SR1_Unblindpad1", 0, 0.3, 1, 1.0); // xlow, ylow, xup, yup
  SR1_Unblindpad1->SetBottomMargin(0); SR1_Unblindpad1->Draw();
  TPad *SR1_Unblindpad2 = new TPad("SR1_Unblindpad2", "SR1_Unblindpad2", 0, 0.0, 1, 0.29);
  SR1_Unblindpad2->SetTopMargin(0); SR1_Unblindpad2->SetBottomMargin(0.3); SR1_Unblindpad2->SetGridy(); SR1_Unblindpad2->Draw();
  // MC vs DATA
  SR1_Unblindpad1->cd();
  // Plot stacked histogram from MC
  MC_hs_SR_m1->Draw("HIST"); MC_hs_SR_m1->SetMaximum(9); MC_hs_SR_m1->GetYaxis()->SetTitle("Events/3.5GeV");
  // Plot MC error
  MC_SR_m1->SetLineColor(2); MC_SR_m1->SetFillColor(2); MC_SR_m1->SetFillStyle(3004); MC_SR_m1->Draw("E2 SAME");
  //Double_t MC_SR_m1_error;
  //Double_t MC_SR_m1_integral = MC_SR_m1->IntegralAndError(1, HM_m_bins, MC_SR_m1_error, "");
  std::cout << "MC SR m1 integral = " << MC_SR_m1_integral << " +/- " << MC_SR_m1_error << std::endl;
  // Overlay data
  DATA_SR_m1->SetFillColor(1); DATA_SR_m1->SetLineColor(1); DATA_SR_m1->SetMarkerStyle(20); DATA_SR_m1->SetMarkerSize(0.6); DATA_SR_m1->Draw("E1 X0 SAME"); txtHeader1->Draw("SAME"); // Draw Error bars
  Double_t DATA_SR_m1_error;
  Double_t DATA_SR_m1_integral = DATA_SR_m1->IntegralAndError(1, HM_m_bins, DATA_SR_m1_error, ""); // "": width
  std::cout << "DATA SR m1 integral = " << DATA_SR_m1_integral << " +/- " << DATA_SR_m1_error << std::endl;
  // Build Legend
  TLegend* SR1_Unblindpad1L = SR1_Unblindpad1->BuildLegend();
  SR1_Unblindpad1L->SetBorderSize(0); SR1_Unblindpad1L->SetFillStyle(0); SR1_Unblindpad1L->SetNColumns(2);
  TList *SR1_Unblindpad1P = SR1_Unblindpad1L->GetListOfPrimitives();
  TIter SR1_Unblindpad1next(SR1_Unblindpad1P);
  TObject *SR1_Unblindpad1obj;
  TLegendEntry *SR1_Unblindpad1li;
  int SR1_Unblindpad1iEntry = 0;
  while ((SR1_Unblindpad1obj = SR1_Unblindpad1next())) {
    SR1_Unblindpad1li = (TLegendEntry*)SR1_Unblindpad1obj;
    SR1_Unblindpad1iEntry++;
    if (SR1_Unblindpad1iEntry==1) SR1_Unblindpad1li->SetLabel("DYToLL (0J)");
    if (SR1_Unblindpad1iEntry==2) SR1_Unblindpad1li->SetLabel("DYToLL (1J)");
    if (SR1_Unblindpad1iEntry==3) SR1_Unblindpad1li->SetLabel("DYToLL (2J)");
    if (SR1_Unblindpad1iEntry==4) SR1_Unblindpad1li->SetLabel("qqToZZTo4L");
    if (SR1_Unblindpad1iEntry==5) SR1_Unblindpad1li->SetLabel("TTJetsToLL");
    if (SR1_Unblindpad1iEntry==6) SR1_Unblindpad1li->SetLabel("ggHToZZTo4L");
    if (SR1_Unblindpad1iEntry==7) SR1_Unblindpad1li->SetLabel("ggToZZTo4mu");
    if (SR1_Unblindpad1iEntry==8) {SR1_Unblindpad1li->SetLabel("MC Error"); SR1_Unblindpad1li->SetOption("f");}
    if (SR1_Unblindpad1iEntry==9) {SR1_Unblindpad1li->SetLabel("Data"); SR1_Unblindpad1li->SetOption("ep");}
  }
  SR1_Unblindpad1->Update(); SR1_Unblindpad1L->SetX1NDC(0.15); SR1_Unblindpad1L->SetX2NDC(0.5); SR1_Unblindpad1L->SetY1NDC(0.5); SR1_Unblindpad1L->SetY2NDC(0.9); SR1_Unblindpad1->Modified();
  gPad->RedrawAxis();
  SR1_Unblind->cd(); SR1_Unblind->Update();
  // Plot pull distribution
  SR1_Unblindpad2->cd();
  // fill pull histogram
  TH1F *pull_SR_m1 = new TH1F("pull_SR_m1","", HM_m_bins, HM_m_min, HM_m_max);
  for (unsigned int iB = 1; iB <= HM_m_bins; iB++) {
    // pull definition: considering data and MC error
    if ( DATA_SR_m1->GetBinContent(iB) != 0 && MC_SR_m1->GetBinContent(iB) != 0 ){
      float pull_SR_m1_iB = ( DATA_SR_m1->GetBinContent(iB) - MC_SR_m1->GetBinContent(iB) ) / sqrt( pow(DATA_SR_m1->GetBinError(iB), 2) + pow(MC_SR_m1->GetBinError(iB), 2) );
      pull_SR_m1->SetBinContent(iB, pull_SR_m1_iB );//iB starts from #1
    }
  }
  pull_SR_m1->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]");
  pull_SR_m1->GetXaxis()->SetTitleSize(15);
  pull_SR_m1->GetXaxis()->SetTitleFont(43);
  pull_SR_m1->GetXaxis()->SetTitleOffset(3.0);
  pull_SR_m1->GetXaxis()->SetLabelSize(15);// labels will be 15 pixels
  pull_SR_m1->GetXaxis()->SetLabelFont(43);// Absolute font size in pixel (precision 3)
  pull_SR_m1->GetYaxis()->SetTitle("Pull");
  pull_SR_m1->GetYaxis()->CenterTitle();
  pull_SR_m1->GetYaxis()->SetTitleSize(15);
  pull_SR_m1->GetYaxis()->SetTitleFont(43);
  pull_SR_m1->GetYaxis()->SetTitleOffset(.9);
  pull_SR_m1->GetYaxis()->SetLabelSize(15);
  pull_SR_m1->GetYaxis()->SetLabelFont(43);
  pull_SR_m1->SetMinimum(-3.5);
  pull_SR_m1->SetMaximum(3.5);
  pull_SR_m1->SetStats(0);
  pull_SR_m1->SetMarkerStyle(20);
  pull_SR_m1->SetMarkerSize(0.6);
  pull_SR_m1->Draw("P");
  SR1_Unblind->Write();
  SR1_Unblind->SaveAs("HighMassShape/SRm1_Unblinded.pdf");

  //==========================================================================
  //=  m2 at SR: data unblinded version (after green light from conveners)   =
  //==========================================================================
  TCanvas *SR2_Unblind = new TCanvas("SR2_Unblind", "SR unblinded m2", 700, 500); SR2_Unblind->Clear();
  TPad *SR2_Unblindpad1 = new TPad("SR2_Unblindpad1", "SR2_Unblindpad1", 0, 0.3, 1, 1.0); // xlow, ylow, xup, yup
  SR2_Unblindpad1->SetBottomMargin(0); SR2_Unblindpad1->Draw();
  TPad *SR2_Unblindpad2 = new TPad("SR2_Unblindpad2", "SR2_Unblindpad2", 0, 0.0, 1, 0.29);
  SR2_Unblindpad2->SetTopMargin(0); SR2_Unblindpad2->SetBottomMargin(0.3); SR2_Unblindpad2->SetGridy(); SR2_Unblindpad2->Draw();
  // MC vs DATA
  SR2_Unblindpad1->cd();
  // Plot stacked histogram from MC
  MC_hs_SR_m2->Draw("HIST"); MC_hs_SR_m2->SetMaximum(9); MC_hs_SR_m2->GetYaxis()->SetTitle("Events/3.5GeV");
  // Plot MC error
  MC_SR_m2->SetLineColor(2); MC_SR_m2->SetFillColor(2); MC_SR_m2->SetFillStyle(3004); MC_SR_m2->Draw("E2 SAME");
  //Double_t MC_SR_m2_error;
  //Double_t MC_SR_m2_integral = MC_SR_m2->IntegralAndError(1, HM_m_bins, MC_SR_m2_error, "");
  std::cout << "MC SR m2 integral = " << MC_SR_m2_integral << " +/- " << MC_SR_m2_error << std::endl;
  // Overlay data
  DATA_SR_m2->SetFillColor(1); DATA_SR_m2->SetLineColor(1); DATA_SR_m2->SetMarkerStyle(20); DATA_SR_m2->SetMarkerSize(0.6); DATA_SR_m2->Draw("E1 X0 SAME"); txtHeader1->Draw("SAME"); // Draw Error bars
  Double_t DATA_SR_m2_error;
  Double_t DATA_SR_m2_integral = DATA_SR_m2->IntegralAndError(1, HM_m_bins, DATA_SR_m2_error, ""); // "": width
  std::cout << "DATA SR m2 integral = " << DATA_SR_m2_integral << " +/- " << DATA_SR_m2_error << std::endl;
  // Build Legend
  TLegend* SR2_Unblindpad1L = SR2_Unblindpad1->BuildLegend();
  SR2_Unblindpad1L->SetBorderSize(0); SR2_Unblindpad1L->SetFillStyle(0); SR2_Unblindpad1L->SetNColumns(2);
  TList *SR2_Unblindpad1P = SR2_Unblindpad1L->GetListOfPrimitives();
  TIter SR2_Unblindpad1next(SR2_Unblindpad1P);
  TObject *SR2_Unblindpad1obj;
  TLegendEntry *SR2_Unblindpad1li;
  int SR2_Unblindpad1iEntry = 0;
  while ((SR2_Unblindpad1obj = SR2_Unblindpad1next())) {
    SR2_Unblindpad1li = (TLegendEntry*)SR2_Unblindpad1obj;
    SR2_Unblindpad1iEntry++;
    if (SR2_Unblindpad1iEntry==1) SR2_Unblindpad1li->SetLabel("DYToLL (0J)");
    if (SR2_Unblindpad1iEntry==2) SR2_Unblindpad1li->SetLabel("DYToLL (1J)");
    if (SR2_Unblindpad1iEntry==3) SR2_Unblindpad1li->SetLabel("DYToLL (2J)");
    if (SR2_Unblindpad1iEntry==4) SR2_Unblindpad1li->SetLabel("qqToZZTo4L");
    if (SR2_Unblindpad1iEntry==5) SR2_Unblindpad1li->SetLabel("TTJetsToLL");
    if (SR2_Unblindpad1iEntry==6) SR2_Unblindpad1li->SetLabel("ggHToZZTo4L");
    if (SR2_Unblindpad1iEntry==7) SR2_Unblindpad1li->SetLabel("ggToZZTo4mu");
    if (SR2_Unblindpad1iEntry==8) {SR2_Unblindpad1li->SetLabel("MC Error"); SR2_Unblindpad1li->SetOption("f");}
    if (SR2_Unblindpad1iEntry==9) {SR2_Unblindpad1li->SetLabel("Data"); SR2_Unblindpad1li->SetOption("ep");}
  }
  SR2_Unblindpad1->Update(); SR2_Unblindpad1L->SetX1NDC(0.15); SR2_Unblindpad1L->SetX2NDC(0.5); SR2_Unblindpad1L->SetY1NDC(0.5); SR2_Unblindpad1L->SetY2NDC(0.9); SR2_Unblindpad1->Modified();
  gPad->RedrawAxis();
  SR2_Unblind->cd(); SR2_Unblind->Update();
  // Plot pull distribution
  SR2_Unblindpad2->cd();
  // fill pull histogram
  TH1F *pull_SR_m2 = new TH1F("pull_SR_m2","", HM_m_bins, HM_m_min, HM_m_max);
  for (unsigned int iB = 1; iB <= HM_m_bins; iB++) {
    // pull definition: considering data and MC error
    if ( DATA_SR_m2->GetBinContent(iB) != 0 && MC_SR_m2->GetBinContent(iB) != 0 ){
      float pull_SR_m2_iB = ( DATA_SR_m2->GetBinContent(iB) - MC_SR_m2->GetBinContent(iB) ) / sqrt( pow(DATA_SR_m2->GetBinError(iB), 2) + pow(MC_SR_m2->GetBinError(iB), 2) );
      pull_SR_m2->SetBinContent(iB, pull_SR_m2_iB );//iB starts from #1
    }
  }
  pull_SR_m2->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]");
  pull_SR_m2->GetXaxis()->SetTitleSize(15);
  pull_SR_m2->GetXaxis()->SetTitleFont(43);
  pull_SR_m2->GetXaxis()->SetTitleOffset(3.0);
  pull_SR_m2->GetXaxis()->SetLabelSize(15);// labels will be 15 pixels
  pull_SR_m2->GetXaxis()->SetLabelFont(43);// Absolute font size in pixel (precision 3)
  pull_SR_m2->GetYaxis()->SetTitle("Pull");
  pull_SR_m2->GetYaxis()->CenterTitle();
  pull_SR_m2->GetYaxis()->SetTitleSize(15);
  pull_SR_m2->GetYaxis()->SetTitleFont(43);
  pull_SR_m2->GetYaxis()->SetTitleOffset(.9);
  pull_SR_m2->GetYaxis()->SetLabelSize(15);
  pull_SR_m2->GetYaxis()->SetLabelFont(43);
  pull_SR_m2->SetMinimum(-3.5);
  pull_SR_m2->SetMaximum(3.5);
  pull_SR_m2->SetStats(0);
  pull_SR_m2->SetMarkerStyle(20);
  pull_SR_m2->SetMarkerSize(0.6);
  pull_SR_m2->Draw("P");
  SR2_Unblind->Write();
  SR2_Unblind->SaveAs("HighMassShape/SRm2_Unblinded.pdf");

  myPlot.Close();

} // End function: void
