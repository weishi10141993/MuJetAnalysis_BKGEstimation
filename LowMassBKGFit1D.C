//*****************************************************************************************************
//* cmsenv                                                                                            *
//* To request more time: sintr -t 480                                                                *
//*                                       Wei Shi @Nov 20, 2019, Rice U.                              *
//*****************************************************************************************************
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

#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooWorkspace.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooClassFactory.h"
#include "RooHistPdf.h"
#include "RooCustomizer.h"
#include "RooMultiVarGaussian.h"
#include "RooTFnBinding.h"
#include "RooArgusBG.h"
#include "RooBernstein.h"
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "macros/tdrStyle.C"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void LowMassBKGFit1D() {
  //Constants
  TString pT_cut  = "17";
  TString eta_cut = "0.9";
  TString iso_cut = "1.5";//isolation cut for leading muon in dimu for Run2 analysis

  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  const double       b_min  = 0.;//for b jets
  const double       b_max  = 10.;
  const unsigned int b_bins = 10;

  //Style
  setTDRStyle();
  TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader("#bf{CMS} #it{Preliminary}    36.734 fb^{-1} (2017 13 TeV)");

  //Output ws
  RooWorkspace* w = new RooWorkspace("w");
  TString Comm = "mkdir -p figures/";
  system( Comm.Data() );

  //Input file
  TChain chain_data_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  TChain chain_data_dimuorphan("cutFlowAnalyzerPXBL4PXFL3/Events_orphan");
  std::ifstream Myfile( "Input_2017CDEF.txt" );
  std::string Line;
  if( !Myfile ) std::cout<<"ERROR opening Myfile."<<std::endl;
  while (std::getline(Myfile, Line)){
    TString Line2(Line);
    if( Line2.Contains("root") ){
      chain_data_dimudimu.Add(Line2.Data());
      chain_data_dimuorphan.Add(Line2.Data());
    }
  }

  //Define RooRealVar
  RooRealVar m1("m1","m_{#mu#mu_{1}}",m_min,m_max,"GeV");
  RooRealVar m1TightBJet("m1TightBJet","m_{#mu#mu_{1}} TightBJet",b_min,b_max,"");
  RooRealVar m1MediumBJet("m1MediumBJet","m_{#mu#mu_{1}} MediumBJet",b_min,b_max,"");
  RooRealVar m1LooseBJet("m1LooseBJet","m_{#mu#mu_{1}} LooseBJet",b_min,b_max,"");
  RooRealVar m2("m2","m_{#mu#mu_{2}}",m_min,m_max,"GeV");
  RooRealVar m2TightBJet("m2TightBJet","m_{#mu#mu_{2}} TightBJet",b_min,b_max,"");
  RooRealVar m2MediumBJet("m2MediumBJet","m_{#mu#mu_{2}} MediumBJet",b_min,b_max,"");
  RooRealVar m2LooseBJet("m2LooseBJet","m_{#mu#mu_{2}} LooseBJet",b_min,b_max,"");

  m1.setBins(m_bins);
  m1TightBJet.setBins(b_bins);
  m1MediumBJet.setBins(b_bins);
  m1LooseBJet.setBins(b_bins);
  m2.setBins(m_bins);
  m2TightBJet.setBins(b_bins);
  m2MediumBJet.setBins(b_bins);
  m2LooseBJet.setBins(b_bins);

  w->import(m1);
  w->import(m1TightBJet);
  w->import(m1MediumBJet);
  w->import(m1LooseBJet);
  w->import(m2);
  w->import(m2TightBJet);
  w->import(m2MediumBJet);
  w->import(m2LooseBJet);

  //**********************************************************************************
  //     Select events for constructing 1D templates, identify high pT muon, m1!=m2
  //**********************************************************************************
  //m1: High pT mu is in orphan dimu
  ostringstream stream_cut_bg_m1_iso;
  stream_cut_bg_m1_iso << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) || ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso = stream_cut_bg_m1_iso.str();
  TTree* tree_dimuorphan_bg_m1 = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso);

  //m2: Orphan mu is high pT
  ostringstream stream_cut_bg_m2_iso;
  stream_cut_bg_m2_iso << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m2_iso = stream_cut_bg_m2_iso.str();
  TTree* tree_dimuorphan_bg_m2 = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso);

  //Setting Names: Control region
  tree_dimuorphan_bg_m1->GetBranch("orph_dimu_mass")->SetName("m1");
  tree_dimuorphan_bg_m1->GetBranch("NPATJetTightB")->SetName("m1TightBJet");
  tree_dimuorphan_bg_m1->GetBranch("NPATJetMediumB")->SetName("m1MediumBJet");
  tree_dimuorphan_bg_m1->GetBranch("NPATJetLooseB")->SetName("m1LooseBJet");

  tree_dimuorphan_bg_m2->GetBranch("orph_dimu_mass")->SetName("m2");
  tree_dimuorphan_bg_m2->GetBranch("NPATJetTightB")->SetName("m2TightBJet");
  tree_dimuorphan_bg_m2->GetBranch("NPATJetMediumB")->SetName("m2MediumBJet");
  tree_dimuorphan_bg_m2->GetBranch("NPATJetLooseB")->SetName("m2LooseBJet");

  //Creating dataset using orphan dimu tree
  RooDataSet* ds_dimuorphan_bg_m1 = new RooDataSet("ds_dimuorphan_bg_m1","ds_dimuorphan_bg_m1", tree_dimuorphan_bg_m1, RooArgSet(m1));
  RooDataSet* ds_dimuorphan_bg_m1TightBJet = new RooDataSet("ds_dimuorphan_bg_m1TightBJet","ds_dimuorphan_bg_m1TightBJet", tree_dimuorphan_bg_m1, RooArgSet(m1TightBJet));
  RooDataSet* ds_dimuorphan_bg_m1MediumBJet = new RooDataSet("ds_dimuorphan_bg_m1MediumBJet","ds_dimuorphan_bg_m1MediumBJet", tree_dimuorphan_bg_m1, RooArgSet(m1MediumBJet));
  RooDataSet* ds_dimuorphan_bg_m1LooseBJet = new RooDataSet("ds_dimuorphan_bg_m1LooseBJet","ds_dimuorphan_bg_m1LooseBJet", tree_dimuorphan_bg_m1, RooArgSet(m1LooseBJet));

  RooDataSet* ds_dimuorphan_bg_m2 = new RooDataSet("ds_dimuorphan_bg_m2","ds_dimuorphan_bg_m2", tree_dimuorphan_bg_m2, RooArgSet(m2));
  RooDataSet* ds_dimuorphan_bg_m2TightBJet = new RooDataSet("ds_dimuorphan_bg_m2TightBJet","ds_dimuorphan_bg_m2TightBJet", tree_dimuorphan_bg_m2, RooArgSet(m2TightBJet));
  RooDataSet* ds_dimuorphan_bg_m2MediumBJet = new RooDataSet("ds_dimuorphan_bg_m2MediumBJet","ds_dimuorphan_bg_m2MediumBJet", tree_dimuorphan_bg_m2, RooArgSet(m2MediumBJet));
  RooDataSet* ds_dimuorphan_bg_m2LooseBJet = new RooDataSet("ds_dimuorphan_bg_m2LooseBJet","ds_dimuorphan_bg_m2LooseBJet", tree_dimuorphan_bg_m2, RooArgSet(m2LooseBJet));

  cout<<"-----Now Printing and importing the Datasets:-----"<<endl;
  ds_dimuorphan_bg_m1->Print("s");
  ds_dimuorphan_bg_m1TightBJet->Print("s");
  ds_dimuorphan_bg_m1MediumBJet->Print("s");
  ds_dimuorphan_bg_m1LooseBJet->Print("s");
  ds_dimuorphan_bg_m2->Print("s");
  ds_dimuorphan_bg_m2TightBJet->Print("s");
  ds_dimuorphan_bg_m2MediumBJet->Print("s");
  ds_dimuorphan_bg_m2LooseBJet->Print("s");

  w->import(*ds_dimuorphan_bg_m1);
  w->import(*ds_dimuorphan_bg_m1TightBJet);
  w->import(*ds_dimuorphan_bg_m1MediumBJet);
  w->import(*ds_dimuorphan_bg_m1LooseBJet);
  w->import(*ds_dimuorphan_bg_m2);
  w->import(*ds_dimuorphan_bg_m2TightBJet);
  w->import(*ds_dimuorphan_bg_m2MediumBJet);
  w->import(*ds_dimuorphan_bg_m2LooseBJet);

  //Draw before fiting
  RooPlot* plotC1 = w->var("m1")->frame(Title("m1 data tempalate NO FIT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC1, DataError(RooAbsData::SumW2), Name("data_m1"));
  float SizeBin1 = plotC1->GetXaxis()->GetBinCenter(3) - plotC1->GetXaxis()->GetBinCenter(2);
  char c_SizeBin1[10];
  snprintf(c_SizeBin1,50,"%.2f",SizeBin1);
  TString Yname1 = "Events / (" + std::string(c_SizeBin1) + "GeV)";
  plotC1->GetYaxis()->SetTitle( Yname1.Data() );
  TCanvas * c1 = new TCanvas("c1");
  c1->cd(); plotC1->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m1.pdf");

  //plot b-jets under m_i template cuts
  RooPlot* plotB1 = w->var("m1TightBJet")->frame(Title("m1: TightBJet"),Bins(b_bins));
  w->data("ds_dimuorphan_bg_m1TightBJet")->plotOn(plotB1, DataError(RooAbsData::SumW2), Name("data_m1TightBJet"));
  plotB1->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m1_TightBJet.png");

  RooPlot* plotB2 = w->var("m1MediumBJet")->frame(Title("m1: MediumBJet"),Bins(b_bins));
  w->data("ds_dimuorphan_bg_m1MediumBJet")->plotOn(plotB2, DataError(RooAbsData::SumW2), Name("data_m1MediumBJet"));
  plotB2->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m1_MediumBJet.png");

  RooPlot* plotB3 = w->var("m1LooseBJet")->frame(Title("m1: LooseBJet"),Bins(b_bins));
  w->data("ds_dimuorphan_bg_m1LooseBJet")->plotOn(plotB3, DataError(RooAbsData::SumW2), Name("data_m1LooseBJet"));
  plotB3->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m1_LooseBJet.png");

  RooPlot* plotC2 = w->var("m2")->frame(Title("m2 data tempalate NO FIT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotC2, DataError(RooAbsData::SumW2), Name("data_m2"));
  plotC2->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m2.pdf");

  RooPlot* plotB4 = w->var("m2TightBJet")->frame(Title("m2: TightBJet"),Bins(b_bins));
  w->data("ds_dimuorphan_bg_m2TightBJet")->plotOn(plotB4, DataError(RooAbsData::SumW2), Name("data_m2TightBJet"));
  plotB4->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m2_TightBJet.png");

  RooPlot* plotB5 = w->var("m2MediumBJet")->frame(Title("m2: MediumBJet"),Bins(b_bins));
  w->data("ds_dimuorphan_bg_m2MediumBJet")->plotOn(plotB5, DataError(RooAbsData::SumW2), Name("data_m2MediumBJet"));
  plotB5->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m2_MediumBJet.png");

  RooPlot* plotB6 = w->var("m2LooseBJet")->frame(Title("m2: LooseBJet"),Bins(b_bins));
  w->data("ds_dimuorphan_bg_m2LooseBJet")->plotOn(plotB6, DataError(RooAbsData::SumW2), Name("data_m2LooseBJet"));
  plotB6->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m2_LooseBJet.png");

  delete c1;

  //****************************************************************************
  //                         Create template for m1
  //****************************************************************************
  cout<<"-----Creating template for m1:-----"<<endl;
  //===========
  //2017 PDF m1
  //===========
  //Initial combianatorial
  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.01, 0.0, 0.3], MmumuC_p[0.05, 0.0, 1.5])");
  //Bulk shape
  w->factory("Bernstein::bgC(m1, {bC06[9.5766e+00, 0.1, 15.], bC16[ 1.0705e-01, 0., 3.], bC26[3.5184e-05, 0., 3.], bC36[1.259, 0., 5.], bC46[2.5370e-03, 0., 3.], bC56[5.7432e-01, 0., 3.], bC66[3.9353e-01, 0.1, 4.]})");
  //Ad hoc gaussian to cover first bump and help other functions
  w->factory("Gaussian::adHocC(m1, adHocC_mass[0.2, 0., 0.6], adHocC_sigma[6.3097e-02, 0.0001, 0.1])");
  //negligible resonances
  w->factory("Gaussian::etaC(m1, 0.54786, 0.007)");
  w->factory("Gaussian::rhoC(m1, 0.78265, 0.009)");
  w->factory("Gaussian::phiC(m1, 1.01946, 0.01)");
  w->factory("CBShape::JpsiC(m1, JpsiC_mean[3.0969, 3.0, 3.35], JpsiC_sigma[0.1, 0.001, 0.3], JpsiC_alpha[1.2, 0.4, 7.0], JpsiC_n[2.0])");
  w->factory("Gaussian::psiC(m1, 3.68609, psiC_sigma[0.031, 0.01, 0.04])");
  // FINAL PDF: Can be used to generate toy dataset: https://root.cern.ch/roofit-20-minutes
  w->factory("SUM::template1D_m1(norm_adHocC[20., 0., 10000.]*adHocC, norm_MmumuC[200., 0., 25000.]*MmumuC, norm_bgC[4400., 1000., 20000.]*bgC, norm_etaC[1.3151e+01, 0., 1000.]*etaC, norm_rhoC[1.0107e+02, 0., 1000.]*rhoC, norm_phiC[9.8640e+01, 0., 1000.]*phiC, norm_JpsiC[8000., 10., 10000.]*JpsiC, norm_psiC[50., 0., 1000.]*psiC)");
  //===============
  //End 2017 PDF m1
  //===============

  RooFitResult *rC = w->pdf("template1D_m1")->fitTo(*(w->data("ds_dimuorphan_bg_m1")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1---------------------"<<endl;
  rC->Print();

  RooPlot* plotC = w->var("m1")->frame(Title("BG template for orphan dimuon high pT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC, DataError(RooAbsData::SumW2), Name("data_m1"));
  w->pdf("template1D_m1")->plotOn(plotC,LineColor(kRed),Precision(0.0001),Name("template1D_m1"));
  float SizeBin = plotC->GetXaxis()->GetBinCenter(3) - plotC->GetXaxis()->GetBinCenter(2);
  char c_SizeBin[10];
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  TString Yname = "Events / (" + std::string(c_SizeBin) + "GeV)";
  plotC->GetYaxis()->SetTitle( Yname.Data() );
  //DATA-Fit ratio
  TCanvas * c_template1D_m1_RooPlot = new TCanvas("c_template1D_m1_RooPlot", "c_template1D_m1_RooPlot",800,800);
  c_template1D_m1_RooPlot->Clear();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  // THE FIT PART
  plotC->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m1_RooPlot->cd(); c_template1D_m1_RooPlot->Update();
  // Data/Fit
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
  TH1F *h_dataFit1 = new TH1F("h_dataFit1","", m_bins, m_min, m_max);
  h_dataFit1->SetLineColor(kBlack); h_dataFit1->SetMarkerStyle(20); h_dataFit1->SetMarkerColor(1); h_dataFit1->SetStats(0);
  h_dataFit1->GetXaxis()->SetTitle("Mass [GeV]");
  h_dataFit1->GetXaxis()->SetTitleSize(20);
  h_dataFit1->GetXaxis()->SetTitleFont(43);
  h_dataFit1->GetXaxis()->SetTitleOffset(3.0);
  h_dataFit1->GetXaxis()->SetLabelSize(15);
  h_dataFit1->GetXaxis()->SetLabelFont(43);//Absolute font size in pixel (precision 3)
  h_dataFit1->GetYaxis()->SetTitle("Data/Fit");
  h_dataFit1->GetYaxis()->SetNdivisions(505);
  h_dataFit1->GetYaxis()->CenterTitle();
  h_dataFit1->GetYaxis()->SetTitleSize(20);
  h_dataFit1->GetYaxis()->SetTitleFont(43);
  h_dataFit1->GetYaxis()->SetTitleOffset(.9);
  h_dataFit1->GetYaxis()->SetLabelSize(15);
  h_dataFit1->GetYaxis()->SetLabelFont(43);//#Absolute font size in pixel (precision 3)
  TH1F *hdata = (TH1F*) ds_dimuorphan_bg_m1->createHistogram("hdata",m1,Binning(m_bins, m_min, m_max));
  TH1F* h_func = new TH1F("h_func","", m_bins, m_min, m_max);
  w->pdf("template1D_m1")->fillHistogram(h_func, m1, hdata->GetEntries());
  for(unsigned int iB=1; iB<m_bins; iB++){
    float ratio = h_func->GetBinContent(iB)/hdata->GetBinContent(iB);
    h_dataFit1->SetBinContent(iB, ratio );
    h_dataFit1->SetBinError(iB, sqrt(hdata->GetBinContent(iB))/hdata->GetBinContent(iB) );
  }
  h_dataFit1->SetMinimum(0.5);
  h_dataFit1->SetMaximum(1.5);
  h_dataFit1->Sumw2();
  h_dataFit1->SetStats(0);
  h_dataFit1->SetMarkerStyle(21);
  h_dataFit1->Draw("ep");
  c_template1D_m1_RooPlot->SaveAs("figures/template1D_m1_RooPlot.pdf");
  c_template1D_m1_RooPlot->SaveAs("figures/template1D_m1_RooPlot.png");
  float chi2_C = plotC->chiSquare(23); //d.o.f = 23

  //****************************************************************************
  //                         Create template for m2
  //****************************************************************************
  cout<<"-----Creating template for m2:-----"<<endl;
  //===========
  //2017 PDF m2
  //===========
  //Initial combianatorial
  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[0.01, 0.0, 0.3], MmumuF_p[0.05, 0.0, 2.])");
  //FIRST KINK
  w->factory("Bernstein::bgF(m2,{bF06[9.9751, 0, 15.], bF16[3.2971e-05, 0., 3.], bF26[1.2361e-07, 0., 3.], bF36[5.1545e-08, 0., 2.], bF46[9.9017e-01, 0., 3.], bF56[3.0607e-01, 0., 3.], bF66[0.5, 0.1, 4.]})");
  // Ad hoc gaussian to cover first bump and help other functions
  w->factory("Gaussian::adHocF(m2,adHocF_mass[0.4, 0.2, 0.6],adHocF_sigma[0.01, 0.001, 0.1])");
  // Resonances
  w->factory("Gaussian::etaF(m2, 0.54786, 0.007)");
  w->factory("Gaussian::rhoF(m2, 0.78265, 0.009)");
  w->factory("Gaussian::phiF(m2, 1.01946, 0.01)");
  w->factory("CBShape::JpsiF(m2, JpsiF_mean[3.0969, 3.0, 3.35], JpsiF_sigma[0.1, 0.01, 0.3], JpsiF_alpha[1.2, 0.4, 10.0], JpsiF_n[2.0])");
  w->factory("Gaussian::psiF(m2, 3.68609, psiF_sigma[0.031, 0.01, 0.04])");
  // FINAL PDF
  w->factory("SUM::template1D_m2(norm_adHocF[150., 0., 500.]*adHocF, norm_MmumuF[10000., 0., 15000.]*MmumuF, norm_bgF[4400., 1000., 20000.]*bgF, norm_etaF[1., 0., 2.]*etaF, norm_rhoF[65., 1., 100.]*rhoF, norm_phiF[110., 1., 1000.]*phiF, norm_JpsiF[6400., 0., 10000.]*JpsiF, norm_psiF[250., 0., 1000.]*psiF)");
  //===============
  //End 2017 PDF m2
  //===============

  RooFitResult *rF = w->pdf("template1D_m2")->fitTo(*(w->data("ds_dimuorphan_bg_m2")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2---------------------"<<endl;
  rF->Print();

  RooPlot* plotF = w->var("m2")->frame(Title("BG template for orphan dimuon no high pT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotF, DataError(RooAbsData::SumW2), Name("data_m2"));
  w->pdf("template1D_m2")->plotOn(plotF,LineColor(kRed),Precision(0.0001),Name("template_m2"));

  SizeBin = plotF->GetXaxis()->GetBinCenter(3) - plotF->GetXaxis()->GetBinCenter(2);
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  Yname = "Events / (" + std::string(c_SizeBin) + "GeV)";
  plotF->GetYaxis()->SetTitle( Yname.Data() );
  //DATA-Fit ratio
  TCanvas * c_template1D_m2_RooPlot = new TCanvas("c_template1D_m2_RooPlot", "c_template1D_m2_RooPlot",800,800);
  c_template1D_m2_RooPlot->Clear();
  TPad *pad1b = new TPad("pad1b", "pad1b", 0, 0.3, 1, 1.0);
  pad1b->SetBottomMargin(0);
  pad1b->Draw();
  pad1b->cd();
  // THE FIT PART
  plotF->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m2_RooPlot->cd(); c_template1D_m2_RooPlot->Update();
  // Data/Fit
  TPad *pad2b = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2b->SetTopMargin(0);
  pad2b->SetBottomMargin(0.35);
  pad2b->SetGridy();
  pad2b->Draw();
  pad2b->cd();       // pad2 becomes the current pad
  TH1F *h_dataFit2 = new TH1F("h_dataFit2","", m_bins, m_min, m_max);
  h_dataFit2->SetLineColor(kBlack); h_dataFit2->SetMarkerStyle(20); h_dataFit2->SetMarkerColor(1); h_dataFit2->SetStats(0);
  h_dataFit2->GetXaxis()->SetTitle("Mass [GeV]");
  h_dataFit2->GetXaxis()->SetTitleSize(20);
  h_dataFit2->GetXaxis()->SetTitleFont(43);
  h_dataFit2->GetXaxis()->SetTitleOffset(3.0);
  h_dataFit2->GetXaxis()->SetLabelSize(15);
  h_dataFit2->GetXaxis()->SetLabelFont(43);//Absolute font size in pixel (precision 3)
  h_dataFit2->GetYaxis()->SetTitle("Data/Fit");
  h_dataFit2->GetYaxis()->SetNdivisions(505);
  h_dataFit2->GetYaxis()->CenterTitle();
  h_dataFit2->GetYaxis()->SetTitleSize(20);
  h_dataFit2->GetYaxis()->SetTitleFont(43);
  h_dataFit2->GetYaxis()->SetTitleOffset(.9);
  h_dataFit2->GetYaxis()->SetLabelSize(15);
  h_dataFit2->GetYaxis()->SetLabelFont(43);//#Absolute font size in pixel (precision 3)
  TH1F *hdata2 = (TH1F*) ds_dimuorphan_bg_m2->createHistogram("hdata2",m2,Binning(m_bins, m_min, m_max));
  TH1F* h_func2 = new TH1F("h_func2","", m_bins, m_min, m_max);
  w->pdf("template1D_m2")->fillHistogram(h_func2, m2, hdata2->GetEntries());
  for(unsigned int iB=1; iB<m_bins; iB++){
    float ratio = h_func2->GetBinContent(iB)/hdata2->GetBinContent(iB);
    h_dataFit2->SetBinContent(iB, ratio );
    h_dataFit2->SetBinError(iB, sqrt(hdata2->GetBinContent(iB))/hdata2->GetBinContent(iB) );
  }
  h_dataFit2->SetMinimum(0.5);  // Define Y ..
  h_dataFit2->SetMaximum(1.5); // .. range
  h_dataFit2->Sumw2();
  h_dataFit2->SetStats(0);      // No statistics on lower plot
  h_dataFit2->SetMarkerStyle(21);
  h_dataFit2->Draw("ep");
  c_template1D_m2_RooPlot->SaveAs("figures/template1D_m2_RooPlot.pdf");
  c_template1D_m2_RooPlot->SaveAs("figures/template1D_m2_RooPlot.png");
  float chi2_F = plotF->chiSquare(23);

  //****************************************************************************
  //                     Create 2D template = m1 x m2
  //****************************************************************************
  cout<<"-----Creating 2D template m1 * m2:-----"<<endl;
  w->factory("PROD::template2D(template1D_m1,template1D_m2)");
  cout<<"1D template m1 fit chi^2/dof: "<<chi2_C<<endl;
  cout<<"1D template m2 fit chi^2/dof: "<<chi2_F<<endl;

  //****************************************************************************
  //                     Extract 1D J/psi template from template m1
  //****************************************************************************
  cout << "Create templates for J/psi from m1" << endl;
  //Get final fit parameter for the J/psi peak in rC and plot the Jpsi peak
  RooRealVar* JpsiC_mean_fitresult = (RooRealVar*) rC->floatParsFinal().find("JpsiC_mean");
  cout << "JpsiC_mean " << JpsiC_mean_fitresult->getVal() << endl;
  RooRealVar* JpsiC_sigma_fitresult = (RooRealVar*) rC->floatParsFinal().find("JpsiC_sigma");
  cout << "JpsiC_sigma " << JpsiC_sigma_fitresult->getVal() << endl;
  RooRealVar* JpsiC_alpha_fitresult = (RooRealVar*) rC->floatParsFinal().find("JpsiC_alpha");
  cout << "JpsiC_alpha " << JpsiC_alpha_fitresult->getVal() << endl;

  RooRealVar Jpsi_m1_mean("Jpsi_m1_mean", "Jpsi_m1_mean", JpsiC_mean_fitresult->getVal());
  RooRealVar Jpsi_m1_sigma("Jpsi_m1_sigma", "Jpsi_m1_sigma", JpsiC_sigma_fitresult->getVal());
  RooRealVar Jpsi_m1_alpha("Jpsi_m1_alpha", "Jpsi_m1_alpha", JpsiC_alpha_fitresult->getVal());
  RooRealVar Jpsi_m1_n("Jpsi_m1_n", "Jpsi_m1_n", 2.0);

  RooCBShape Jpsi_m1("Jpsi_m1", "Jpsi_m1", m1, Jpsi_m1_mean, Jpsi_m1_sigma, Jpsi_m1_alpha, Jpsi_m1_n);
  w->import(Jpsi_m1);

  RooPlot* plot_Jpsi_m1 = w->var("m1")->frame(Title("J/psi template m1"),Bins(m_bins));
  w->pdf("Jpsi_m1")->plotOn(plot_Jpsi_m1,LineColor(kRed),Precision(0.0001),Name("plot_Jpsi_m1"));

  TCanvas * c_template1D_Jpsi_m1_RooPlot = new TCanvas("c_template1D_Jpsi_m1_RooPlot", "c_template1D_Jpsi_m1_RooPlot");
  c_template1D_Jpsi_m1_RooPlot->cd();
  SizeBin = plot_Jpsi_m1->GetXaxis()->GetBinCenter(3) - plot_Jpsi_m1->GetXaxis()->GetBinCenter(2);
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  Yname = "Events / (" + std::string(c_SizeBin) + "GeV)";
  plot_Jpsi_m1->GetYaxis()->SetTitle( Yname.Data() );
  plot_Jpsi_m1->Draw();
  txtHeader->Draw();
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1_RooPlot.pdf");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1_RooPlot.png");

  //****************************************************************************
  //                       Extract 1D J/psi template from template m2
  //****************************************************************************
  cout << "Create templates for J/psi from m2" << endl;
  RooRealVar* JpsiF_mean_fitresult = (RooRealVar*) rF->floatParsFinal().find("JpsiF_mean");
  cout << "JpsiF_mean " << JpsiF_mean_fitresult->getVal() << endl;
  RooRealVar* JpsiF_sigma_fitresult = (RooRealVar*) rF->floatParsFinal().find("JpsiF_sigma");
  cout << "JpsiF_sigma " << JpsiF_sigma_fitresult->getVal() << endl;
  RooRealVar* JpsiF_alpha_fitresult = (RooRealVar*) rF->floatParsFinal().find("JpsiF_alpha");
  cout << "JpsiF_alpha " << JpsiF_alpha_fitresult->getVal() << endl;

  RooRealVar Jpsi_m2_mean("Jpsi_m2_mean", "Jpsi_m2_mean", JpsiF_mean_fitresult->getVal());
  RooRealVar Jpsi_m2_sigma("Jpsi_m2_sigma", "Jpsi_m2_sigma", JpsiF_sigma_fitresult->getVal());
  RooRealVar Jpsi_m2_alpha("Jpsi_m2_alpha", "Jpsi_m2_alpha", JpsiF_alpha_fitresult->getVal());
  RooRealVar Jpsi_m2_n("Jpsi_m2_n", "Jpsi_m2_n", 2.0);

  RooCBShape Jpsi_m2("Jpsi_m2", "Jpsi_m2", m2, Jpsi_m2_mean, Jpsi_m2_sigma, Jpsi_m2_alpha, Jpsi_m2_n);
  w->import(Jpsi_m2);//normalized to 1

  RooPlot* plot_Jpsi_m2 = w->var("m2")->frame(Title("J/psi template m2"),Bins(m_bins));
  w->pdf("Jpsi_m2")->plotOn(plot_Jpsi_m2,LineColor(kRed),Precision(0.0001),Name("plot_Jpsi_m2"));

  TCanvas * c_template1D_Jpsi_m2_RooPlot = new TCanvas("c_template1D_Jpsi_m2_RooPlot", "c_template1D_Jpsi_m2_RooPlot");
  c_template1D_Jpsi_m2_RooPlot->cd();
  SizeBin = plot_Jpsi_m2->GetXaxis()->GetBinCenter(3) - plot_Jpsi_m2->GetXaxis()->GetBinCenter(2);
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  Yname = "Events / (" + std::string(c_SizeBin) + "GeV)";
  plot_Jpsi_m2->GetYaxis()->SetTitle( Yname.Data() );
  plot_Jpsi_m2->Draw();
  txtHeader->Draw();
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2_RooPlot.pdf");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2_RooPlot.png");

  //****************************************************************************
  //                     Create 2D template (m1 x m2) for J/psi
  //****************************************************************************
  cout << "Create 2D template (m1 * m2) for J/psi" << endl;
  w->factory("PROD::Jpsi_2D(Jpsi_m1,Jpsi_m2)");

  //**************************************************************************************
  //     First tryout to construct a generic 1D template, m1=m2,
  //     Don't distinguish high pT muon in orphan dimu or not
  //**************************************************************************************
  //Placeholder

  //****************************************************************************
  //           For later use in LowMassBKGPlot2D.C: datasets of 2 dimu events
  //****************************************************************************
  cout << "Create trees on signal events" << endl;

  //To be used for validate the method using 2 dimu data events at CR (removed iso for large stats)
  ostringstream stream_cut_control_offDiagonal;
  stream_cut_control_offDiagonal << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Abs(massC-massF) >= 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_min << " && massC < " << m_max << " && massF > " << m_min << " && massF < " << m_max;
  TString cut_control_offDiagonal = stream_cut_control_offDiagonal.str();
  TTree* tree_dimudimu_control_offDiagonal_2D       = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  
  //To be used for scatter plot 2 dimu events @ CR later in LowMassBKGPlot2D.C
  ostringstream stream_cut_control_Iso_offDiagonal;
  stream_cut_control_Iso_offDiagonal << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) >= 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_min << " && massC < " << m_max << " && massF > " << m_min << " && massF < " << m_max;
  TString cut_control_Iso_offDiagonal = stream_cut_control_Iso_offDiagonal.str();
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D       = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);

  //To be used to unblind the SR (scatter plot) later in LowMassBKGPlot2D.C (difference to above: mass cut <, not >)
  ostringstream stream_cut_signal;
  stream_cut_signal << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) < 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_min << " && massC < " << m_max << " && massF > " << m_min << " && massF < " << m_max;
  TString cut_signal = stream_cut_signal.str();
  TTree* tree_dimudimu_signal_2D = chain_data_dimudimu.CopyTree(cut_signal);

  //Setting Names
  tree_dimudimu_control_offDiagonal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_offDiagonal_2D->GetBranch("massF")->SetName("m2");
  tree_dimudimu_control_Iso_offDiagonal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_Iso_offDiagonal_2D->GetBranch("massF")->SetName("m2");
  tree_dimudimu_signal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_signal_2D->GetBranch("massF")->SetName("m2");

  //Creating 2 dimu dataset
  RooDataSet* ds_dimudimu_control_offDiagonal_2D = new RooDataSet("ds_dimudimu_control_offDiagonal_2D","ds_dimudimu_control_offDiagonal_2D", tree_dimudimu_control_offDiagonal_2D, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D","ds_dimudimu_control_Iso_offDiagonal_2D", tree_dimudimu_control_Iso_offDiagonal_2D, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_signal_2D = new RooDataSet("ds_dimudimu_signal_2D","ds_dimudimu_signal_2D", tree_dimudimu_signal_2D, RooArgSet(m1,m2));

  ds_dimudimu_control_offDiagonal_2D->Print("s");
  ds_dimudimu_control_Iso_offDiagonal_2D->Print("s");
  ds_dimudimu_signal_2D->Print("s");

  w->import(*ds_dimudimu_control_offDiagonal_2D);
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D);
  w->import(*ds_dimudimu_signal_2D);

  //****************************************************************************
  //                           Save to Workspace
  //****************************************************************************
  cout<<"Save to workspace"<<endl;
  w->writeToFile("ws_FINAL.root");
}
