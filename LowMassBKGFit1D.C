//******************************************************************************************************************
//* cmsenv                                                                                                         *
//* To request more time and memory: sintr -t 480 -m 102400                                                        *
//*                                       Wei Shi @Nov 20, 2019, Rice U.                                           *
//******************************************************************************************************************
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
#include "Constants.h"
#include "Config.h"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void LowMassBKGFit1D() {

  //Configure inputs for year
  BKG_cfg::ConfigureInput(year);
  setTDRStyle();

  TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader(header);

  //Output ws
  RooWorkspace* w = new RooWorkspace("w");
  TString Comm = "mkdir -p figures/";
  system( Comm.Data() );

  //Input file
  TChain chain_data_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  TChain chain_data_dimuorphan("cutFlowAnalyzerPXBL4PXFL3/Events_orphan");
  std::ifstream Myfile(inputFile1);
  std::string Line;
  if( !Myfile ) std::cout<<"ERROR opening Myfile."<<std::endl;
  while (std::getline(Myfile, Line)){
    TString Line2(Line);
    if( Line2.Contains("root") ){
      chain_data_dimudimu.Add(Line2.Data());
      chain_data_dimuorphan.Add(Line2.Data());
    }
  }

  //Define RooRealVar in workspace
  RooRealVar m1("m1", "m_{#mu#mu_{1}}", m_min, m_max, "GeV");
  RooRealVar m2("m2", "m_{#mu#mu_{2}}", m_min, m_max, "GeV");
  m1.setBins(m_bins);
  m2.setBins(m_bins);

  //For below/above J/psi
  RooRealVar m1_below_Jpsi("m1_below_Jpsi", "m_{#mu#mu_{1}}", m_min, m_Jpsi_dn, "GeV");
  RooRealVar m1_above_Jpsi("m1_above_Jpsi", "m_{#mu#mu_{1}}", m_Jpsi_up, m_max, "GeV");
  RooRealVar m2_below_Jpsi("m2_below_Jpsi", "m_{#mu#mu_{2}}", m_min, m_Jpsi_dn, "GeV");
  RooRealVar m2_above_Jpsi("m2_above_Jpsi", "m_{#mu#mu_{2}}", m_Jpsi_up, m_max, "GeV");
  m1_below_Jpsi.setBins(m_bins_below_Jpsi);
  m1_above_Jpsi.setBins(m_bins_above_Jpsi);
  m2_below_Jpsi.setBins(m_bins_below_Jpsi);
  m2_above_Jpsi.setBins(m_bins_above_Jpsi);

  w->import(m1);
  w->import(m2);
  w->import(m1_below_Jpsi);
  w->import(m1_above_Jpsi);
  w->import(m2_below_Jpsi);
  w->import(m2_above_Jpsi);

  //**********************************************************************************
  //     Select events for constructing 1D templates, identify high pT muon, m1!=m2
  //**********************************************************************************
  //m1: High pT mu is in orphan dimu
  ostringstream stream_cut_bg_m1_iso;
  stream_cut_bg_m1_iso << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) || ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso = stream_cut_bg_m1_iso.str();
  //Print selections to check
  std::cout << "m1 selctions (low mass): " << cut_bg_m1_iso.Data() << std::endl;
  TTree* tree_dimuorphan_bg_m1 = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso);

  //m1_below_Jpsi: Same as above except below J/psi
  ostringstream stream_cut_bg_m1_iso_below_Jpsi;
  stream_cut_bg_m1_iso_below_Jpsi << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) || ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_Jpsi_dn;
  TString cut_bg_m1_iso_below_Jpsi = stream_cut_bg_m1_iso_below_Jpsi.str();
  std::cout << "m1 selctions (low mass below J/psi): " << cut_bg_m1_iso_below_Jpsi.Data() << std::endl;
  TTree* tree_dimuorphan_bg_m1_below_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso_below_Jpsi);

  //m1_above_Jpsi: Same as above except above J/psi
  ostringstream stream_cut_bg_m1_iso_above_Jpsi;
  stream_cut_bg_m1_iso_above_Jpsi << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) || ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_Jpsi_up << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso_above_Jpsi = stream_cut_bg_m1_iso_above_Jpsi.str();
  std::cout << "m1 selctions (low mass above J/psi): " << cut_bg_m1_iso_above_Jpsi.Data() << std::endl;
  TTree* tree_dimuorphan_bg_m1_above_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso_above_Jpsi);

  //m2: Orphan mu is high pT
  ostringstream stream_cut_bg_m2_iso;
  stream_cut_bg_m2_iso << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m2_iso = stream_cut_bg_m2_iso.str();
  std::cout << "m2 selctions (low mass): " << cut_bg_m2_iso.Data() << std::endl;
  TTree* tree_dimuorphan_bg_m2 = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso);

  //m2: Same as above except below J/psi
  ostringstream stream_cut_bg_m2_iso_below_Jpsi;
  stream_cut_bg_m2_iso_below_Jpsi << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_Jpsi_dn;
  TString cut_bg_m2_iso_below_Jpsi = stream_cut_bg_m2_iso_below_Jpsi.str();
  std::cout << "m2 selctions (low mass below J/psi): " << cut_bg_m2_iso_below_Jpsi.Data() << std::endl;
  TTree* tree_dimuorphan_bg_m2_below_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso_below_Jpsi);

  //m2: Same as above except above J/psi
  ostringstream stream_cut_bg_m2_iso_above_Jpsi;
  stream_cut_bg_m2_iso_above_Jpsi << "orph_dimu_Mu0_isoTk0p3 < " << iso_cut << " && orph_dimu_Mu0_isoTk0p3 >= 0 && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_AllTrackerMu && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_Jpsi_up << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m2_iso_above_Jpsi = stream_cut_bg_m2_iso_above_Jpsi.str();
  std::cout << "m2 selctions (low mass above J/psi): " << cut_bg_m2_iso_above_Jpsi.Data() << std::endl;
  TTree* tree_dimuorphan_bg_m2_above_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso_above_Jpsi);

  //Setting Names: Control region
  tree_dimuorphan_bg_m1->GetBranch("orph_dimu_mass")->SetName("m1");
  tree_dimuorphan_bg_m2->GetBranch("orph_dimu_mass")->SetName("m2");

  tree_dimuorphan_bg_m1_below_Jpsi->GetBranch("orph_dimu_mass")->SetName("m1_below_Jpsi");
  tree_dimuorphan_bg_m1_above_Jpsi->GetBranch("orph_dimu_mass")->SetName("m1_above_Jpsi");

  tree_dimuorphan_bg_m2_below_Jpsi->GetBranch("orph_dimu_mass")->SetName("m2_below_Jpsi");
  tree_dimuorphan_bg_m2_above_Jpsi->GetBranch("orph_dimu_mass")->SetName("m2_above_Jpsi");

  //Constructor of a data set from (part of) an ROOT TTRee:
  //The dimensions of the data set are defined by the 'vars' RooArgSet. For each dimension specified, the TTree must have a branch with the same name.
  RooDataSet* ds_dimuorphan_bg_m1 = new RooDataSet("ds_dimuorphan_bg_m1", "ds_dimuorphan_bg_m1", tree_dimuorphan_bg_m1, RooArgSet(m1));
  RooDataSet* ds_dimuorphan_bg_m2 = new RooDataSet("ds_dimuorphan_bg_m2", "ds_dimuorphan_bg_m2", tree_dimuorphan_bg_m2, RooArgSet(m2));

  RooDataSet* ds_dimuorphan_bg_m1_below_Jpsi = new RooDataSet("ds_dimuorphan_bg_m1_below_Jpsi", "ds_dimuorphan_bg_m1_below_Jpsi", tree_dimuorphan_bg_m1_below_Jpsi, RooArgSet(m1_below_Jpsi));
  RooDataSet* ds_dimuorphan_bg_m1_above_Jpsi = new RooDataSet("ds_dimuorphan_bg_m1_above_Jpsi", "ds_dimuorphan_bg_m1_above_Jpsi", tree_dimuorphan_bg_m1_above_Jpsi, RooArgSet(m1_above_Jpsi));

  RooDataSet* ds_dimuorphan_bg_m2_below_Jpsi = new RooDataSet("ds_dimuorphan_bg_m2_below_Jpsi", "ds_dimuorphan_bg_m2_below_Jpsi", tree_dimuorphan_bg_m2_below_Jpsi, RooArgSet(m2_below_Jpsi));
  RooDataSet* ds_dimuorphan_bg_m2_above_Jpsi = new RooDataSet("ds_dimuorphan_bg_m2_above_Jpsi", "ds_dimuorphan_bg_m2_above_Jpsi", tree_dimuorphan_bg_m2_above_Jpsi, RooArgSet(m2_above_Jpsi));

  cout<<"-----Now Printing and importing the Datasets:-----"<<endl;
  ds_dimuorphan_bg_m1->Print("s");
  ds_dimuorphan_bg_m2->Print("s");
  ds_dimuorphan_bg_m1_below_Jpsi->Print("s");
  ds_dimuorphan_bg_m1_above_Jpsi->Print("s");
  ds_dimuorphan_bg_m2_below_Jpsi->Print("s");
  ds_dimuorphan_bg_m2_above_Jpsi->Print("s");

  w->import(*ds_dimuorphan_bg_m1);
  w->import(*ds_dimuorphan_bg_m2);
  w->import(*ds_dimuorphan_bg_m1_below_Jpsi);
  w->import(*ds_dimuorphan_bg_m2_below_Jpsi);
  w->import(*ds_dimuorphan_bg_m1_above_Jpsi);
  w->import(*ds_dimuorphan_bg_m2_above_Jpsi);

  //Draw before fiting
  RooPlot* plotC1 = w->var("m1")->frame(Title("m1 data tempalate NO FIT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC1, DataError(RooAbsData::SumW2), Name("data_m1"));
  plotC1->GetYaxis()->SetTitle("Events/0.04GeV");
  TCanvas * c1 = new TCanvas("c1");
  c1->cd(); plotC1->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1.pdf");
  c1->SaveAs("figures/dimuorphan_m1.png");
  c1->SaveAs("figures/dimuorphan_m1.root");

  RooPlot* plotC2 = w->var("m2")->frame(Title("m2 data tempalate NO FIT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotC2, DataError(RooAbsData::SumW2), Name("data_m2"));
  plotC2->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC2->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2.pdf");
  c1->SaveAs("figures/dimuorphan_m2.png");
  c1->SaveAs("figures/dimuorphan_m2.root");

  RooPlot* plotC1_below_Jpsi = w->var("m1_below_Jpsi")->frame(Title("m1 data tempalate NO FIT below Jpsi"),Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m1_below_Jpsi")->plotOn(plotC1_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_below_Jpsi"));
  plotC1_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC1_below_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1_below_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m1_below_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m1_below_Jpsi.root");

  RooPlot* plotC1_above_Jpsi = w->var("m1_above_Jpsi")->frame(Title("m1 data tempalate NO FIT above Jpsi"),Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m1_above_Jpsi")->plotOn(plotC1_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_above_Jpsi"));
  plotC1_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC1_above_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1_above_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m1_above_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m1_above_Jpsi.root");

  RooPlot* plotC2_below_Jpsi = w->var("m2_below_Jpsi")->frame(Title("m2 data tempalate NO FIT below Jpsi"),Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m2_below_Jpsi")->plotOn(plotC2_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_below_Jpsi"));
  plotC2_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC2_below_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2_below_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m2_below_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m2_below_Jpsi.root");

  RooPlot* plotC2_above_Jpsi = w->var("m2_above_Jpsi")->frame(Title("m2 data tempalate NO FIT above Jpsi"),Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m2_above_Jpsi")->plotOn(plotC2_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_above_Jpsi"));
  plotC2_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC2_above_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2_above_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m2_above_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m2_above_Jpsi.root");

  delete c1;

  //****************************************************************************
  //                         Create template for m1
  //****************************************************************************
  cout<<"-----Creating template for m1:-----"<<endl;
  //===========
  //2017 PDF m1
  //===========
  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.01, 0., 0.3], MmumuC_p[0.05, 0., 2.5])");
  w->factory("Bernstein::bgC(m1, {bC06[9.5766e+00, 0.1, 25.], bC16[ 1.0705e-01, -5., 3.], bC26[3.5184e-05, -5., 3.], bC36[1.259, -3., 5.], bC46[2.5370e-03, -3., 3.], bC56[5.7432e-01, -3., 3.], bC66[3.9353e-01, -3., 4.]})");
  w->factory("Gaussian::adHocC(m1, adHocC_mass[0.2, 0., 0.6], adHocC_sigma[6.3097e-02, 0., 0.5])");
  w->factory("Gaussian::etaC(m1, 0.54786, 0.007)");
  w->factory("Gaussian::rhoC(m1, 0.78265, 0.009)");
  w->factory("Gaussian::phiC(m1, 1.01946, 0.01)");
  w->factory("CBShape::JpsiC(m1, JpsiC_mean[3.0969, 3.0, 3.35], JpsiC_sigma[0.1, 0., 0.3], JpsiC_alpha[1.2, 0.4, 7.0], JpsiC_n[2.0])");
  w->factory("Gaussian::psiC(m1, 3.68609, psiC_sigma[0.031, 0.01, 0.04])");

  w->factory("SUM::template1D_m1(norm_adHocC[20., 0., 10000.]*adHocC, norm_MmumuC[200., 0., 35000.]*MmumuC, norm_bgC[4400., 1000., 30000.]*bgC, norm_etaC[1.3151e+01, 0., 1000.]*etaC, norm_rhoC[1.0107e+02, 0., 1000.]*rhoC, norm_phiC[9.8640e+01, 0., 1000.]*phiC, norm_JpsiC[8000., 10., 20000.]*JpsiC, norm_psiC[50., 0., 1000.]*psiC)");
  RooFitResult *rC = w->pdf("template1D_m1")->fitTo(*(w->data("ds_dimuorphan_bg_m1")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1---------------------"<<endl;
  rC->Print();

  RooPlot* plotC = w->var("m1")->frame(Title("1D template for orphan dimuon high pT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC, DataError(RooAbsData::SumW2), Name("data_m1"));
  w->pdf("template1D_m1")->plotOn(plotC,LineColor(kRed),Precision(0.0001),Name("template1D_m1"));
  plotC->GetYaxis()->SetTitle("Events/0.04GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1 = new TCanvas("c_template1D_m1", "c_template1D_m1",800,800);
  c_template1D_m1->Clear();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); pad1->Draw(); pad1->cd();
  plotC->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m1->cd(); c_template1D_m1->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2->SetTopMargin(0); pad2->SetBottomMargin(0.35); pad2->SetGridy(); pad2->Draw(); pad2->cd();
  TH1F *h_dataFit1 = new TH1F("h_dataFit1","", m_bins, m_min, m_max);
  h_dataFit1->SetLineColor(kBlack); h_dataFit1->SetMarkerStyle(20); h_dataFit1->SetMarkerColor(1); h_dataFit1->SetStats(0);
  h_dataFit1->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1->GetXaxis()->SetTitleSize(20); h_dataFit1->GetXaxis()->SetTitleFont(43); h_dataFit1->GetXaxis()->SetTitleOffset(3.0); h_dataFit1->GetXaxis()->SetLabelSize(15); h_dataFit1->GetXaxis()->SetLabelFont(43);
  h_dataFit1->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1->GetYaxis()->SetNdivisions(505); h_dataFit1->GetYaxis()->CenterTitle(); h_dataFit1->GetYaxis()->SetTitleSize(20); h_dataFit1->GetYaxis()->SetTitleFont(43); h_dataFit1->GetYaxis()->SetTitleOffset(.9); h_dataFit1->GetYaxis()->SetLabelSize(15); h_dataFit1->GetYaxis()->SetLabelFont(43);
  TH1F *hdata = (TH1F*) ds_dimuorphan_bg_m1->createHistogram("hdata",m1,Binning(m_bins, m_min, m_max));
  TH1F* h_func = new TH1F("h_func","", m_bins, m_min, m_max);
  w->pdf("template1D_m1")->fillHistogram(h_func, m1, hdata->GetEntries());
  for(unsigned int iB=1; iB<m_bins; iB++){
    float ratio = h_func->GetBinContent(iB)/hdata->GetBinContent(iB);
    h_dataFit1->SetBinContent(iB, ratio );
    h_dataFit1->SetBinError(iB, ratio/sqrt(hdata->GetBinContent(iB)) );
  }
  h_dataFit1->SetMinimum(0.5); h_dataFit1->SetMaximum(1.5); h_dataFit1->Sumw2(); h_dataFit1->SetStats(0); h_dataFit1->SetMarkerStyle(21); h_dataFit1->Draw("ep");
  c_template1D_m1->SaveAs("figures/template1D_m1.pdf");
  c_template1D_m1->SaveAs("figures/template1D_m1.png");
  c_template1D_m1->SaveAs("figures/template1D_m1.root");
  float chi2_C = plotC->chiSquare(23); //d.o.f = 23, Refer to: https://root.cern.ch/doc/master/RooPlot_8h_source.html
  cout<<"------------------ End m1 ---------------------"<<endl;

  //=========================
  //2017 PDF m1 (below J/psi)
  //=========================
  w->factory("EXPR::MmumuC_below_Jpsi('m1_below_Jpsi*pow( (m1_below_Jpsi/m_below_Jpsi)*(m1_below_Jpsi/m_below_Jpsi) - 1.0, MmumuC_p_below_Jpsi )*exp( -MmumuC_c_below_Jpsi*( (m1_below_Jpsi/m_below_Jpsi)*(m1_below_Jpsi/m_below_Jpsi) - 1.0 ) )', m1_below_Jpsi, m_below_Jpsi[0.2113], MmumuC_c_below_Jpsi[0.1, -0.3, 0.3], MmumuC_p_below_Jpsi[0.1, 0.0, 1.5])");
  w->factory("Bernstein::bgC_below_Jpsi(m1_below_Jpsi, {bC06_below_Jpsi[1., 0.1, 15.], bC16_below_Jpsi[1.3, 0., 5.], bC26_below_Jpsi[0.1, -10., 3.], bC36_below_Jpsi[0.1, -3., 5.], bC46_below_Jpsi[1.5, -3., 3.], bC56_below_Jpsi[0.5, 0., 3.], bC66_below_Jpsi[0.1, -1., 4.]})");
  w->factory("Gaussian::adHocC_below_Jpsi(m1_below_Jpsi, adHocC_mass_below_Jpsi[0.3, 0., 0.6], adHocC_sigma_below_Jpsi[0.07, 0., 0.1])");
  w->factory("Gaussian::etaC_below_Jpsi(m1_below_Jpsi, 0.54786, 0.007)");
  w->factory("Gaussian::rhoC_below_Jpsi(m1_below_Jpsi, 0.78265, 0.009)");
  w->factory("Gaussian::phiC_below_Jpsi(m1_below_Jpsi, 1.01946, 0.01)");

  w->factory("SUM::template1D_m1_below_Jpsi(norm_adHocC_below_Jpsi[1000., 0., 10000.]*adHocC_below_Jpsi, norm_MmumuC_below_Jpsi[14000., 0., 35000.]*MmumuC_below_Jpsi, norm_bgC_below_Jpsi[6000., 1000., 25000.]*bgC_below_Jpsi, norm_etaC_below_Jpsi[0.08, 0., 1000.]*etaC_below_Jpsi, norm_rhoC_below_Jpsi[40., 0., 1000.]*rhoC_below_Jpsi, norm_phiC_below_Jpsi[58., 0., 1000.]*phiC_below_Jpsi)");
  RooFitResult *rC_below_Jpsi = w->pdf("template1D_m1_below_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m1_below_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1 (below J/psi)---------------------"<<endl;
  rC_below_Jpsi->Print();

  RooPlot* plotC_below_Jpsi = w->var("m1_below_Jpsi")->frame(Title("1D template for orphan dimuon high pT below Jpsi"),Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m1_below_Jpsi")->plotOn(plotC_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_below_Jpsi"));
  w->pdf("template1D_m1_below_Jpsi")->plotOn(plotC_below_Jpsi, LineColor(kRed), Precision(0.0001), Name("template1D_m1_below_Jpsi")); // By default only fitted range is shown
  plotC_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1_below_Jpsi = new TCanvas("c_template1D_m1_below_Jpsi", "c_template1D_m1_below_Jpsi",800,800);
  c_template1D_m1_below_Jpsi->Clear();
  TPad *pad1_below_Jpsi = new TPad("pad1_below_Jpsi", "pad1_below_Jpsi", 0, 0.3, 1, 1.0);
  pad1_below_Jpsi->SetBottomMargin(0); pad1_below_Jpsi->Draw(); pad1_below_Jpsi->cd();
  plotC_below_Jpsi->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m1_below_Jpsi->cd(); c_template1D_m1_below_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2_below_Jpsi = new TPad("pad2_below_Jpsi", "pad2_below_Jpsi", 0, 0.0, 1, 0.29);
  pad2_below_Jpsi->SetTopMargin(0); pad2_below_Jpsi->SetBottomMargin(0.35); pad2_below_Jpsi->SetGridy(); pad2_below_Jpsi->Draw(); pad2_below_Jpsi->cd();
  TH1F *h_dataFit1_below_Jpsi = new TH1F("h_dataFit1_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  h_dataFit1_below_Jpsi->SetLineColor(kBlack); h_dataFit1_below_Jpsi->SetMarkerStyle(20); h_dataFit1_below_Jpsi->SetMarkerColor(1); h_dataFit1_below_Jpsi->SetStats(0);
  h_dataFit1_below_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1_below_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit1_below_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit1_below_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit1_below_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit1_below_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit1_below_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1_below_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit1_below_Jpsi->GetYaxis()->CenterTitle(); h_dataFit1_below_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit1_below_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit1_below_Jpsi->GetYaxis()->SetTitleOffset(.9); h_dataFit1_below_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit1_below_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata_below_Jpsi = (TH1F*) ds_dimuorphan_bg_m1_below_Jpsi->createHistogram("hdata_below_Jpsi", m1_below_Jpsi, Binning(m_bins_below_Jpsi, m_min, m_Jpsi_dn));
  TH1F* h_func_below_Jpsi = new TH1F("h_func_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  w->pdf("template1D_m1_below_Jpsi")->fillHistogram(h_func_below_Jpsi, m1_below_Jpsi, hdata_below_Jpsi->GetEntries());
  for(unsigned int iB=1; iB<m_bins_below_Jpsi; iB++){
    float ratio = h_func_below_Jpsi->GetBinContent(iB)/hdata_below_Jpsi->GetBinContent(iB);
    h_dataFit1_below_Jpsi->SetBinContent(iB, ratio );
    h_dataFit1_below_Jpsi->SetBinError(iB, ratio/sqrt(hdata_below_Jpsi->GetBinContent(iB)) );
  }
  h_dataFit1_below_Jpsi->SetMinimum(0.5); h_dataFit1_below_Jpsi->SetMaximum(1.5); h_dataFit1_below_Jpsi->Sumw2(); h_dataFit1_below_Jpsi->SetStats(0); h_dataFit1_below_Jpsi->SetMarkerStyle(21); h_dataFit1_below_Jpsi->Draw("ep");
  c_template1D_m1_below_Jpsi->SaveAs("figures/template1D_m1_below_Jpsi.pdf");
  c_template1D_m1_below_Jpsi->SaveAs("figures/template1D_m1_below_Jpsi.png");
  c_template1D_m1_below_Jpsi->SaveAs("figures/template1D_m1_below_Jpsi.root");
  float chi2_C_below_Jpsi = plotC_below_Jpsi->chiSquare(17);
  cout<<"------------------ End m1 (below Jpsi) ---------------------"<<endl;

  //=========================
  //2017 PDF m1 (above J/psi)
  //=========================
  w->factory("Bernstein::bgC_above_Jpsi(m1_above_Jpsi, {bC06_above_Jpsi[4.6, 0.1, 15.], bC16_above_Jpsi[-3., -5., 3.], bC26_above_Jpsi[2.9, -3., 5.], bC36_above_Jpsi[0.1, -2., 5.], bC46_above_Jpsi[-0.1, -3., 3.], bC56_above_Jpsi[0.6, 0., 3.], bC66_above_Jpsi[0.6, 0.1, 4.]})");
  w->factory("Gaussian::psiC_above_Jpsi(m1_above_Jpsi, 3.68609, psiC_sigma_above_Jpsi[0.03, 0.01, 0.04])");

  w->factory("SUM::template1D_m1_above_Jpsi(norm_bgC_above_Jpsi[3500., 1000., 20000.]*bgC_above_Jpsi, norm_psiC_above_Jpsi[400., 0., 1000.]*psiC_above_Jpsi)");
  RooFitResult *rC_above_Jpsi = w->pdf("template1D_m1_above_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m1_above_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1 (above J/psi)---------------------"<<endl;
  rC_above_Jpsi->Print();

  RooPlot* plotC_above_Jpsi = w->var("m1_above_Jpsi")->frame(Title("1D template for orphan dimuon high pT above Jpsi"),Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m1_above_Jpsi")->plotOn(plotC_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_above_Jpsi"));
  w->pdf("template1D_m1_above_Jpsi")->plotOn(plotC_above_Jpsi, LineColor(kRed), Precision(0.0001), Name("template1D_m1_above_Jpsi"));
  plotC_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1_above_Jpsi = new TCanvas("c_template1D_m1_above_Jpsi", "c_template1D_m1_above_Jpsi",800,800);
  c_template1D_m1_above_Jpsi->Clear();
  TPad *pad1_above_Jpsi = new TPad("pad1_above_Jpsi", "pad1_above_Jpsi", 0, 0.3, 1, 1.0);
  pad1_above_Jpsi->SetBottomMargin(0); pad1_above_Jpsi->Draw(); pad1_above_Jpsi->cd();
  plotC_above_Jpsi->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m1_above_Jpsi->cd(); c_template1D_m1_above_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2_above_Jpsi = new TPad("pad2_above_Jpsi", "pad2_above_Jpsi", 0, 0.0, 1, 0.29);
  pad2_above_Jpsi->SetTopMargin(0); pad2_above_Jpsi->SetBottomMargin(0.35); pad2_above_Jpsi->SetGridy(); pad2_above_Jpsi->Draw(); pad2_above_Jpsi->cd();
  TH1F *h_dataFit1_above_Jpsi = new TH1F("h_dataFit1_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  h_dataFit1_above_Jpsi->SetLineColor(kBlack); h_dataFit1_above_Jpsi->SetMarkerStyle(20); h_dataFit1_above_Jpsi->SetMarkerColor(1); h_dataFit1_above_Jpsi->SetStats(0);
  h_dataFit1_above_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1_above_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit1_above_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit1_above_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit1_above_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit1_above_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit1_above_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1_above_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit1_above_Jpsi->GetYaxis()->CenterTitle(); h_dataFit1_above_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit1_above_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit1_above_Jpsi->GetYaxis()->SetTitleOffset(.9); h_dataFit1_above_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit1_above_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata_above_Jpsi = (TH1F*) ds_dimuorphan_bg_m1_above_Jpsi->createHistogram("hdata_above_Jpsi", m1_above_Jpsi, Binning(m_bins_above_Jpsi, m_Jpsi_up, m_max));
  TH1F* h_func_above_Jpsi = new TH1F("h_func_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  w->pdf("template1D_m1_above_Jpsi")->fillHistogram(h_func_above_Jpsi, m1_above_Jpsi, hdata_above_Jpsi->GetEntries());
  for(unsigned int iB=1; iB<m_bins_above_Jpsi; iB++){
    float ratio = h_func_above_Jpsi->GetBinContent(iB)/hdata_above_Jpsi->GetBinContent(iB);
    h_dataFit1_above_Jpsi->SetBinContent(iB, ratio );
    h_dataFit1_above_Jpsi->SetBinError(iB, ratio/sqrt(hdata_above_Jpsi->GetBinContent(iB)) );
  }
  h_dataFit1_above_Jpsi->SetMinimum(0.5); h_dataFit1_above_Jpsi->SetMaximum(1.5); h_dataFit1_above_Jpsi->Sumw2(); h_dataFit1_above_Jpsi->SetStats(0); h_dataFit1_above_Jpsi->SetMarkerStyle(21); h_dataFit1_above_Jpsi->Draw("ep");
  c_template1D_m1_above_Jpsi->SaveAs("figures/template1D_m1_above_Jpsi.pdf");
  c_template1D_m1_above_Jpsi->SaveAs("figures/template1D_m1_above_Jpsi.png");
  c_template1D_m1_above_Jpsi->SaveAs("figures/template1D_m1_above_Jpsi.root");
  float chi2_C_above_Jpsi = plotC_above_Jpsi->chiSquare(10);
  cout<<"------------------ End m1 (above Jpsi) ---------------------"<<endl;

  //****************************************************************************
  //                         Create template for m2
  //****************************************************************************
  cout<<"-----Creating template for m2:-----"<<endl;
  //===========
  //2017 PDF m2
  //===========
  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[0.01, 0., 0.3], MmumuF_p[0.05, 0., 4.])");
  w->factory("Bernstein::bgF(m2, {bF06[9.9751, 0, 30.], bF16[3.2971e-05, -6., 5.], bF26[1.2361e-07, -5., 5.], bF36[5.1545e-08, -3., 10.], bF46[9.9017e-01, -3., 5.], bF56[3.0607e-01, -3., 3.], bF66[0.5, 0.1, 4.]})");
  w->factory("Gaussian::adHocF(m2, adHocF_mass[0.4, 0.1, 0.6], adHocF_sigma[0.01, 0.001, 0.1])");
  w->factory("Gaussian::etaF(m2, 0.54786, 0.007)");
  w->factory("Gaussian::rhoF(m2, 0.78265, 0.009)");
  w->factory("Gaussian::phiF(m2, 1.01946, 0.01)");
  w->factory("CBShape::JpsiF(m2, JpsiF_mean[3.0969, 3.0, 3.35], JpsiF_sigma[0.1, 0., 0.3], JpsiF_alpha[1.2, 0.4, 10.0], JpsiF_n[2.0])");
  w->factory("Gaussian::psiF(m2, 3.68609, psiF_sigma[0.031, 0.01, 0.05])");

  w->factory("SUM::template1D_m2(norm_adHocF[150., 0., 500.]*adHocF, norm_MmumuF[10000., 0., 35000.]*MmumuF, norm_bgF[4400., 1000., 25000.]*bgF, norm_etaF[1., 0., 2.]*etaF, norm_rhoF[65., 1., 100.]*rhoF, norm_phiF[110., 1., 1000.]*phiF, norm_JpsiF[6400., 0., 20000.]*JpsiF, norm_psiF[250., 0., 1000.]*psiF)");
  RooFitResult *rF = w->pdf("template1D_m2")->fitTo(*(w->data("ds_dimuorphan_bg_m2")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2---------------------"<<endl;
  rF->Print();

  RooPlot* plotF = w->var("m2")->frame(Title("BG template for orphan dimuon no high pT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotF, DataError(RooAbsData::SumW2), Name("data_m2"));
  w->pdf("template1D_m2")->plotOn(plotF,LineColor(kRed),Precision(0.0001),Name("template1D_m2"));
  plotF->GetYaxis()->SetTitle("Events/0.04GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2 = new TCanvas("c_template1D_m2", "c_template1D_m2", 800, 800);
  c_template1D_m2->Clear();
  TPad *pad1b = new TPad("pad1b", "pad1b", 0, 0.3, 1, 1.0);
  pad1b->SetBottomMargin(0); pad1b->Draw(); pad1b->cd();
  plotF->Draw("same"); txtHeader->Draw("same");
  c_template1D_m2->cd(); c_template1D_m2->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2b->SetTopMargin(0); pad2b->SetBottomMargin(0.35); pad2b->SetGridy(); pad2b->Draw(); pad2b->cd();
  TH1F *h_dataFit2 = new TH1F("h_dataFit2","", m_bins, m_min, m_max);
  h_dataFit2->SetLineColor(kBlack); h_dataFit2->SetMarkerStyle(20); h_dataFit2->SetMarkerColor(1); h_dataFit2->SetStats(0);
  h_dataFit2->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2->GetXaxis()->SetTitleSize(20); h_dataFit2->GetXaxis()->SetTitleFont(43); h_dataFit2->GetXaxis()->SetTitleOffset(3.0); h_dataFit2->GetXaxis()->SetLabelSize(15); h_dataFit2->GetXaxis()->SetLabelFont(43);
  h_dataFit2->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2->GetYaxis()->SetNdivisions(505); h_dataFit2->GetYaxis()->CenterTitle(); h_dataFit2->GetYaxis()->SetTitleSize(20); h_dataFit2->GetYaxis()->SetTitleFont(43); h_dataFit2->GetYaxis()->SetTitleOffset(.9); h_dataFit2->GetYaxis()->SetLabelSize(15); h_dataFit2->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2 = (TH1F*) ds_dimuorphan_bg_m2->createHistogram("hdata2",m2,Binning(m_bins, m_min, m_max));
  TH1F* h_func2 = new TH1F("h_func2","", m_bins, m_min, m_max);
  w->pdf("template1D_m2")->fillHistogram(h_func2, m2, hdata2->GetEntries());
  for(unsigned int iB=1; iB<m_bins; iB++){
    float ratio = h_func2->GetBinContent(iB)/hdata2->GetBinContent(iB);
    h_dataFit2->SetBinContent(iB, ratio );
    h_dataFit2->SetBinError(iB, ratio/sqrt(hdata2->GetBinContent(iB)) );
  }
  h_dataFit2->SetMinimum(0.5); h_dataFit2->SetMaximum(1.5); h_dataFit2->Sumw2(); h_dataFit2->SetStats(0); h_dataFit2->SetMarkerStyle(21); h_dataFit2->Draw("ep");
  c_template1D_m2->SaveAs("figures/template1D_m2.pdf");
  c_template1D_m2->SaveAs("figures/template1D_m2.png");
  c_template1D_m2->SaveAs("figures/template1D_m2.root");
  float chi2_F = plotF->chiSquare(23);
  cout<<"------------------ End m2 ---------------------"<<endl;

  //=========================
  //2017 PDF m2 (below J/psi)
  //=========================
  w->factory("EXPR::MmumuF_below_Jpsi('m2_below_Jpsi*pow( (m2_below_Jpsi/m_below_Jpsi)*(m2_below_Jpsi/m_below_Jpsi) - 1.0, MmumuF_p_below_Jpsi )*exp( -MmumuF_c_below_Jpsi*( (m2_below_Jpsi/m_below_Jpsi)*(m2_below_Jpsi/m_below_Jpsi) - 1.0 ) )', m2_below_Jpsi, m_below_Jpsi[0.2113], MmumuF_c_below_Jpsi[0.1, -0.3, 0.3], MmumuF_p_below_Jpsi[0.05, 0.0, 2.])");
  w->factory("Bernstein::bgF_below_Jpsi(m2_below_Jpsi, {bF06_below_Jpsi[0.7, 0, 15.], bF16_below_Jpsi[1., 0., 3.], bF26_below_Jpsi[0.1, -3., 3.], bF36_below_Jpsi[0.3, 0., 2.], bF46_below_Jpsi[3., 0., 5.], bF56_below_Jpsi[0.5, 0., 3.], bF66_below_Jpsi[0.1, -1, 4.]})");
  w->factory("Gaussian::adHocF_below_Jpsi(m2_below_Jpsi, adHocF_mass_below_Jpsi[0.3, 0.2, 0.6], adHocF_sigma_below_Jpsi[0.06, 0.001, 0.3])");
  w->factory("Gaussian::etaF_below_Jpsi(m2_below_Jpsi, 0.54786, 0.007)");
  w->factory("Gaussian::rhoF_below_Jpsi(m2_below_Jpsi, 0.78265, 0.009)");
  w->factory("Gaussian::phiF_below_Jpsi(m2_below_Jpsi, 1.01946, 0.01)");

  w->factory("SUM::template1D_m2_below_Jpsi(norm_adHocF_below_Jpsi[400., 0., 700.]*adHocF_below_Jpsi, norm_MmumuF_below_Jpsi[10000., 0., 15000.]*MmumuF_below_Jpsi, norm_bgF_below_Jpsi[5000., 1000., 20000.]*bgF_below_Jpsi, norm_etaF_below_Jpsi[2., 0., 4.]*etaF_below_Jpsi, norm_rhoF_below_Jpsi[17., 1., 100.]*rhoF_below_Jpsi, norm_phiF_below_Jpsi[32., 1., 1000.]*phiF_below_Jpsi)");
  RooFitResult *rF_below_Jpsi = w->pdf("template1D_m2_below_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m2_below_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2 (below J/psi)---------------------"<<endl;
  rF_below_Jpsi->Print();

  RooPlot* plotF_below_Jpsi = w->var("m2_below_Jpsi")->frame(Title("BG template for orphan dimuon no high pT below Jpsi"),Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m2_below_Jpsi")->plotOn(plotF_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_below_Jpsi"));
  w->pdf("template1D_m2_below_Jpsi")->plotOn(plotF_below_Jpsi,LineColor(kRed),Precision(0.0001),Name("template1D_m2_below_Jpsi"));
  plotF_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2_below_Jpsi = new TCanvas("c_template1D_m2_below_Jpsi", "c_template1D_m2_below_Jpsi", 800, 800);
  c_template1D_m2_below_Jpsi->Clear();
  TPad *pad1b_below_Jpsi = new TPad("pad1b_below_Jpsi", "pad1b_below_Jpsi", 0, 0.3, 1, 1.0);
  pad1b_below_Jpsi->SetBottomMargin(0); pad1b_below_Jpsi->Draw(); pad1b_below_Jpsi->cd();
  plotF_below_Jpsi->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m2_below_Jpsi->cd(); c_template1D_m2_below_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b_below_Jpsi = new TPad("pad2b_below_Jpsi", "pad2b_below_Jpsi", 0, 0.0, 1, 0.29);
  pad2b_below_Jpsi->SetTopMargin(0); pad2b_below_Jpsi->SetBottomMargin(0.35); pad2b_below_Jpsi->SetGridy(); pad2b_below_Jpsi->Draw(); pad2b_below_Jpsi->cd();
  TH1F *h_dataFit2_below_Jpsi = new TH1F("h_dataFit2_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  h_dataFit2_below_Jpsi->SetLineColor(kBlack); h_dataFit2_below_Jpsi->SetMarkerStyle(20); h_dataFit2_below_Jpsi->SetMarkerColor(1); h_dataFit2_below_Jpsi->SetStats(0);
  h_dataFit2_below_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2_below_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit2_below_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit2_below_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit2_below_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit2_below_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit2_below_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2_below_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit2_below_Jpsi->GetYaxis()->CenterTitle(); h_dataFit2_below_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit2_below_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit2_below_Jpsi->GetYaxis()->SetTitleOffset(.9); h_dataFit2_below_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit2_below_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2_below_Jpsi = (TH1F*) ds_dimuorphan_bg_m2_below_Jpsi->createHistogram("hdata2_below_Jpsi", m2_below_Jpsi, Binning(m_bins_below_Jpsi, m_min, m_Jpsi_dn));
  TH1F* h_func2_below_Jpsi = new TH1F("h_func2_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  w->pdf("template1D_m2_below_Jpsi")->fillHistogram(h_func2_below_Jpsi, m2_below_Jpsi, hdata2_below_Jpsi->GetEntries());
  for(unsigned int iB=1; iB<m_bins_below_Jpsi; iB++){
    float ratio = h_func2_below_Jpsi->GetBinContent(iB)/hdata2_below_Jpsi->GetBinContent(iB);
    h_dataFit2_below_Jpsi->SetBinContent(iB, ratio );
    h_dataFit2_below_Jpsi->SetBinError(iB, ratio/sqrt(hdata2_below_Jpsi->GetBinContent(iB)) );
  }
  h_dataFit2_below_Jpsi->SetMinimum(0.5); h_dataFit2_below_Jpsi->SetMaximum(1.5); h_dataFit2_below_Jpsi->Sumw2(); h_dataFit2_below_Jpsi->SetStats(0); h_dataFit2_below_Jpsi->SetMarkerStyle(21); h_dataFit2_below_Jpsi->Draw("ep");
  c_template1D_m2_below_Jpsi->SaveAs("figures/template1D_m2_below_Jpsi.pdf");
  c_template1D_m2_below_Jpsi->SaveAs("figures/template1D_m2_below_Jpsi.png");
  c_template1D_m2_below_Jpsi->SaveAs("figures/template1D_m2_below_Jpsi.root");
  float chi2_F_below_Jpsi = plotF_below_Jpsi->chiSquare(17);
  cout<<"------------------ End m2 (below Jpsi) ---------------------"<<endl;

  //=========================
  //2017 PDF m2 (above J/psi)
  //=========================
  w->factory("Bernstein::bgF_above_Jpsi(m2_above_Jpsi, {bF06_above_Jpsi[3.5, 0, 15.], bF16_above_Jpsi[-3., -4., 3.], bF26_above_Jpsi[2.4, -5., 5.], bF36_above_Jpsi[0.1, -5., 5.], bF46_above_Jpsi[0.1, -3., 3.], bF56_above_Jpsi[0.5, -3., 3.], bF66_above_Jpsi[0.4, -3, 4.]})");
  w->factory("Gaussian::psiF_above_Jpsi(m2_above_Jpsi, 3.68609, psiF_sigma_above_Jpsi[0.031, 0.01, 0.1])");

  w->factory("SUM::template1D_m2_above_Jpsi(norm_bgF_above_Jpsi[1950., 500., 20000.]*bgF_above_Jpsi, norm_psiF_above_Jpsi[220., 0., 1000.]*psiF_above_Jpsi)");
  RooFitResult *rF_above_Jpsi = w->pdf("template1D_m2_above_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m2_above_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2 (above J/psi)---------------------"<<endl;
  rF_above_Jpsi->Print();

  RooPlot* plotF_above_Jpsi = w->var("m2_above_Jpsi")->frame(Title("BG template for orphan dimuon no high pT above Jpsi"),Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m2_above_Jpsi")->plotOn(plotF_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_above_Jpsi"));
  w->pdf("template1D_m2_above_Jpsi")->plotOn(plotF_above_Jpsi,LineColor(kRed),Precision(0.0001),Name("template1D_m2_above_Jpsi"));
  plotF_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2_above_Jpsi = new TCanvas("c_template1D_m2_above_Jpsi", "c_template1D_m2_above_Jpsi", 800, 800);
  c_template1D_m2_above_Jpsi->Clear();
  TPad *pad1b_above_Jpsi = new TPad("pad1b_above_Jpsi", "pad1b_above_Jpsi", 0, 0.3, 1, 1.0);
  pad1b_above_Jpsi->SetBottomMargin(0); pad1b_above_Jpsi->Draw(); pad1b_above_Jpsi->cd();
  plotF_above_Jpsi->Draw("same");
  txtHeader->Draw("same");
  c_template1D_m2_above_Jpsi->cd(); c_template1D_m2_above_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b_above_Jpsi = new TPad("pad2b_above_Jpsi", "pad2b_above_Jpsi", 0, 0.0, 1, 0.29);
  pad2b_above_Jpsi->SetTopMargin(0); pad2b_above_Jpsi->SetBottomMargin(0.35); pad2b_above_Jpsi->SetGridy(); pad2b_above_Jpsi->Draw(); pad2b_above_Jpsi->cd();
  TH1F *h_dataFit2_above_Jpsi = new TH1F("h_dataFit2_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  h_dataFit2_above_Jpsi->SetLineColor(kBlack); h_dataFit2_above_Jpsi->SetMarkerStyle(20); h_dataFit2_above_Jpsi->SetMarkerColor(1); h_dataFit2_above_Jpsi->SetStats(0);
  h_dataFit2_above_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2_above_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit2_above_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit2_above_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit2_above_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit2_above_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit2_above_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2_above_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit2_above_Jpsi->GetYaxis()->CenterTitle(); h_dataFit2_above_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit2_above_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit2_above_Jpsi->GetYaxis()->SetTitleOffset(.9); h_dataFit2_above_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit2_above_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2_above_Jpsi = (TH1F*) ds_dimuorphan_bg_m2_above_Jpsi->createHistogram("hdata2_above_Jpsi", m2_above_Jpsi, Binning(m_bins_above_Jpsi, m_Jpsi_up, m_max));
  TH1F* h_func2_above_Jpsi = new TH1F("h_func2_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  w->pdf("template1D_m2_above_Jpsi")->fillHistogram(h_func2_above_Jpsi, m2_above_Jpsi, hdata2_above_Jpsi->GetEntries());
  for(unsigned int iB=1; iB<m_bins_above_Jpsi; iB++){
    float ratio = h_func2_above_Jpsi->GetBinContent(iB)/hdata2_above_Jpsi->GetBinContent(iB);
    h_dataFit2_above_Jpsi->SetBinContent(iB, ratio );
    h_dataFit2_above_Jpsi->SetBinError(iB, ratio/sqrt(hdata2_above_Jpsi->GetBinContent(iB)) );
  }
  h_dataFit2_above_Jpsi->SetMinimum(0.5); h_dataFit2_above_Jpsi->SetMaximum(1.5); h_dataFit2_above_Jpsi->Sumw2(); h_dataFit2_above_Jpsi->SetStats(0); h_dataFit2_above_Jpsi->SetMarkerStyle(21); h_dataFit2_above_Jpsi->Draw("ep");
  c_template1D_m2_above_Jpsi->SaveAs("figures/template1D_m2_above_Jpsi.pdf");
  c_template1D_m2_above_Jpsi->SaveAs("figures/template1D_m2_above_Jpsi.png");
  c_template1D_m2_above_Jpsi->SaveAs("figures/template1D_m2_above_Jpsi.root");
  float chi2_F_above_Jpsi = plotF_above_Jpsi->chiSquare(10);
  cout<<"------------------ End m2 (above Jpsi) ---------------------"<<endl;

  //****************************************************************************
  //                     Create 2D template = m1 x m2
  //****************************************************************************
  cout << "-----Creating 2D template m1 * m2:-----" << endl;
  cout << "1D template m1 fit chi^2/dof: " << chi2_C << endl;
  cout << "1D template m2 fit chi^2/dof: " << chi2_F << endl;
  w->factory("PROD::template2D(template1D_m1, template1D_m2)");//Assume: two 1D p.d.f.s are not correlated

  cout << "-----Creating 2D templates m1 * m2 (exclude J/psi: below + above):-----" << endl;
  cout << "1D template m1 (below J/psi) fit chi^2/dof: " << chi2_C_below_Jpsi << endl;
  cout << "1D template m1 (above J/psi) fit chi^2/dof: " << chi2_C_above_Jpsi << endl;
  cout << "1D template m2 (below J/psi) fit chi^2/dof: " << chi2_F_below_Jpsi << endl;
  cout << "1D template m2 (above J/psi) fit chi^2/dof: " << chi2_F_above_Jpsi << endl;
  //These 2D pdf will be used to evaluate the contributions only below and only above Jpsi
  w->factory("PROD::template2D_below_Jpsi(template1D_m1_below_Jpsi, template1D_m2_below_Jpsi)");
  w->factory("PROD::template2D_above_Jpsi(template1D_m1_above_Jpsi, template1D_m2_above_Jpsi)");

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
  plot_Jpsi_m1->GetYaxis()->SetTitle("Events/0.04GeV");
  plot_Jpsi_m1->Draw(); txtHeader->Draw();
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1.pdf");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1.png");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1.root");

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
  plot_Jpsi_m2->GetYaxis()->SetTitle("Events/0.04GeV");
  plot_Jpsi_m2->Draw(); txtHeader->Draw();
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2.pdf");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2.png");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2.root");

  //****************************************************************************
  //                     Create 2D template (m1 x m2) for J/psi
  //****************************************************************************
  cout << "Create 2D template (m1 * m2) for J/psi" << endl;
  w->factory("PROD::Jpsi_2D(Jpsi_m1,Jpsi_m2)");

  //****************************************************************************
  //           For later use in LowMassBKGPlot2D.C: datasets of 2 dimu events
  //****************************************************************************
  cout << "Create trees on signal events" << endl;

  //To be used for scatter plot 2 dimu events @ CR later in LowMassBKGPlot2D.C
  ostringstream stream_cut_control_Iso_offDiagonal;
  stream_cut_control_Iso_offDiagonal << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) >= 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_min << " && massC < " << m_max << " && massF > " << m_min << " && massF < " << m_max;
  TString cut_control_Iso_offDiagonal = stream_cut_control_Iso_offDiagonal.str();
  std::cout << "2-dimu CR selctions (low mass): " << cut_control_Iso_offDiagonal.Data() << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);

  //Same as above but without J/psi
  ostringstream stream_cut_control_Iso_offDiagonal_no_Jpsi;
  stream_cut_control_Iso_offDiagonal_no_Jpsi << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) >= 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && ( (massC > " << m_min << " && massC < " << m_Jpsi_dn << ") || (massC > " << m_Jpsi_up << " && massC < "<< m_max << ") ) && ( (massF > " << m_min << " && massF < " << m_Jpsi_dn << ") || (massF > "<< m_Jpsi_up << " && massF < " << m_max <<") )";
  TString cut_control_Iso_offDiagonal_no_Jpsi = stream_cut_control_Iso_offDiagonal_no_Jpsi.str();
  std::cout << "2-dimu CR selctions (low mass no J/psi): " << cut_control_Iso_offDiagonal_no_Jpsi.Data() << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_no_Jpsi);

  //To be used to unblind the SR (scatter plot) later in LowMassBKGPlot2D.C (difference to above: mass cut <, not >)
  ostringstream stream_cut_signal;
  stream_cut_signal << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) < 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_min << " && massC < " << m_max << " && massF > " << m_min << " && massF < " << m_max;
  TString cut_signal = stream_cut_signal.str();
  std::cout << "2-dimu SR selctions (low mass): " << cut_signal.Data() << std::endl;
  TTree* tree_dimudimu_signal_2D = chain_data_dimudimu.CopyTree(cut_signal);

  //Same as above but without J/psi
  ostringstream stream_cut_signal_no_Jpsi;
  stream_cut_signal_no_Jpsi << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) < 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && ( (massC > " << m_min << " && massC < " << m_Jpsi_dn << ") || (massC > " << m_Jpsi_up << " && massC < "<< m_max << ") ) && ( (massF > " << m_min << " && massF < " << m_Jpsi_dn << ") || (massF > "<< m_Jpsi_up << " && massF < " << m_max <<") )";
  TString cut_signal_no_Jpsi = stream_cut_signal_no_Jpsi.str();
  std::cout << "2-dimu SR selctions (low mass no J/psi): " << cut_signal_no_Jpsi.Data() << std::endl;
  TTree* tree_dimudimu_signal_2D_no_Jpsi = chain_data_dimudimu.CopyTree(cut_signal_no_Jpsi);

  //Below Jpsi only
  ostringstream stream_cut_control_Iso_offDiagonal_below_Jpsi;
  stream_cut_control_Iso_offDiagonal_below_Jpsi << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) >= 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_min << " && massC < " << m_Jpsi_dn << " && massF > " << m_min << " && massF < " << m_Jpsi_dn;
  TString cut_control_Iso_offDiagonal_below_Jpsi = stream_cut_control_Iso_offDiagonal_below_Jpsi.str();
  std::cout << "2-dimu CR selctions (low mass below J/psi only): " << cut_control_Iso_offDiagonal_below_Jpsi.Data() << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_below_Jpsi);

  //Above Jpsi only
  ostringstream stream_cut_control_Iso_offDiagonal_above_Jpsi;
  stream_cut_control_Iso_offDiagonal_above_Jpsi << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && diMuonCMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonCMu0_IsoTk0p3_FittedVtx >= 0 && diMuonFMu0_IsoTk0p3_FittedVtx < " << iso_cut << " && diMuonFMu0_IsoTk0p3_FittedVtx >= 0 && TMath::Abs(massC-massF) >= 3*(0.003044 + 0.007025*(massC+massF)/2.0 + 0.000053*(massC+massF)*(massC+massF)/4.0) && massC > " << m_Jpsi_up << " && massC < " << m_max << " && massF > " << m_Jpsi_up << " && massF < " << m_max;
  TString cut_control_Iso_offDiagonal_above_Jpsi = stream_cut_control_Iso_offDiagonal_above_Jpsi.str();
  std::cout << "2-dimu CR selctions (low mass below J/psi only): " << cut_control_Iso_offDiagonal_above_Jpsi.Data() << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_above_Jpsi);

  //Setting Names
  tree_dimudimu_control_Iso_offDiagonal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_Iso_offDiagonal_2D->GetBranch("massF")->SetName("m2");
  tree_dimudimu_signal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_signal_2D->GetBranch("massF")->SetName("m2");
  //No J/psi version
  tree_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->GetBranch("massF")->SetName("m2");
  tree_dimudimu_signal_2D_no_Jpsi->GetBranch("massC")->SetName("m1");
  tree_dimudimu_signal_2D_no_Jpsi->GetBranch("massF")->SetName("m2");
  //Below Jpsi only
  tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->GetBranch("massC")->SetName("m1_below_Jpsi");
  tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->GetBranch("massF")->SetName("m2_below_Jpsi");
  //Above Jpsi only
  tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->GetBranch("massC")->SetName("m1_above_Jpsi");
  tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->GetBranch("massF")->SetName("m2_above_Jpsi");

  //Creating 2 dimu dataset
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D","ds_dimudimu_control_Iso_offDiagonal_2D", tree_dimudimu_control_Iso_offDiagonal_2D, RooArgSet(m1,m2));
  //RooDataSet* ds_dimudimu_control_offDiagonal_2D = new RooDataSet("ds_dimudimu_control_offDiagonal_2D","ds_dimudimu_control_offDiagonal_2D", tree_dimudimu_control_offDiagonal_2D, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_signal_2D = new RooDataSet("ds_dimudimu_signal_2D","ds_dimudimu_signal_2D", tree_dimudimu_signal_2D, RooArgSet(m1,m2));
  //No J/psi version
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi","ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi", tree_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi, RooArgSet(m1,m2));
  //RooDataSet* ds_dimudimu_control_offDiagonal_2D_no_Jpsi = new RooDataSet("ds_dimudimu_control_offDiagonal_2D_no_Jpsi","ds_dimudimu_control_offDiagonal_2D_no_Jpsi", tree_dimudimu_control_offDiagonal_2D_no_Jpsi, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_signal_2D_no_Jpsi = new RooDataSet("ds_dimudimu_signal_2D_no_Jpsi", "ds_dimudimu_signal_2D_no_Jpsi", tree_dimudimu_signal_2D_no_Jpsi, RooArgSet(m1,m2));
  //Below Jpsi only
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi","ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi", tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi, RooArgSet(m1_below_Jpsi,m2_below_Jpsi));
  //Above Jpsi only
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi","ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi", tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi, RooArgSet(m1_above_Jpsi,m2_above_Jpsi));

  ds_dimudimu_control_Iso_offDiagonal_2D->Print("s");
  ds_dimudimu_signal_2D->Print("s");
  //No J/psi version
  ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->Print("s");
  ds_dimudimu_signal_2D_no_Jpsi->Print("s");
  //Below Jpsi only
  ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->Print("s");
  //Above Jpsi only
  ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->Print("s");

  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D);
  w->import(*ds_dimudimu_signal_2D);
  //No J/psi version
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi);
  w->import(*ds_dimudimu_signal_2D_no_Jpsi);
  //Below Jpsi only
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi);
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi);

  //****************************************************************************
  //                           Save to Workspace
  //****************************************************************************
  cout<<"Save to workspace"<<endl;
  w->writeToFile(outFileLM);
}
