//=========================================================================
//= cmsenv                                                                =
//= Run it as: root -l -b -q LowMassBKGFit1D18.C                          =
//=          Wei Shi @Nov 20, 2019, Rice U.                               =
//=========================================================================
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
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "TAxis.h"

#include "macros/tdrStyle.C"
#include "Constants.h"
#include "Config.h"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void LowMassBKGFit1D18() {

  //Configure inputs for year
  BKG_cfg::ConfigureInput(year);
  setTDRStyle();

  TLegend *txtHeaderfit = new TLegend(0.14, 0.935, 0.98, 1);
  txtHeaderfit->SetFillColor(kWhite);
  txtHeaderfit->SetFillStyle(0);
  txtHeaderfit->SetBorderSize(0);
  txtHeaderfit->SetTextFont(42);
  txtHeaderfit->SetTextSize(0.045);
  txtHeaderfit->SetTextAlign(22);
  txtHeaderfit->SetHeader(headerfit);

  TLegend *txtHeader = new TLegend(0.14, 0.935, 0.98, 1);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader(headertpl);

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

  // Define RooRealVar in workspace
  RooRealVar m1("m1", "m_{#mu#mu_{1}}", m_min, m_max, "GeV");
  RooRealVar m2("m2", "m_{#mu#mu_{2}}", m_min, m_max, "GeV");
  m1.setBins(m_bins);
  m2.setBins(m_bins);

  // For below/above J/psi
  RooRealVar m1_below_Jpsi("m1_below_Jpsi", "m_{#mu#mu_{1}}", m_min, m_Jpsi_dn, "GeV");
  RooRealVar m1_above_Jpsi("m1_above_Jpsi", "m_{#mu#mu_{1}}", m_Jpsi_up, m_max, "GeV");
  RooRealVar m2_below_Jpsi("m2_below_Jpsi", "m_{#mu#mu_{2}}", m_min, m_Jpsi_dn, "GeV");
  RooRealVar m2_above_Jpsi("m2_above_Jpsi", "m_{#mu#mu_{2}}", m_Jpsi_up, m_max, "GeV");
  // For test: use with caution
  RooRealVar m1_above_Upsilon("m1_above_Upsilon", "m_{#mu#mu_{1}}", m_Upsilon_up, m_highmax, "GeV");
  RooRealVar m2_above_Upsilon("m2_above_Upsilon", "m_{#mu#mu_{2}}", m_Upsilon_up, m_highmax, "GeV");
  // all range
  RooRealVar m1_all("m1_all", "m_{#mu#mu_{1}}", m_min, m_highmax, "GeV");
  RooRealVar m2_all("m2_all", "m_{#mu#mu_{2}}", m_min, m_highmax, "GeV");

  m1_below_Jpsi.setBins(m_bins_below_Jpsi);
  m1_above_Jpsi.setBins(m_bins_above_Jpsi);
  m2_below_Jpsi.setBins(m_bins_below_Jpsi);
  m2_above_Jpsi.setBins(m_bins_above_Jpsi);
  // For test: use with caution
  m1_above_Upsilon.setBins(m_bins_above_Upsilon);
  m2_above_Upsilon.setBins(m_bins_above_Upsilon);
  // all range
  m1_all.setBins(m_bins_all);
  m2_all.setBins(m_bins_all);

  w->import(m1);
  w->import(m2);
  w->import(m1_below_Jpsi);
  w->import(m1_above_Jpsi);
  w->import(m2_below_Jpsi);
  w->import(m2_above_Jpsi);
  // For test: use with caution
  w->import(m1_above_Upsilon);
  w->import(m2_above_Upsilon);
  // all mass range
  w->import(m1_all);
  w->import(m2_all);

  //=*********************************************************************************
  //     Select events for constructing 1D templates, identify high pT muon, m1!=m2
  //=*********************************************************************************
  //m1: Orphan associated dimu both mu have high pT mu 24GeV and |eta|<2
  ostringstream stream_cut_bg_m1_iso;
  stream_cut_bg_m1_iso << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) && ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso = stream_cut_bg_m1_iso.str();
  //Print selections to check
  std::cout << "m1 selctions (low mass): " << cut_bg_m1_iso.Data() << std::endl;
  std::cout << "                         " << std::endl;
  TTree* tree_dimuorphan_bg_m1 = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso);

  //m1_below_Jpsi: Same as above except below J/psi
  ostringstream stream_cut_bg_m1_iso_below_Jpsi;
  stream_cut_bg_m1_iso_below_Jpsi << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) && ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_Jpsi_dn;
  TString cut_bg_m1_iso_below_Jpsi = stream_cut_bg_m1_iso_below_Jpsi.str();
  std::cout << "m1 selctions (low mass below J/psi): " << cut_bg_m1_iso_below_Jpsi.Data() << std::endl;
  std::cout << "                                     " << std::endl;
  TTree* tree_dimuorphan_bg_m1_below_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso_below_Jpsi);

  //m1_above_Jpsi: Same as above except above J/psi
  ostringstream stream_cut_bg_m1_iso_above_Jpsi;
  stream_cut_bg_m1_iso_above_Jpsi << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) && ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_Jpsi_up << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso_above_Jpsi = stream_cut_bg_m1_iso_above_Jpsi.str();
  std::cout << "m1 selctions (low mass above J/psi): " << cut_bg_m1_iso_above_Jpsi.Data() << std::endl;
  std::cout << "                                     " << std::endl;
  TTree* tree_dimuorphan_bg_m1_above_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso_above_Jpsi);

  //m2: Orphan mu has high pT mu 24GeV and |eta|<2
  ostringstream stream_cut_bg_m2_iso;
  stream_cut_bg_m2_iso << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m2_iso = stream_cut_bg_m2_iso.str();
  std::cout << "m2 selctions (low mass): " << cut_bg_m2_iso.Data() << std::endl;
  std::cout << "                         " << std::endl;
  TTree* tree_dimuorphan_bg_m2 = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso);

  //m2: Same as above except below J/psi
  ostringstream stream_cut_bg_m2_iso_below_Jpsi;
  stream_cut_bg_m2_iso_below_Jpsi << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_min << " && orph_dimu_mass < " << m_Jpsi_dn;
  TString cut_bg_m2_iso_below_Jpsi = stream_cut_bg_m2_iso_below_Jpsi.str();
  std::cout << "m2 selctions (low mass below J/psi): " << cut_bg_m2_iso_below_Jpsi.Data() << std::endl;
  std::cout << "                                     " << std::endl;
  TTree* tree_dimuorphan_bg_m2_below_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso_below_Jpsi);

  //m2: Same as above except above J/psi
  ostringstream stream_cut_bg_m2_iso_above_Jpsi;
  stream_cut_bg_m2_iso_above_Jpsi << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_Jpsi_up << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m2_iso_above_Jpsi = stream_cut_bg_m2_iso_above_Jpsi.str();
  std::cout << "m2 selctions (low mass above J/psi): " << cut_bg_m2_iso_above_Jpsi.Data() << std::endl;
  std::cout << "                                     " << std::endl;
  TTree* tree_dimuorphan_bg_m2_above_Jpsi = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso_above_Jpsi);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!   Note: some high mass cuts are hard to implement in this data-driven method with ONLY one dimuon,
  //!!!         so use this bkg template (above Upsilon) with caution! Result could be wrong!!!
  //!!!         For example, the DY veto cut can't be applied because there is no alternative pair from orphan-dimu
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //m1_above_Upsilon: Similar to above except above Upsilon + some high mass cuts like DY BKG veto and SA BKG veto)
  ostringstream stream_cut_bg_m1_iso_above_Upsilon;
  stream_cut_bg_m1_iso_above_Upsilon << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && ( orph_dimu_nSAMu == 0 || ( orph_dimu_nSAMu == 1 && orph_dimu_Lxy > 0.1 && orph_matched_segs_Mu0>=2 && orph_matched_segs_Mu1>=2 ) ) && (  ( orph_PtMu0 > " << pT_cut << " && TMath::Abs(orph_EtaMu0) < " << eta_cut << " ) && ( orph_PtMu1 > " << pT_cut << " && TMath::Abs(orph_EtaMu1) < " << eta_cut << " )  ) && orph_dimu_mass > " << m_Upsilon_up << " && orph_dimu_mass < " << m_highmax;
  TString cut_bg_m1_iso_above_Upsilon = stream_cut_bg_m1_iso_above_Upsilon.str();
  std::cout << "m1 selctions (Above Upsilon): " << cut_bg_m1_iso_above_Upsilon.Data() << std::endl;
  std::cout << "                                       " << std::endl;
  TTree* tree_dimuorphan_bg_m1_above_Upsilon = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso_above_Upsilon);

  //m2_above_Upsilon: Similar to above except above Upsilon + some high mass cuts like DY BKG veto and SA BKG veto)
  ostringstream stream_cut_bg_m2_iso_above_Upsilon;
  stream_cut_bg_m2_iso_above_Upsilon << "orph_dimu_isoTk < " << iso_cut << " && ( orph_dimu_Mu0_hitpix_Phase1 == 1 || orph_dimu_Mu1_hitpix_Phase1 == 1 ) && orph_isSignalHLTFired && orph_isVertexOK && orph_passOffLineSelPtEta && orph_dimu_nSAMu <= 1 && orph_dimu_prob > 0.2*(1 - orph_dimu_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(orph_dimu_dR)) + 109.83*pow(sqrt(orph_dimu_dR), 2) - 92.7445*pow(sqrt(orph_dimu_dR), 3) + 36.8351*pow(sqrt(orph_dimu_dR), 4) )*pow(fabs(orph_dimu_Lxy/10.0), 2.0) ) && ( orph_dimu_nSAMu == 0 || ( orph_dimu_nSAMu == 1 && orph_dimu_Lxy > 0.1 && orph_matched_segs_Mu0>=2 && orph_matched_segs_Mu1>=2 ) ) && orph_PtOrph > " << pT_cut << " && TMath::Abs(orph_EtaOrph) < " << eta_cut << " && orph_dimu_mass > " << m_Upsilon_up << " && orph_dimu_mass < " << m_highmax;
  TString cut_bg_m2_iso_above_Upsilon = stream_cut_bg_m2_iso_above_Upsilon.str();
  std::cout << "m2 selctions (Above Upsilon): " << cut_bg_m2_iso_above_Upsilon.Data() << std::endl;
  std::cout << "                                       " << std::endl;
  TTree* tree_dimuorphan_bg_m2_above_Upsilon = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso_above_Upsilon);

  //Setting Names: Control region
  tree_dimuorphan_bg_m1->GetBranch("orph_dimu_mass")->SetName("m1");
  tree_dimuorphan_bg_m2->GetBranch("orph_dimu_mass")->SetName("m2");
  tree_dimuorphan_bg_m1_below_Jpsi->GetBranch("orph_dimu_mass")->SetName("m1_below_Jpsi");
  tree_dimuorphan_bg_m1_above_Jpsi->GetBranch("orph_dimu_mass")->SetName("m1_above_Jpsi");
  tree_dimuorphan_bg_m2_below_Jpsi->GetBranch("orph_dimu_mass")->SetName("m2_below_Jpsi");
  tree_dimuorphan_bg_m2_above_Jpsi->GetBranch("orph_dimu_mass")->SetName("m2_above_Jpsi");
  //For test: use with caution
  tree_dimuorphan_bg_m1_above_Upsilon->GetBranch("orph_dimu_mass")->SetName("m1_above_Upsilon");
  tree_dimuorphan_bg_m2_above_Upsilon->GetBranch("orph_dimu_mass")->SetName("m2_above_Upsilon");

  //Constructor of a data set from (part of) an ROOT TTRee:
  //The dimensions of the data set are defined by the 'vars' RooArgSet. For each dimension specified, the TTree must have a branch with the same name.
  RooDataSet* ds_dimuorphan_bg_m1 = new RooDataSet("ds_dimuorphan_bg_m1", "ds_dimuorphan_bg_m1", tree_dimuorphan_bg_m1, RooArgSet(m1));
  RooDataSet* ds_dimuorphan_bg_m2 = new RooDataSet("ds_dimuorphan_bg_m2", "ds_dimuorphan_bg_m2", tree_dimuorphan_bg_m2, RooArgSet(m2));
  RooDataSet* ds_dimuorphan_bg_m1_below_Jpsi = new RooDataSet("ds_dimuorphan_bg_m1_below_Jpsi", "ds_dimuorphan_bg_m1_below_Jpsi", tree_dimuorphan_bg_m1_below_Jpsi, RooArgSet(m1_below_Jpsi));
  RooDataSet* ds_dimuorphan_bg_m1_above_Jpsi = new RooDataSet("ds_dimuorphan_bg_m1_above_Jpsi", "ds_dimuorphan_bg_m1_above_Jpsi", tree_dimuorphan_bg_m1_above_Jpsi, RooArgSet(m1_above_Jpsi));
  RooDataSet* ds_dimuorphan_bg_m2_below_Jpsi = new RooDataSet("ds_dimuorphan_bg_m2_below_Jpsi", "ds_dimuorphan_bg_m2_below_Jpsi", tree_dimuorphan_bg_m2_below_Jpsi, RooArgSet(m2_below_Jpsi));
  RooDataSet* ds_dimuorphan_bg_m2_above_Jpsi = new RooDataSet("ds_dimuorphan_bg_m2_above_Jpsi", "ds_dimuorphan_bg_m2_above_Jpsi", tree_dimuorphan_bg_m2_above_Jpsi, RooArgSet(m2_above_Jpsi));
  //For test: use with caution
  RooDataSet* ds_dimuorphan_bg_m1_above_Upsilon = new RooDataSet("ds_dimuorphan_bg_m1_above_Upsilon", "ds_dimuorphan_bg_m1_above_Upsilon", tree_dimuorphan_bg_m1_above_Upsilon, RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimuorphan_bg_m2_above_Upsilon = new RooDataSet("ds_dimuorphan_bg_m2_above_Upsilon", "ds_dimuorphan_bg_m2_above_Upsilon", tree_dimuorphan_bg_m2_above_Upsilon, RooArgSet(m2_above_Upsilon));

  cout<<"-----Now Printing and importing the Datasets:-----"<<endl;
  ds_dimuorphan_bg_m1->Print("s");
  ds_dimuorphan_bg_m2->Print("s");
  ds_dimuorphan_bg_m1_below_Jpsi->Print("s");
  ds_dimuorphan_bg_m1_above_Jpsi->Print("s");
  ds_dimuorphan_bg_m2_below_Jpsi->Print("s");
  ds_dimuorphan_bg_m2_above_Jpsi->Print("s");
  //For test: use with caution
  ds_dimuorphan_bg_m1_above_Upsilon->Print("s");
  ds_dimuorphan_bg_m2_above_Upsilon->Print("s");

  w->import(*ds_dimuorphan_bg_m1);
  w->import(*ds_dimuorphan_bg_m2);
  w->import(*ds_dimuorphan_bg_m1_below_Jpsi);
  w->import(*ds_dimuorphan_bg_m2_below_Jpsi);
  w->import(*ds_dimuorphan_bg_m1_above_Jpsi);
  w->import(*ds_dimuorphan_bg_m2_above_Jpsi);
  //For test: use with caution
  w->import(*ds_dimuorphan_bg_m1_above_Upsilon);
  w->import(*ds_dimuorphan_bg_m2_above_Upsilon);

  //Draw before fiting
  RooPlot* plotC1 = w->var("m1")->frame(Title("m1 data tempalate NO FIT"), Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC1, DataError(RooAbsData::SumW2), Name("data_m1"));
  plotC1->GetYaxis()->SetTitle("Events/0.04GeV");
  TCanvas * c1 = new TCanvas("c1");
  c1->cd(); plotC1->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1.pdf");
  c1->SaveAs("figures/dimuorphan_m1.png");
  c1->SaveAs("figures/dimuorphan_m1.root");

  RooPlot* plotC2 = w->var("m2")->frame(Title("m2 data tempalate NO FIT"), Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotC2, DataError(RooAbsData::SumW2), Name("data_m2"));
  plotC2->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC2->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2.pdf");
  c1->SaveAs("figures/dimuorphan_m2.png");
  c1->SaveAs("figures/dimuorphan_m2.root");

  RooPlot* plotC1_below_Jpsi = w->var("m1_below_Jpsi")->frame(Title("m1 data tempalate NO FIT below Jpsi"), Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m1_below_Jpsi")->plotOn(plotC1_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_below_Jpsi"));
  plotC1_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC1_below_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1_below_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m1_below_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m1_below_Jpsi.root");

  RooPlot* plotC1_above_Jpsi = w->var("m1_above_Jpsi")->frame(Title("m1 data tempalate NO FIT above Jpsi"), Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m1_above_Jpsi")->plotOn(plotC1_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_above_Jpsi"));
  plotC1_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC1_above_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1_above_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m1_above_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m1_above_Jpsi.root");

  RooPlot* plotC2_below_Jpsi = w->var("m2_below_Jpsi")->frame(Title("m2 data tempalate NO FIT below Jpsi"), Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m2_below_Jpsi")->plotOn(plotC2_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_below_Jpsi"));
  plotC2_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC2_below_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2_below_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m2_below_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m2_below_Jpsi.root");

  RooPlot* plotC2_above_Jpsi = w->var("m2_above_Jpsi")->frame(Title("m2 data tempalate NO FIT above Jpsi"), Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m2_above_Jpsi")->plotOn(plotC2_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_above_Jpsi"));
  plotC2_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV");
  plotC2_above_Jpsi->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2_above_Jpsi.pdf");
  c1->SaveAs("figures/dimuorphan_m2_above_Jpsi.png");
  c1->SaveAs("figures/dimuorphan_m2_above_Jpsi.root");

  //For test: use with caution
  RooPlot* plotC1_above_Upsilon = w->var("m1_above_Upsilon")->frame(Title("m1 data tempalate NO FIT above Upsilon"), Bins(m_bins_above_Upsilon));
  w->data("ds_dimuorphan_bg_m1_above_Upsilon")->plotOn(plotC1_above_Upsilon, DataError(RooAbsData::SumW2), Name("data_m1_above_Upsilon"));
  plotC1_above_Upsilon->GetYaxis()->SetTitle("Events/0.5GeV");
  plotC1_above_Upsilon->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m1_above_Upsilon.pdf");
  c1->SaveAs("figures/dimuorphan_m1_above_Upsilon.png");
  c1->SaveAs("figures/dimuorphan_m1_above_Upsilon.root");
  //For test: use with caution
  RooPlot* plotC2_above_Upsilon = w->var("m2_above_Upsilon")->frame(Title("m2 data tempalate NO FIT above Upsilon"), Bins(m_bins_above_Upsilon));
  w->data("ds_dimuorphan_bg_m2_above_Upsilon")->plotOn(plotC2_above_Upsilon, DataError(RooAbsData::SumW2), Name("data_m2_above_Upsilon"));
  plotC2_above_Upsilon->GetYaxis()->SetTitle("Events/0.5GeV");
  plotC2_above_Upsilon->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/dimuorphan_m2_above_Upsilon.pdf");
  c1->SaveAs("figures/dimuorphan_m2_above_Upsilon.png");
  c1->SaveAs("figures/dimuorphan_m2_above_Upsilon.root");

  delete c1;

  //=***************************************************************************
  //                         Create template for m1
  //=***************************************************************************
  cout<<"-----Creating template for m1:-----"<<endl;
  //===========
  //2018 PDF m1
  //===========
  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.01, -1.0, 0.3], MmumuC_p[0.05, -1.0, 2.5])");
  w->factory("Bernstein::bgC(m1, {bC06[4.0, 0.1, 25.], bC16[-2.3, -5., 3.], bC26[-0.01, -5., 3.], bC36[1.2, -3., 5.], bC46[0.54, -3., 3.], bC56[0.9, -3., 3.], bC66[0.7, -3., 4.]})");
  w->factory("Gaussian::adHocC(m1, adHocC_mass[0.3, 0., 0.6], adHocC_sigma[0.1, 0., 0.5])");
  w->factory("Gaussian::etaC(m1, 0.54786, 0.007)");
  w->factory("Gaussian::rhoC(m1, 0.78265, 0.009)");
  w->factory("Gaussian::phiC(m1, 1.01946, 0.01)");
  w->factory("CBShape::JpsiC(m1, JpsiC_mean[3.0969, 3.0, 3.35], JpsiC_sigma[0.1, 0., 0.3], JpsiC_alpha[1.5, 0.4, 7.0], JpsiC_n[2.0])");
  w->factory("Gaussian::psiC(m1, 3.68609, psiC_sigma[0.031, 0.01, 0.04])");
  w->factory("SUM::template1D_m1(norm_adHocC[1500., 0., 10000.]*adHocC, norm_MmumuC[11000., 0., 35000.]*MmumuC, norm_bgC[10000., 1000., 30000.]*bgC, norm_etaC[0.01, 0., 1000.]*etaC, norm_rhoC[150, 0., 1000.]*rhoC, norm_phiC[110, 0., 1000.]*phiC, norm_JpsiC[5000., 10., 20000.]*JpsiC, norm_psiC[200., 0., 1000.]*psiC)");

  RooFitResult *rC = w->pdf("template1D_m1")->fitTo(*(w->data("ds_dimuorphan_bg_m1")), Extended(1), Save(), SumW2Error(kTRUE));//this is fit to unbinned data (data is not binned into a histogram), extended means extended parameter: number of events
  cout<<"------------------RooFitResult for m1---------------------"<<endl;
  rC->Print();

  RooPlot* plotC = w->var("m1")->frame(Title("1D template for orphan dimuon high pT"), Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC, DataError(RooAbsData::SumW2), Name("data_m1"));//data is binned now: RooHist
  w->pdf("template1D_m1")->plotOn(plotC, LineColor(kRed), Precision(0.0001), Name("template1D_m1"));//template: RooCurve
  plotC->GetYaxis()->SetTitle("Events/0.04GeV"); plotC->GetYaxis()->SetTitleOffset(1.1);

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1 = new TCanvas("c_template1D_m1", "c_template1D_m1", 800, 800);
  c_template1D_m1->Clear();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); pad1->Draw(); pad1->cd();
  plotC->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m1->cd(); c_template1D_m1->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2->SetTopMargin(0); pad2->SetBottomMargin(0.35); pad2->SetGridy(); pad2->Draw(); pad2->cd();
  TH1F *h_dataFit1 = new TH1F("h_dataFit1","", m_bins, m_min, m_max);
  h_dataFit1->SetLineColor(kBlack); h_dataFit1->SetMarkerStyle(20); h_dataFit1->SetMarkerColor(1); h_dataFit1->SetStats(0);
  h_dataFit1->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1->GetXaxis()->SetTitleSize(20); h_dataFit1->GetXaxis()->SetTitleFont(43); h_dataFit1->GetXaxis()->SetTitleOffset(3.0); h_dataFit1->GetXaxis()->SetLabelSize(15); h_dataFit1->GetXaxis()->SetLabelFont(43);
  h_dataFit1->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1->GetYaxis()->SetNdivisions(505); h_dataFit1->GetYaxis()->CenterTitle(); h_dataFit1->GetYaxis()->SetTitleSize(20); h_dataFit1->GetYaxis()->SetTitleFont(43); h_dataFit1->GetYaxis()->SetTitleOffset(1.2); h_dataFit1->GetYaxis()->SetLabelSize(15); h_dataFit1->GetYaxis()->SetLabelFont(43);
  TH1F *hdata = (TH1F*) ds_dimuorphan_bg_m1->createHistogram("hdata", m1, Binning(m_bins, m_min, m_max));
  TH1F* h_func = new TH1F("h_func","", m_bins, m_min, m_max);
  w->pdf("template1D_m1")->fillHistogram(h_func, m1, hdata->GetEntries());
  for (unsigned int iB=1; iB < m_bins; iB++){
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
  //2018 PDF m1 (below J/psi)
  //=========================
  w->factory("EXPR::MmumuC_below_Jpsi('m1_below_Jpsi*pow( (m1_below_Jpsi/m_below_Jpsi)*(m1_below_Jpsi/m_below_Jpsi) - 1.0, MmumuC_p_below_Jpsi )*exp( -MmumuC_c_below_Jpsi*( (m1_below_Jpsi/m_below_Jpsi)*(m1_below_Jpsi/m_below_Jpsi) - 1.0 ) )', m1_below_Jpsi, m_below_Jpsi[0.2113], MmumuC_c_below_Jpsi[0.1, -0.3, 0.3], MmumuC_p_below_Jpsi[0.1, -1.0, 1.5])");
  w->factory("Bernstein::bgC_below_Jpsi(m1_below_Jpsi, {bC06_below_Jpsi[1., 0.1, 15.], bC16_below_Jpsi[1.5, 0., 5.], bC26_below_Jpsi[0.1, -10., 3.], bC36_below_Jpsi[-0.1, -3., 5.], bC46_below_Jpsi[1.5, -3., 3.], bC56_below_Jpsi[0.5, 0., 3.], bC66_below_Jpsi[0.1, -1., 4.]})");
  w->factory("Gaussian::adHocC_below_Jpsi(m1_below_Jpsi, adHocC_mass_below_Jpsi[0.3, 0., 0.6], adHocC_sigma_below_Jpsi[0.07, 0., 0.3])");
  w->factory("Gaussian::etaC_below_Jpsi(m1_below_Jpsi, 0.54786, 0.007)");
  w->factory("Gaussian::rhoC_below_Jpsi(m1_below_Jpsi, 0.78265, 0.009)");
  w->factory("Gaussian::phiC_below_Jpsi(m1_below_Jpsi, 1.01946, 0.01)");
  w->factory("SUM::template1D_m1_below_Jpsi(norm_adHocC_below_Jpsi[2000., 0., 10000.]*adHocC_below_Jpsi, norm_MmumuC_below_Jpsi[7000., 0., 35000.]*MmumuC_below_Jpsi, norm_bgC_below_Jpsi[5000., 1000., 25000.]*bgC_below_Jpsi, norm_etaC_below_Jpsi[10, 0., 1000.]*etaC_below_Jpsi, norm_rhoC_below_Jpsi[133., 0., 1000.]*rhoC_below_Jpsi, norm_phiC_below_Jpsi[100., 0., 1000.]*phiC_below_Jpsi)");

  RooFitResult *rC_below_Jpsi = w->pdf("template1D_m1_below_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m1_below_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1 (below J/psi)---------------------"<<endl;
  rC_below_Jpsi->Print();

  RooPlot* plotC_below_Jpsi = w->var("m1_below_Jpsi")->frame(Title("1D template for orphan dimuon high pT below Jpsi"), Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m1_below_Jpsi")->plotOn(plotC_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_below_Jpsi"));
  w->pdf("template1D_m1_below_Jpsi")->plotOn(plotC_below_Jpsi, LineColor(kRed), Precision(0.0001), Name("template1D_m1_below_Jpsi")); // By default only fitted range is shown
  plotC_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV"); plotC_below_Jpsi->GetYaxis()->SetTitleOffset(1.1);

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1_below_Jpsi = new TCanvas("c_template1D_m1_below_Jpsi", "c_template1D_m1_below_Jpsi", 800, 800);
  c_template1D_m1_below_Jpsi->Clear();
  TPad *pad1_below_Jpsi = new TPad("pad1_below_Jpsi", "pad1_below_Jpsi", 0, 0.3, 1, 1.0);
  pad1_below_Jpsi->SetBottomMargin(0); pad1_below_Jpsi->Draw(); pad1_below_Jpsi->cd();
  plotC_below_Jpsi->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m1_below_Jpsi->cd(); c_template1D_m1_below_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2_below_Jpsi = new TPad("pad2_below_Jpsi", "pad2_below_Jpsi", 0, 0.0, 1, 0.29);
  pad2_below_Jpsi->SetTopMargin(0); pad2_below_Jpsi->SetBottomMargin(0.35); pad2_below_Jpsi->SetGridy(); pad2_below_Jpsi->Draw(); pad2_below_Jpsi->cd();
  TH1F *h_dataFit1_below_Jpsi = new TH1F("h_dataFit1_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  h_dataFit1_below_Jpsi->SetLineColor(kBlack); h_dataFit1_below_Jpsi->SetMarkerStyle(20); h_dataFit1_below_Jpsi->SetMarkerColor(1); h_dataFit1_below_Jpsi->SetStats(0);
  h_dataFit1_below_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1_below_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit1_below_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit1_below_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit1_below_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit1_below_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit1_below_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1_below_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit1_below_Jpsi->GetYaxis()->CenterTitle(); h_dataFit1_below_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit1_below_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit1_below_Jpsi->GetYaxis()->SetTitleOffset(1.2); h_dataFit1_below_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit1_below_Jpsi->GetYaxis()->SetLabelFont(43);
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
  //2018 PDF m1 (above J/psi)
  //=========================
  w->factory("Bernstein::bgC_above_Jpsi(m1_above_Jpsi, {bC06_above_Jpsi[2.0, 0.1, 15.], bC16_above_Jpsi[0.1., -5., 3.], bC26_above_Jpsi[2.5, -3., 5.], bC36_above_Jpsi[0.1, -3., 5.], bC46_above_Jpsi[-0.1, -3., 5.], bC56_above_Jpsi[0.3, 0., 3.], bC66_above_Jpsi[1, 0.1, 4.]})");
  w->factory("Gaussian::psiC_above_Jpsi(m1_above_Jpsi, 3.68609, psiC_sigma_above_Jpsi[0.03, 0.01, 0.05])");
  w->factory("SUM::template1D_m1_above_Jpsi(norm_bgC_above_Jpsi[3500., 1000., 20000.]*bgC_above_Jpsi, norm_psiC_above_Jpsi[300., 0., 1000.]*psiC_above_Jpsi)");

  RooFitResult *rC_above_Jpsi = w->pdf("template1D_m1_above_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m1_above_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1 (above J/psi)---------------------"<<endl;
  rC_above_Jpsi->Print();

  RooPlot* plotC_above_Jpsi = w->var("m1_above_Jpsi")->frame(Title("1D template for orphan dimuon high pT above Jpsi"), Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m1_above_Jpsi")->plotOn(plotC_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m1_above_Jpsi"));
  w->pdf("template1D_m1_above_Jpsi")->plotOn(plotC_above_Jpsi, LineColor(kRed), Precision(0.0001), Name("template1D_m1_above_Jpsi"));
  plotC_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV"); plotC_above_Jpsi->GetYaxis()->SetTitleOffset(1.1);

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1_above_Jpsi = new TCanvas("c_template1D_m1_above_Jpsi", "c_template1D_m1_above_Jpsi", 800, 800);
  c_template1D_m1_above_Jpsi->Clear();
  TPad *pad1_above_Jpsi = new TPad("pad1_above_Jpsi", "pad1_above_Jpsi", 0, 0.3, 1, 1.0);
  pad1_above_Jpsi->SetBottomMargin(0); pad1_above_Jpsi->Draw(); pad1_above_Jpsi->cd();
  plotC_above_Jpsi->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m1_above_Jpsi->cd(); c_template1D_m1_above_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2_above_Jpsi = new TPad("pad2_above_Jpsi", "pad2_above_Jpsi", 0, 0.0, 1, 0.29);
  pad2_above_Jpsi->SetTopMargin(0); pad2_above_Jpsi->SetBottomMargin(0.35); pad2_above_Jpsi->SetGridy(); pad2_above_Jpsi->Draw(); pad2_above_Jpsi->cd();
  TH1F *h_dataFit1_above_Jpsi = new TH1F("h_dataFit1_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  h_dataFit1_above_Jpsi->SetLineColor(kBlack); h_dataFit1_above_Jpsi->SetMarkerStyle(20); h_dataFit1_above_Jpsi->SetMarkerColor(1); h_dataFit1_above_Jpsi->SetStats(0);
  h_dataFit1_above_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1_above_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit1_above_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit1_above_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit1_above_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit1_above_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit1_above_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1_above_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit1_above_Jpsi->GetYaxis()->CenterTitle(); h_dataFit1_above_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit1_above_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit1_above_Jpsi->GetYaxis()->SetTitleOffset(1.2); h_dataFit1_above_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit1_above_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata_above_Jpsi = (TH1F*) ds_dimuorphan_bg_m1_above_Jpsi->createHistogram("hdata_above_Jpsi", m1_above_Jpsi, Binning(m_bins_above_Jpsi, m_Jpsi_up, m_max));
  TH1F* h_func_above_Jpsi = new TH1F("h_func_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  w->pdf("template1D_m1_above_Jpsi")->fillHistogram(h_func_above_Jpsi, m1_above_Jpsi, hdata_above_Jpsi->GetEntries());
  for (unsigned int iB=1; iB < m_bins_above_Jpsi; iB++) {
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

  //===========================
  //2018 PDF m1 (above Upsilon)
  //===========================
  w->factory("Bernstein::bgC_above_Upsilon(m1_above_Upsilon, {bC06_above_Upsilon[2.0, 0.1, 15.], bC16_above_Upsilon[0.1., -5., 3.], bC26_above_Upsilon[2.5, -3., 5.], bC36_above_Upsilon[0.1, -3., 5.], bC46_above_Upsilon[-0.1, -3., 5.], bC56_above_Upsilon[0.3, 0., 3.], bC66_above_Upsilon[1, 0.1, 4.]})");
  w->factory("SUM::template1D_m1_above_Upsilon(norm_bgC_above_Upsilon[20000., 1000., 30000.]*bgC_above_Upsilon)");

  RooFitResult *rC_above_Upsilon = w->pdf("template1D_m1_above_Upsilon")->fitTo(*(w->data("ds_dimuorphan_bg_m1_above_Upsilon")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m1 (above Upsilon)---------------------"<<endl;
  rC_above_Upsilon->Print();

  RooPlot* plotC_above_Upsilon = w->var("m1_above_Upsilon")->frame(Title("1D template for orphan dimuon high pT above Upsilon"), Bins(m_bins_above_Upsilon));
  w->data("ds_dimuorphan_bg_m1_above_Upsilon")->plotOn(plotC_above_Upsilon, DataError(RooAbsData::SumW2), Name("data_m1_above_Upsilon"));
  w->pdf("template1D_m1_above_Upsilon")->plotOn(plotC_above_Upsilon, LineColor(kRed), Precision(0.0001), Name("template1D_m1_above_Upsilon"));
  plotC_above_Upsilon->GetYaxis()->SetTitle("Events/0.5GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m1_above_Upsilon = new TCanvas("c_template1D_m1_above_Upsilon", "c_template1D_m1_above_Upsilon", 800, 800);
  c_template1D_m1_above_Upsilon->Clear();
  TPad *pad1_above_Upsilon = new TPad("pad1_above_Upsilon", "pad1_above_Upsilon", 0, 0.3, 1, 1.0);
  pad1_above_Upsilon->SetBottomMargin(0); pad1_above_Upsilon->Draw(); pad1_above_Upsilon->cd();
  plotC_above_Upsilon->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m1_above_Upsilon->cd(); c_template1D_m1_above_Upsilon->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2_above_Upsilon = new TPad("pad2_above_Upsilon", "pad2_above_Upsilon", 0, 0.0, 1, 0.29);
  pad2_above_Upsilon->SetTopMargin(0); pad2_above_Upsilon->SetBottomMargin(0.35); pad2_above_Upsilon->SetGridy(); pad2_above_Upsilon->Draw(); pad2_above_Upsilon->cd();
  TH1F *h_dataFit1_above_Upsilon = new TH1F("h_dataFit1_above_Upsilon","", m_bins_above_Upsilon, m_Upsilon_up, m_highmax);
  h_dataFit1_above_Upsilon->SetLineColor(kBlack); h_dataFit1_above_Upsilon->SetMarkerStyle(20); h_dataFit1_above_Upsilon->SetMarkerColor(1); h_dataFit1_above_Upsilon->SetStats(0);
  h_dataFit1_above_Upsilon->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit1_above_Upsilon->GetXaxis()->SetTitleSize(20); h_dataFit1_above_Upsilon->GetXaxis()->SetTitleFont(43); h_dataFit1_above_Upsilon->GetXaxis()->SetTitleOffset(3.0); h_dataFit1_above_Upsilon->GetXaxis()->SetLabelSize(15); h_dataFit1_above_Upsilon->GetXaxis()->SetLabelFont(43);
  h_dataFit1_above_Upsilon->GetYaxis()->SetTitle("Fit/Data"); h_dataFit1_above_Upsilon->GetYaxis()->SetNdivisions(505); h_dataFit1_above_Upsilon->GetYaxis()->CenterTitle(); h_dataFit1_above_Upsilon->GetYaxis()->SetTitleSize(20); h_dataFit1_above_Upsilon->GetYaxis()->SetTitleFont(43); h_dataFit1_above_Upsilon->GetYaxis()->SetTitleOffset(1.2); h_dataFit1_above_Upsilon->GetYaxis()->SetLabelSize(15); h_dataFit1_above_Upsilon->GetYaxis()->SetLabelFont(43);
  TH1F *hdata_above_Upsilon = (TH1F*) ds_dimuorphan_bg_m1_above_Upsilon->createHistogram("hdata_above_Upsilon", m1_above_Upsilon, Binning(m_bins_above_Upsilon, m_Upsilon_up, m_highmax));
  TH1F* h_func_above_Upsilon = new TH1F("h_func_above_Upsilon", "", m_bins_above_Upsilon, m_Upsilon_up, m_highmax);
  w->pdf("template1D_m1_above_Upsilon")->fillHistogram(h_func_above_Upsilon, m1_above_Upsilon, hdata_above_Upsilon->GetEntries());
  for (unsigned int iB=1; iB < m_bins_above_Upsilon; iB++) {
    float ratio = h_func_above_Upsilon->GetBinContent(iB)/hdata_above_Upsilon->GetBinContent(iB);
    h_dataFit1_above_Upsilon->SetBinContent(iB, ratio );
    h_dataFit1_above_Upsilon->SetBinError(iB, ratio/sqrt(hdata_above_Upsilon->GetBinContent(iB)) );
  }
  h_dataFit1_above_Upsilon->SetMinimum(0.5); h_dataFit1_above_Upsilon->SetMaximum(1.5); h_dataFit1_above_Upsilon->Sumw2(); h_dataFit1_above_Upsilon->SetStats(0); h_dataFit1_above_Upsilon->SetMarkerStyle(21); h_dataFit1_above_Upsilon->Draw("ep");
  c_template1D_m1_above_Upsilon->SaveAs("figures/template1D_m1_above_Upsilon.pdf");
  c_template1D_m1_above_Upsilon->SaveAs("figures/template1D_m1_above_Upsilon.png");
  c_template1D_m1_above_Upsilon->SaveAs("figures/template1D_m1_above_Upsilon.root");
  float chi2_C_above_Upsilon = plotC_above_Upsilon->chiSquare(10);
  cout<<"------------------ End m1 (above Upsilon) ---------------------"<<endl;

  //=***************************************************************************
  //                         Create template for m2
  //=***************************************************************************
  cout<<"-----Creating template for m2:-----"<<endl;
  //===========
  //2018 PDF m2
  //===========
  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[0.01, -0.3, 0.3], MmumuF_p[0.05, -1.0, 4.])");
  w->factory("Bernstein::bgF(m2, {bF06[3, 0, 30.], bF16[0.9, -6., 5.], bF26[0.1, -5., 5.], bF36[0.5, -3., 10.], bF46[0.1, -3., 5.], bF56[0.1, -3., 3.], bF66[0.5, 0.1, 4.]})");
  w->factory("Gaussian::adHocF(m2, adHocF_mass[0.4, 0.1, 0.6], adHocF_sigma[0.01, 0.001, 0.3])");
  w->factory("Gaussian::etaF(m2, 0.54786, 0.007)");
  w->factory("Gaussian::rhoF(m2, 0.78265, 0.009)");
  w->factory("Gaussian::phiF(m2, 1.01946, 0.01)");
  w->factory("CBShape::JpsiF(m2, JpsiF_mean[3.0969, 3.0, 3.35], JpsiF_sigma[0.1, 0., 0.3], JpsiF_alpha[1.2, 0.4, 10.0], JpsiF_n[2.0])");
  w->factory("Gaussian::psiF(m2, 3.68609, psiF_sigma[0.031, 0.01, 0.05])");
  w->factory("SUM::template1D_m2(norm_adHocF[500., 0., 2000.]*adHocF, norm_MmumuF[10000., 0., 35000.]*MmumuF, norm_bgF[20000., 1000., 50000.]*bgF, norm_etaF[3., 0., 10.]*etaF, norm_rhoF[400., 1., 1000.]*rhoF, norm_phiF[110., 1., 1000.]*phiF, norm_JpsiF[6400., 0., 20000.]*JpsiF, norm_psiF[250., 0., 1000.]*psiF)");

  RooFitResult *rF = w->pdf("template1D_m2")->fitTo(*(w->data("ds_dimuorphan_bg_m2")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2---------------------"<<endl;
  rF->Print();

  RooPlot* plotF = w->var("m2")->frame(Title("BG template for orphan dimuon no high pT"), Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotF, DataError(RooAbsData::SumW2), Name("data_m2"));
  w->pdf("template1D_m2")->plotOn(plotF,LineColor(kRed),Precision(0.0001),Name("template1D_m2"));
  plotF->GetYaxis()->SetTitle("Events/0.04GeV"); plotF->GetYaxis()->SetTitleOffset(1.1);

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2 = new TCanvas("c_template1D_m2", "c_template1D_m2", 800, 800);
  c_template1D_m2->Clear();
  TPad *pad1b = new TPad("pad1b", "pad1b", 0, 0.3, 1, 1.0);
  pad1b->SetBottomMargin(0); pad1b->Draw(); pad1b->cd();
  plotF->Draw("same"); txtHeaderfit->Draw("same");
  c_template1D_m2->cd(); c_template1D_m2->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b = new TPad("pad2", "pad2", 0, 0.0, 1, 0.29);
  pad2b->SetTopMargin(0); pad2b->SetBottomMargin(0.35); pad2b->SetGridy(); pad2b->Draw(); pad2b->cd();
  TH1F *h_dataFit2 = new TH1F("h_dataFit2","", m_bins, m_min, m_max);
  h_dataFit2->SetLineColor(kBlack); h_dataFit2->SetMarkerStyle(20); h_dataFit2->SetMarkerColor(1); h_dataFit2->SetStats(0);
  h_dataFit2->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2->GetXaxis()->SetTitleSize(20); h_dataFit2->GetXaxis()->SetTitleFont(43); h_dataFit2->GetXaxis()->SetTitleOffset(3.0); h_dataFit2->GetXaxis()->SetLabelSize(15); h_dataFit2->GetXaxis()->SetLabelFont(43);
  h_dataFit2->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2->GetYaxis()->SetNdivisions(505); h_dataFit2->GetYaxis()->CenterTitle(); h_dataFit2->GetYaxis()->SetTitleSize(20); h_dataFit2->GetYaxis()->SetTitleFont(43); h_dataFit2->GetYaxis()->SetTitleOffset(1.2); h_dataFit2->GetYaxis()->SetLabelSize(15); h_dataFit2->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2 = (TH1F*) ds_dimuorphan_bg_m2->createHistogram("hdata2", m2, Binning(m_bins, m_min, m_max));
  TH1F* h_func2 = new TH1F("h_func2","", m_bins, m_min, m_max);
  w->pdf("template1D_m2")->fillHistogram(h_func2, m2, hdata2->GetEntries());
  for (unsigned int iB = 1; iB < m_bins; iB++) {
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
  //2018 PDF m2 (below J/psi)
  //=========================
  w->factory("EXPR::MmumuF_below_Jpsi('m2_below_Jpsi*pow( (m2_below_Jpsi/m_below_Jpsi)*(m2_below_Jpsi/m_below_Jpsi) - 1.0, MmumuF_p_below_Jpsi )*exp( -MmumuF_c_below_Jpsi*( (m2_below_Jpsi/m_below_Jpsi)*(m2_below_Jpsi/m_below_Jpsi) - 1.0 ) )', m2_below_Jpsi, m_below_Jpsi[0.2113], MmumuF_c_below_Jpsi[0.01, -0.3, 0.3], MmumuF_p_below_Jpsi[0.1, 0.0, 2.0])");
  w->factory("Bernstein::bgF_below_Jpsi(m2_below_Jpsi, {bF06_below_Jpsi[1.5, 0, 15.], bF16_below_Jpsi[1.5, 0., 3.], bF26_below_Jpsi[-0.2, -3., 3.], bF36_below_Jpsi[1.0, 0., 3.], bF46_below_Jpsi[0.01, -2.0, 5.], bF56_below_Jpsi[0.5, 0., 3.], bF66_below_Jpsi[0.01, -1, 4.]})");
  w->factory("Gaussian::adHocF_below_Jpsi(m2_below_Jpsi, adHocF_mass_below_Jpsi[0.3, 0.2, 0.6], adHocF_sigma_below_Jpsi[0.05, 0.001, 0.3])");
  w->factory("Gaussian::etaF_below_Jpsi(m2_below_Jpsi, 0.54786, 0.007)");
  w->factory("Gaussian::rhoF_below_Jpsi(m2_below_Jpsi, 0.78265, 0.009)");
  w->factory("Gaussian::phiF_below_Jpsi(m2_below_Jpsi, 1.01946, 0.01)");
  w->factory("SUM::template1D_m2_below_Jpsi(norm_adHocF_below_Jpsi[1500., 0., 2000.]*adHocF_below_Jpsi, norm_MmumuF_below_Jpsi[15000., 0., 20000.]*MmumuF_below_Jpsi, norm_bgF_below_Jpsi[7000., 1000., 20000.]*bgF_below_Jpsi, norm_etaF_below_Jpsi[10., 0., 20.]*etaF_below_Jpsi, norm_rhoF_below_Jpsi[250., 1., 1000.]*rhoF_below_Jpsi, norm_phiF_below_Jpsi[250., 1., 1000.]*phiF_below_Jpsi)");

  RooFitResult *rF_below_Jpsi = w->pdf("template1D_m2_below_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m2_below_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2 (below J/psi)---------------------"<<endl;
  rF_below_Jpsi->Print();

  RooPlot* plotF_below_Jpsi = w->var("m2_below_Jpsi")->frame(Title("BG template for orphan dimuon no high pT below Jpsi"), Bins(m_bins_below_Jpsi));
  w->data("ds_dimuorphan_bg_m2_below_Jpsi")->plotOn(plotF_below_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_below_Jpsi"));
  w->pdf("template1D_m2_below_Jpsi")->plotOn(plotF_below_Jpsi,LineColor(kRed),Precision(0.0001),Name("template1D_m2_below_Jpsi"));
  plotF_below_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV"); plotF_below_Jpsi->GetYaxis()->SetTitleOffset(1.1);

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2_below_Jpsi = new TCanvas("c_template1D_m2_below_Jpsi", "c_template1D_m2_below_Jpsi", 800, 800);
  c_template1D_m2_below_Jpsi->Clear();
  TPad *pad1b_below_Jpsi = new TPad("pad1b_below_Jpsi", "pad1b_below_Jpsi", 0, 0.3, 1, 1.0);
  pad1b_below_Jpsi->SetBottomMargin(0); pad1b_below_Jpsi->Draw(); pad1b_below_Jpsi->cd();
  plotF_below_Jpsi->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m2_below_Jpsi->cd(); c_template1D_m2_below_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b_below_Jpsi = new TPad("pad2b_below_Jpsi", "pad2b_below_Jpsi", 0, 0.0, 1, 0.29);
  pad2b_below_Jpsi->SetTopMargin(0); pad2b_below_Jpsi->SetBottomMargin(0.35); pad2b_below_Jpsi->SetGridy(); pad2b_below_Jpsi->Draw(); pad2b_below_Jpsi->cd();
  TH1F *h_dataFit2_below_Jpsi = new TH1F("h_dataFit2_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  h_dataFit2_below_Jpsi->SetLineColor(kBlack); h_dataFit2_below_Jpsi->SetMarkerStyle(20); h_dataFit2_below_Jpsi->SetMarkerColor(1); h_dataFit2_below_Jpsi->SetStats(0);
  h_dataFit2_below_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2_below_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit2_below_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit2_below_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit2_below_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit2_below_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit2_below_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2_below_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit2_below_Jpsi->GetYaxis()->CenterTitle(); h_dataFit2_below_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit2_below_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit2_below_Jpsi->GetYaxis()->SetTitleOffset(1.2); h_dataFit2_below_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit2_below_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2_below_Jpsi = (TH1F*) ds_dimuorphan_bg_m2_below_Jpsi->createHistogram("hdata2_below_Jpsi", m2_below_Jpsi, Binning(m_bins_below_Jpsi, m_min, m_Jpsi_dn));
  TH1F* h_func2_below_Jpsi = new TH1F("h_func2_below_Jpsi","", m_bins_below_Jpsi, m_min, m_Jpsi_dn);
  w->pdf("template1D_m2_below_Jpsi")->fillHistogram(h_func2_below_Jpsi, m2_below_Jpsi, hdata2_below_Jpsi->GetEntries());
  for (unsigned int iB = 1; iB < m_bins_below_Jpsi; iB++) {
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
  //2018 PDF m2 (above J/psi)
  //=========================
  w->factory("Bernstein::bgF_above_Jpsi(m2_above_Jpsi, {bF06_above_Jpsi[0.8, -2, 15.], bF16_above_Jpsi[-0.7, -4., 3.], bF26_above_Jpsi[3.0, -5., 5.], bF36_above_Jpsi[-2.5, -5., 5.], bF46_above_Jpsi[3, -3., 5.], bF56_above_Jpsi[-0.2, -3., 3.], bF66_above_Jpsi[1.0, -3, 4.]})");
  w->factory("Gaussian::psiF_above_Jpsi(m2_above_Jpsi, 3.68609, psiF_sigma_above_Jpsi[0.031, 0.01, 0.5])");
  w->factory("SUM::template1D_m2_above_Jpsi(norm_bgF_above_Jpsi[25000, 500., 50000.]*bgF_above_Jpsi, norm_psiF_above_Jpsi[500., 0., 3000.]*psiF_above_Jpsi)");

  RooFitResult *rF_above_Jpsi = w->pdf("template1D_m2_above_Jpsi")->fitTo(*(w->data("ds_dimuorphan_bg_m2_above_Jpsi")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2 (above J/psi)---------------------"<<endl;
  rF_above_Jpsi->Print();

  RooPlot* plotF_above_Jpsi = w->var("m2_above_Jpsi")->frame(Title("BG template for orphan dimuon no high pT above Jpsi"), Bins(m_bins_above_Jpsi));
  w->data("ds_dimuorphan_bg_m2_above_Jpsi")->plotOn(plotF_above_Jpsi, DataError(RooAbsData::SumW2), Name("data_m2_above_Jpsi"));
  w->pdf("template1D_m2_above_Jpsi")->plotOn(plotF_above_Jpsi,LineColor(kRed),Precision(0.0001),Name("template1D_m2_above_Jpsi"));
  plotF_above_Jpsi->GetYaxis()->SetTitle("Events/0.04GeV"); plotF_above_Jpsi->GetYaxis()->SetTitleOffset(1.1);

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2_above_Jpsi = new TCanvas("c_template1D_m2_above_Jpsi", "c_template1D_m2_above_Jpsi", 800, 800);
  c_template1D_m2_above_Jpsi->Clear();
  TPad *pad1b_above_Jpsi = new TPad("pad1b_above_Jpsi", "pad1b_above_Jpsi", 0, 0.3, 1, 1.0);
  pad1b_above_Jpsi->SetBottomMargin(0); pad1b_above_Jpsi->Draw(); pad1b_above_Jpsi->cd();
  plotF_above_Jpsi->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m2_above_Jpsi->cd(); c_template1D_m2_above_Jpsi->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b_above_Jpsi = new TPad("pad2b_above_Jpsi", "pad2b_above_Jpsi", 0, 0.0, 1, 0.29);
  pad2b_above_Jpsi->SetTopMargin(0); pad2b_above_Jpsi->SetBottomMargin(0.35); pad2b_above_Jpsi->SetGridy(); pad2b_above_Jpsi->Draw(); pad2b_above_Jpsi->cd();
  TH1F *h_dataFit2_above_Jpsi = new TH1F("h_dataFit2_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  h_dataFit2_above_Jpsi->SetLineColor(kBlack); h_dataFit2_above_Jpsi->SetMarkerStyle(20); h_dataFit2_above_Jpsi->SetMarkerColor(1); h_dataFit2_above_Jpsi->SetStats(0);
  h_dataFit2_above_Jpsi->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2_above_Jpsi->GetXaxis()->SetTitleSize(20); h_dataFit2_above_Jpsi->GetXaxis()->SetTitleFont(43); h_dataFit2_above_Jpsi->GetXaxis()->SetTitleOffset(3.0); h_dataFit2_above_Jpsi->GetXaxis()->SetLabelSize(15); h_dataFit2_above_Jpsi->GetXaxis()->SetLabelFont(43);
  h_dataFit2_above_Jpsi->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2_above_Jpsi->GetYaxis()->SetNdivisions(505); h_dataFit2_above_Jpsi->GetYaxis()->CenterTitle(); h_dataFit2_above_Jpsi->GetYaxis()->SetTitleSize(20); h_dataFit2_above_Jpsi->GetYaxis()->SetTitleFont(43); h_dataFit2_above_Jpsi->GetYaxis()->SetTitleOffset(1.2); h_dataFit2_above_Jpsi->GetYaxis()->SetLabelSize(15); h_dataFit2_above_Jpsi->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2_above_Jpsi = (TH1F*) ds_dimuorphan_bg_m2_above_Jpsi->createHistogram("hdata2_above_Jpsi", m2_above_Jpsi, Binning(m_bins_above_Jpsi, m_Jpsi_up, m_max));
  TH1F* h_func2_above_Jpsi = new TH1F("h_func2_above_Jpsi","", m_bins_above_Jpsi, m_Jpsi_up, m_max);
  w->pdf("template1D_m2_above_Jpsi")->fillHistogram(h_func2_above_Jpsi, m2_above_Jpsi, hdata2_above_Jpsi->GetEntries());
  for (unsigned int iB = 1; iB < m_bins_above_Jpsi; iB++) {
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

  //===========================
  //2018 PDF m2 (above Upsilon)
  //===========================
  w->factory("Bernstein::bgF_above_Upsilon(m2_above_Upsilon, {bF06_above_Upsilon[0.8, -2, 15.], bF16_above_Upsilon[-0.7, -4., 3.], bF26_above_Upsilon[3.0, -5., 5.], bF36_above_Upsilon[-2.5, -5., 5.], bF46_above_Upsilon[3, -3., 5.], bF56_above_Upsilon[-0.2, -3., 3.], bF66_above_Upsilon[1.0, -3, 4.]})");
  w->factory("SUM::template1D_m2_above_Upsilon(norm_bgF_above_Upsilon[50000, 500., 100000.]*bgF_above_Upsilon)");

  RooFitResult *rF_above_Upsilon = w->pdf("template1D_m2_above_Upsilon")->fitTo(*(w->data("ds_dimuorphan_bg_m2_above_Upsilon")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooFitResult for m2 (above Upsilon)---------------------"<<endl;
  rF_above_Upsilon->Print();

  RooPlot* plotF_above_Upsilon = w->var("m2_above_Upsilon")->frame(Title("BG template for orphan dimuon no high pT above Upsilon"), Bins(m_bins_above_Upsilon));
  w->data("ds_dimuorphan_bg_m2_above_Upsilon")->plotOn(plotF_above_Upsilon, DataError(RooAbsData::SumW2), Name("data_m2_above_Upsilon"));
  w->pdf("template1D_m2_above_Upsilon")->plotOn(plotF_above_Upsilon,LineColor(kRed),Precision(0.0001),Name("template1D_m2_above_Upsilon"));
  plotF_above_Upsilon->GetYaxis()->SetTitle("Events/0.5GeV");

  // Upper pad: fit overlay data
  TCanvas * c_template1D_m2_above_Upsilon = new TCanvas("c_template1D_m2_above_Upsilon", "c_template1D_m2_above_Upsilon", 800, 800);
  c_template1D_m2_above_Upsilon->Clear();
  TPad *pad1b_above_Upsilon = new TPad("pad1b_above_Upsilon", "pad1b_above_Upsilon", 0, 0.3, 1, 1.0);
  pad1b_above_Upsilon->SetBottomMargin(0); pad1b_above_Upsilon->Draw(); pad1b_above_Upsilon->cd();
  plotF_above_Upsilon->Draw("same");
  txtHeaderfit->Draw("same");
  c_template1D_m2_above_Upsilon->cd(); c_template1D_m2_above_Upsilon->Update();

  // Lower pad: Fit/Data ratio
  TPad *pad2b_above_Upsilon = new TPad("pad2b_above_Upsilon", "pad2b_above_Upsilon", 0, 0.0, 1, 0.29);
  pad2b_above_Upsilon->SetTopMargin(0); pad2b_above_Upsilon->SetBottomMargin(0.35); pad2b_above_Upsilon->SetGridy(); pad2b_above_Upsilon->Draw(); pad2b_above_Upsilon->cd();
  TH1F *h_dataFit2_above_Upsilon = new TH1F("h_dataFit2_above_Upsilon","", m_bins_above_Upsilon, m_Upsilon_up, m_highmax);
  h_dataFit2_above_Upsilon->SetLineColor(kBlack); h_dataFit2_above_Upsilon->SetMarkerStyle(20); h_dataFit2_above_Upsilon->SetMarkerColor(1); h_dataFit2_above_Upsilon->SetStats(0);
  h_dataFit2_above_Upsilon->GetXaxis()->SetTitle("Mass [GeV]"); h_dataFit2_above_Upsilon->GetXaxis()->SetTitleSize(20); h_dataFit2_above_Upsilon->GetXaxis()->SetTitleFont(43); h_dataFit2_above_Upsilon->GetXaxis()->SetTitleOffset(3.0); h_dataFit2_above_Upsilon->GetXaxis()->SetLabelSize(15); h_dataFit2_above_Upsilon->GetXaxis()->SetLabelFont(43);
  h_dataFit2_above_Upsilon->GetYaxis()->SetTitle("Fit/Data"); h_dataFit2_above_Upsilon->GetYaxis()->SetNdivisions(505); h_dataFit2_above_Upsilon->GetYaxis()->CenterTitle(); h_dataFit2_above_Upsilon->GetYaxis()->SetTitleSize(20); h_dataFit2_above_Upsilon->GetYaxis()->SetTitleFont(43); h_dataFit2_above_Upsilon->GetYaxis()->SetTitleOffset(1.2); h_dataFit2_above_Upsilon->GetYaxis()->SetLabelSize(15); h_dataFit2_above_Upsilon->GetYaxis()->SetLabelFont(43);
  TH1F *hdata2_above_Upsilon = (TH1F*) ds_dimuorphan_bg_m2_above_Upsilon->createHistogram("hdata2_above_Upsilon", m2_above_Upsilon, Binning(m_bins_above_Upsilon, m_Upsilon_up, m_highmax));
  TH1F* h_func2_above_Upsilon = new TH1F("h_func2_above_Upsilon","", m_bins_above_Upsilon, m_Upsilon_up, m_highmax);
  w->pdf("template1D_m2_above_Upsilon")->fillHistogram(h_func2_above_Upsilon, m2_above_Upsilon, hdata2_above_Upsilon->GetEntries());
  for (unsigned int iB=1; iB < m_bins_above_Upsilon; iB++) {
    float ratio = h_func2_above_Upsilon->GetBinContent(iB)/hdata2_above_Upsilon->GetBinContent(iB);
    h_dataFit2_above_Upsilon->SetBinContent(iB, ratio );
    h_dataFit2_above_Upsilon->SetBinError(iB, ratio/sqrt(hdata2_above_Upsilon->GetBinContent(iB)) );
  }
  h_dataFit2_above_Upsilon->SetMinimum(0.5); h_dataFit2_above_Upsilon->SetMaximum(1.5); h_dataFit2_above_Upsilon->Sumw2(); h_dataFit2_above_Upsilon->SetStats(0); h_dataFit2_above_Upsilon->SetMarkerStyle(21); h_dataFit2_above_Upsilon->Draw("ep");
  c_template1D_m2_above_Upsilon->SaveAs("figures/template1D_m2_above_Upsilon.pdf");
  c_template1D_m2_above_Upsilon->SaveAs("figures/template1D_m2_above_Upsilon.png");
  c_template1D_m2_above_Upsilon->SaveAs("figures/template1D_m2_above_Upsilon.root");
  float chi2_F_above_Upsilon = plotF_above_Upsilon->chiSquare(10);
  cout<<"------------------ End m2 (above Upsilon) ---------------------"<<endl;

  //=***************************************************************************
  //                     Create 2D template = m1 x m2
  //=***************************************************************************
  cout << "-----Creating 2D template m1 * m2:-----" << endl;
  cout << "1D template m1 fit chi^2/dof: " << chi2_C << endl;
  cout << "1D template m2 fit chi^2/dof: " << chi2_F << endl;
  w->factory("PROD::template2D(template1D_m1, template1D_m2)");//Assume: two 1D p.d.f.s are not correlated

  cout << "-----Creating 2D templates m1 * m2 (exclude J/psi: below + above):-----" << endl;
  cout << "1D template m1 (below J/psi) fit chi^2/dof: " << chi2_C_below_Jpsi << endl;
  cout << "1D template m1 (above J/psi) fit chi^2/dof: " << chi2_C_above_Jpsi << endl;
  cout << "1D template m2 (below J/psi) fit chi^2/dof: " << chi2_F_below_Jpsi << endl;
  cout << "1D template m2 (above J/psi) fit chi^2/dof: " << chi2_F_above_Jpsi << endl;

  cout << "1D template m1 (above Upsilon) fit chi^2/dof: " << chi2_C_above_Upsilon << endl;
  cout << "1D template m2 (above Upsilon) fit chi^2/dof: " << chi2_F_above_Upsilon << endl;
  //These 2D pdf will be used to evaluate the contributions only below and only above Jpsi
  w->factory("PROD::template2D_below_Jpsi(template1D_m1_below_Jpsi, template1D_m2_below_Jpsi)");
  w->factory("PROD::template2D_above_Jpsi(template1D_m1_above_Jpsi, template1D_m2_above_Jpsi)");
  w->factory("PROD::template2D_above_Upsilon(template1D_m1_above_Upsilon, template1D_m2_above_Upsilon)");

  //=***************************************************************************
  //                     Extract 1D J/psi template from template m1
  //=***************************************************************************
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

  RooPlot* plot_Jpsi_m1 = w->var("m1")->frame(Title("J/psi template m1"), Bins(m_bins));
  w->pdf("Jpsi_m1")->plotOn(plot_Jpsi_m1, LineColor(kRed), Precision(0.0001), Name("plot_Jpsi_m1"));

  TCanvas * c_template1D_Jpsi_m1_RooPlot = new TCanvas("c_template1D_Jpsi_m1_RooPlot", "c_template1D_Jpsi_m1_RooPlot");
  c_template1D_Jpsi_m1_RooPlot->cd();
  plot_Jpsi_m1->GetYaxis()->SetTitle("Events/0.04GeV");
  plot_Jpsi_m1->Draw(); txtHeader->Draw();
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1.pdf");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1.png");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1.root");

  //=***************************************************************************
  //                       Extract 1D J/psi template from template m2
  //=***************************************************************************
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

  RooPlot* plot_Jpsi_m2 = w->var("m2")->frame(Title("J/psi template m2"), Bins(m_bins));
  w->pdf("Jpsi_m2")->plotOn(plot_Jpsi_m2, LineColor(kRed), Precision(0.0001), Name("plot_Jpsi_m2"));

  TCanvas * c_template1D_Jpsi_m2_RooPlot = new TCanvas("c_template1D_Jpsi_m2_RooPlot", "c_template1D_Jpsi_m2_RooPlot");
  c_template1D_Jpsi_m2_RooPlot->cd();
  plot_Jpsi_m2->GetYaxis()->SetTitle("Events/0.04GeV");
  plot_Jpsi_m2->Draw(); txtHeader->Draw();
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2.pdf");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2.png");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2.root");

  //=***************************************************************************
  //                     Create 2D template (m1 x m2) for J/psi
  //=***************************************************************************
  cout << "Create 2D template (m1 * m2) for J/psi" << endl;
  w->factory("PROD::Jpsi_2D(Jpsi_m1,Jpsi_m2)");

  //=***************************************************************************
  //           For later use in LowMassBKGPlot2D.C: datasets of 2 dimu events
  //=***************************************************************************
  //To be used for scatter plot 2 dimu events @ CR later in LowMassBKGPlot2D.C
  cout<<"                               " <<endl;
  cout << "Create trees on signal events" << endl;

  //===================
  //= Below Jpsi only =
  //===================
  //HERE: No need to apply cuts associated to higher mass signals above 11 GeV, such as DY cut and SAmu bkg cut
  // CR for bkg estimate and validate
  ostringstream stream_cut_control_Iso_offDiagonal_below_Jpsi;
  stream_cut_control_Iso_offDiagonal_below_Jpsi << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) >= BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_min << " && diMuonC_FittedVtx_m < " << m_Jpsi_dn << " && diMuonF_FittedVtx_m > " << m_min << " && diMuonF_FittedVtx_m < " << m_Jpsi_dn;
  TString cut_control_Iso_offDiagonal_below_Jpsi = stream_cut_control_Iso_offDiagonal_below_Jpsi.str();
  std::cout << "2-dimu CR selctions (low mass below J/psi only): " << cut_control_Iso_offDiagonal_below_Jpsi.Data() << std::endl;
  std::cout << "                                                 " << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_below_Jpsi);

  // SR for unblinding
  ostringstream stream_cut_signal_below_Jpsi;
  stream_cut_signal_below_Jpsi << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_min << " && diMuonC_FittedVtx_m < " << m_Jpsi_dn << " && diMuonF_FittedVtx_m > " << m_min << " && diMuonF_FittedVtx_m < " << m_Jpsi_dn;
  TString cut_signal_below_Jpsi = stream_cut_signal_below_Jpsi.str();
  std::cout << "2-dimu SR selctions (low mass below J/psi only): " << cut_signal_below_Jpsi.Data() << std::endl;
  std::cout << "                                                 " << std::endl;
  TTree* tree_dimudimu_signal_2D_below_Jpsi = chain_data_dimudimu.CopyTree(cut_signal_below_Jpsi);

  //===================
  //= Above Jpsi only =
  //===================
  //HERE: No need to apply cuts associated to higher mass signals above 11 GeV, such as DY cut and SAmu bkg cut
  // CR for bkg estimate and validate
  ostringstream stream_cut_control_Iso_offDiagonal_above_Jpsi;
  stream_cut_control_Iso_offDiagonal_above_Jpsi << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) >= BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Jpsi_up << " && diMuonC_FittedVtx_m < " << m_max << " && diMuonF_FittedVtx_m > " << m_Jpsi_up << " && diMuonF_FittedVtx_m < " << m_max;
  TString cut_control_Iso_offDiagonal_above_Jpsi = stream_cut_control_Iso_offDiagonal_above_Jpsi.str();
  std::cout << "2-dimu CR selctions (low mass above J/psi only): " << cut_control_Iso_offDiagonal_above_Jpsi.Data() << std::endl;
  std::cout << "                                                 " << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_above_Jpsi);

  // SR for unblinding
  ostringstream stream_cut_signal_above_Jpsi;
  stream_cut_signal_above_Jpsi << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Jpsi_up << " && diMuonC_FittedVtx_m < " << m_max << " && diMuonF_FittedVtx_m > " << m_Jpsi_up << " && diMuonF_FittedVtx_m < " << m_max;
  TString cut_signal_above_Jpsi = stream_cut_signal_above_Jpsi.str();
  std::cout << "2-dimu SR selctions (low mass above J/psi only): " << cut_signal_above_Jpsi.Data() << std::endl;
  std::cout << "                                                 " << std::endl;
  TTree* tree_dimudimu_signal_2D_above_Jpsi = chain_data_dimudimu.CopyTree(cut_signal_above_Jpsi);

  //===========================
  //= Above Upsion only: data =
  //===========================
  //HERE: Need to apply all cuts associated to higher mass signals above 11 GeV, such as DY cut and SAmu bkg cut
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!   Note: the data-driven template from above is not accurate,
  //!!!         so this signal dataset should be used just as a simple reference when you use it later in 2D plot,
  //!!!         you should NOT use signal dataset below to infer bkg events at signal region based on the 2D template
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // CR for bkg estimate and validate (For test: use with caution!!!)
  ostringstream stream_cut_control_Iso_offDiagonal_above_Upsilon;
  stream_cut_control_Iso_offDiagonal_above_Upsilon << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) && ( (dimuC_Mu0_SA==1 && muJetC_Mu0_matched_segs>=2) || (dimuC_Mu1_SA==1 && muJetC_Mu1_matched_segs>=2) || (dimuF_Mu0_SA==1 && muJetF_Mu0_matched_segs>=2) || (dimuF_Mu1_SA==1 && muJetF_Mu1_matched_segs>=2) ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) >= BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Upsilon_up << " && diMuonC_FittedVtx_m < " << m_highmax << " && diMuonF_FittedVtx_m > " << m_Upsilon_up << " && diMuonF_FittedVtx_m < " << m_highmax;
  TString cut_control_Iso_offDiagonal_above_Upsilon = stream_cut_control_Iso_offDiagonal_above_Upsilon.str();
  std::cout << "2-dimu CR selctions (high mass above Upsilon only): " << cut_control_Iso_offDiagonal_above_Upsilon.Data() << std::endl;
  std::cout << "                                                    " << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_above_Upsilon);

  // SR for unblinding (For test: use with caution!!!)
  ostringstream stream_cut_signal_above_Upsilon;
  stream_cut_signal_above_Upsilon << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) && ( (dimuC_Mu0_SA==1 && muJetC_Mu0_matched_segs>=2) || (dimuC_Mu1_SA==1 && muJetC_Mu1_matched_segs>=2) || (dimuF_Mu0_SA==1 && muJetF_Mu0_matched_segs>=2) || (dimuF_Mu1_SA==1 && muJetF_Mu1_matched_segs>=2) ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Upsilon_up << " && diMuonC_FittedVtx_m < " << m_highmax << " && diMuonF_FittedVtx_m > " << m_Upsilon_up << " && diMuonF_FittedVtx_m < " << m_highmax;
  TString cut_signal_above_Upsilon = stream_cut_signal_above_Upsilon.str();
  std::cout << "2-dimu SR selctions (high mass above Upsilon only): " << cut_signal_above_Upsilon.Data() << std::endl;
  std::cout << "                                                    " << std::endl;
  TTree* tree_dimudimu_signal_2D_above_Upsilon = chain_data_dimudimu.CopyTree(cut_signal_above_Upsilon);

  //======================================
  //= all mass range (0.21-60 GeV): data =
  //======================================
  // CR
  ostringstream stream_cut_control_Iso_offDiagonal_all; // Add cuts from three regions
  stream_cut_control_Iso_offDiagonal_all << "( is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) && ( (dimuC_Mu0_SA==1 && muJetC_Mu0_matched_segs>=2) || (dimuC_Mu1_SA==1 && muJetC_Mu1_matched_segs>=2) || (dimuF_Mu0_SA==1 && muJetF_Mu0_matched_segs>=2) || (dimuF_Mu1_SA==1 && muJetF_Mu1_matched_segs>=2) ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) >= BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Upsilon_up << " && diMuonC_FittedVtx_m < " << m_highmax << " && diMuonF_FittedVtx_m > " << m_Upsilon_up << " && diMuonF_FittedVtx_m < " << m_highmax << ") || ( is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) >= BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Jpsi_up << " && diMuonC_FittedVtx_m < " << m_max << " && diMuonF_FittedVtx_m > " << m_Jpsi_up << " && diMuonF_FittedVtx_m < " << m_max << ") || ( is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) >= BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_min << " && diMuonC_FittedVtx_m < " << m_Jpsi_dn << " && diMuonF_FittedVtx_m > " << m_min << " && diMuonF_FittedVtx_m < " << m_Jpsi_dn << ")";
  TString cut_control_Iso_offDiagonal_all = stream_cut_control_Iso_offDiagonal_all.str();
  std::cout << "2-dimu CR selctions (ALL mass range): " << cut_control_Iso_offDiagonal_all.Data() << std::endl;
  std::cout << "                                                    " << std::endl;
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D_all = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal_all);

  // SR
  ostringstream stream_cut_signal_all; // Add cuts from three regions
  stream_cut_signal_all << "( is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) && ( (dimuC_Mu0_SA==1 && muJetC_Mu0_matched_segs>=2) || (dimuC_Mu1_SA==1 && muJetC_Mu1_matched_segs>=2) || (dimuF_Mu0_SA==1 && muJetF_Mu0_matched_segs>=2) || (dimuF_Mu1_SA==1 && muJetF_Mu1_matched_segs>=2) ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Upsilon_up << " && diMuonC_FittedVtx_m < " << m_highmax << " && diMuonF_FittedVtx_m > " << m_Upsilon_up << " && diMuonF_FittedVtx_m < " << m_highmax << ") || ( is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Jpsi_up << " && diMuonC_FittedVtx_m < " << m_max << " && diMuonF_FittedVtx_m > " << m_Jpsi_up << " && diMuonF_FittedVtx_m < " << m_max << ") || ( is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_min << " && diMuonC_FittedVtx_m < " << m_Jpsi_dn << " && diMuonF_FittedVtx_m > " << m_min << " && diMuonF_FittedVtx_m < " << m_Jpsi_dn << ")";
  TString cut_signal_all = stream_cut_signal_all.str();
  std::cout << "2-dimu SR selctions (high mass above Upsilon only): " << cut_signal_all.Data() << std::endl;
  std::cout << "                                                    " << std::endl;
  TTree* tree_dimudimu_signal_2D_all = chain_data_dimudimu.CopyTree(cut_signal_all);

  // Below Jpsi only
  tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_below_Jpsi");
  tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_below_Jpsi");
  tree_dimudimu_signal_2D_below_Jpsi->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_below_Jpsi");
  tree_dimudimu_signal_2D_below_Jpsi->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_below_Jpsi");
  // Above Jpsi only
  tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Jpsi");
  tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Jpsi");
  tree_dimudimu_signal_2D_above_Jpsi->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Jpsi");
  tree_dimudimu_signal_2D_above_Jpsi->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Jpsi");
  // Above Upsilon only (for test: use with caution!!!)
  tree_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  tree_dimudimu_signal_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_signal_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  // all range
  tree_dimudimu_control_Iso_offDiagonal_2D_all->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_all");
  tree_dimudimu_control_Iso_offDiagonal_2D_all->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_all");
  tree_dimudimu_signal_2D_all->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_all");
  tree_dimudimu_signal_2D_all->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_all");

  // Below Jpsi only version
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi", "ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi", tree_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi, RooArgSet(m1_below_Jpsi, m2_below_Jpsi));
  RooDataSet* ds_dimudimu_signal_2D_below_Jpsi = new RooDataSet("ds_dimudimu_signal_2D_below_Jpsi", "ds_dimudimu_signal_2D_below_Jpsi", tree_dimudimu_signal_2D_below_Jpsi, RooArgSet(m1_below_Jpsi, m2_below_Jpsi));
  // Above Jpsi only version
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi", "ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi", tree_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi, RooArgSet(m1_above_Jpsi, m2_above_Jpsi));
  RooDataSet* ds_dimudimu_signal_2D_above_Jpsi = new RooDataSet("ds_dimudimu_signal_2D_above_Jpsi", "ds_dimudimu_signal_2D_above_Jpsi", tree_dimudimu_signal_2D_above_Jpsi, RooArgSet(m1_above_Jpsi, m2_above_Jpsi));
  // Above Upsilon only version (for test: use with caution!!!)
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon", "ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon", tree_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon, RooArgSet(m1_above_Upsilon, m2_above_Upsilon));
  RooDataSet* ds_dimudimu_signal_2D_above_Upsilon = new RooDataSet("ds_dimudimu_signal_2D_above_Upsilon", "ds_dimudimu_signal_2D_above_Upsilon", tree_dimudimu_signal_2D_above_Upsilon, RooArgSet(m1_above_Upsilon, m2_above_Upsilon));
  // For all range 0.21-60 GeV (sanity check after unblinding)
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D_all = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D_all", "ds_dimudimu_control_Iso_offDiagonal_2D_all", tree_dimudimu_control_Iso_offDiagonal_2D_all, RooArgSet(m1_all, m2_all));
  RooDataSet* ds_dimudimu_signal_2D_all = new RooDataSet("ds_dimudimu_signal_2D_all", "ds_dimudimu_signal_2D_all", tree_dimudimu_signal_2D_all, RooArgSet(m1_all, m2_all));

  // Below Jpsi only
  ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->Print("s");
  ds_dimudimu_signal_2D_below_Jpsi->Print("s");
  // Above Jpsi only
  ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->Print("s");
  ds_dimudimu_signal_2D_above_Jpsi->Print("s");
  // Above Upsilon only (for test: use with caution!!!)
  ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->Print("s");
  ds_dimudimu_signal_2D_above_Upsilon->Print("s");
  // all range
  ds_dimudimu_control_Iso_offDiagonal_2D_all->Print("s");
  ds_dimudimu_signal_2D_all->Print("s");

  // Below Jpsi only
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi);
  w->import(*ds_dimudimu_signal_2D_below_Jpsi);
  // Above Jpsi only
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi);
  w->import(*ds_dimudimu_signal_2D_above_Jpsi);
  // Above Upsilon only (for test: use with caution!!!)
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon);
  w->import(*ds_dimudimu_signal_2D_above_Upsilon);
  // all range
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D_all);
  w->import(*ds_dimudimu_signal_2D_all);

  //=============================
  //= Above Upsion only: BKG MC =
  //=============================
  //---------------------------------------------------------------------------------------------------------------------
  //Below we use kernel density estimation to get smooth background shape at SR using unbinned MC event (scaled to data):
  //Kernal choices are Gaussian and Crystal Ball functions. Use the maximum difference in the shape as systematic unc.
  //Note:
  //1. Contributions from each MC needs to be plotted separately using their own scale factors to data.
  //2. Bandwidth mimics signal resolution; so increase the bandwidth at higher mass as the resolution.
  //---------------------------------------------------------------------------------------------------------------------
  // SR cut
  ostringstream stream_cut_bkgmc_above_Upsilon;
  stream_cut_bkgmc_above_Upsilon << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) && ( (dimuC_Mu0_SA==1 && muJetC_Mu0_matched_segs>=2) || (dimuC_Mu1_SA==1 && muJetC_Mu1_matched_segs>=2) || (dimuF_Mu0_SA==1 && muJetF_Mu0_matched_segs>=2) || (dimuF_Mu1_SA==1 && muJetF_Mu1_matched_segs>=2) ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && diMuonC_IsoTk_FittedVtx < " << iso_cut << " && diMuonF_IsoTk_FittedVtx < " << iso_cut << " && TMath::Abs(diMuonC_FittedVtx_m-diMuonF_FittedVtx_m) < BKG_cfg::My_MassWindow(diMuonC_FittedVtx_m, diMuonF_FittedVtx_m) && diMuonC_FittedVtx_m > " << m_Upsilon_up << " && diMuonC_FittedVtx_m < " << m_highmax << " && diMuonF_FittedVtx_m > " << m_Upsilon_up << " && diMuonF_FittedVtx_m < " << m_highmax;
  TString cut_bkgmc_above_Upsilon = stream_cut_bkgmc_above_Upsilon.str();
  std::cout << "MC SR 2-dimu events (high mass above Upsilon only): " << cut_bkgmc_above_Upsilon.Data() << std::endl;
  std::cout << "                                                    " << std::endl;

  //MC Ntuples input file from each SM process: // DY0J: no contribution, 0 entries, not included
  //Include: DY1J, DY2J, qqToZZTo4L, TTJets_DiLept, ggHToZZTo4L, ggToZZTo4mu
  TChain chain_DY1J_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream MyDY1Jfile(DY1JFile);//this is BKG MC Ntuple from Config.h
  std::string MyDY1JLine;
  if( !MyDY1Jfile ) std::cout<<"ERROR opening DY1J file."<<std::endl;
  while (std::getline(MyDY1Jfile, MyDY1JLine)){
    TString DY1JLine(MyDY1JLine);
    if( DY1JLine.Contains("root") ) chain_DY1J_dimudimu.Add(DY1JLine.Data());
  }

  TChain chain_DY2J_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream MyDY2Jfile(DY2JFile);
  std::string MyDY2JLine;
  if( !MyDY2Jfile ) std::cout<<"ERROR opening DY2J file."<<std::endl;
  while (std::getline(MyDY2Jfile, MyDY2JLine)){
    TString DY2JLine(MyDY2JLine);
    if( DY2JLine.Contains("root") ) chain_DY2J_dimudimu.Add(DY2JLine.Data());
  }

  TChain chain_qqToZZ_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream MyqqToZZfile(qqToZZFile);
  std::string MyqqToZZLine;
  if( !MyqqToZZfile ) std::cout<<"ERROR opening qqToZZ file."<<std::endl;
  while (std::getline(MyqqToZZfile, MyqqToZZLine)){
    TString qqToZZLine(MyqqToZZLine);
    if( qqToZZLine.Contains("root") ) chain_qqToZZ_dimudimu.Add(qqToZZLine.Data());
  }

  TChain chain_TTJets_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream MyTTJetsfile(TTJetsFile);
  std::string MyTTJetsLine;
  if( !MyTTJetsfile ) std::cout<<"ERROR opening TTJets file."<<std::endl;
  while (std::getline(MyTTJetsfile, MyTTJetsLine)){
    TString TTJetsLine(MyTTJetsLine);
    if( TTJetsLine.Contains("root") ) chain_TTJets_dimudimu.Add(TTJetsLine.Data());
  }

  TChain chain_ggHToZZ_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream MyggHToZZfile(ggHToZZFile);
  std::string MyggHToZZLine;
  if( !MyggHToZZfile ) std::cout<<"ERROR opening ggHToZZ file."<<std::endl;
  while (std::getline(MyggHToZZfile, MyggHToZZLine)){
    TString ggHToZZLine(MyggHToZZLine);
    if( ggHToZZLine.Contains("root") ) chain_ggHToZZ_dimudimu.Add(ggHToZZLine.Data());
  }

  TChain chain_ggToZZ_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream MyggToZZfile(ggToZZFile);
  std::string MyggToZZLine;
  if( !MyggToZZfile ) std::cout<<"ERROR opening ggToZZ file."<<std::endl;
  while (std::getline(MyggToZZfile, MyggToZZLine)){
    TString ggToZZLine(MyggToZZLine);
    if( ggToZZLine.Contains("root") ) chain_ggToZZ_dimudimu.Add(ggToZZLine.Data());
  }

  TTree* tree_dimudimu_DY1J_2D_above_Upsilon    = chain_DY1J_dimudimu.CopyTree(cut_bkgmc_above_Upsilon);
  TTree* tree_dimudimu_DY2J_2D_above_Upsilon    = chain_DY2J_dimudimu.CopyTree(cut_bkgmc_above_Upsilon);
  TTree* tree_dimudimu_qqToZZ_2D_above_Upsilon  = chain_qqToZZ_dimudimu.CopyTree(cut_bkgmc_above_Upsilon);
  TTree* tree_dimudimu_TTJets_2D_above_Upsilon  = chain_TTJets_dimudimu.CopyTree(cut_bkgmc_above_Upsilon);
  TTree* tree_dimudimu_ggHToZZ_2D_above_Upsilon = chain_ggHToZZ_dimudimu.CopyTree(cut_bkgmc_above_Upsilon);
  TTree* tree_dimudimu_ggToZZ_2D_above_Upsilon  = chain_ggToZZ_dimudimu.CopyTree(cut_bkgmc_above_Upsilon);

  tree_dimudimu_DY1J_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_DY1J_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  tree_dimudimu_DY2J_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_DY2J_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  tree_dimudimu_qqToZZ_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_qqToZZ_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  tree_dimudimu_TTJets_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_TTJets_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  tree_dimudimu_ggHToZZ_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_ggHToZZ_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");
  tree_dimudimu_ggToZZ_2D_above_Upsilon->GetBranch("diMuonC_FittedVtx_m")->SetName("m1_above_Upsilon");
  tree_dimudimu_ggToZZ_2D_above_Upsilon->GetBranch("diMuonF_FittedVtx_m")->SetName("m2_above_Upsilon");

  // 1-D dataset
  RooDataSet* ds_dimudimu_DY1J_m1_above_Upsilon    = new RooDataSet("ds_dimudimu_DY1J_m1_above_Upsilon",    "ds_dimudimu_DY1J_m1_above_Upsilon",    tree_dimudimu_DY1J_2D_above_Upsilon,    RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimudimu_DY1J_m2_above_Upsilon    = new RooDataSet("ds_dimudimu_DY1J_m2_above_Upsilon",    "ds_dimudimu_DY1J_m2_above_Upsilon",    tree_dimudimu_DY1J_2D_above_Upsilon,    RooArgSet(m2_above_Upsilon));
  RooDataSet* ds_dimudimu_DY2J_m1_above_Upsilon    = new RooDataSet("ds_dimudimu_DY2J_m1_above_Upsilon",    "ds_dimudimu_DY2J_m1_above_Upsilon",    tree_dimudimu_DY2J_2D_above_Upsilon,    RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimudimu_DY2J_m2_above_Upsilon    = new RooDataSet("ds_dimudimu_DY2J_m2_above_Upsilon",    "ds_dimudimu_DY2J_m2_above_Upsilon",    tree_dimudimu_DY2J_2D_above_Upsilon,    RooArgSet(m2_above_Upsilon));
  RooDataSet* ds_dimudimu_qqToZZ_m1_above_Upsilon  = new RooDataSet("ds_dimudimu_qqToZZ_m1_above_Upsilon",  "ds_dimudimu_qqToZZ_m1_above_Upsilon",  tree_dimudimu_qqToZZ_2D_above_Upsilon,  RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimudimu_qqToZZ_m2_above_Upsilon  = new RooDataSet("ds_dimudimu_qqToZZ_m2_above_Upsilon",  "ds_dimudimu_qqToZZ_m2_above_Upsilon",  tree_dimudimu_qqToZZ_2D_above_Upsilon,  RooArgSet(m2_above_Upsilon));
  RooDataSet* ds_dimudimu_TTJets_m1_above_Upsilon  = new RooDataSet("ds_dimudimu_TTJets_m1_above_Upsilon",  "ds_dimudimu_TTJets_m1_above_Upsilon",  tree_dimudimu_TTJets_2D_above_Upsilon,  RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimudimu_TTJets_m2_above_Upsilon  = new RooDataSet("ds_dimudimu_TTJets_m2_above_Upsilon",  "ds_dimudimu_TTJets_m2_above_Upsilon",  tree_dimudimu_TTJets_2D_above_Upsilon,  RooArgSet(m2_above_Upsilon));
  RooDataSet* ds_dimudimu_ggHToZZ_m1_above_Upsilon = new RooDataSet("ds_dimudimu_ggHToZZ_m1_above_Upsilon", "ds_dimudimu_ggHToZZ_m1_above_Upsilon", tree_dimudimu_ggHToZZ_2D_above_Upsilon, RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimudimu_ggHToZZ_m2_above_Upsilon = new RooDataSet("ds_dimudimu_ggHToZZ_m2_above_Upsilon", "ds_dimudimu_ggHToZZ_m2_above_Upsilon", tree_dimudimu_ggHToZZ_2D_above_Upsilon, RooArgSet(m2_above_Upsilon));
  RooDataSet* ds_dimudimu_ggToZZ_m1_above_Upsilon  = new RooDataSet("ds_dimudimu_ggToZZ_m1_above_Upsilon",  "ds_dimudimu_ggToZZ_m1_above_Upsilon",  tree_dimudimu_ggToZZ_2D_above_Upsilon,  RooArgSet(m1_above_Upsilon));
  RooDataSet* ds_dimudimu_ggToZZ_m2_above_Upsilon  = new RooDataSet("ds_dimudimu_ggToZZ_m2_above_Upsilon",  "ds_dimudimu_ggToZZ_m2_above_Upsilon",  tree_dimudimu_ggToZZ_2D_above_Upsilon,  RooArgSet(m2_above_Upsilon));

  ds_dimudimu_DY1J_m1_above_Upsilon->Print("s");
  ds_dimudimu_DY1J_m2_above_Upsilon->Print("s");
  ds_dimudimu_DY2J_m1_above_Upsilon->Print("s");
  ds_dimudimu_DY2J_m2_above_Upsilon->Print("s");
  ds_dimudimu_qqToZZ_m1_above_Upsilon->Print("s");
  ds_dimudimu_qqToZZ_m2_above_Upsilon->Print("s");
  ds_dimudimu_TTJets_m1_above_Upsilon->Print("s");
  ds_dimudimu_TTJets_m2_above_Upsilon->Print("s");
  ds_dimudimu_ggHToZZ_m1_above_Upsilon->Print("s");
  ds_dimudimu_ggHToZZ_m2_above_Upsilon->Print("s");
  ds_dimudimu_ggToZZ_m1_above_Upsilon->Print("s");
  ds_dimudimu_ggToZZ_m2_above_Upsilon->Print("s");
  w->import(*ds_dimudimu_DY1J_m1_above_Upsilon);
  w->import(*ds_dimudimu_DY1J_m2_above_Upsilon);
  w->import(*ds_dimudimu_DY2J_m1_above_Upsilon);
  w->import(*ds_dimudimu_DY2J_m2_above_Upsilon);
  w->import(*ds_dimudimu_qqToZZ_m1_above_Upsilon);
  w->import(*ds_dimudimu_qqToZZ_m2_above_Upsilon);
  w->import(*ds_dimudimu_TTJets_m1_above_Upsilon);
  w->import(*ds_dimudimu_TTJets_m2_above_Upsilon);
  w->import(*ds_dimudimu_ggHToZZ_m1_above_Upsilon);
  w->import(*ds_dimudimu_ggHToZZ_m2_above_Upsilon);
  w->import(*ds_dimudimu_ggToZZ_m1_above_Upsilon);
  w->import(*ds_dimudimu_ggToZZ_m2_above_Upsilon);

  //===========
  // 1-D pdf m1
  //===========
  //ROOT example on kernel estimation: https://root.cern.ch/doc/v608/rf707__kernelestimation_8C_source.html
  RooKeysPdf DY1J_pdf_m1("DY1J_pdf_m1",       "DY1J_pdf_m1",    m1_above_Upsilon, *ds_dimudimu_DY1J_m1_above_Upsilon, RooKeysPdf::NoMirror);//Adaptive Gaussian kernel
  RooKeysPdf DY1J_pdf_m1_dn("DY1J_pdf_m1_dn", "DY1J_pdf_m1_dn", m1_above_Upsilon, *ds_dimudimu_DY1J_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);// narrow width
  RooKeysPdf DY1J_pdf_m1_up("DY1J_pdf_m1_up", "DY1J_pdf_m1_up", m1_above_Upsilon, *ds_dimudimu_DY1J_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);// wider width

  RooKeysPdf DY2J_pdf_m1("DY2J_pdf_m1",       "DY2J_pdf_m1",    m1_above_Upsilon, *ds_dimudimu_DY2J_m1_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf DY2J_pdf_m1_dn("DY2J_pdf_m1_dn", "DY2J_pdf_m1_dn", m1_above_Upsilon, *ds_dimudimu_DY2J_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf DY2J_pdf_m1_up("DY2J_pdf_m1_up", "DY2J_pdf_m1_up", m1_above_Upsilon, *ds_dimudimu_DY2J_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf qqToZZ_pdf_m1("qqToZZ_pdf_m1",       "qqToZZ_pdf_m1",    m1_above_Upsilon, *ds_dimudimu_qqToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf qqToZZ_pdf_m1_dn("qqToZZ_pdf_m1_dn", "qqToZZ_pdf_m1_dn", m1_above_Upsilon, *ds_dimudimu_qqToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf qqToZZ_pdf_m1_up("qqToZZ_pdf_m1_up", "qqToZZ_pdf_m1_up", m1_above_Upsilon, *ds_dimudimu_qqToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf TTJets_pdf_m1("TTJets_pdf_m1",       "TTJets_pdf_m1",    m1_above_Upsilon, *ds_dimudimu_TTJets_m1_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf TTJets_pdf_m1_dn("TTJets_pdf_m1_dn", "TTJets_pdf_m1_dn", m1_above_Upsilon, *ds_dimudimu_TTJets_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf TTJets_pdf_m1_up("TTJets_pdf_m1_up", "TTJets_pdf_m1_up", m1_above_Upsilon, *ds_dimudimu_TTJets_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf ggHToZZ_pdf_m1("ggHToZZ_pdf_m1",       "ggHToZZ_pdf_m1",    m1_above_Upsilon, *ds_dimudimu_ggHToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf ggHToZZ_pdf_m1_dn("ggHToZZ_pdf_m1_dn", "ggHToZZ_pdf_m1_dn", m1_above_Upsilon, *ds_dimudimu_ggHToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf ggHToZZ_pdf_m1_up("ggHToZZ_pdf_m1_up", "ggHToZZ_pdf_m1_up", m1_above_Upsilon, *ds_dimudimu_ggHToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf ggToZZ_pdf_m1("ggToZZ_pdf_m1",       "ggToZZ_pdf_m1",    m1_above_Upsilon, *ds_dimudimu_ggToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf ggToZZ_pdf_m1_dn("ggToZZ_pdf_m1_dn", "ggToZZ_pdf_m1_dn", m1_above_Upsilon, *ds_dimudimu_ggToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf ggToZZ_pdf_m1_up("ggToZZ_pdf_m1_up", "ggToZZ_pdf_m1_up", m1_above_Upsilon, *ds_dimudimu_ggToZZ_m1_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  w->import(DY1J_pdf_m1);
  w->import(DY1J_pdf_m1_dn);
  w->import(DY1J_pdf_m1_up);
  w->import(DY2J_pdf_m1);
  w->import(DY2J_pdf_m1_dn);
  w->import(DY2J_pdf_m1_up);
  w->import(qqToZZ_pdf_m1);
  w->import(qqToZZ_pdf_m1_dn);
  w->import(qqToZZ_pdf_m1_up);
  w->import(TTJets_pdf_m1);
  w->import(TTJets_pdf_m1_dn);
  w->import(TTJets_pdf_m1_up);
  w->import(ggHToZZ_pdf_m1);
  w->import(ggHToZZ_pdf_m1_dn);
  w->import(ggHToZZ_pdf_m1_up);
  w->import(ggToZZ_pdf_m1);
  w->import(ggToZZ_pdf_m1_dn);
  w->import(ggToZZ_pdf_m1_up);

  // Overlay pdf on MC events (NOT DATA!) in binned frame
  RooPlot* DY1Jframe1    = m1_above_Upsilon.frame(Title("DY1J MC + PDF: SR m1 binned"),    Bins(HM_m_bins));//same binning as MC template for reference
  RooPlot* DY2Jframe1    = m1_above_Upsilon.frame(Title("DY2J MC + PDF: SR m1 binned"),    Bins(HM_m_bins));
  RooPlot* qqToZZframe1  = m1_above_Upsilon.frame(Title("qqToZZ MC + PDF: SR m1 binned"),  Bins(HM_m_bins));
  RooPlot* TTJetsframe1  = m1_above_Upsilon.frame(Title("TTJets MC + PDF: SR m1 binned"),  Bins(HM_m_bins));
  RooPlot* ggHToZZframe1 = m1_above_Upsilon.frame(Title("ggHToZZ MC + PDF: SR m1 binned"), Bins(HM_m_bins));
  RooPlot* ggToZZframe1  = m1_above_Upsilon.frame(Title("ggToZZ MC + PDF: SR m1 binned"),  Bins(HM_m_bins));

  ds_dimudimu_DY1J_m1_above_Upsilon->plotOn(DY1Jframe1);       DY1J_pdf_m1.plotOn(DY1Jframe1);       DY1J_pdf_m1_dn.plotOn(DY1Jframe1, LineStyle(kDashed), LineColor(kRed));       DY1J_pdf_m1_up.plotOn(DY1Jframe1, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_DY2J_m1_above_Upsilon->plotOn(DY2Jframe1);       DY2J_pdf_m1.plotOn(DY2Jframe1);       DY2J_pdf_m1_dn.plotOn(DY2Jframe1, LineStyle(kDashed), LineColor(kRed));       DY2J_pdf_m1_up.plotOn(DY2Jframe1, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_qqToZZ_m1_above_Upsilon->plotOn(qqToZZframe1);   qqToZZ_pdf_m1.plotOn(qqToZZframe1);   qqToZZ_pdf_m1_dn.plotOn(qqToZZframe1, LineStyle(kDashed), LineColor(kRed));   qqToZZ_pdf_m1_up.plotOn(qqToZZframe1, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_TTJets_m1_above_Upsilon->plotOn(TTJetsframe1);   TTJets_pdf_m1.plotOn(TTJetsframe1);   TTJets_pdf_m1_dn.plotOn(TTJetsframe1, LineStyle(kDashed), LineColor(kRed));   TTJets_pdf_m1_up.plotOn(TTJetsframe1, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_ggHToZZ_m1_above_Upsilon->plotOn(ggHToZZframe1); ggHToZZ_pdf_m1.plotOn(ggHToZZframe1); ggHToZZ_pdf_m1_dn.plotOn(ggHToZZframe1, LineStyle(kDashed), LineColor(kRed)); ggHToZZ_pdf_m1_up.plotOn(ggHToZZframe1, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_ggToZZ_m1_above_Upsilon->plotOn(ggToZZframe1);   ggToZZ_pdf_m1.plotOn(ggToZZframe1);   ggToZZ_pdf_m1_dn.plotOn(ggToZZframe1, LineStyle(kDashed), LineColor(kRed));   ggToZZ_pdf_m1_up.plotOn(ggToZZframe1, LineStyle(kDashed), LineColor(kMagenta));

  // unbinned, smooth pdf
  RooPlot* DY1Jframe1nobin    = m1_above_Upsilon.frame(Title("DY1J PDF: SR m1 unbinned"));
  RooPlot* DY2Jframe1nobin    = m1_above_Upsilon.frame(Title("DY2J PDF: SR m1 unbinned"));
  RooPlot* qqToZZframe1nobin  = m1_above_Upsilon.frame(Title("qqToZZ PDF: SR m1 unbinned"));
  RooPlot* TTJetsframe1nobin  = m1_above_Upsilon.frame(Title("TTJets PDF: SR m1 unbinned"));
  RooPlot* ggHToZZframe1nobin = m1_above_Upsilon.frame(Title("ggHToZZ PDF: SR m1 unbinned"));
  RooPlot* ggToZZframe1nobin  = m1_above_Upsilon.frame(Title("ggToZZ PDF: SR m1 unbinned"));

  DY1J_pdf_m1.plotOn(DY1Jframe1nobin);       DY1J_pdf_m1_dn.plotOn(DY1Jframe1nobin, LineStyle(kDashed), LineColor(kRed));       DY1J_pdf_m1_up.plotOn(DY1Jframe1nobin, LineStyle(kDashed), LineColor(kMagenta));
  DY2J_pdf_m1.plotOn(DY2Jframe1nobin);       DY2J_pdf_m1_dn.plotOn(DY2Jframe1nobin, LineStyle(kDashed), LineColor(kRed));       DY2J_pdf_m1_up.plotOn(DY2Jframe1nobin, LineStyle(kDashed), LineColor(kMagenta));
  qqToZZ_pdf_m1.plotOn(qqToZZframe1nobin);   qqToZZ_pdf_m1_dn.plotOn(qqToZZframe1nobin, LineStyle(kDashed), LineColor(kRed));   qqToZZ_pdf_m1_up.plotOn(qqToZZframe1nobin, LineStyle(kDashed), LineColor(kMagenta));
  TTJets_pdf_m1.plotOn(TTJetsframe1nobin);   TTJets_pdf_m1_dn.plotOn(TTJetsframe1nobin, LineStyle(kDashed), LineColor(kRed));   TTJets_pdf_m1_up.plotOn(TTJetsframe1nobin, LineStyle(kDashed), LineColor(kMagenta));
  ggHToZZ_pdf_m1.plotOn(ggHToZZframe1nobin); ggHToZZ_pdf_m1_dn.plotOn(ggHToZZframe1nobin, LineStyle(kDashed), LineColor(kRed)); ggHToZZ_pdf_m1_up.plotOn(ggHToZZframe1nobin, LineStyle(kDashed), LineColor(kMagenta));
  ggToZZ_pdf_m1.plotOn(ggToZZframe1nobin);   ggToZZ_pdf_m1_dn.plotOn(ggToZZframe1nobin, LineStyle(kDashed), LineColor(kRed));   ggToZZ_pdf_m1_up.plotOn(ggToZZframe1nobin, LineStyle(kDashed), LineColor(kMagenta));

  //===========
  // 1-D pdf m2
  //===========
  RooKeysPdf DY1J_pdf_m2("DY1J_pdf_m2", "DY1J_pdf_m2",          m2_above_Upsilon, *ds_dimudimu_DY1J_m2_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf DY1J_pdf_m2_dn("DY1J_pdf_m2_dn", "DY1J_pdf_m2_dn", m2_above_Upsilon, *ds_dimudimu_DY1J_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf DY1J_pdf_m2_up("DY1J_pdf_m2_up", "DY1J_pdf_m2_up", m2_above_Upsilon, *ds_dimudimu_DY1J_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf DY2J_pdf_m2("DY2J_pdf_m2",       "DY2J_pdf_m2",    m2_above_Upsilon, *ds_dimudimu_DY2J_m2_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf DY2J_pdf_m2_dn("DY2J_pdf_m2_dn", "DY2J_pdf_m2_dn", m2_above_Upsilon, *ds_dimudimu_DY2J_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf DY2J_pdf_m2_up("DY2J_pdf_m2_up", "DY2J_pdf_m2_up", m2_above_Upsilon, *ds_dimudimu_DY2J_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf qqToZZ_pdf_m2("qqToZZ_pdf_m2",       "qqToZZ_pdf_m2",    m2_above_Upsilon, *ds_dimudimu_qqToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf qqToZZ_pdf_m2_dn("qqToZZ_pdf_m2_dn", "qqToZZ_pdf_m2_dn", m2_above_Upsilon, *ds_dimudimu_qqToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf qqToZZ_pdf_m2_up("qqToZZ_pdf_m2_up", "qqToZZ_pdf_m2_up", m2_above_Upsilon, *ds_dimudimu_qqToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf TTJets_pdf_m2("TTJets_pdf_m2",       "TTJets_pdf_m2",    m2_above_Upsilon, *ds_dimudimu_TTJets_m2_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf TTJets_pdf_m2_dn("TTJets_pdf_m2_dn", "TTJets_pdf_m2_dn", m2_above_Upsilon, *ds_dimudimu_TTJets_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf TTJets_pdf_m2_up("TTJets_pdf_m2_up", "TTJets_pdf_m2_up", m2_above_Upsilon, *ds_dimudimu_TTJets_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf ggHToZZ_pdf_m2("ggHToZZ_pdf_m2",       "ggHToZZ_pdf_m2",    m2_above_Upsilon, *ds_dimudimu_ggHToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf ggHToZZ_pdf_m2_dn("ggHToZZ_pdf_m2_dn", "ggHToZZ_pdf_m2_dn", m2_above_Upsilon, *ds_dimudimu_ggHToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf ggHToZZ_pdf_m2_up("ggHToZZ_pdf_m2_up", "ggHToZZ_pdf_m2_up", m2_above_Upsilon, *ds_dimudimu_ggHToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  RooKeysPdf ggToZZ_pdf_m2("ggToZZ_pdf_m2",       "ggToZZ_pdf_m2",    m2_above_Upsilon, *ds_dimudimu_ggToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror);
  RooKeysPdf ggToZZ_pdf_m2_dn("ggToZZ_pdf_m2_dn", "ggToZZ_pdf_m2_dn", m2_above_Upsilon, *ds_dimudimu_ggToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_dn);
  RooKeysPdf ggToZZ_pdf_m2_up("ggToZZ_pdf_m2_up", "ggToZZ_pdf_m2_up", m2_above_Upsilon, *ds_dimudimu_ggToZZ_m2_above_Upsilon, RooKeysPdf::NoMirror, width_scale_up);

  w->import(DY1J_pdf_m2);
  w->import(DY1J_pdf_m2_dn);
  w->import(DY1J_pdf_m2_up);
  w->import(DY2J_pdf_m2);
  w->import(DY2J_pdf_m2_dn);
  w->import(DY2J_pdf_m2_up);
  w->import(qqToZZ_pdf_m2);
  w->import(qqToZZ_pdf_m2_dn);
  w->import(qqToZZ_pdf_m2_up);
  w->import(TTJets_pdf_m2);
  w->import(TTJets_pdf_m2_dn);
  w->import(TTJets_pdf_m2_up);
  w->import(ggHToZZ_pdf_m2);
  w->import(ggHToZZ_pdf_m2_dn);
  w->import(ggHToZZ_pdf_m2_up);
  w->import(ggToZZ_pdf_m2);
  w->import(ggToZZ_pdf_m2_dn);
  w->import(ggToZZ_pdf_m2_up);

  // Overlay pdf on MC events (NOT DATA!) in binned frame
  RooPlot* DY1Jframe2    = m2_above_Upsilon.frame(Title("DY1J MC + PDF: SR m2 binned"), Bins(HM_m_bins));
  RooPlot* DY2Jframe2    = m2_above_Upsilon.frame(Title("DY2J MC + PDF: SR m2 binned"),    Bins(HM_m_bins));
  RooPlot* qqToZZframe2  = m2_above_Upsilon.frame(Title("qqToZZ MC + PDF: SR m2 binned"),  Bins(HM_m_bins));
  RooPlot* TTJetsframe2  = m2_above_Upsilon.frame(Title("TTJets MC + PDF: SR m2 binned"),  Bins(HM_m_bins));
  RooPlot* ggHToZZframe2 = m2_above_Upsilon.frame(Title("ggHToZZ MC + PDF: SR m2 binned"), Bins(HM_m_bins));
  RooPlot* ggToZZframe2  = m2_above_Upsilon.frame(Title("ggToZZ MC + PDF: SR m2 binned"),  Bins(HM_m_bins));

  ds_dimudimu_DY1J_m2_above_Upsilon->plotOn(DY1Jframe2);       DY1J_pdf_m2.plotOn(DY1Jframe2);       DY1J_pdf_m2_dn.plotOn(DY1Jframe2, LineStyle(kDashed), LineColor(kRed));       DY1J_pdf_m2_up.plotOn(DY1Jframe2, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_DY2J_m2_above_Upsilon->plotOn(DY2Jframe2);       DY2J_pdf_m2.plotOn(DY2Jframe2);       DY2J_pdf_m2_dn.plotOn(DY2Jframe2, LineStyle(kDashed), LineColor(kRed));       DY2J_pdf_m2_up.plotOn(DY2Jframe2, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_qqToZZ_m2_above_Upsilon->plotOn(qqToZZframe2);   qqToZZ_pdf_m2.plotOn(qqToZZframe2);   qqToZZ_pdf_m2_dn.plotOn(qqToZZframe2, LineStyle(kDashed), LineColor(kRed));   qqToZZ_pdf_m2_up.plotOn(qqToZZframe2, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_TTJets_m2_above_Upsilon->plotOn(TTJetsframe2);   TTJets_pdf_m2.plotOn(TTJetsframe2);   TTJets_pdf_m2_dn.plotOn(TTJetsframe2, LineStyle(kDashed), LineColor(kRed));   TTJets_pdf_m2_up.plotOn(TTJetsframe2, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_ggHToZZ_m2_above_Upsilon->plotOn(ggHToZZframe2); ggHToZZ_pdf_m2.plotOn(ggHToZZframe2); ggHToZZ_pdf_m2_dn.plotOn(ggHToZZframe2, LineStyle(kDashed), LineColor(kRed)); ggHToZZ_pdf_m2_up.plotOn(ggHToZZframe2, LineStyle(kDashed), LineColor(kMagenta));
  ds_dimudimu_ggToZZ_m2_above_Upsilon->plotOn(ggToZZframe2);   ggToZZ_pdf_m2.plotOn(ggToZZframe2);   ggToZZ_pdf_m2_dn.plotOn(ggToZZframe2, LineStyle(kDashed), LineColor(kRed));   ggToZZ_pdf_m2_up.plotOn(ggToZZframe2, LineStyle(kDashed), LineColor(kMagenta));

  // unbinned, smooth pdf
  RooPlot* DY1Jframe2nobin    = m2_above_Upsilon.frame(Title("DY1J PDF: SR m2 unbinned"));
  RooPlot* DY2Jframe2nobin    = m2_above_Upsilon.frame(Title("DY2J PDF: SR m2 unbinned"));
  RooPlot* qqToZZframe2nobin  = m2_above_Upsilon.frame(Title("qqToZZ PDF: SR m2 unbinned"));
  RooPlot* TTJetsframe2nobin  = m2_above_Upsilon.frame(Title("TTJets PDF: SR m2 unbinned"));
  RooPlot* ggHToZZframe2nobin = m2_above_Upsilon.frame(Title("ggHToZZ PDF: SR m2 unbinned"));
  RooPlot* ggToZZframe2nobin  = m2_above_Upsilon.frame(Title("ggToZZ PDF: SR m2 unbinned"));

  DY1J_pdf_m2.plotOn(DY1Jframe2nobin);       DY1J_pdf_m2_dn.plotOn(DY1Jframe2nobin, LineStyle(kDashed), LineColor(kRed));       DY1J_pdf_m2_up.plotOn(DY1Jframe2nobin, LineStyle(kDashed), LineColor(kMagenta));
  DY2J_pdf_m2.plotOn(DY2Jframe2nobin);       DY2J_pdf_m2_dn.plotOn(DY2Jframe2nobin, LineStyle(kDashed), LineColor(kRed));       DY2J_pdf_m2_up.plotOn(DY2Jframe2nobin, LineStyle(kDashed), LineColor(kMagenta));
  qqToZZ_pdf_m2.plotOn(qqToZZframe2nobin);   qqToZZ_pdf_m2_dn.plotOn(qqToZZframe2nobin, LineStyle(kDashed), LineColor(kRed));   qqToZZ_pdf_m2_up.plotOn(qqToZZframe2nobin, LineStyle(kDashed), LineColor(kMagenta));
  TTJets_pdf_m2.plotOn(TTJetsframe2nobin);   TTJets_pdf_m2_dn.plotOn(TTJetsframe2nobin, LineStyle(kDashed), LineColor(kRed));   TTJets_pdf_m2_up.plotOn(TTJetsframe2nobin, LineStyle(kDashed), LineColor(kMagenta));
  ggHToZZ_pdf_m2.plotOn(ggHToZZframe2nobin); ggHToZZ_pdf_m2_dn.plotOn(ggHToZZframe2nobin, LineStyle(kDashed), LineColor(kRed)); ggHToZZ_pdf_m2_up.plotOn(ggHToZZframe2nobin, LineStyle(kDashed), LineColor(kMagenta));
  ggToZZ_pdf_m2.plotOn(ggToZZframe2nobin);   ggToZZ_pdf_m2_dn.plotOn(ggToZZframe2nobin, LineStyle(kDashed), LineColor(kRed));   ggToZZ_pdf_m2_up.plotOn(ggToZZframe2nobin, LineStyle(kDashed), LineColor(kMagenta));

  //Draw
  TCanvas* DY1JPDF = new TCanvas("DY1JPDF", "DY1JPDF", 800, 800);
  DY1JPDF->Divide(2, 2);
  DY1JPDF->cd(1); gPad->SetLeftMargin(0.15); DY1Jframe1->GetYaxis()->SetTitleOffset(1.2);      DY1Jframe1->Draw();
  DY1JPDF->cd(2); gPad->SetLeftMargin(0.2);  DY1Jframe1nobin->GetYaxis()->SetTitleOffset(1.5); DY1Jframe1nobin->Draw();
  DY1JPDF->cd(3); gPad->SetLeftMargin(0.15); DY1Jframe2->GetYaxis()->SetTitleOffset(1.2);      DY1Jframe2->Draw();
  DY1JPDF->cd(4); gPad->SetLeftMargin(0.2);  DY1Jframe2nobin->GetYaxis()->SetTitleOffset(1.5); DY1Jframe2nobin->Draw();
  DY1JPDF->SaveAs("HighMassShape/DY1JPDF.pdf");
  DY1JPDF->SaveAs("HighMassShape/DY1JPDF.root");

  TCanvas* DY2JPDF = new TCanvas("DY2JPDF", "DY2JPDF", 800, 800);
  DY2JPDF->Divide(2, 2);
  DY2JPDF->cd(1); gPad->SetLeftMargin(0.15); DY2Jframe1->GetYaxis()->SetTitleOffset(1.2);      DY2Jframe1->Draw();
  DY2JPDF->cd(2); gPad->SetLeftMargin(0.2);  DY2Jframe1nobin->GetYaxis()->SetTitleOffset(1.5); DY2Jframe1nobin->Draw();
  DY2JPDF->cd(3); gPad->SetLeftMargin(0.15); DY2Jframe2->GetYaxis()->SetTitleOffset(1.2);      DY2Jframe2->Draw();
  DY2JPDF->cd(4); gPad->SetLeftMargin(0.2);  DY2Jframe2nobin->GetYaxis()->SetTitleOffset(1.5); DY2Jframe2nobin->Draw();
  DY2JPDF->SaveAs("HighMassShape/DY2JPDF.pdf");
  DY2JPDF->SaveAs("HighMassShape/DY2JPDF.root");

  TCanvas* qqToZZPDF = new TCanvas("qqToZZPDF", "qqToZZPDF", 800, 800);
  qqToZZPDF->Divide(2, 2);
  qqToZZPDF->cd(1); gPad->SetLeftMargin(0.15); qqToZZframe1->GetYaxis()->SetTitleOffset(1.2);      qqToZZframe1->Draw();
  qqToZZPDF->cd(2); gPad->SetLeftMargin(0.2);  qqToZZframe1nobin->GetYaxis()->SetTitleOffset(1.5); qqToZZframe1nobin->Draw();
  qqToZZPDF->cd(3); gPad->SetLeftMargin(0.15); qqToZZframe2->GetYaxis()->SetTitleOffset(1.2);      qqToZZframe2->Draw();
  qqToZZPDF->cd(4); gPad->SetLeftMargin(0.2);  qqToZZframe2nobin->GetYaxis()->SetTitleOffset(1.5); qqToZZframe2nobin->Draw();
  qqToZZPDF->SaveAs("HighMassShape/qqToZZPDF.pdf");
  qqToZZPDF->SaveAs("HighMassShape/qqToZZPDF.root");

  TCanvas* TTJetsPDF = new TCanvas("TTJetsPDF", "TTJetsPDF", 800, 800);
  TTJetsPDF->Divide(2, 2);
  TTJetsPDF->cd(1); gPad->SetLeftMargin(0.15); TTJetsframe1->GetYaxis()->SetTitleOffset(1.2);      TTJetsframe1->Draw();
  TTJetsPDF->cd(2); gPad->SetLeftMargin(0.2);  TTJetsframe1nobin->GetYaxis()->SetTitleOffset(1.5); TTJetsframe1nobin->Draw();
  TTJetsPDF->cd(3); gPad->SetLeftMargin(0.15); TTJetsframe2->GetYaxis()->SetTitleOffset(1.2);      TTJetsframe2->Draw();
  TTJetsPDF->cd(4); gPad->SetLeftMargin(0.2);  TTJetsframe2nobin->GetYaxis()->SetTitleOffset(1.5); TTJetsframe2nobin->Draw();
  TTJetsPDF->SaveAs("HighMassShape/TTJetsPDF.pdf");
  TTJetsPDF->SaveAs("HighMassShape/TTJetsPDF.root");

  TCanvas* ggHToZZPDF = new TCanvas("ggHToZZPDF", "ggHToZZPDF", 800, 800);
  ggHToZZPDF->Divide(2, 2);
  ggHToZZPDF->cd(1); gPad->SetLeftMargin(0.15); ggHToZZframe1->GetYaxis()->SetTitleOffset(1.2);      ggHToZZframe1->Draw();
  ggHToZZPDF->cd(2); gPad->SetLeftMargin(0.2);  ggHToZZframe1nobin->GetYaxis()->SetTitleOffset(1.5); ggHToZZframe1nobin->Draw();
  ggHToZZPDF->cd(3); gPad->SetLeftMargin(0.15); ggHToZZframe2->GetYaxis()->SetTitleOffset(1.2);      ggHToZZframe2->Draw();
  ggHToZZPDF->cd(4); gPad->SetLeftMargin(0.2);  ggHToZZframe2nobin->GetYaxis()->SetTitleOffset(1.5); ggHToZZframe2nobin->Draw();
  ggHToZZPDF->SaveAs("HighMassShape/ggHToZZPDF.pdf");
  ggHToZZPDF->SaveAs("HighMassShape/ggHToZZPDF.root");

  TCanvas* ggToZZPDF = new TCanvas("ggToZZPDF", "ggToZZPDF", 800, 800);
  ggToZZPDF->Divide(2, 2);
  ggToZZPDF->cd(1); gPad->SetLeftMargin(0.15); ggToZZframe1->GetYaxis()->SetTitleOffset(1.2);      ggToZZframe1->Draw();
  ggToZZPDF->cd(2); gPad->SetLeftMargin(0.2);  ggToZZframe1nobin->GetYaxis()->SetTitleOffset(1.5); ggToZZframe1nobin->Draw();
  ggToZZPDF->cd(3); gPad->SetLeftMargin(0.15); ggToZZframe2->GetYaxis()->SetTitleOffset(1.2);      ggToZZframe2->Draw();
  ggToZZPDF->cd(4); gPad->SetLeftMargin(0.2);  ggToZZframe2nobin->GetYaxis()->SetTitleOffset(1.5); ggToZZframe2nobin->Draw();
  ggToZZPDF->SaveAs("HighMassShape/ggToZZPDF.pdf");
  ggToZZPDF->SaveAs("HighMassShape/ggToZZPDF.root");

  //==============================================================================================
  //Sum PDF based on Entries * scale factor: // DY0J: no contribution, 0 entries, not included
  //see worksheet "Bootstrapping":
  //  https://docs.google.com/spreadsheets/d/10qXKHSw9QLrpPm2s0pASCbB27k6kVyEtTEwYi_oTra0/edit?usp=sharing
  //Include 6 processes: DY1J, DY2J, qqToZZTo4L, TTJets_DiLept, ggHToZZTo4L, ggToZZTo4mu
  //          Entries   *   Scale factor     =   PDF weight +/- 1sigma error
  // DY1J       2           6.2974242E-01        1.2594848      0.8905902711
  // DY2J       6           3.5720940E-01        2.1432564      0.8749807577
  // qqToZZ     1006        4.1624890E-03        4.1874639      0.1320237575
  // TTJets     25          1.0998854E-01        2.7497135      0.5499427205
  // ggHToZZ    1331        7.6245783E-04        1.0148314      0.0278166519
  // ggToZZ     3053        1.1110732E-04        0.3392106      0.006139119188
  //
  // The entries are also varied +/-1 sigma and then propagate to PDF weight error in the end column
  // Four different weights are used to evaluate the unc:
  //   PDF nominal weight +1 sigma
  //   PDF nominal weight -1 sigma
  //   Braid I  (nominal wgt +/- 1 sigma, start from +)
  //   Braid II (nominal wgt -/+ 1 sigma, start from -)
  // m1 and m2 use the same weight set
  //==============================================================================================

  // default nominal
  w->factory("SUM::HighMassFit2018_m1(1.2594848*DY1J_pdf_m1, 2.1432564*DY2J_pdf_m1, 4.1874639*qqToZZ_pdf_m1, 2.7497135*TTJets_pdf_m1, 1.0148314*ggHToZZ_pdf_m1, 0.3392106*ggToZZ_pdf_m1)");
  w->factory("SUM::HighMassFit2018_m2(1.2594848*DY1J_pdf_m2, 2.1432564*DY2J_pdf_m2, 4.1874639*qqToZZ_pdf_m2, 2.7497135*TTJets_pdf_m2, 1.0148314*ggHToZZ_pdf_m2, 0.3392106*ggToZZ_pdf_m2)");

  // bandwidth vary dn
  w->factory("SUM::HighMassFit2018_m1_dn(1.2594848*DY1J_pdf_m1_dn, 2.1432564*DY2J_pdf_m1_dn, 4.1874639*qqToZZ_pdf_m1_dn, 2.7497135*TTJets_pdf_m1_dn, 1.0148314*ggHToZZ_pdf_m1_dn, 0.3392106*ggToZZ_pdf_m1_dn)");
  w->factory("SUM::HighMassFit2018_m2_dn(1.2594848*DY1J_pdf_m2_dn, 2.1432564*DY2J_pdf_m2_dn, 4.1874639*qqToZZ_pdf_m2_dn, 2.7497135*TTJets_pdf_m2_dn, 1.0148314*ggHToZZ_pdf_m2_dn, 0.3392106*ggToZZ_pdf_m2_dn)");
  // bandwidth vary up
  w->factory("SUM::HighMassFit2018_m1_up(1.2594848*DY1J_pdf_m1_up, 2.1432564*DY2J_pdf_m1_up, 4.1874639*qqToZZ_pdf_m1_up, 2.7497135*TTJets_pdf_m1_up, 1.0148314*ggHToZZ_pdf_m1_up, 0.3392106*ggToZZ_pdf_m1_up)");
  w->factory("SUM::HighMassFit2018_m2_up(1.2594848*DY1J_pdf_m2_up, 2.1432564*DY2J_pdf_m2_up, 4.1874639*qqToZZ_pdf_m2_up, 2.7497135*TTJets_pdf_m2_up, 1.0148314*ggHToZZ_pdf_m2_up, 0.3392106*ggToZZ_pdf_m2_up)");

  // PDF nominal weight + 1 sigma
  w->factory("SUM::HighMassFit2018_m1_weight_up(2.1500751*DY1J_pdf_m1, 3.0182371*DY2J_pdf_m1, 4.3194877*qqToZZ_pdf_m1, 3.2996563*TTJets_pdf_m1, 1.0426480*ggHToZZ_pdf_m1, 0.3453498*ggToZZ_pdf_m1)");
  w->factory("SUM::HighMassFit2018_m2_weight_up(2.1500751*DY1J_pdf_m2, 3.0182371*DY2J_pdf_m2, 4.3194877*qqToZZ_pdf_m2, 3.2996563*TTJets_pdf_m2, 1.0426480*ggHToZZ_pdf_m2, 0.3453498*ggToZZ_pdf_m2)");
  // PDF nominal weight - 1 sigma
  w->factory("SUM::HighMassFit2018_m1_weight_dn(0.3688946*DY1J_pdf_m1, 1.2682756*DY2J_pdf_m1, 4.0554402*qqToZZ_pdf_m1, 2.1997709*TTJets_pdf_m1, 0.9870147*ggHToZZ_pdf_m1, 0.3330715*ggToZZ_pdf_m1)");
  w->factory("SUM::HighMassFit2018_m2_weight_dn(0.3688946*DY1J_pdf_m2, 1.2682756*DY2J_pdf_m2, 4.0554402*qqToZZ_pdf_m2, 2.1997709*TTJets_pdf_m2, 0.9870147*ggHToZZ_pdf_m2, 0.3330715*ggToZZ_pdf_m2)");
  // PDF Braid I (nominal wgt +/- 1 sigma, start from +)
  w->factory("SUM::HighMassFit2018_m1_weight_braidI(2.1500751*DY1J_pdf_m1, 1.2682756*DY2J_pdf_m1, 4.3194877*qqToZZ_pdf_m1, 2.1997709*TTJets_pdf_m1, 1.0426480*ggHToZZ_pdf_m1, 0.3330715*ggToZZ_pdf_m1)");
  w->factory("SUM::HighMassFit2018_m2_weight_braidI(2.1500751*DY1J_pdf_m2, 1.2682756*DY2J_pdf_m2, 4.3194877*qqToZZ_pdf_m2, 2.1997709*TTJets_pdf_m2, 1.0426480*ggHToZZ_pdf_m2, 0.3330715*ggToZZ_pdf_m2)");
  // PDF Braid II (nominal wgt -/+ 1 sigma, start from -)
  w->factory("SUM::HighMassFit2018_m1_weight_braidII(0.3688946*DY1J_pdf_m1, 3.0182371*DY2J_pdf_m1, 4.0554402*qqToZZ_pdf_m1, 3.2996563*TTJets_pdf_m1, 0.9870147*ggHToZZ_pdf_m1, 0.3453498*ggToZZ_pdf_m1)");
  w->factory("SUM::HighMassFit2018_m2_weight_braidII(0.3688946*DY1J_pdf_m2, 3.0182371*DY2J_pdf_m2, 4.0554402*qqToZZ_pdf_m2, 3.2996563*TTJets_pdf_m2, 0.9870147*ggHToZZ_pdf_m2, 0.3453498*ggToZZ_pdf_m2)");

  RooPlot* plotShapem1 = w->var("m1_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m1)"));
  RooPlot* plotShapem2 = w->var("m2_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m2)"));
  RooPlot* plotShapem1dn = w->var("m1_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m1 dn)"));
  RooPlot* plotShapem2dn = w->var("m2_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m2 dn)"));
  RooPlot* plotShapem1up = w->var("m1_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m1 up)"));
  RooPlot* plotShapem2up = w->var("m2_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m2 up)"));
  RooPlot* plotShapem1All = w->var("m1_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m1 dn + up)"));
  RooPlot* plotShapem2All = w->var("m2_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m2 dn + up)"));
  RooPlot* plotShapem1PDFWeight = w->var("m1_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m1 PDF weight)"));
  RooPlot* plotShapem2PDFWeight = w->var("m2_above_Upsilon")->frame(Title("Background PDF above 11 GeV (m2 PDF weight)"));

  // plot individually
  w->pdf("HighMassFit2018_m1")->plotOn(plotShapem1, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m1"));
  w->pdf("HighMassFit2018_m2")->plotOn(plotShapem2, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m2"));
  w->pdf("HighMassFit2018_m1_dn")->plotOn(plotShapem1dn, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m1_dn"));
  w->pdf("HighMassFit2018_m2_dn")->plotOn(plotShapem2dn, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m2_dn"));
  w->pdf("HighMassFit2018_m1_up")->plotOn(plotShapem1up, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m1_up"));
  w->pdf("HighMassFit2018_m2_up")->plotOn(plotShapem2up, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m2_up"));

  //overlay
  w->pdf("HighMassFit2018_m1")->plotOn(plotShapem1All, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m1"));
  w->pdf("HighMassFit2018_m1_dn")->plotOn(plotShapem1All, LineColor(kRed), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m1_dn"));
  w->pdf("HighMassFit2018_m1_up")->plotOn(plotShapem1All, LineColor(kMagenta), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m1_up"));

  w->pdf("HighMassFit2018_m2")->plotOn(plotShapem2All, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m2"));
  w->pdf("HighMassFit2018_m2_dn")->plotOn(plotShapem2All, LineColor(kRed), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m2_dn"));
  w->pdf("HighMassFit2018_m2_up")->plotOn(plotShapem2All, LineColor(kMagenta), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m2_up"));

  w->pdf("HighMassFit2018_m1")->plotOn(plotShapem1PDFWeight, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m1"));
  w->pdf("HighMassFit2018_m1_weight_up")->plotOn(plotShapem1PDFWeight, LineColor(kRed), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m1_weight_up"));
  w->pdf("HighMassFit2018_m1_weight_dn")->plotOn(plotShapem1PDFWeight, LineColor(kMagenta), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m1_weight_dn"));
  w->pdf("HighMassFit2018_m1_weight_braidI")->plotOn(plotShapem1PDFWeight, LineColor(kGreen), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m1_weight_braidI"));
  w->pdf("HighMassFit2018_m1_weight_braidII")->plotOn(plotShapem1PDFWeight, LineColor(kCyan), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m1_weight_braidII"));

  w->pdf("HighMassFit2018_m2")->plotOn(plotShapem2PDFWeight, LineColor(kBlue), Precision(0.0001), Name("HighMassFit2018_m2"));
  w->pdf("HighMassFit2018_m2_weight_up")->plotOn(plotShapem2PDFWeight, LineColor(kRed), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m2_weight_up"));
  w->pdf("HighMassFit2018_m2_weight_dn")->plotOn(plotShapem2PDFWeight, LineColor(kMagenta), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m2_weight_dn"));
  w->pdf("HighMassFit2018_m2_weight_braidI")->plotOn(plotShapem2PDFWeight, LineColor(kGreen), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m2_weight_braidI"));
  w->pdf("HighMassFit2018_m2_weight_braidII")->plotOn(plotShapem2PDFWeight, LineColor(kCyan), LineStyle(kDashed), Precision(0.0001), Name("HighMassFit2018_m2_weight_braidII"));

  TCanvas * c_shape1D_above_Upsilon = new TCanvas("c_shape1D_above_Upsilon", "c_shape1D_above_Upsilon", 800, 800);
  c_shape1D_above_Upsilon->cd();
  plotShapem1->Draw();   plotShapem1->GetYaxis()->SetTitle("p.d.f.");   c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon.pdf");
  plotShapem2->Draw();   plotShapem2->GetYaxis()->SetTitle("p.d.f.");   c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon.pdf");
  plotShapem1dn->Draw(); plotShapem1dn->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_dn.pdf");
  plotShapem2dn->Draw(); plotShapem2dn->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_dn.pdf");
  plotShapem1up->Draw(); plotShapem1up->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_up.pdf");
  plotShapem2up->Draw(); plotShapem2up->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_up.pdf");
  plotShapem1All->Draw(); plotShapem1All->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_All.pdf");
  plotShapem2All->Draw(); plotShapem2All->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_All.pdf");
  plotShapem1PDFWeight->Draw(); plotShapem1PDFWeight->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_PDFWeight.pdf");
  plotShapem2PDFWeight->Draw(); plotShapem2PDFWeight->GetYaxis()->SetTitle("p.d.f."); c_shape1D_above_Upsilon->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_PDFWeight.pdf");

  // Draw systematic unc. for m1: nom + scale bandwidth up/dn
  TCanvas * c_m1_above_Upsilon_sys = new TCanvas("c_m1_above_Upsilon_sys", "c_m1_above_Upsilon_sys", 800, 800);
  c_m1_above_Upsilon_sys->Clear();
  // Upper pad: pdfs: nom, up, dn
  TPad *padup_m1_above_Upsilon_sys = new TPad("padup_m1_above_Upsilon_sys", "padup_m1_above_Upsilon_sys", 0, 0.3, 1, 1.0);
  padup_m1_above_Upsilon_sys->SetBottomMargin(0); padup_m1_above_Upsilon_sys->Draw(); padup_m1_above_Upsilon_sys->cd();
  plotShapem1All->Draw(); plotShapem1All->GetYaxis()->SetTitle("Background PDF");
  c_m1_above_Upsilon_sys->cd(); c_m1_above_Upsilon_sys->Update();
  // Lower pad: (up/dn-nom)/nom; quote largest abs as sys. unc.
  TPad *paddn_m1_above_Upsilon_sys = new TPad("paddn_m1_above_Upsilon_sys", "paddn_m1_above_Upsilon_sys", 0, 0.0, 1, 0.29);
  paddn_m1_above_Upsilon_sys->SetTopMargin(0); paddn_m1_above_Upsilon_sys->SetBottomMargin(0.35); paddn_m1_above_Upsilon_sys->SetGridy(); paddn_m1_above_Upsilon_sys->Draw(); paddn_m1_above_Upsilon_sys->cd();
  TH1F *h_sys_m1_above_Upsilon_dn = new TH1F("h_sys_m1_above_Upsilon_dn", "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m1_above_Upsilon_up = new TH1F("h_sys_m1_above_Upsilon_up", "", 500, m_Upsilon_up, m_highmax);
  h_sys_m1_above_Upsilon_up->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]"); h_sys_m1_above_Upsilon_up->GetXaxis()->SetTitleSize(20); h_sys_m1_above_Upsilon_up->GetXaxis()->SetTitleFont(43); h_sys_m1_above_Upsilon_up->GetXaxis()->SetTitleOffset(3.0); h_sys_m1_above_Upsilon_up->GetXaxis()->SetLabelSize(15); h_sys_m1_above_Upsilon_up->GetXaxis()->SetLabelFont(43);
  h_sys_m1_above_Upsilon_up->GetYaxis()->SetTitle("Sys. unc."); h_sys_m1_above_Upsilon_up->GetYaxis()->SetNdivisions(505); h_sys_m1_above_Upsilon_up->GetYaxis()->CenterTitle(); h_sys_m1_above_Upsilon_up->GetYaxis()->SetTitleSize(20); h_sys_m1_above_Upsilon_up->GetYaxis()->SetTitleFont(43); h_sys_m1_above_Upsilon_up->GetYaxis()->SetTitleOffset(1.4); h_sys_m1_above_Upsilon_up->GetYaxis()->SetLabelSize(15); h_sys_m1_above_Upsilon_up->GetYaxis()->SetLabelFont(43);
  //histogram pdf: nom, dn, up
  TH1F* h_m1_above_Upsilon = new TH1F("h_m1_above_Upsilon", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m1_above_Upsilon_dn = new TH1F("h_m1_above_Upsilon_dn", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m1_above_Upsilon_up = new TH1F("h_m1_above_Upsilon_up", "", 500, m_Upsilon_up, m_highmax);
  w->pdf("HighMassFit2018_m1")->fillHistogram(h_m1_above_Upsilon, m1_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m1_dn")->fillHistogram(h_m1_above_Upsilon_dn, m1_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m1_up")->fillHistogram(h_m1_above_Upsilon_up, m1_above_Upsilon, 5000000);
  for (unsigned int iB=1; iB < 500; iB++) {
    float unc_dn = ( h_m1_above_Upsilon_dn->GetBinContent(iB) - h_m1_above_Upsilon->GetBinContent(iB) )*1.0/h_m1_above_Upsilon->GetBinContent(iB);
    float unc_up = ( h_m1_above_Upsilon_up->GetBinContent(iB) - h_m1_above_Upsilon->GetBinContent(iB) )*1.0/h_m1_above_Upsilon->GetBinContent(iB);
    h_sys_m1_above_Upsilon_dn->SetBinContent(iB, unc_dn );
    h_sys_m1_above_Upsilon_up->SetBinContent(iB, unc_up );
  }
  h_sys_m1_above_Upsilon_up->SetMinimum(shape_unc_plot_min); h_sys_m1_above_Upsilon_up->SetMaximum(shape_unc_plot_max); h_sys_m1_above_Upsilon_up->SetMarkerStyle(21); h_sys_m1_above_Upsilon_up->SetMarkerSize(0.2); h_sys_m1_above_Upsilon_up->SetMarkerColor(kMagenta); h_sys_m1_above_Upsilon_up->SetStats(0); h_sys_m1_above_Upsilon_up->Draw("p");
  h_sys_m1_above_Upsilon_dn->SetMinimum(shape_unc_plot_min); h_sys_m1_above_Upsilon_dn->SetMaximum(shape_unc_plot_max); h_sys_m1_above_Upsilon_dn->SetMarkerStyle(21); h_sys_m1_above_Upsilon_dn->SetMarkerSize(0.2); h_sys_m1_above_Upsilon_dn->SetMarkerColor(kRed); h_sys_m1_above_Upsilon_dn->SetStats(0); h_sys_m1_above_Upsilon_dn->Draw("p same");
  c_m1_above_Upsilon_sys->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_All_and_SysUnc.pdf");
  c_m1_above_Upsilon_sys->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_All_and_SysUnc.root");

  // Draw systematic unc. for m2: nom + scale bandwidth up/dn
  TCanvas * c_m2_above_Upsilon_sys = new TCanvas("c_m2_above_Upsilon_sys", "c_m2_above_Upsilon_sys", 800, 800);
  c_m2_above_Upsilon_sys->Clear();
  // Upper pad: pdfs: nom, up, dn
  TPad *padup_m2_above_Upsilon_sys = new TPad("padup_m2_above_Upsilon_sys", "padup_m2_above_Upsilon_sys", 0, 0.3, 1, 1.0);
  padup_m2_above_Upsilon_sys->SetBottomMargin(0); padup_m2_above_Upsilon_sys->Draw(); padup_m2_above_Upsilon_sys->cd();
  plotShapem2All->Draw(); plotShapem2All->GetYaxis()->SetTitle("Background PDF");
  c_m2_above_Upsilon_sys->cd(); c_m2_above_Upsilon_sys->Update();
  // Lower pad: (up/dn-nom)/nom; quote largest abs as sys. unc.
  TPad *paddn_m2_above_Upsilon_sys = new TPad("paddn_m2_above_Upsilon_sys", "paddn_m2_above_Upsilon_sys", 0, 0.0, 1, 0.29);
  paddn_m2_above_Upsilon_sys->SetTopMargin(0); paddn_m2_above_Upsilon_sys->SetBottomMargin(0.35); paddn_m2_above_Upsilon_sys->SetGridy(); paddn_m2_above_Upsilon_sys->Draw(); paddn_m2_above_Upsilon_sys->cd();
  TH1F *h_sys_m2_above_Upsilon_dn = new TH1F("h_sys_m2_above_Upsilon_dn", "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m2_above_Upsilon_up = new TH1F("h_sys_m2_above_Upsilon_up", "", 500, m_Upsilon_up, m_highmax);
  h_sys_m2_above_Upsilon_up->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]"); h_sys_m2_above_Upsilon_up->GetXaxis()->SetTitleSize(20); h_sys_m2_above_Upsilon_up->GetXaxis()->SetTitleFont(43); h_sys_m2_above_Upsilon_up->GetXaxis()->SetTitleOffset(3.0); h_sys_m2_above_Upsilon_up->GetXaxis()->SetLabelSize(15); h_sys_m2_above_Upsilon_up->GetXaxis()->SetLabelFont(43);
  h_sys_m2_above_Upsilon_up->GetYaxis()->SetTitle("Sys. unc."); h_sys_m2_above_Upsilon_up->GetYaxis()->SetNdivisions(505); h_sys_m2_above_Upsilon_up->GetYaxis()->CenterTitle(); h_sys_m2_above_Upsilon_up->GetYaxis()->SetTitleSize(20); h_sys_m2_above_Upsilon_up->GetYaxis()->SetTitleFont(43); h_sys_m2_above_Upsilon_up->GetYaxis()->SetTitleOffset(1.4); h_sys_m2_above_Upsilon_up->GetYaxis()->SetLabelSize(15); h_sys_m2_above_Upsilon_up->GetYaxis()->SetLabelFont(43);
  //histogram pdf: nom, dn, up
  TH1F* h_m2_above_Upsilon = new TH1F("h_m2_above_Upsilon","", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m2_above_Upsilon_dn = new TH1F("h_m2_above_Upsilon_dn", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m2_above_Upsilon_up = new TH1F("h_m2_above_Upsilon_up", "", 500, m_Upsilon_up, m_highmax);
  w->pdf("HighMassFit2018_m2")->fillHistogram(h_m2_above_Upsilon, m2_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m2_dn")->fillHistogram(h_m2_above_Upsilon_dn, m2_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m2_up")->fillHistogram(h_m2_above_Upsilon_up, m2_above_Upsilon, 5000000);
  for (unsigned int iB=1; iB < 500; iB++) {
    float unc_dn = ( h_m2_above_Upsilon_dn->GetBinContent(iB) - h_m2_above_Upsilon->GetBinContent(iB) )*1.0/h_m2_above_Upsilon->GetBinContent(iB);
    float unc_up = ( h_m2_above_Upsilon_up->GetBinContent(iB) - h_m2_above_Upsilon->GetBinContent(iB) )*1.0/h_m2_above_Upsilon->GetBinContent(iB);
    h_sys_m2_above_Upsilon_dn->SetBinContent(iB, unc_dn );
    h_sys_m2_above_Upsilon_up->SetBinContent(iB, unc_up );
  }
  h_sys_m2_above_Upsilon_up->SetMinimum(shape_unc_plot_min); h_sys_m2_above_Upsilon_up->SetMaximum(shape_unc_plot_max); h_sys_m2_above_Upsilon_up->SetMarkerStyle(21); h_sys_m2_above_Upsilon_up->SetMarkerSize(0.2); h_sys_m2_above_Upsilon_up->SetMarkerColor(kMagenta); h_sys_m2_above_Upsilon_up->SetStats(0); h_sys_m2_above_Upsilon_up->Draw("p");
  h_sys_m2_above_Upsilon_dn->SetMinimum(shape_unc_plot_min); h_sys_m2_above_Upsilon_dn->SetMaximum(shape_unc_plot_max); h_sys_m2_above_Upsilon_dn->SetMarkerStyle(21); h_sys_m2_above_Upsilon_dn->SetMarkerSize(0.2); h_sys_m2_above_Upsilon_dn->SetMarkerColor(kRed); h_sys_m2_above_Upsilon_dn->SetStats(0); h_sys_m2_above_Upsilon_dn->Draw("p same");
  c_m2_above_Upsilon_sys->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_All_and_SysUnc.pdf");
  c_m2_above_Upsilon_sys->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_All_and_SysUnc.root");

  // Draw systematic unc. for m1: MC PDF weight vary within +/-1sigma
  TCanvas * c_m1_above_Upsilon_sys_PDFWeight = new TCanvas("c_m1_above_Upsilon_sys_PDFWeight", "c_m1_above_Upsilon_sys_PDFWeight", 800, 800);
  c_m1_above_Upsilon_sys_PDFWeight->Clear();
  // Upper pad: pdfs with nom + 4 categories weight
  TPad *padup_m1_above_Upsilon_sys_PDFWeight = new TPad("padup_m1_above_Upsilon_sys_PDFWeight", "padup_m1_above_Upsilon_sys_PDFWeight", 0, 0.3, 1, 1.0);
  padup_m1_above_Upsilon_sys_PDFWeight->SetBottomMargin(0); padup_m1_above_Upsilon_sys_PDFWeight->Draw(); padup_m1_above_Upsilon_sys_PDFWeight->cd();
  plotShapem1PDFWeight->Draw(); plotShapem1PDFWeight->GetYaxis()->SetTitle("Background PDF");
  c_m1_above_Upsilon_sys_PDFWeight->cd(); c_m1_above_Upsilon_sys_PDFWeight->Update();
  // Lower pad: (pdf weight varied-nom)/nom; quote largest abs as sys. unc.
  TPad *paddn_m1_above_Upsilon_sys_PDFWeight = new TPad("paddn_m1_above_Upsilon_sys_PDFWeight", "paddn_m1_above_Upsilon_sys_PDFWeight", 0, 0.0, 1, 0.29);
  paddn_m1_above_Upsilon_sys_PDFWeight->SetTopMargin(0); paddn_m1_above_Upsilon_sys_PDFWeight->SetBottomMargin(0.35); paddn_m1_above_Upsilon_sys_PDFWeight->SetGridy(); paddn_m1_above_Upsilon_sys_PDFWeight->Draw(); paddn_m1_above_Upsilon_sys_PDFWeight->cd();
  TH1F *h_sys_m1_above_Upsilon_weight_up      = new TH1F("h_sys_m1_above_Upsilon_weight_up",      "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m1_above_Upsilon_weight_dn      = new TH1F("h_sys_m1_above_Upsilon_weight_dn",      "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m1_above_Upsilon_weight_braidI  = new TH1F("h_sys_m1_above_Upsilon_weight_braidI",  "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m1_above_Upsilon_weight_braidII = new TH1F("h_sys_m1_above_Upsilon_weight_braidII", "", 500, m_Upsilon_up, m_highmax);
  h_sys_m1_above_Upsilon_weight_up->GetXaxis()->SetTitle("m_{#mu#mu1} [GeV]"); h_sys_m1_above_Upsilon_weight_up->GetXaxis()->SetTitleSize(20); h_sys_m1_above_Upsilon_weight_up->GetXaxis()->SetTitleFont(43); h_sys_m1_above_Upsilon_weight_up->GetXaxis()->SetTitleOffset(3.0); h_sys_m1_above_Upsilon_weight_up->GetXaxis()->SetLabelSize(15); h_sys_m1_above_Upsilon_weight_up->GetXaxis()->SetLabelFont(43);
  h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetTitle("Sys. unc."); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetNdivisions(505); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->CenterTitle(); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetTitleSize(20); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetTitleFont(43); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetTitleOffset(1.4); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetLabelSize(15); h_sys_m1_above_Upsilon_weight_up->GetYaxis()->SetLabelFont(43);
  //histogram pdf: nom + four weights varied
  TH1F* h_m1_above_Upsilon_wgt = new TH1F("h_m1_above_Upsilon_wgt", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m1_above_Upsilon_wgt_up = new TH1F("h_m1_above_Upsilon_wgt_up", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m1_above_Upsilon_wgt_dn = new TH1F("h_m1_above_Upsilon_wgt_dn", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m1_above_Upsilon_wgt_braidi = new TH1F("h_m1_above_Upsilon_wgt_braidi", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m1_above_Upsilon_wgt_braidii = new TH1F("h_m1_above_Upsilon_wgt_braidii", "", 500, m_Upsilon_up, m_highmax);
  w->pdf("HighMassFit2018_m1")->fillHistogram(h_m1_above_Upsilon_wgt, m1_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m1_weight_up")->fillHistogram(h_m1_above_Upsilon_wgt_up, m1_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m1_weight_dn")->fillHistogram(h_m1_above_Upsilon_wgt_dn, m1_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m1_weight_braidI")->fillHistogram(h_m1_above_Upsilon_wgt_braidi, m1_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m1_weight_braidII")->fillHistogram(h_m1_above_Upsilon_wgt_braidii, m1_above_Upsilon, 5000000);
  for (unsigned int iB=1; iB < 500; iB++) {
    float unc_wgt_up = ( h_m1_above_Upsilon_wgt_up->GetBinContent(iB) - h_m1_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m1_above_Upsilon_wgt->GetBinContent(iB);
    float unc_wgt_dn = ( h_m1_above_Upsilon_wgt_dn->GetBinContent(iB) - h_m1_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m1_above_Upsilon_wgt->GetBinContent(iB);
    float unc_wgt_braidi = ( h_m1_above_Upsilon_wgt_braidi->GetBinContent(iB) - h_m1_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m1_above_Upsilon_wgt->GetBinContent(iB);
    float unc_wgt_braidii = ( h_m1_above_Upsilon_wgt_braidii->GetBinContent(iB) - h_m1_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m1_above_Upsilon_wgt->GetBinContent(iB);
    h_sys_m1_above_Upsilon_weight_up->SetBinContent(iB, unc_wgt_up );
    h_sys_m1_above_Upsilon_weight_dn->SetBinContent(iB, unc_wgt_dn );
    h_sys_m1_above_Upsilon_weight_braidI->SetBinContent(iB, unc_wgt_braidi );
    h_sys_m1_above_Upsilon_weight_braidII->SetBinContent(iB, unc_wgt_braidii );
  }
  h_sys_m1_above_Upsilon_weight_up->SetMinimum(shape_unc_plot_min); h_sys_m1_above_Upsilon_weight_up->SetMaximum(shape_unc_plot_max); h_sys_m1_above_Upsilon_weight_up->SetMarkerStyle(21); h_sys_m1_above_Upsilon_weight_up->SetMarkerSize(0.2); h_sys_m1_above_Upsilon_weight_up->SetMarkerColor(kRed); h_sys_m1_above_Upsilon_weight_up->SetStats(0); h_sys_m1_above_Upsilon_weight_up->Draw("p");
  h_sys_m1_above_Upsilon_weight_dn->SetMinimum(shape_unc_plot_min); h_sys_m1_above_Upsilon_weight_dn->SetMaximum(shape_unc_plot_max); h_sys_m1_above_Upsilon_weight_dn->SetMarkerStyle(21); h_sys_m1_above_Upsilon_weight_dn->SetMarkerSize(0.2); h_sys_m1_above_Upsilon_weight_dn->SetMarkerColor(kMagenta); h_sys_m1_above_Upsilon_weight_dn->SetStats(0); h_sys_m1_above_Upsilon_weight_dn->Draw("p same");
  h_sys_m1_above_Upsilon_weight_braidI->SetMinimum(shape_unc_plot_min); h_sys_m1_above_Upsilon_weight_braidI->SetMaximum(shape_unc_plot_max); h_sys_m1_above_Upsilon_weight_braidI->SetMarkerStyle(21); h_sys_m1_above_Upsilon_weight_braidI->SetMarkerSize(0.2); h_sys_m1_above_Upsilon_weight_braidI->SetMarkerColor(kGreen); h_sys_m1_above_Upsilon_weight_braidI->SetStats(0); h_sys_m1_above_Upsilon_weight_braidI->Draw("p same");
  h_sys_m1_above_Upsilon_weight_braidII->SetMinimum(shape_unc_plot_min); h_sys_m1_above_Upsilon_weight_braidII->SetMaximum(shape_unc_plot_max); h_sys_m1_above_Upsilon_weight_braidII->SetMarkerStyle(21); h_sys_m1_above_Upsilon_weight_braidII->SetMarkerSize(0.2); h_sys_m1_above_Upsilon_weight_braidII->SetMarkerColor(kCyan); h_sys_m1_above_Upsilon_weight_braidII->SetStats(0); h_sys_m1_above_Upsilon_weight_braidII->Draw("p same");
  c_m1_above_Upsilon_sys_PDFWeight->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_PDFWeight_and_SysUnc.pdf");
  c_m1_above_Upsilon_sys_PDFWeight->SaveAs("HighMassShape/BKG_Shape1D_m1_above_Upsilon_PDFWeight_and_SysUnc.root");

  // Draw systematic unc. for m2: MC PDF weight vary within +/-1sigma
  TCanvas * c_m2_above_Upsilon_sys_PDFWeight = new TCanvas("c_m2_above_Upsilon_sys_PDFWeight", "c_m2_above_Upsilon_sys_PDFWeight", 800, 800);
  c_m2_above_Upsilon_sys_PDFWeight->Clear();
  // Upper pad: pdfs with nom + 4 categories weight
  TPad *padup_m2_above_Upsilon_sys_PDFWeight = new TPad("padup_m2_above_Upsilon_sys_PDFWeight", "padup_m2_above_Upsilon_sys_PDFWeight", 0, 0.3, 1, 1.0);
  padup_m2_above_Upsilon_sys_PDFWeight->SetBottomMargin(0); padup_m2_above_Upsilon_sys_PDFWeight->Draw(); padup_m2_above_Upsilon_sys_PDFWeight->cd();
  plotShapem2PDFWeight->Draw(); plotShapem2PDFWeight->GetYaxis()->SetTitle("Background PDF");
  c_m2_above_Upsilon_sys_PDFWeight->cd(); c_m2_above_Upsilon_sys_PDFWeight->Update();
  // Lower pad: (pdf weight varied-nom)/nom; quote largest abs as sys. unc.
  TPad *paddn_m2_above_Upsilon_sys_PDFWeight = new TPad("paddn_m2_above_Upsilon_sys_PDFWeight", "paddn_m2_above_Upsilon_sys_PDFWeight", 0, 0.0, 1, 0.29);
  paddn_m2_above_Upsilon_sys_PDFWeight->SetTopMargin(0); paddn_m2_above_Upsilon_sys_PDFWeight->SetBottomMargin(0.35); paddn_m2_above_Upsilon_sys_PDFWeight->SetGridy(); paddn_m2_above_Upsilon_sys_PDFWeight->Draw(); paddn_m2_above_Upsilon_sys_PDFWeight->cd();
  TH1F *h_sys_m2_above_Upsilon_weight_up      = new TH1F("h_sys_m2_above_Upsilon_weight_up",      "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m2_above_Upsilon_weight_dn      = new TH1F("h_sys_m2_above_Upsilon_weight_dn",      "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m2_above_Upsilon_weight_braidI  = new TH1F("h_sys_m2_above_Upsilon_weight_braidI",  "", 500, m_Upsilon_up, m_highmax);
  TH1F *h_sys_m2_above_Upsilon_weight_braidII = new TH1F("h_sys_m2_above_Upsilon_weight_braidII", "", 500, m_Upsilon_up, m_highmax);
  h_sys_m2_above_Upsilon_weight_up->GetXaxis()->SetTitle("m_{#mu#mu2} [GeV]"); h_sys_m2_above_Upsilon_weight_up->GetXaxis()->SetTitleSize(20); h_sys_m2_above_Upsilon_weight_up->GetXaxis()->SetTitleFont(43); h_sys_m2_above_Upsilon_weight_up->GetXaxis()->SetTitleOffset(3.0); h_sys_m2_above_Upsilon_weight_up->GetXaxis()->SetLabelSize(15); h_sys_m2_above_Upsilon_weight_up->GetXaxis()->SetLabelFont(43);
  h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetTitle("Sys. unc."); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetNdivisions(505); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->CenterTitle(); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetTitleSize(20); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetTitleFont(43); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetTitleOffset(1.4); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetLabelSize(15); h_sys_m2_above_Upsilon_weight_up->GetYaxis()->SetLabelFont(43);
  //histogram pdf: nom + four weights varied
  TH1F* h_m2_above_Upsilon_wgt = new TH1F("h_m2_above_Upsilon_wgt", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m2_above_Upsilon_wgt_up = new TH1F("h_m2_above_Upsilon_wgt_up", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m2_above_Upsilon_wgt_dn = new TH1F("h_m2_above_Upsilon_wgt_dn", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m2_above_Upsilon_wgt_braidi = new TH1F("h_m2_above_Upsilon_wgt_braidi", "", 500, m_Upsilon_up, m_highmax);
  TH1F* h_m2_above_Upsilon_wgt_braidii = new TH1F("h_m2_above_Upsilon_wgt_braidii", "", 500, m_Upsilon_up, m_highmax);
  w->pdf("HighMassFit2018_m2")->fillHistogram(h_m2_above_Upsilon_wgt, m2_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m2_weight_up")->fillHistogram(h_m2_above_Upsilon_wgt_up, m2_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m2_weight_dn")->fillHistogram(h_m2_above_Upsilon_wgt_dn, m2_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m2_weight_braidI")->fillHistogram(h_m2_above_Upsilon_wgt_braidi, m2_above_Upsilon, 5000000);
  w->pdf("HighMassFit2018_m2_weight_braidII")->fillHistogram(h_m2_above_Upsilon_wgt_braidii, m2_above_Upsilon, 5000000);
  for (unsigned int iB=1; iB < 500; iB++) {
    float unc_wgt_up = ( h_m2_above_Upsilon_wgt_up->GetBinContent(iB) - h_m2_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m2_above_Upsilon_wgt->GetBinContent(iB);
    float unc_wgt_dn = ( h_m2_above_Upsilon_wgt_dn->GetBinContent(iB) - h_m2_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m2_above_Upsilon_wgt->GetBinContent(iB);
    float unc_wgt_braidi = ( h_m2_above_Upsilon_wgt_braidi->GetBinContent(iB) - h_m2_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m2_above_Upsilon_wgt->GetBinContent(iB);
    float unc_wgt_braidii = ( h_m2_above_Upsilon_wgt_braidii->GetBinContent(iB) - h_m2_above_Upsilon_wgt->GetBinContent(iB) )*1.0/h_m2_above_Upsilon_wgt->GetBinContent(iB);
    h_sys_m2_above_Upsilon_weight_up->SetBinContent(iB, unc_wgt_up );
    h_sys_m2_above_Upsilon_weight_dn->SetBinContent(iB, unc_wgt_dn );
    h_sys_m2_above_Upsilon_weight_braidI->SetBinContent(iB, unc_wgt_braidi );
    h_sys_m2_above_Upsilon_weight_braidII->SetBinContent(iB, unc_wgt_braidii );
  }
  h_sys_m2_above_Upsilon_weight_up->SetMinimum(shape_unc_plot_min); h_sys_m2_above_Upsilon_weight_up->SetMaximum(shape_unc_plot_max); h_sys_m2_above_Upsilon_weight_up->SetMarkerStyle(21); h_sys_m2_above_Upsilon_weight_up->SetMarkerSize(0.2); h_sys_m2_above_Upsilon_weight_up->SetMarkerColor(kRed); h_sys_m2_above_Upsilon_weight_up->SetStats(0); h_sys_m2_above_Upsilon_weight_up->Draw("p");
  h_sys_m2_above_Upsilon_weight_dn->SetMinimum(shape_unc_plot_min); h_sys_m2_above_Upsilon_weight_dn->SetMaximum(shape_unc_plot_max); h_sys_m2_above_Upsilon_weight_dn->SetMarkerStyle(21); h_sys_m2_above_Upsilon_weight_dn->SetMarkerSize(0.2); h_sys_m2_above_Upsilon_weight_dn->SetMarkerColor(kMagenta); h_sys_m2_above_Upsilon_weight_dn->SetStats(0); h_sys_m2_above_Upsilon_weight_dn->Draw("p same");
  h_sys_m2_above_Upsilon_weight_braidI->SetMinimum(shape_unc_plot_min); h_sys_m2_above_Upsilon_weight_braidI->SetMaximum(shape_unc_plot_max); h_sys_m2_above_Upsilon_weight_braidI->SetMarkerStyle(21); h_sys_m2_above_Upsilon_weight_braidI->SetMarkerSize(0.2); h_sys_m2_above_Upsilon_weight_braidI->SetMarkerColor(kGreen); h_sys_m2_above_Upsilon_weight_braidI->SetStats(0); h_sys_m2_above_Upsilon_weight_braidI->Draw("p same");
  h_sys_m2_above_Upsilon_weight_braidII->SetMinimum(shape_unc_plot_min); h_sys_m2_above_Upsilon_weight_braidII->SetMaximum(shape_unc_plot_max); h_sys_m2_above_Upsilon_weight_braidII->SetMarkerStyle(21); h_sys_m2_above_Upsilon_weight_braidII->SetMarkerSize(0.2); h_sys_m2_above_Upsilon_weight_braidII->SetMarkerColor(kCyan); h_sys_m2_above_Upsilon_weight_braidII->SetStats(0); h_sys_m2_above_Upsilon_weight_braidII->Draw("p same");
  c_m2_above_Upsilon_sys_PDFWeight->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_PDFWeight_and_SysUnc.pdf");
  c_m2_above_Upsilon_sys_PDFWeight->SaveAs("HighMassShape/BKG_Shape1D_m2_above_Upsilon_PDFWeight_and_SysUnc.root");

  w->factory("PROD::HighMassBKG(HighMassFit2018_m1, HighMassFit2018_m2)");
  // scale bandwidth up and down
  w->factory("PROD::HighMassBKG_dn(HighMassFit2018_m1_dn, HighMassFit2018_m2_dn)");
  w->factory("PROD::HighMassBKG_up(HighMassFit2018_m1_up, HighMassFit2018_m2_up)");
  // vary PDF weight +/1 sigma
  w->factory("PROD::HighMassBKG_weight_up(HighMassFit2018_m1_weight_up, HighMassFit2018_m2_weight_up)");
  w->factory("PROD::HighMassBKG_weight_dn(HighMassFit2018_m1_weight_dn, HighMassFit2018_m2_weight_dn)");
  w->factory("PROD::HighMassBKG_weight_braidI(HighMassFit2018_m1_weight_braidI, HighMassFit2018_m2_weight_braidI)");
  w->factory("PROD::HighMassBKG_weight_braidII(HighMassFit2018_m1_weight_braidII, HighMassFit2018_m2_weight_braidII)");

  // assume a flat background pdf for above 11 GeV
  RooPolynomial flatm1("flatm1", "flatm1", m1_above_Upsilon);
  RooPolynomial flatm2("flatm2", "flatm2", m2_above_Upsilon);
  w->import(flatm1);
  w->import(flatm2);
  w->factory("PROD::flat2D(flatm1,flatm2)");

  //=***************************************************************************
  //                           Save to Workspace
  //=***************************************************************************
  cout<<"Save to workspace"<<endl;
  w->writeToFile(outFileLM);
}
