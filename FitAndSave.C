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

void FitAndSave() {
  //Parameteres
  bool useTrig=true;
  TString iso_cut= "2";
  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  //Style
  setTDRStyle();
  TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader("CMS Prelim. 2017C-F  #sqrt{s} = 13 TeV   L_{int} = 36.734 fb^{-1}");
  //Output ws
  RooWorkspace* w = new RooWorkspace("w");
  TString Comm = "mkdir -p figures/";
  system( Comm.Data() );
  //Input file
  TChain chain_data_dimudimu("cutFlowAnalyzerPXBL3PXFL2/Events");
  TChain chain_data_dimuorphan("cutFlowAnalyzerPXBL3PXFL2/Events_orphan");
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
  RooRealVar m1("m1","m_{#mu#mu_{1}}",m_min,m_max,"GeV/#it{c}^{2}");
  RooRealVar m2("m2","m_{#mu#mu_{2}}",m_min,m_max,"GeV/#it{c}^{2}");
  m1.setBins(m_bins);
  m2.setBins(m_bins);
  w->import(m1);
  w->import(m2);
  //Selection bb Control Region
  ostringstream stream_cut_bg_m1_iso;
  if(useTrig) stream_cut_bg_m1_iso << "orph_dimu_isoTk < " << iso_cut << " && orph_dimu_isoTk >= 0 && containstrig2 > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  else        stream_cut_bg_m1_iso << "orph_dimu_isoTk < " << iso_cut << " && orph_dimu_isoTk >= 0 && ((orph_PtMu0>17 && TMath::Abs(orph_EtaMu0<0.9)) || (orph_PtMu1>17 && TMath::Abs(orph_EtaMu1<0.9))) && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  ostringstream stream_cut_bg_m2_iso;
  if(useTrig) stream_cut_bg_m2_iso << "orph_dimu_isoTk < " << iso_cut << " && orph_dimu_isoTk >= 0 && containstrig > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  else        stream_cut_bg_m2_iso << "orph_dimu_isoTk < " << iso_cut << " && orph_dimu_isoTk >= 0 && orph_PtOrph>17 && TMath::Abs(orph_EtaOrph)<0.9 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso = stream_cut_bg_m1_iso.str();
  TString cut_bg_m2_iso = stream_cut_bg_m2_iso.str();
  //Selection Signal
  TString cut_diagonal                = "abs(massC-massF) <= (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 9. && massF > 0.25 && massF < 9.";
  TString cut_control_offDiagonal     = "abs(massC-massF) > (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 9. && massF > 0.25 && massF < 9.";
  TString cut_control_Iso_offDiagonal = "isoC_1mm >= 0 && isoC_1mm < 2. && isoF_1mm >= 0 && isoF_1mm < 2. && abs(massC-massF) > (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 9. && massF > 0.25 && massF < 9.";
  TString cut_signal                  = "isoC_1mm>=0 && isoC_1mm<2. && isoF_1mm>=0 && isoF_1mm<2. && abs(massC-massF) <= (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 9. && massF > 0.25 && massF < 9.";
  //TString cut_control_nonIso          = "isoC_1mm > 2. && isoC_1mm < 8. && isoF_1mm > 2. && isoF_1mm < 8. && massC > 0.25 && massC < 9. && massF > 0.25 && massF < 9.";
  //TTree: bb Control Region
  TTree* tree_dimuorphan_bg_m1                          = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso);
  TTree* tree_dimuorphan_bg_m2                          = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso);
  //TTree: not used in this macro, to be used in other macros (PlotSignal_and_Background and PlotBBbar)
  TTree* tree_dimudimu_diagonal_2D                      = chain_data_dimudimu.CopyTree(cut_diagonal);
  TTree* tree_dimudimu_diagonal_1D_massC                = chain_data_dimudimu.CopyTree(cut_diagonal);
  TTree* tree_dimudimu_diagonal_1D_massF                = chain_data_dimudimu.CopyTree(cut_diagonal);
  TTree* tree_dimudimu_control_offDiagonal_2D           = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  TTree* tree_dimudimu_control_offDiagonal_1D_massC     = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  TTree* tree_dimudimu_control_offDiagonal_1D_massF     = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D       = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);
  TTree* tree_dimudimu_control_Iso_offDiagonal_1D_massC = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);
  TTree* tree_dimudimu_control_Iso_offDiagonal_1D_massF = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);
  TTree* tree_dimudimu_signal_2D                        = chain_data_dimudimu.CopyTree(cut_signal);
  //TTree* tree_dimudimu_control_nonIso                   = chain_data_dimudimu.CopyTree(cut_control_nonIso);

  //cout<<"------OffDiagonal SCAN------"<<endl;
  //tree_dimudimu_control_Iso_offDiagonal_2D->Scan("massC:massF:isoC_1mm:isoF_1mm:run:event:lumi");
  //cout<<"------Signal SCAN------"<<endl;
  //tree_dimudimu_signal_2D->Scan("massC:massF:run:event:lumi:isoC_1mm:isoF_1mm");

  //Setting Names bb Control region
  tree_dimuorphan_bg_m1->GetBranch("orph_dimu_mass")->SetName("m1");
  tree_dimuorphan_bg_m2->GetBranch("orph_dimu_mass")->SetName("m2");
  //Setting Names bb Signal
  tree_dimudimu_diagonal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_diagonal_2D->GetBranch("massF")->SetName("m2");
  tree_dimudimu_diagonal_1D_massC->GetBranch("massC")->SetName("m1");
  tree_dimudimu_diagonal_1D_massF->GetBranch("massF")->SetName("m1");
  tree_dimudimu_control_offDiagonal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_offDiagonal_2D->GetBranch("massF")->SetName("m2");
  tree_dimudimu_control_offDiagonal_1D_massC->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_offDiagonal_1D_massF->GetBranch("massF")->SetName("m1");
  tree_dimudimu_control_Iso_offDiagonal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_Iso_offDiagonal_2D->GetBranch("massF")->SetName("m2");
  tree_dimudimu_control_Iso_offDiagonal_1D_massC->GetBranch("massC")->SetName("m1");
  tree_dimudimu_control_Iso_offDiagonal_1D_massF->GetBranch("massF")->SetName("m1");
  tree_dimudimu_signal_2D->GetBranch("massC")->SetName("m1");
  tree_dimudimu_signal_2D->GetBranch("massF")->SetName("m2");
  //Creating dataset bb Control Region
  RooDataSet* ds_dimuorphan_bg_m1 = new RooDataSet("ds_dimuorphan_bg_m1","ds_dimuorphan_bg_m1", tree_dimuorphan_bg_m1, RooArgSet(m1));
  RooDataSet* ds_dimuorphan_bg_m2 = new RooDataSet("ds_dimuorphan_bg_m2","ds_dimuorphan_bg_m2", tree_dimuorphan_bg_m2, RooArgSet(m2));
  ////Creating dataset Signal
  RooDataSet* ds_dimudimu_diagonal_2D = new RooDataSet("ds_dimudimu_diagonal_2D","ds_dimudimu_diagonal_2D", tree_dimudimu_diagonal_2D, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_diagonal_1D_massC = new RooDataSet("ds_dimudimu_diagonal_1D_massC","ds_dimudimu_diagonal_1D_massC", tree_dimudimu_diagonal_1D_massC, RooArgSet(m1));
  RooDataSet* ds_dimudimu_diagonal_1D_massF = new RooDataSet("ds_dimudimu_diagonal_1D_massF","ds_dimudimu_diagonal_1D_massF", tree_dimudimu_diagonal_1D_massF, RooArgSet(m1));
  RooDataSet* ds_dimudimu_diagonal_1D = new RooDataSet("ds_dimudimu_diagonal_1D","ds_dimudimu_diagonal_1D", tree_dimudimu_diagonal_1D_massC, RooArgSet(m1));
  ds_dimudimu_diagonal_1D->append(*ds_dimudimu_diagonal_1D_massF);
  RooDataSet* ds_dimudimu_control_offDiagonal_2D = new RooDataSet("ds_dimudimu_control_offDiagonal_2D","ds_dimudimu_control_offDiagonal_2D", tree_dimudimu_control_offDiagonal_2D, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_control_offDiagonal_1D_massC = new RooDataSet("ds_dimudimu_control_offDiagonal_1D_massC","ds_dimudimu_control_offDiagonal_1D_massC", tree_dimudimu_control_offDiagonal_1D_massC, RooArgSet(m1));
  RooDataSet* ds_dimudimu_control_offDiagonal_1D_massF = new RooDataSet("ds_dimudimu_control_offDiagonal_1D_massF","ds_dimudimu_control_offDiagonal_1D_massF", tree_dimudimu_control_offDiagonal_1D_massF, RooArgSet(m1));
  RooDataSet* ds_dimudimu_control_offDiagonal_1D = new RooDataSet("ds_dimudimu_control_offDiagonal_1D","ds_dimudimu_control_offDiagonal_1D", tree_dimudimu_control_offDiagonal_1D_massC, RooArgSet(m1));
  ds_dimudimu_control_offDiagonal_1D->append(*ds_dimudimu_control_offDiagonal_1D_massF);
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_2D = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_2D","ds_dimudimu_control_Iso_offDiagonal_2D", tree_dimudimu_control_Iso_offDiagonal_2D, RooArgSet(m1,m2));
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_1D_massC = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_1D_massC","ds_dimudimu_control_Iso_offDiagonal_1D_massC", tree_dimudimu_control_Iso_offDiagonal_1D_massC, RooArgSet(m1));
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_1D_massF = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_1D_massF","ds_dimudimu_control_Iso_offDiagonal_1D_massF", tree_dimudimu_control_Iso_offDiagonal_1D_massF, RooArgSet(m1));
  RooDataSet* ds_dimudimu_control_Iso_offDiagonal_1D = new RooDataSet("ds_dimudimu_control_Iso_offDiagonal_1D","ds_dimudimu_control_Iso_offDiagonal_1D", tree_dimudimu_control_Iso_offDiagonal_1D_massC, RooArgSet(m1));
  ds_dimudimu_control_Iso_offDiagonal_1D->append(*ds_dimudimu_control_Iso_offDiagonal_1D_massF);
  RooDataSet* ds_dimudimu_signal_2D = new RooDataSet("ds_dimudimu_signal_2D","ds_dimudimu_signal_2D", tree_dimudimu_signal_2D, RooArgSet(m1,m2));

  cout<<"-----Now Printing and importing the Datasets:-----"<<endl;
  ds_dimuorphan_bg_m1->Print("s");
  ds_dimuorphan_bg_m2->Print("s");
  ds_dimudimu_diagonal_2D->Print("s");
  ds_dimudimu_diagonal_1D_massC->Print("s");
  ds_dimudimu_diagonal_1D_massF->Print("s");
  ds_dimudimu_diagonal_1D->Print("s");
  ds_dimudimu_control_offDiagonal_2D->Print("s");
  ds_dimudimu_control_offDiagonal_1D_massC->Print("s");
  ds_dimudimu_control_offDiagonal_1D_massF->Print("s");
  ds_dimudimu_control_offDiagonal_1D->Print("s");
  ds_dimudimu_control_Iso_offDiagonal_2D->Print("s");
  ds_dimudimu_control_Iso_offDiagonal_1D_massC->Print("s");
  ds_dimudimu_control_Iso_offDiagonal_1D_massF->Print("s");
  ds_dimudimu_control_Iso_offDiagonal_1D->Print("s");
  ds_dimudimu_signal_2D->Print("s");
  w->import(*ds_dimuorphan_bg_m1);
  w->import(*ds_dimuorphan_bg_m2);
  w->import(*ds_dimudimu_diagonal_2D);
  w->import(*ds_dimudimu_diagonal_1D_massC);
  w->import(*ds_dimudimu_diagonal_1D_massF);
  w->import(*ds_dimudimu_diagonal_1D);
  w->import(*ds_dimudimu_control_offDiagonal_2D);
  w->import(*ds_dimudimu_control_offDiagonal_1D_massC);
  w->import(*ds_dimudimu_control_offDiagonal_1D_massF);
  w->import(*ds_dimudimu_control_offDiagonal_1D);
  w->import(*ds_dimudimu_control_Iso_offDiagonal_2D);
  w->import(*ds_dimudimu_control_Iso_offDiagonal_1D_massC);
  w->import(*ds_dimudimu_control_Iso_offDiagonal_1D_massF);
  w->import(*ds_dimudimu_control_Iso_offDiagonal_1D);
  w->import(*ds_dimudimu_signal_2D);

  //Draw before fiting
  RooPlot* plotC1 = w->var("m1")->frame(Title("bb Tempale no FIT"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC1, DataError(RooAbsData::SumW2), Name("data_m1"));
  float SizeBin1 = plotC1->GetXaxis()->GetBinCenter(3) - plotC1->GetXaxis()->GetBinCenter(2);
  char c_SizeBin1[10];
  snprintf(c_SizeBin1,50,"%.2f",SizeBin1);
  TString Yname1 = "Events / (" + std::string(c_SizeBin1) + "[GeV])";
  plotC1->GetYaxis()->SetTitle( Yname1.Data() );
  TCanvas * c1 = new TCanvas("c1");
  c1->cd(); plotC1->Draw(); txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m1.pdf");
  RooPlot* plotC2 = w->var("m2")->frame(Title("bb data"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotC2, DataError(RooAbsData::SumW2), Name("data_m2"));
  plotC2->Draw();
  txtHeader->Draw();
  c1->SaveAs("figures/h_dimuorphan_bg_m2.pdf");
  delete c1;

  //****************************************************************************
  //                         Create template for m1
  //****************************************************************************
  cout<<"-----Creating templates:-----"<<endl;
  /*
  //===========
  //2016 fit m1
  //===========
  //Initial combianatorial
  //also very important, before fix to 0.5
  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.01, 0.0, 0.3], MmumuC_p[0.05, 0.0, 1.5])");
  //FIRST KINK
  w->factory("Bernstein::bgC(m1,{bC06[9.5766e+00,0.1,15.], bC16[ 1.0705e-01,0.,3.], bC26[3.5184e-05,0.,3.], bC36[1.259,0.,5.], bC46[2.5370e-03,0.,3.], bC56[5.7432e-01,0.,3.], bC66[3.9353e-01,0.1,4.]})");
  // Resonances
  w->factory("Gaussian::etaC(m1,0.548,0.033)");
  w->factory("Gaussian::rhoC(m1,0.782,0.036)");
  w->factory("Gaussian::phiC(m1,1.019,0.027)");
  w->factory("Gaussian::psiC(m1,3.7,psiC_sigma[0.033,0.001,0.1])");
  // Ad hoc gaussian to cover first bump and help other functions
  w->factory("Gaussian::adHocC(m1,adHocC_mass[0.2,0.,0.6],adHocC_sigma[6.3097e-02,0.0001,0.1])");
  // Visible Resonances
  w->factory("CBShape::JpsiC(m1, JpsiC_mean[3.12,3.0,3.35], JpsiC_sigma[0.1,0.001,0.3], JpsiC_alpha[1.2,0.4,7.0], JpsiC_n[2.0])");
  w->factory("Gaussian::Up1C(m1,Up1C_mean[9.43,9.39,9.500], Up1C_sigma[0.101,0.01,0.20])");
  w->factory("Gaussian::Up2C(m1,Up2C_mean[10.0,9.90,10.20], Up2C_sigma[0.060,0.01,0.15])");
  w->factory("Gaussian::Up3C(m1,Up3C_mean[10.45,10.2,10.7], Up3C_sigma[0.050,0.01,0.15])");
  // FINAL PDF
  w->factory("SUM::template1D_m1(norm_adHocC[20., 0., 10000.]*adHocC,norm_MmumuC[200., 0., 25000.]*MmumuC, norm_bgC[4400.,100.,20000.]*bgC, norm_etaC[1.3151e+01,0.,1000.]*etaC, norm_rhoC[1.0107e+02,0.,1000.]*rhoC, norm_phiC[9.8640e+01,0.,1000.]*phiC, norm_JpsiC[50.,10.,6000.]*JpsiC, norm_psiC[50.,0.,1000.]*psiC)");
  */

  //===========
  //2017 fit m1
  //===========
  //Initial combianatorial
  //also very important, before fix to 0.5
  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.01, 0.0, 0.3], MmumuC_p[0.05, 0.0, 1.5])");
  //FIRST KINK
  w->factory("Bernstein::bgC(m1,{bC06[9.5766e+00,0.1,15.], bC16[ 1.0705e-01,0.,3.], bC26[3.5184e-05,0.,3.], bC36[1.259,0.,5.], bC46[2.5370e-03,0.,3.], bC56[5.7432e-01,0.,3.], bC66[3.9353e-01,0.1,4.]})");
  // Resonances
  w->factory("Gaussian::etaC(m1,0.548,0.033)");
  w->factory("Gaussian::rhoC(m1,0.782,0.036)");
  w->factory("Gaussian::phiC(m1,1.019,0.027)");
  w->factory("Gaussian::psiC(m1,3.7,psiC_sigma[0.033,0.001,0.1])");
  // Ad hoc gaussian to cover first bump and help other functions
  w->factory("Gaussian::adHocC(m1,adHocC_mass[0.2,0.,0.6],adHocC_sigma[6.3097e-02,0.0001,0.1])");
  // Visible Resonances
  w->factory("CBShape::JpsiC(m1, JpsiC_mean[3.12,3.0,3.35], JpsiC_sigma[0.1,0.001,0.3], JpsiC_alpha[1.2,0.4,7.0], JpsiC_n[2.0])");
  w->factory("Gaussian::Up1C(m1,Up1C_mean[9.43,9.39,9.500], Up1C_sigma[0.101,0.01,0.20])");
  w->factory("Gaussian::Up2C(m1,Up2C_mean[10.0,9.90,10.20], Up2C_sigma[0.060,0.01,0.15])");
  w->factory("Gaussian::Up3C(m1,Up3C_mean[10.45,10.2,10.7], Up3C_sigma[0.050,0.01,0.15])");
  // FINAL PDF
  w->factory("SUM::template1D_m1(norm_adHocC[20., 0., 10000.]*adHocC,norm_MmumuC[200., 0., 25000.]*MmumuC, norm_bgC[4400.,100.,20000.]*bgC, norm_etaC[1.3151e+01,0.,1000.]*etaC, norm_rhoC[1.0107e+02,0.,1000.]*rhoC, norm_phiC[9.8640e+01,0.,1000.]*phiC, norm_JpsiC[50.,10.,6000.]*JpsiC, norm_psiC[50.,0.,1000.]*psiC)");

  RooFitResult *rC = w->pdf("template1D_m1")->fitTo(*(w->data("ds_dimuorphan_bg_m1")), Extended(1), Save(), SumW2Error(kTRUE));
  cout<<"------------------RooPrintable 1---------------------"<<endl;
  rC->Print();

  RooPlot* plotC = w->var("m1")->frame(Title("BG template for trigger dimuon"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC, DataError(RooAbsData::SumW2), Name("data_m1"));
  w->pdf("template1D_m1")->plotOn(plotC,LineColor(kRed),Precision(0.0001),Name("template1D_m1"));
  float SizeBin = plotC->GetXaxis()->GetBinCenter(3) - plotC->GetXaxis()->GetBinCenter(2);
  char c_SizeBin[10];
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  TString Yname = "Events / (" + std::string(c_SizeBin) + "[GeV])";
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
  h_dataFit1->SetMinimum(0.5);  // Define Y ..
  h_dataFit1->SetMaximum(1.5); // .. range
  h_dataFit1->Sumw2();
  h_dataFit1->SetStats(0);      // No statistics on lower plot
  h_dataFit1->SetMarkerStyle(21);
  h_dataFit1->Draw("ep");
  c_template1D_m1_RooPlot->SaveAs("figures/template1D_m1_RooPlot.pdf");
  c_template1D_m1_RooPlot->SaveAs("figures/template1D_m1_RooPlot.png");
  float chi2_C = plotC->chiSquare(23);

  //****************************************************************************
  //                         Create template for m2
  //****************************************************************************
  /*
  //===========
  //2016 fit m2
  //===========
  //Initial combianatorial
  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[0.01, 0.0, 0.3], MmumuF_p[0.05, 0.0, 2.])");
  //FIRST KINK
  w->factory("Bernstein::bgF(m2,{bF06[9.9751,0,15.], bF16[3.2971e-05,0.,3.], bF26[ 1.2361e-07,0.,3.], bF36[5.1545e-08,0.,2.], bF46[9.9017e-01,0.,3.], bF56[3.0607e-01,0.,3.], bF66[0.5,0.1,4.]})");
  // Resonances
  w->factory("Gaussian::etaF(m2,0.548,0.02)");
  w->factory("Gaussian::rhoF(m2,0.782,0.025)");
  w->factory("Gaussian::phiF(m2,1.019,0.02)");
  w->factory("Gaussian::psiF(m2,3.675,psiF_sigma[0.033,0.001,0.1])");
  // Ad hoc gaussian to cover first bump and help other functions
  w->factory("Gaussian::adHocF(m2,adHocF_mass[0.4,0.2,0.6],adHocF_sigma[0.01,0.001,0.1])");
  // Visible Resonances
  w->factory("CBShape::JpsiF(m2, JpsiF_mean[3.12,3.0,3.35], JpsiF_sigma[0.1,0.02,0.3], JpsiF_alpha[1.2,0.4,10.0], JpsiF_n[2.0])");
  w->factory("Gaussian::Up1F(m2,Up1F_mean[9.43,9.39,9.500], Up1F_sigma[0.101,0.01,0.20])");
  w->factory("Gaussian::Up2F(m2,Up2F_mean[10.0,9.90,10.20], Up2F_sigma[0.060,0.01,0.15])");
  w->factory("Gaussian::Up3F(m2,Up3F_mean[10.35,10.2,10.5], Up3F_sigma[0.050,0.01,0.10])");
  // FINAL PDF
  w->factory("SUM::template1D_m2(norm_adHocF[10., 0., 500.]*adHocF,norm_MmumuF[200., 0., 10000.]*MmumuF, norm_bgF[4400.,0.,8000.]*bgF, norm_etaF[3.3894e-05,0.,1.]*etaF, norm_rhoF[9.,1.,100.]*rhoF, norm_phiF[100.,1.,1000.]*phiF, norm_JpsiF[50.,0.,2500.]*JpsiF, norm_psiF[50.,0.,1000.]*psiF)");
  */

  //===========
  //2017 fit m2
  //===========
  //Initial combianatorial
  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[0.01, 0.0, 0.3], MmumuF_p[0.05, 0.0, 2.])");
  //FIRST KINK
  w->factory("Bernstein::bgF(m2,{bF06[9.9751,0,15.], bF16[3.2971e-05,0.,3.], bF26[ 1.2361e-07,0.,3.], bF36[5.1545e-08,0.,2.], bF46[9.9017e-01,0.,3.], bF56[3.0607e-01,0.,3.], bF66[0.5,0.1,4.]})");
  // Resonances
  w->factory("Gaussian::etaF(m2,0.548,0.02)");
  w->factory("Gaussian::rhoF(m2,0.782,0.025)");
  w->factory("Gaussian::phiF(m2,1.019,0.02)");
  w->factory("Gaussian::psiF(m2,3.675,psiF_sigma[0.033,0.001,0.1])");
  // Ad hoc gaussian to cover first bump and help other functions
  w->factory("Gaussian::adHocF(m2,adHocF_mass[0.4,0.2,0.6],adHocF_sigma[0.01,0.001,0.1])");
  // Visible Resonances
  w->factory("CBShape::JpsiF(m2, JpsiF_mean[3.12,3.0,3.35], JpsiF_sigma[0.1,0.02,0.3], JpsiF_alpha[1.2,0.4,10.0], JpsiF_n[2.0])");
  w->factory("Gaussian::Up1F(m2,Up1F_mean[9.43,9.39,9.500], Up1F_sigma[0.101,0.01,0.20])");
  w->factory("Gaussian::Up2F(m2,Up2F_mean[10.0,9.90,10.20], Up2F_sigma[0.060,0.01,0.15])");
  w->factory("Gaussian::Up3F(m2,Up3F_mean[10.35,10.2,10.5], Up3F_sigma[0.050,0.01,0.10])");
  // FINAL PDF
  w->factory("SUM::template1D_m2(norm_adHocF[10., 0., 500.]*adHocF,norm_MmumuF[200., 0., 10000.]*MmumuF, norm_bgF[4400.,0.,8000.]*bgF, norm_etaF[3.3894e-05,0.,1.]*etaF, norm_rhoF[9.,1.,100.]*rhoF, norm_phiF[100.,1.,1000.]*phiF, norm_JpsiF[50.,0.,2500.]*JpsiF, norm_psiF[50.,0.,1000.]*psiF)");

  RooFitResult *rF = w->pdf("template1D_m2")->fitTo(*(w->data("ds_dimuorphan_bg_m2")), Extended(1), Save(), SumW2Error(kTRUE));
  rF->Print();

  RooPlot* plotF = w->var("m2")->frame(Title("BG template for non-trigger dimuon"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotF, DataError(RooAbsData::SumW2), Name("data_m2"));
  w->pdf("template1D_m2")->plotOn(plotF,LineColor(kRed),Precision(0.0001),Name("template_m2"));

  SizeBin = plotF->GetXaxis()->GetBinCenter(3) - plotF->GetXaxis()->GetBinCenter(2);
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  Yname = "Events / (" + std::string(c_SizeBin) + "[GeV])";
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
  w->factory("PROD::template2D(template1D_m1,template1D_m2)");

  //****************************************************************************
  //                       Extract 1D J/psi template from template m1
  //****************************************************************************
  cout << "Create templates for J/psi" << endl;

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
  Yname = "Events / (" + std::string(c_SizeBin) + "[GeV])";
  plot_Jpsi_m1->GetYaxis()->SetTitle( Yname.Data() );
  plot_Jpsi_m1->Draw();
  txtHeader->Draw();
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1_RooPlot.pdf");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("figures/template1D_Jpsi_m1_RooPlot.png");

  //****************************************************************************
  //                       Extract 1D J/psi template from template m2
  //****************************************************************************
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
  w->import(Jpsi_m2);

  RooPlot* plot_Jpsi_m2 = w->var("m2")->frame(Title("J/psi template m2"),Bins(m_bins));
  w->pdf("Jpsi_m2")->plotOn(plot_Jpsi_m2,LineColor(kRed),Precision(0.0001),Name("plot_Jpsi_m2"));

  TCanvas * c_template1D_Jpsi_m2_RooPlot = new TCanvas("c_template1D_Jpsi_m2_RooPlot", "c_template1D_Jpsi_m2_RooPlot");
  c_template1D_Jpsi_m2_RooPlot->cd();
  SizeBin = plot_Jpsi_m2->GetXaxis()->GetBinCenter(3) - plot_Jpsi_m2->GetXaxis()->GetBinCenter(2);
  snprintf(c_SizeBin,50,"%.2f",SizeBin);
  Yname = "Events / (" + std::string(c_SizeBin) + "[GeV])";
  plot_Jpsi_m2->GetYaxis()->SetTitle( Yname.Data() );
  plot_Jpsi_m2->Draw();
  txtHeader->Draw();
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2_RooPlot.pdf");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("figures/template1D_Jpsi_m2_RooPlot.png");

  //****************************************************************************
  //                     Create 2D template (m1 x m2) for J/psi
  //****************************************************************************
  w->factory("PROD::Jpsi_2D(Jpsi_m1,Jpsi_m2)");

  //****************************************************************************
  //                           Save to Workspace
  //****************************************************************************
  w->writeToFile("ws_FINAL.root");
  cout<<"template1D_m1_RooPlot has "<<chi2_C<<endl;
  cout<<"template1D_m2_RooPlot has "<<chi2_F<<endl;
}
