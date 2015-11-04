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
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "tdrStyle.C"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void FitAndSave() {

  setTDRStyle();

  TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader("CMS Prelim. 2015  #sqrt{s} = 13 TeV   L_{int} = ~ pb^{-1}");

  RooWorkspace* w = new RooWorkspace("w");

//  ifstream fin( "FitNtuple_Run2015D2_DoubleMu.txt" );
//  string line;
//  getline( fin, line );
//  stringstream ss( line ); int nLine; ss >> nLine;
//  vector< string > lines( nLine );
//  int i = 0;
//  cout<<"nLine "<<nLine<<endl;
//  while( i < nLine && getline( fin, line ) ) {
//      lines[i] = line;
//      cout<<"line) "<<line<<endl;
//      ++i;
//  }

  TString file_name = "FitNtuple_Run2015D3_DoubleMu.root";
  TChain chain_data_dimudimu("cutFlowAnalyzer/Events");
  TChain chain_data_dimuorphan("cutFlowAnalyzer/Events_orphan");
  chain_data_dimudimu.Add(file_name.Data());
  chain_data_dimuorphan.Add(file_name.Data());

  ////Load Ntuples
  //string line[9999];
  //int nLine = 0;
  //std::ifstream myfile( "FitNtuple_Run2015D2_DoubleMu.txt" );
  //for(int i = 0; myfile.good(); i++) nLine++;
  //for(int i = 0; i<nLine-1; i++) {
  //  getline(myfile, line[i]);
  //  cout<<i<<") "<<line[i]<<endl;    
  //  chain_data_dimudimu.Add( line[i].c_str() );
  //  chain_data_dimuorphan.Add( line[i].c_str() );
  //}

  const double       m_min  = 0.2113;
  const double       m_max  = 3.5536;
  const unsigned int m_bins = 66;

  //Efficiencies studies
  float Pass_OffLine                 = chain_data_dimuorphan.Draw("orph_passOffLineSel>>hist","orph_passOffLineSel>0","goff");
  float Pass_OffLine_pt              = chain_data_dimuorphan.Draw("orph_passOffLineSelPt>>hist","orph_passOffLineSelPt>0","goff");
  float Pass_OffLine_pt1788          = chain_data_dimuorphan.Draw("orph_passOffLineSelPt1788>>hist","orph_passOffLineSelPt1788>0","goff");
  float Pass_FiredTrig               = chain_data_dimuorphan.Draw("orph_FiredTrig>>hist","orph_FiredTrig>0","goff");
  float Pass_FiredTrig_pt            = chain_data_dimuorphan.Draw("orph_FiredTrig_pt>>hist","orph_FiredTrig_pt>0","goff");
  float Pass_FiredTrig_ptColl        = chain_data_dimuorphan.Draw("orph_FiredTrig_ptColl>>hist","orph_FiredTrig_ptColl>0","goff");

  float Pass_FiredTrig_notOff        = chain_data_dimuorphan.Draw("orph_FiredTrig>>hist","orph_FiredTrig>0 && orph_passOffLineSelPt1788<1","goff");
  float Pass_Offl_noFiredTrig        = chain_data_dimuorphan.Draw("orph_FiredTrig>>hist","orph_FiredTrig<1 && orph_passOffLineSelPt1788>0","goff");

  float Pass_OffLineFiredTrig        = chain_data_dimuorphan.Draw("orph_FiredTrig>>hist","orph_FiredTrig>0 && orph_passOffLineSelPt1788>0","goff");
  float Pass_OffLineFiredTrig_pt     = chain_data_dimuorphan.Draw("orph_FiredTrig_pt>>hist","orph_FiredTrig_pt>0 && orph_passOffLineSelPt1788>0","goff");
  float Pass_OffLineFiredTrig_ptColl = chain_data_dimuorphan.Draw("orph_FiredTrig_ptColl>>hist","orph_FiredTrig_ptColl>0 && orph_passOffLineSelPt1788>0","goff");

  RooRealVar m1("m1","m_{#mu#mu_{1}}",m_min,m_max,"GeV/#it{c}^{2}");
  RooRealVar m2("m2","m_{#mu#mu_{2}}",m_min,m_max,"GeV/#it{c}^{2}");
  m1.setBins(m_bins);
  m2.setBins(m_bins);
  w->import(m1);
  w->import(m2);

  ostringstream stream_cut_bg_m1_iso;
  //stream_cut_bg_m1_iso << "orph_dimu_isoTk < 2. && orph_dimu_isoTk >= 0 && containstrig2 > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  //stream_cut_bg_m1_iso << "orph_dimu_isoTk >= 0 && containstrig2 > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  stream_cut_bg_m1_iso << "orph_dimu_isoTk < 100. && orph_dimu_isoTk >= 0 && containstrig2 > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;

  ostringstream stream_cut_bg_m2_iso;
  //stream_cut_bg_m2_iso << "orph_dimu_isoTk < 2. && orph_dimu_isoTk >= 0 && containstrig > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  //stream_cut_bg_m2_iso << "orph_dimu_isoTk > 0 && containstrig > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  stream_cut_bg_m2_iso << "orph_dimu_isoTk < 100. && orph_dimu_isoTk >= 0 && containstrig > 0 && orph_dimu_mass > "<< m_min << " && orph_dimu_mass < " << m_max;
  TString cut_bg_m1_iso = stream_cut_bg_m1_iso.str();
  TString cut_bg_m2_iso = stream_cut_bg_m2_iso.str();
  TString cut_diagonal                = "abs(massC-massF) <= (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 3.55 && massF > 0.25 && massF < 3.55";
  TString cut_control_offDiagonal     = "abs(massC-massF) > (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 3.55 && massF > 0.25 && massF < 3.55";
  TString cut_control_Iso_offDiagonal = "isoC_1mm >= 0 && isoC_1mm < 2. && isoF_1mm >= 0 && isoF_1mm < 2. && abs(massC-massF) > (0.13 + 0.065*(massC+massF)/2.) && massC > 0.25 && massC < 3.55 && massF > 0.25 && massF < 3.55";
  TString cut_control_nonIso          = "isoC_1mm > 2. && isoC_1mm < 8. && isoF_1mm > 2. && isoF_1mm < 8. && massC > 0.25 && massC < 3.55 && massF > 0.25 && massF < 3.55";

  TString cut_signal                  = "isoC_1mm>=0 && isoC_1mm<2. && isoF_1mm>=0 && isoF_1mm<2. && abs(massC-massF) <= (0.13 + 0.065*(massC+massF)/2.)";

  TTree* tree_dimuorphan_bg_m1        = chain_data_dimuorphan.CopyTree(cut_bg_m1_iso);
  TTree* tree_dimuorphan_bg_m2        = chain_data_dimuorphan.CopyTree(cut_bg_m2_iso);

  TTree* tree_dimudimu_diagonal_2D                      = chain_data_dimudimu.CopyTree(cut_diagonal);
  TTree* tree_dimudimu_diagonal_1D_massC                = chain_data_dimudimu.CopyTree(cut_diagonal);
  TTree* tree_dimudimu_diagonal_1D_massF                = chain_data_dimudimu.CopyTree(cut_diagonal);
  TTree* tree_dimudimu_control_offDiagonal_2D           = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  TTree* tree_dimudimu_control_offDiagonal_1D_massC     = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  TTree* tree_dimudimu_control_offDiagonal_1D_massF     = chain_data_dimudimu.CopyTree(cut_control_offDiagonal);
  TTree* tree_dimudimu_control_Iso_offDiagonal_2D       = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);
  cout<<"------OffDiagonal SCAN------"<<endl;
  tree_dimudimu_control_Iso_offDiagonal_2D->Scan("massC:massF");
  TTree* tree_dimudimu_control_Iso_offDiagonal_1D_massC = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);
  TTree* tree_dimudimu_control_Iso_offDiagonal_1D_massF = chain_data_dimudimu.CopyTree(cut_control_Iso_offDiagonal);
  TTree* tree_dimudimu_control_nonIso                   = chain_data_dimudimu.CopyTree(cut_control_nonIso);

  TTree* tree_dimudimu_signal_2D  = chain_data_dimudimu.CopyTree(cut_signal);
  //cout<<"------Signal SCAN------"<<endl;
  //tree_dimudimu_signal_2D->Scan("massC:massF");

  tree_dimuorphan_bg_m1->GetBranch("orph_dimu_mass")->SetName("m1");
  tree_dimuorphan_bg_m2->GetBranch("orph_dimu_mass")->SetName("m2");
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

  RooDataSet* ds_dimuorphan_bg_m1 = new RooDataSet("ds_dimuorphan_bg_m1","ds_dimuorphan_bg_m1", tree_dimuorphan_bg_m1, RooArgSet(m1));
  RooDataSet* ds_dimuorphan_bg_m2 = new RooDataSet("ds_dimuorphan_bg_m2","ds_dimuorphan_bg_m2", tree_dimuorphan_bg_m2, RooArgSet(m2));

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

  cout<<"-----Now Printing the Datasets:-----"<<endl;
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

  //****************************************************************************
  //                         Create template for m1                             
  //****************************************************************************
  cout<<"-----Creating templates:-----"<<endl;
  //  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.5, 0.0, 2.0], MmumuC_p[0.5])");
  w->factory("EXPR::MmumuC('m1*pow( (m1/m)*(m1/m) - 1.0, MmumuC_p )*exp( -MmumuC_c*( (m1/m)*(m1/m) - 1.0 ) )',m1, m[0.2113], MmumuC_c[0.5, 0.0, 20.0], MmumuC_p[0.5])");
  //  w->factory("Bernstein::bgC(m1,{bC06[0.0], bC16[0.1,0.,3.], bC26[1.,0.,3.], bC36[2.,0.,3.], bC46[0.2,0.,3.], bC56[0.5,0.,3.], bC66[0.1,0.,3.]})");
  w->factory("Bernstein::bgC(m1,{bC06[0.0], bC16[0.1,0.,3.], bC26[1.,0.,3.], bC36[2.,0.,3.], bC46[0.2,0.,3.], bC56[0.5,0.,3.], bC66[0.1,0.,3.]})");
  w->factory("Gaussian::etaC(m1,0.548,0.030)");
  w->factory("Gaussian::rhoC(m1,0.782,0.031)");
  w->factory("Gaussian::phiC(m1,1.019,0.033)");
  w->factory("CBShape::JpsiC(m1, JpsiC_mean[3.097,3.0,3.2], JpsiC_sigma[0.028,0.01,0.06], JpsiC_alpha[1.8,1.0,3.0], JpsiC_n[2.0])");

  //w->factory("SUM::template1D_m1(norm_MmumuC[1000., 0., 50000.]*MmumuC, norm_bgC[4400.,1000.,8000.]*bgC, norm_etaC[150.,0.,800.]*etaC, norm_rhoC[300.,0.,800.]*rhoC, norm_phiC[500.,0.,800.]*phiC, norm_JpsiC[1800.,900.,2700.]*JpsiC)");
  w->factory("SUM::template1D_m1(norm_MmumuC[1000., 0., 50000.]*MmumuC, norm_bgC[4400.,1000.,8000.]*bgC, norm_etaC[150.,0.,800.]*etaC, norm_rhoC[300.,0.,800.]*rhoC, norm_phiC[500.,0.,800.]*phiC, norm_JpsiC[1000.,300.,2700.]*JpsiC)");

  RooFitResult *rC = w->pdf("template1D_m1")->fitTo(*(w->data("ds_dimuorphan_bg_m1")), Extended(1), Save(), SumW2Error(kTRUE));
  rC->Print();

  RooPlot* plotC = w->var("m1")->frame(Title("BG template for trigger dimuon"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m1")->plotOn(plotC, DataError(RooAbsData::SumW2), Name("data_m1"));
  w->pdf("template1D_m1")->plotOn(plotC,LineColor(kRed),Precision(0.0001),Name("template_m1"));
  //  w->pdf("template1D_m1")->paramOn(plotC,Layout(0.2,0.7,0.9),ShowConstants(1));
  //  plotC->getAttText()->SetTextSize(0.02);

  TCanvas * c_template1D_m1_RooPlot = new TCanvas("c_template1D_m1_RooPlot", "c_template1D_m1_RooPlot");
  c_template1D_m1_RooPlot->cd();
  plotC->Draw();
  txtHeader->Draw();
  c_template1D_m1_RooPlot->SaveAs("template1D_m1_RooPlot.pdf");
  c_template1D_m1_RooPlot->SaveAs("template1D_m1_RooPlot.png");

  //****************************************************************************
  //                         Create template for m2                             
  //****************************************************************************

  //  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[1.0, 0.5, 1.0], MmumuF_p[1.0, 0.5, 1.0])");
  w->factory("EXPR::MmumuF('m2*pow( (m2/m)*(m2/m) - 1.0, MmumuF_p )*exp( -MmumuF_c*( (m2/m)*(m2/m) - 1.0 ) )',m2, m[0.2113], MmumuF_c[1.0, 0.5, 10.0], MmumuF_p[0.5])");
  //  w->factory("Bernstein::bgF(m2,{bF06[0.0],bF16[0.1,0.,3.],bF26[1.,0.,3.],bF36[2.,0.,3.],bF46[0.2,0.,3.],bF56[0.5,0.,3.],bF66[0.1,0.,3.]})");
  w->factory("Bernstein::bgF(m2,{bF06[0.0],bF16[0.1,0.,3.],bF26[1.,-1.,3.],bF36[2.,0.,3.],bF46[0.2,0.,3.],bF56[0.5,0.,3.],bF66[0.1,-1.,3.]})");
  w->factory("Gaussian::etaF(m2,0.548,0.033)");
  w->factory("Gaussian::rhoF(m2,0.782,0.036)");
  w->factory("Gaussian::phiF(m2,1.019,0.039)");
  w->factory("CBShape::JpsiF(m2, JpsiF_mean[3.097,3.0,3.2], JpsiF_sigma[0.028,0.01,0.06], JpsiF_alpha[1.8,1.0,3.0], JpsiF_n[2.0])");

  //w->factory("SUM::template1D_m2(norm_MmumuF[100., 0., 5000.]*MmumuF, norm_bgF[2400.,1000.,5000.]*bgF,norm_etaF[15.,0.0,100.]*etaF,norm_rhoF[15.,0.,100.]*rhoF,norm_phiF[15.,10.,100.]*phiF,norm_JpsiF[1200.,300.,3000.]*JpsiF)");
  //m1: w->factory("SUM::template1D_m1(norm_MmumuC[1000., 0., 50000.]*MmumuC, norm_bgC[4400.,1000.,8000.]*bgC, norm_etaC[150.,0.,800.]*etaC, norm_rhoC[300.,0.,800.]*rhoC, norm_phiC[500.,0.,800.]*phiC, norm_JpsiC[1800.,900.,2700.]*JpsiC)");
    w->factory("SUM::template1D_m2(norm_MmumuF[1000., 0., 50000.]*MmumuF, norm_bgF[700.,700.,8000.]*bgF, norm_etaF[15.,0.,800.]*etaF, norm_rhoF[300.,0.,800.]*rhoF, norm_phiF[500.,0.,800.]*phiF, norm_JpsiF[700.,200.,2700.]*JpsiF)");

  RooFitResult *rF = w->pdf("template1D_m2")->fitTo(*(w->data("ds_dimuorphan_bg_m2")), Extended(1), Save(), SumW2Error(kTRUE));
  rF->Print();

  RooPlot* plotF = w->var("m2")->frame(Title("BG template for non-trigger dimuon"),Bins(m_bins));
  w->data("ds_dimuorphan_bg_m2")->plotOn(plotF, DataError(RooAbsData::SumW2), Name("data_m2"));
  w->pdf("template1D_m2")->plotOn(plotF,LineColor(kRed),Precision(0.0001),Name("template_m2"));
  //  w->pdf("template1D_m2")->paramOn(plotF,Layout(0.2,0.7,0.9),ShowConstants(1));
  //  plotF->getAttText()->SetTextSize(0.02);

  TCanvas * c_template1D_m2_RooPlot = new TCanvas("c_template1D_m2_RooPlot", "c_template1D_m2_RooPlot");
  c_template1D_m2_RooPlot->cd();
  plotF->Draw();
  txtHeader->Draw();
  c_template1D_m2_RooPlot->SaveAs("template1D_m2_RooPlot.pdf");
  c_template1D_m2_RooPlot->SaveAs("template1D_m2_RooPlot.png");

  //****************************************************************************
  //                     Create 2D template = m1 x m2                           
  //****************************************************************************

  w->factory("PROD::template2D(template1D_m1,template1D_m2)");

  //****************************************************************************
  //                       Create 1D templates for J/psi                        
  //****************************************************************************

  cout << "Create templates for J/psi" << endl;

  //  template for m1

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
  plot_Jpsi_m1->Draw();
  txtHeader->Draw();
  c_template1D_Jpsi_m1_RooPlot->SaveAs("template1D_Jpsi_m1_RooPlot.pdf");
  c_template1D_Jpsi_m1_RooPlot->SaveAs("template1D_Jpsi_m1_RooPlot.png");

  //  template for m2

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
  plot_Jpsi_m2->Draw();
  txtHeader->Draw();
  c_template1D_Jpsi_m2_RooPlot->SaveAs("template1D_Jpsi_m2_RooPlot.pdf");
  c_template1D_Jpsi_m2_RooPlot->SaveAs("template1D_Jpsi_m2_RooPlot.png");

  //****************************************************************************
  //                     Create 2D template (m1 x m2) for J/psi                 
  //****************************************************************************

  w->factory("PROD::Jpsi_2D(Jpsi_m1,Jpsi_m2)");

  //****************************************************************************
  //                           Save to Workspace                                
  //****************************************************************************

  w->writeToFile("ws.root");

  cout<<"Efficiencies on bb: "<<endl;
  cout<<"--#Ev. Pass OffLine: "<<Pass_OffLine<<" but only "<<Pass_OffLine_pt<<" has Pt>17 and only "<<Pass_OffLine_pt1788<<" have pt 17 8 8"<<endl;
  cout<<" Times you have fired but no selection: "<<Pass_FiredTrig_notOff<<endl;
  cout<<" Times you have not fired but passed selection: "<<Pass_Offl_noFiredTrig<<endl;
  cout<<"--#Ev. FiredTrig: "<<Pass_FiredTrig<<" Eff. to OffLine: "<<Pass_FiredTrig/Pass_OffLine_pt<<". I have "<<Pass_OffLineFiredTrig<<" #Ev. FiredTrigger && Offline. "<<Pass_OffLineFiredTrig/Pass_OffLine_pt<<endl;
  cout<<"--#Ev. FiredTrig+Pt: "<<Pass_FiredTrig_pt<<" Eff. to OffLine: "<<Pass_FiredTrig_pt/Pass_OffLine_pt<<". I have "<<Pass_OffLineFiredTrig_pt<<" #Ev. FiredTrigger+Pt && Offline. "<<Pass_OffLineFiredTrig_ptColl/Pass_OffLine_pt<<endl;
  cout<<"--#Ev. FiredTrig+PT+Coll: "<<Pass_FiredTrig_ptColl<<" Eff. to OffLine: "<<Pass_FiredTrig_ptColl/Pass_OffLine_pt<<". I have "<<Pass_OffLineFiredTrig_ptColl<<" #Ev. FiredTrigger+pt+Coll && Offline. "<<Pass_OffLineFiredTrig_ptColl/Pass_OffLine_pt<<endl;
}
