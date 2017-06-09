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
#include "TColor.h"

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
#include "RooAddPdf.h"
#include "macros/tdrStyle.C"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void PlotBBbar() {

  setTDRStyle();

  TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader("CMS Prelim. 2016  #sqrt{s} = 13 TeV   L_{int} = 35.9 fb^{-1}");

  TFile* file = new TFile("ws_FINAL.root");
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  // Resolution 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

  //****************************************************************************
  //                           Draw template for m1                             
  //****************************************************************************

  TH1D* h1_ds_dimuorphan_bg_m1 = (TH1D*)w->data("ds_dimuorphan_bg_m1")->createHistogram("m1",m_bins);
  h1_ds_dimuorphan_bg_m1->GetYaxis()->SetRangeUser(0.,200.);
  h1_ds_dimuorphan_bg_m1->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_ds_dimuorphan_bg_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");

  TH1D* h1_template1D_m1 = (TH1D*)w->pdf("template1D_m1")->createHistogram("m1",1000);
  h1_template1D_m1->Scale( h1_ds_dimuorphan_bg_m1->Integral("width") / h1_template1D_m1->Integral("width") );
  h1_template1D_m1->SetLineWidth(3);
  h1_template1D_m1->SetLineColor(kRed);

  TCanvas * c_m1 = new TCanvas("c_m1", "c_m1");
  c_m1->cd();
  h1_ds_dimuorphan_bg_m1->Draw("E1");
  h1_template1D_m1->Draw("Lsame");
  txtHeader->Draw();
  c_m1->SaveAs("figures/template1D_m1.pdf");
  c_m1->SaveAs("figures/template1D_m1.png");

  //****************************************************************************
  //                           Draw template for m2                             
  //****************************************************************************

  TH1D* h1_ds_dimuorphan_bg_m2 = (TH1D*)w->data("ds_dimuorphan_bg_m2")->createHistogram("m2",m_bins);
  h1_ds_dimuorphan_bg_m2->GetYaxis()->SetRangeUser(0.,200.);
  h1_ds_dimuorphan_bg_m2->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_ds_dimuorphan_bg_m2->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");

  TH1D* h1_template1D_m2 = (TH1D*)w->pdf("template1D_m2")->createHistogram("m2",1000);
  h1_template1D_m2->Scale( h1_ds_dimuorphan_bg_m2->Integral("width") / h1_template1D_m2->Integral("width") );
  h1_template1D_m2->SetLineWidth(3);
  h1_template1D_m2->SetLineColor(kRed);

  TCanvas * c_m2 = new TCanvas("c_m2", "c_m2");
  c_m2->cd();
  h1_ds_dimuorphan_bg_m2->Draw("E1");
  h1_template1D_m2->Draw("Lsame");
  txtHeader->Draw();
  c_m2->SaveAs("figures/template1D_m2.pdf");
  c_m2->SaveAs("figures/template1D_m2.png");

  //****************************************************************************
  //            Draw 2D template m1 x m2 and data in off-diagonal region        
  //****************************************************************************

  TH2D* h2_Template2D = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",2.0*m_bins,2.0*m_bins);
  cout << "2D template integral: " << h2_Template2D->Integral() << std::endl;

  TH2D* h2_Template2D_diagonal    = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);
  TH2D* h2_Template2D_offDiagonal = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);

  for(int i=1;i<=1000;i++) {
    for(int j=1;j<=1000;j++) {
	double m_1 = h2_Template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
	double m_2 = h2_Template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
	if ( fabs(m_1 - m_2) < 5.*(0.026 + 0.013*(m_1 + m_2)/2.) ) {
	  h2_Template2D_offDiagonal->SetBinContent(i,j,0.);
	} else {
	  h2_Template2D_diagonal->SetBinContent(i,j,0.);
	}
    }
  }
  cout << "Template2D_offDiagonal integral: " << h2_Template2D_offDiagonal->Integral() << std::endl;
  cout << "Template2D_diagonal integral:    " << h2_Template2D_diagonal->Integral() << std::endl;

  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2",1000,1000);
  cout << "Data dimudimu_control_Iso_offDiagonal integral:" << h2_dimudimu_control_Iso_offDiagonal_2D->Integral() << std::endl;

  h2_Template2D->Scale(h2_dimudimu_control_Iso_offDiagonal_2D->Integral()/h2_Template2D->Integral()*(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral());
  cout << "Scaled 2D template integral: " << h2_Template2D->Integral() << std::endl;

  TCanvas * c_template2D_m1_vs_m2 = new TCanvas("c_template2D_m1_vs_m2", "c_template2D_m1_vs_m2",0,1320,1004,928);
  c_template2D_m1_vs_m2->SetCanvasSize(1000,900);
  c_template2D_m1_vs_m2->SetLeftMargin(0.126);
  c_template2D_m1_vs_m2->SetRightMargin(0.154);
  c_template2D_m1_vs_m2->cd();
  c_template2D_m1_vs_m2->SetLogz();

  h2_Template2D->GetXaxis()->SetTitle("");
  h2_Template2D->GetYaxis()->SetTitle("");
  h2_Template2D->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_Template2D->GetZaxis()->CenterTitle(true);
  h2_Template2D->GetZaxis()->SetLabelFont(42);
  h2_Template2D->GetZaxis()->SetLabelOffset(-0.005);
  h2_Template2D->GetZaxis()->SetLabelSize(0.044);
  h2_Template2D->GetZaxis()->SetTitleSize(0.044);
  h2_Template2D->GetZaxis()->SetTitleOffset(1.2);
  h2_Template2D->GetZaxis()->SetTitleFont(42);

  //gStyle->SetPalette(52); //Grey Scale
  const Int_t NCont = 99;
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  Int_t MyPalette[NCont];
  for ( int i=0; i<NCont; i++ ) MyPalette[i] = FI+i;
  gStyle->SetPalette(NCont, MyPalette);
  h2_Template2D->SetContour(NCont);

  h2_Template2D->Draw("Cont4 Colz");

  // This required to draw scatter plot without LogZ 
  TPad* pad = new TPad("pad", "pad",0,0,1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.126);
  pad->SetRightMargin(0.154);

  pad->SetFillColor(0);
  pad->SetFillStyle(4000);
  pad->SetBorderMode(0);
  pad->SetBorderSize(2);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetFrameFillStyle(0);
  pad->SetFrameBorderMode(0);
  pad->SetFrameFillStyle(0);
  pad->SetFrameBorderMode(0);

  h2_dimudimu_control_Iso_offDiagonal_2D->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
  h2_dimudimu_control_Iso_offDiagonal_2D->GetXaxis()->SetTitleOffset(0.93);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
  h2_dimudimu_control_Iso_offDiagonal_2D->GetYaxis()->SetTitleOffset(0.95);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->CenterTitle(true);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetLabelFont(42);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetLabelOffset(-0.005);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetLabelSize(0.044);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetTitleSize(0.044);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetTitleOffset(1.2);
  h2_dimudimu_control_Iso_offDiagonal_2D->GetZaxis()->SetTitleFont(42);
  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerSize(3.0);
  h2_dimudimu_control_Iso_offDiagonal_2D->Draw("");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_tmp = new TH2D( *h2_dimudimu_control_Iso_offDiagonal_2D);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerSize(2.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->Draw("same");

  double diagonal_x1 = ( (1.0+kB/2.0)*m_min + kA )/( 1.0 - kB/2.0 );
  double diagonal_x2 = ( (1.0-kB/2.0)*m_max - kA )/( 1.0 + kB/2.0 );
  std::cout << "diagonal_x1 " << diagonal_x1 << std::endl;
  std::cout << "diagonal_x2 " << diagonal_x2 << std::endl;

  TLine *line1 = new TLine(m_min, diagonal_x1, diagonal_x2, m_max);
  line1->SetLineColor(0);
  line1->SetLineStyle(9);
  line1->SetLineWidth(2);
  line1->Draw();
  TLine *line2 = new TLine(diagonal_x1,m_min,m_max,diagonal_x2);
  line2->SetLineColor(0);
  line2->SetLineStyle(9);
  line2->SetLineWidth(2);
  line2->Draw();

  txtHeader->Draw();

  c_template2D_m1_vs_m2->SaveAs("figures/template2D_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures/template2D_m1_vs_m2.png");

  TCanvas * c_template2D_m1_vs_m2_3DLog = new TCanvas("c_template2D_m1_vs_m2_3DLog", "c_template2D_m1_vs_m2_3DLog");
  c_template2D_m1_vs_m2_3DLog->cd();
  c_template2D_m1_vs_m2_3DLog->SetLogz();
  TH2D *h2_Template2D_3DLog = new TH2D( *h2_Template2D );
  h2_Template2D_3DLog->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
  h2_Template2D_3DLog->GetXaxis()->SetTitleOffset(1.5);
  h2_Template2D_3DLog->GetXaxis()->CenterTitle(1);
  h2_Template2D_3DLog->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
  h2_Template2D_3DLog->GetYaxis()->SetTitleOffset(1.7);
  h2_Template2D_3DLog->GetYaxis()->CenterTitle(1);
  h2_Template2D_3DLog->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_Template2D_3DLog->GetZaxis()->CenterTitle(true);
  h2_Template2D_3DLog->GetZaxis()->SetLabelFont(42);
  //  h2_Template2D_3DLog->GetZaxis()->SetLabelOffset(-0.005);
  h2_Template2D_3DLog->GetZaxis()->SetLabelSize(0.044);
  h2_Template2D_3DLog->GetZaxis()->SetTitleSize(0.044);
  h2_Template2D_3DLog->GetZaxis()->SetTitleOffset(1.8);
  h2_Template2D_3DLog->GetZaxis()->SetTitleFont(42);
  h2_Template2D_3DLog->GetZaxis()->SetRangeUser(0.00004,0.1);
  h2_Template2D_3DLog->Draw("Surf1 FB");
  txtHeader->Draw();

  c_template2D_m1_vs_m2_3DLog->SaveAs("figures/template2D_m1_vs_m2_3D.pdf");
  c_template2D_m1_vs_m2_3DLog->SaveAs("figures/template2D_m1_vs_m2_3D.png");

  //****************************************************************************
  //        Create 2D template = m1 x m2 without J/psi region                   
  //****************************************************************************

  TH2D* h2_Template2D_NoJpsi = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);
  TH2D* h2_Template2D_Jpsi   = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);

  for(int i=1;i<=1000;i++) {
    for(int j=1;j<=1000;j++) {
	double m_1 = h2_Template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
	double m_2 = h2_Template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
	if ( m_1 > 2.95 && m_1 < 3.2 ) {
	  h2_Template2D_NoJpsi->SetBinContent(i,j,0.);
	} else {
	  h2_Template2D_Jpsi->SetBinContent(i,j,0.);
	}
    }
  }
  cout << "h2_Template2D_NoJpsi integral: " << h2_Template2D_NoJpsi->Integral() << std::endl;
  cout << "h2_Template2D_Jpsi integral:   " << h2_Template2D_Jpsi->Integral() << std::endl;

  //****************************************************************************
  //              Control region = off diagonal region                          
  //                             m1                                             
  //****************************************************************************

  // NO isolation requirement applied

  TH1D *h1_control_offDiagonal_massC_data = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_1D_massC")->createHistogram("m1",m_bins);
  h1_control_offDiagonal_massC_data->SetStats(0);
  h1_control_offDiagonal_massC_data->SetMarkerStyle(20);
  h1_control_offDiagonal_massC_data->GetXaxis()->SetTitle("m_{#mu#mu_{1}} (GeV/#it{c}^{2})");
  h1_control_offDiagonal_massC_data->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_control_offDiagonal_massC_data->GetYaxis()->SetRangeUser(0.,280.);

  TH1D *h1_control_offDiagonal_massC_template = new TH1D( *h2_Template2D_offDiagonal->ProjectionX() );
  h1_control_offDiagonal_massC_template->Scale( h1_control_offDiagonal_massC_data->Integral("width") / h1_control_offDiagonal_massC_template->Integral("width") );

  h1_control_offDiagonal_massC_template->SetLineColor(kRed);
  h1_control_offDiagonal_massC_template->SetLineWidth(3);

  TCanvas * c_control_offDiagonal_massC = new TCanvas("c_control_offDiagonal_massC", "c_control_offDiagonal_massC");
  c_control_offDiagonal_massC->cd();
  h1_control_offDiagonal_massC_data->Draw("e1");
  h1_control_offDiagonal_massC_template->Draw("same");
  txtHeader->Draw();
  c_control_offDiagonal_massC->SaveAs("figures/control_offDiagonal_m1.pdf");
  c_control_offDiagonal_massC->SaveAs("figures/control_offDiagonal_m1.png");

  // Isolation requirement applied

  TH1D *h1_control_Iso_offDiagonal_massC_data = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_1D_massC")->createHistogram("m1",m_bins);
  h1_control_Iso_offDiagonal_massC_data->SetStats(0);
  h1_control_Iso_offDiagonal_massC_data->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massC_data->GetXaxis()->SetTitle("m_{#mu#mu_{1}} (GeV/#it{c}^{2})");
  h1_control_Iso_offDiagonal_massC_data->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_control_Iso_offDiagonal_massC_data->GetYaxis()->SetRangeUser(0.,3.);

  TH1D *h1_control_Iso_offDiagonal_massC_template = new TH1D( *h2_Template2D_offDiagonal->ProjectionX() );

  h1_control_Iso_offDiagonal_massC_template->Scale( h1_control_Iso_offDiagonal_massC_data->Integral("width") / h1_control_Iso_offDiagonal_massC_template->Integral("width") );

  h1_control_Iso_offDiagonal_massC_template->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massC_template->SetLineWidth(3);

  TCanvas * c_control_Iso_offDiagonal_massC = new TCanvas("c_control_Iso_offDiagonal_massC", "c_control_Iso_offDiagonal_massC");
  c_control_Iso_offDiagonal_massC->cd();
  h1_control_Iso_offDiagonal_massC_data->Draw("e1");
  h1_control_Iso_offDiagonal_massC_template->Draw("same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massC->SaveAs("control_Iso_offDiagonal_m1.pdf");
  c_control_Iso_offDiagonal_massC->SaveAs("control_Iso_offDiagonal_m1.png");

  //****************************************************************************
  //               Control region = off diagonal region                         
  //                             m2                                             
  //****************************************************************************

  // NO isolation requirement applied

  TH1D *h1_control_offDiagonal_massF_data = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_1D_massF")->createHistogram("m1",m_bins);
  h1_control_offDiagonal_massF_data->SetStats(0);
  h1_control_offDiagonal_massF_data->SetMarkerStyle(20);
  h1_control_offDiagonal_massF_data->GetXaxis()->SetTitle("m_{#mu#mu_{2}} (GeV/#it{c}^{2})");
  h1_control_offDiagonal_massF_data->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_control_offDiagonal_massF_data->GetYaxis()->SetRangeUser(0.,200.);

  TH1D *h1_control_offDiagonal_massF_template = new TH1D( *h2_Template2D_offDiagonal->ProjectionY() );
  h1_control_offDiagonal_massF_template->Scale( h1_control_offDiagonal_massF_data->Integral("width") / h1_control_offDiagonal_massF_template->Integral("width") );

  h1_control_offDiagonal_massF_template->SetLineColor(kRed);
  h1_control_offDiagonal_massF_template->SetLineWidth(3);

  TCanvas * c_control_offDiagonal_massF = new TCanvas("c_control_offDiagonal_massF", "c_control_offDiagonal_massF");
  c_control_offDiagonal_massF->cd();
  h1_control_offDiagonal_massF_data->Draw("e1");
  h1_control_offDiagonal_massF_template->Draw("same");
  txtHeader->Draw();
  c_control_offDiagonal_massF->SaveAs("figures/control_offDiagonal_m2.pdf");
  c_control_offDiagonal_massF->SaveAs("figures/control_offDiagonal_m2.png");

  // Isolation requirement applied

  TH1D *h1_control_Iso_offDiagonal_massF_data = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_1D_massF")->createHistogram("m1",m_bins);
  h1_control_Iso_offDiagonal_massF_data->SetStats(0);
  h1_control_Iso_offDiagonal_massF_data->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massF_data->GetXaxis()->SetTitle("m_{#mu#mu_{2}} (GeV/#it{c}^{2})");
  h1_control_Iso_offDiagonal_massF_data->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_control_Iso_offDiagonal_massF_data->GetYaxis()->SetRangeUser(0.,3.);

  TH1D *h1_control_Iso_offDiagonal_massF_template = new TH1D( *h2_Template2D_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_massF_template->Scale( h1_control_Iso_offDiagonal_massF_data->Integral("width") / h1_control_Iso_offDiagonal_massF_template->Integral("width") );

  h1_control_Iso_offDiagonal_massF_template->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massF_template->SetLineWidth(3);

  TCanvas * c_control_Iso_offDiagonal_massF = new TCanvas("c_control_Iso_offDiagonal_massF", "c_control_Iso_offDiagonal_massF");
  c_control_Iso_offDiagonal_massF->cd();
  h1_control_Iso_offDiagonal_massF_data->Draw("e1");
  h1_control_Iso_offDiagonal_massF_template->Draw("same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massF->SaveAs("figures/control_Iso_offDiagonal_m2.pdf");
  c_control_Iso_offDiagonal_massF->SaveAs("figures/control_Iso_offDiagonal_m2.png");

  //****************************************************************************
  //                    Control region = off diagonal region                    
  //                             m1 and m2 combined                             
  //****************************************************************************

  // NO isolation requirement applied

  TH1D *h1_control_offDiagonal_data = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_1D")->createHistogram("m1",m_bins);
  h1_control_offDiagonal_data->SetStats(0);
  h1_control_offDiagonal_data->SetMarkerStyle(20);
  h1_control_offDiagonal_data->GetXaxis()->SetTitle("m_{#mu#mu_{i}} (i=1,2) [GeV/#it{c}^{2}]");
  h1_control_offDiagonal_data->GetYaxis()->SetTitle("Events #times 2 / (0.05 GeV/#it{c}^{2})");
  h1_control_offDiagonal_data->GetYaxis()->SetRangeUser(0.,240.);

  TH1D *h1_control_offDiagonal_template = new TH1D( *h2_Template2D_offDiagonal->ProjectionX() );
  TH1D *h1_control_offDiagonal_template_Y = new TH1D( *h2_Template2D_offDiagonal->ProjectionY() );
  h1_control_offDiagonal_template->Add(h1_control_offDiagonal_template_Y);
  h1_control_offDiagonal_template->Scale( h1_control_offDiagonal_data->Integral("width") / h1_control_offDiagonal_template->Integral("width") );

  h1_control_offDiagonal_template->SetLineColor(kRed);
  h1_control_offDiagonal_template->SetLineWidth(3);

  TCanvas * c_control_offDiagonal = new TCanvas("c_control_offDiagonal", "c_control_offDiagonal");
  c_control_offDiagonal->cd();
  h1_control_offDiagonal_data->Draw("e1");
  h1_control_offDiagonal_template->Draw("same");
  txtHeader->Draw();
  c_control_offDiagonal->SaveAs("figures/control_offDiagonal.pdf");
  c_control_offDiagonal->SaveAs("figures/control_offDiagonal.png");

  // Isolation requirement applied

  TH1D *h1_control_Iso_offDiagonal_data = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_1D")->createHistogram("m1",m_bins);
  h1_control_Iso_offDiagonal_data->SetStats(0);
  h1_control_Iso_offDiagonal_data->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_data->GetXaxis()->SetTitle("m_{#mu#mu_{i}} (i=1,2) [GeV/#it{c}^{2}]");
  h1_control_Iso_offDiagonal_data->GetYaxis()->SetTitle("Events #times 2 / (0.05 GeV/#it{c}^{2})");
  h1_control_Iso_offDiagonal_data->GetYaxis()->SetRangeUser(0.,3.);

  TH1D *h1_control_Iso_offDiagonal_template = new TH1D( *h2_Template2D_offDiagonal->ProjectionX() );
  TH1D *h1_control_Iso_offDiagonal_templateY = new TH1D( *h2_Template2D_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_template->Add( h1_control_Iso_offDiagonal_templateY );
  h1_control_Iso_offDiagonal_template->Scale( h1_control_Iso_offDiagonal_data->Integral("width") / h1_control_Iso_offDiagonal_template->Integral("width") );

  h1_control_Iso_offDiagonal_template->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_template->SetLineWidth(3);

  TCanvas * c_control_Iso_offDiagonal = new TCanvas("c_control_Iso_offDiagonal", "c_control_Iso_offDiagonal");
  c_control_Iso_offDiagonal->cd();
  h1_control_Iso_offDiagonal_data->Draw("e1");
  h1_control_Iso_offDiagonal_template->Draw("same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal->SaveAs("figures/control_Iso_offDiagonal.pdf");
  c_control_Iso_offDiagonal->SaveAs("figures/control_Iso_offDiagonal.png");

  //****************************************************************************
  //       Control region = off diagonal region, m1 and m2 combined             
  //                            J/psi vicinity                                  
  //****************************************************************************

  // NO isolation requirement applied

  TH1D *h1_control_offDiagonal_Jpsi = new TH1D("h1_control_offDiagonal_Jpsi", "h1_control_offDiagonal_Jpsi", 5, 2.95, 3.2);
  h1_control_offDiagonal_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{i}} (i=1,2) [GeV/#it{c}^{2}]");
  h1_control_offDiagonal_Jpsi->GetYaxis()->SetTitle("Events #times 2 / (0.05 GeV/#it{c}^{2})");
  h1_control_offDiagonal_Jpsi->GetYaxis()->SetRangeUser(0.,500.);

  TCanvas * c_control_offDiagonal_Jpsi = new TCanvas("c_control_offDiagonal_Jpsi", "c_control_offDiagonal_Jpsi");
  c_control_offDiagonal_Jpsi->cd();
  h1_control_offDiagonal_Jpsi->Draw();
  h1_control_offDiagonal_data->Draw("e1same");
  h1_control_offDiagonal_template->Draw("Csame");
  txtHeader->Draw();
  c_control_offDiagonal_Jpsi->SaveAs("figures/control_offDiagonal_Jpsi.pdf");
  c_control_offDiagonal_Jpsi->SaveAs("figures/control_offDiagonal_Jpsi.png");

  // Isolation requirement applied

  TH1D *h1_control_Iso_offDiagonal_Jpsi = new TH1D("h1_control_Iso_offDiagonal_Jpsi", "h1_control_Iso_offDiagonal_Jpsi", 5, 2.95, 3.2);
  h1_control_Iso_offDiagonal_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{i}} (i=1,2) [GeV/#it{c}^{2}]");
  h1_control_Iso_offDiagonal_Jpsi->GetYaxis()->SetTitle("Events #times 2 / (0.05 GeV/#it{c}^{2})");
  h1_control_Iso_offDiagonal_Jpsi->GetYaxis()->SetRangeUser(0.,3.);

  TCanvas * c_control_Iso_offDiagonal_Jpsi = new TCanvas("c_control_Iso_offDiagonal_Jpsi", "c_control_Iso_offDiagonal_Jpsi");
  c_control_Iso_offDiagonal_Jpsi->cd();
  h1_control_Iso_offDiagonal_Jpsi->Draw();
  h1_control_Iso_offDiagonal_data->Draw("e1same");
  h1_control_Iso_offDiagonal_template->Draw("Csame");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_Jpsi->SaveAs("figures/control_Iso_offDiagonal_Jpsi.pdf");
  c_control_Iso_offDiagonal_Jpsi->SaveAs("figures/control_Iso_offDiagonal_Jpsi.png");

  w->writeToFile("ws.root");
}
