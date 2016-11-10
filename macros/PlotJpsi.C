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

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void PlotJpsi() {
  
	setTDRStyle();
	
	TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
//  txtHeader->SetHeader("CMS Simulation");
//  txtHeader->SetHeader("CMS Simulation #sqrt{s} = 8 TeV");
//  txtHeader->SetHeader("CMS Prelim. 2011   #sqrt{s} = 7 TeV   L_{int} = 5.3 fb^{-1}");
//  txtHeader->SetHeader("CMS 2011   #sqrt{s} = 7 TeV   L_{int} = 5.3 fb^{-1}");
  txtHeader->SetHeader("CMS Prelim. 2012  #sqrt{s} = 8 TeV   L_{int} = 20.7 fb^{-1}");
//  txtHeader->SetHeader("CMS 2012   #sqrt{s} = 8 TeV   L_{int} = 20.65 fb^{-1}");
//  txtHeader->Draw();
	
	RooWorkspace* w = new RooWorkspace("w");
	
	TFile* file = new TFile("ws.root");
	RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
	const double       m_max  = 3.5536;
	const unsigned int m_bins = 66;
  
  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

//  double nEvents_Jpsi = 2.0;
  double nEvents_Jpsi = 0.058;

  //****************************************************************************
  //                           Draw template for m1                             
  //****************************************************************************

  TH1D* h1_Jpsi_m1_fakeData = (TH1D*)w->pdf("Jpsi_m1")->createHistogram("m1",m_bins);
  TH1D* h1_Jpsi_m1 = (TH1D*)w->pdf("Jpsi_m1")->createHistogram("m1",1000);
  h1_Jpsi_m1->Scale( nEvents_Jpsi * h1_Jpsi_m1_fakeData->Integral("width") / h1_Jpsi_m1->Integral("width") );
  double maxY = 1.2*h1_Jpsi_m1->GetMaximum();
  cout << "maxY = " << maxY << endl;
  h1_Jpsi_m1->GetYaxis()->SetRangeUser(0.,maxY);
  h1_Jpsi_m1->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_Jpsi_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
  h1_Jpsi_m1->SetLineWidth(3);
  h1_Jpsi_m1->SetLineColor(kRed);
  
  TCanvas * c_m1 = new TCanvas("c_m1", "c_m1");
	c_m1->cd();
	h1_Jpsi_m1->Draw("L");
  txtHeader->Draw();
  c_m1->SaveAs("template1D_Jpsi_m1.pdf");
  c_m1->SaveAs("template1D_Jpsi_m1.png");

  //****************************************************************************
  //                           Draw template for m2                             
  //****************************************************************************
  
  TH1D* h1_Jpsi_m2_fakeData = (TH1D*)w->pdf("Jpsi_m2")->createHistogram("m2",m_bins);
  TH1D* h1_Jpsi_m2 = (TH1D*)w->pdf("Jpsi_m2")->createHistogram("m2",1000);
  h1_Jpsi_m2->Scale( nEvents_Jpsi * h1_Jpsi_m2_fakeData->Integral("width") / h1_Jpsi_m2->Integral("width") );
  h1_Jpsi_m2->GetYaxis()->SetRangeUser(0.,maxY);
  h1_Jpsi_m2->GetYaxis()->SetTitle("Events / (0.05 GeV/#it{c}^{2})");
  h1_Jpsi_m2->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
  h1_Jpsi_m2->SetLineWidth(3);
  h1_Jpsi_m2->SetLineColor(kRed);
  
  TCanvas * c_m2 = new TCanvas("c_m2", "c_m2");
	c_m2->cd();
	h1_Jpsi_m2->Draw("L");
  txtHeader->Draw();
  c_m2->SaveAs("template1D_Jpsi_m2.pdf");
  c_m2->SaveAs("template1D_Jpsi_m2.png");
  
	//****************************************************************************
  //                         Draw 2D template m1 x m2                           
  //****************************************************************************
  
  TH2D* h2_Jpsi_2D = (TH2D*)w->pdf("Jpsi_2D")->createHistogram("m1,m2",2.0*m_bins,2.0*m_bins);
  cout << "Jpsi 2D template integral: " << h2_Jpsi_2D->Integral() << std::endl;
  
  TH2D* h2_Jpsi_2D_diagonal    = (TH2D*)w->pdf("Jpsi_2D")->createHistogram("m1,m2",1000,1000);
  TH2D* h2_Jpsi_2D_offDiagonal = (TH2D*)w->pdf("Jpsi_2D")->createHistogram("m1,m2",1000,1000);
  
	for(int i=1;i<=1000;i++) {
	  for(int j=1;j<=1000;j++) {
	    double m_1 = h2_Jpsi_2D_diagonal->GetXaxis()->GetBinCenter(i);
	    double m_2 = h2_Jpsi_2D_diagonal->GetYaxis()->GetBinCenter(j);
	    if ( fabs(m_1 - m_2) < 5.*(0.026 + 0.013*(m_1 + m_2)/2.) ) {
	      h2_Jpsi_2D_diagonal->SetBinContent(i,j,0.);
      } else {
        h2_Jpsi_2D_offDiagonal->SetBinContent(i,j,0.);
      }
    }
	}
	cout << "Jpsi 2D template offDiagonal integral: " << h2_Jpsi_2D_offDiagonal->Integral() << std::endl;
	cout << "Jpsi 2D template diagonal integral:    " << h2_Jpsi_2D_diagonal->Integral() << std::endl;

	h2_Jpsi_2D->Scale(nEvents_Jpsi);
	cout << "Scaled 2D template integral: " << h2_Jpsi_2D->Integral() << std::endl;

  TCanvas * c_template2D_m1_vs_m2 = new TCanvas("c_template2D_m1_vs_m2", "c_template2D_m1_vs_m2",0,1320,1004,928);
  c_template2D_m1_vs_m2->SetCanvasSize(1000,900);
  c_template2D_m1_vs_m2->SetLeftMargin(0.126);
  c_template2D_m1_vs_m2->SetRightMargin(0.154);
	c_template2D_m1_vs_m2->cd();
	c_template2D_m1_vs_m2->SetLogz();
	
	h2_Jpsi_2D->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h2_Jpsi_2D->GetXaxis()->SetTitleOffset(0.93);
	h2_Jpsi_2D->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h2_Jpsi_2D->GetYaxis()->SetTitleOffset(0.95);
  h2_Jpsi_2D->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_Jpsi_2D->GetZaxis()->CenterTitle(true);
  h2_Jpsi_2D->GetZaxis()->SetLabelFont(42);
  h2_Jpsi_2D->GetZaxis()->SetLabelOffset(-0.005);
  h2_Jpsi_2D->GetZaxis()->SetLabelSize(0.044);
  h2_Jpsi_2D->GetZaxis()->SetTitleSize(0.044);
  h2_Jpsi_2D->GetZaxis()->SetTitleOffset(1.2);
  h2_Jpsi_2D->GetZaxis()->SetTitleFont(42);
	
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
	h2_Jpsi_2D->SetContour(NCont);
	
  h2_Jpsi_2D->Draw("Cont4 Colz");
  
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
  
  TH2D *h2_Jpsi_2D_Fake = new TH2D("h2_Jpsi_2D_Fake", "h2_Jpsi_2D_Fake", 2.0*m_bins, m_min, m_max, 2.0*m_bins, m_min, m_max);
  h2_Jpsi_2D_Fake->Draw();
  
  double diagonal_x1 = ( (1.0+kB/2.0)*m_min + kA )/( 1.0 - kB/2.0 );
  double diagonal_x2 = ( (1.0-kB/2.0)*m_max - kA )/( 1.0 + kB/2.0 );
  std::cout << "diagonal_x1 " << diagonal_x1 << std::endl;
  std::cout << "diagonal_x2 " << diagonal_x2 << std::endl;
  
  TLine *line1 = new TLine(m_min, diagonal_x1, diagonal_x2, m_max);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(9);
  line1->SetLineWidth(2);
  line1->Draw();
  line2 = new TLine(diagonal_x1,m_min,m_max,diagonal_x2);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(9);
  line2->SetLineWidth(2);
  line2->Draw();
  
  txtHeader->Draw();
  
	c_template2D_m1_vs_m2->SaveAs("template2D_Jpsi_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("template2D_Jpsi_m1_vs_m2.png");
  
  TCanvas * c_template2D_m1_vs_m2_3DLog = new TCanvas("c_template2D_m1_vs_m2_3DLog", "c_template2D_m1_vs_m2_3DLog");
	c_template2D_m1_vs_m2_3DLog->cd();
	c_template2D_m1_vs_m2_3DLog->SetLogz();
	TH2D *h2_Jpsi_3DLog = new TH2D( *h2_Jpsi_2D );
	h2_Jpsi_3DLog->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h2_Jpsi_3DLog->GetXaxis()->SetTitleOffset(1.5);
	h2_Jpsi_3DLog->GetXaxis()->CenterTitle(1);
	h2_Jpsi_3DLog->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h2_Jpsi_3DLog->GetYaxis()->SetTitleOffset(1.7);
	h2_Jpsi_3DLog->GetYaxis()->CenterTitle(1);
  h2_Jpsi_3DLog->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_Jpsi_3DLog->GetZaxis()->CenterTitle(true);
  h2_Jpsi_3DLog->GetZaxis()->SetLabelFont(42);
//  h2_Jpsi_3DLog->GetZaxis()->SetLabelOffset(-0.005);
  h2_Jpsi_3DLog->GetZaxis()->SetLabelSize(0.044);
  h2_Jpsi_3DLog->GetZaxis()->SetTitleSize(0.044);
  h2_Jpsi_3DLog->GetZaxis()->SetTitleOffset(1.8);
  h2_Jpsi_3DLog->GetZaxis()->SetTitleFont(42);
	h2_Jpsi_3DLog->GetZaxis()->SetRangeUser(0.00004,0.1);
  h2_Jpsi_3DLog->Draw("Surf1 FB");
  txtHeader->Draw();
  
	c_template2D_m1_vs_m2_3DLog->SaveAs("template2D_Jpsi_m1_vs_m2_3D.pdf");
  c_template2D_m1_vs_m2_3DLog->SaveAs("template2D_Jpsi_m1_vs_m2_3D.png");

}
