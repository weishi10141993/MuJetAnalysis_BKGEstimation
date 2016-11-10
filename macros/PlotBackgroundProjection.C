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

void PlotBackgroundProjection() {
  
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
  txtHeader->SetHeader("CMS Prelim. 2012  #sqrt{s} = 8 TeV   L_{int} = 20.65 fb^{-1}");
//  txtHeader->SetHeader("CMS 2012   #sqrt{s} = 8 TeV   L_{int} = 20.65 fb^{-1}");
//  txtHeader->Draw();
	
	RooWorkspace* w = new RooWorkspace("w");
	
	TFile* file = new TFile("ws.root");
	RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
	const double       m_max  = 3.5536;
	const unsigned int m_bins = 66;
	const unsigned int m_fine_bins = 4000; // use 4000 for final plot
  
  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

//  double nEvents_Jpsi = 2.0;
  double nEvents_Jpsi = 0.058;
  
	//****************************************************************************
  //                         Draw 2D template m1 x m2                           
  //****************************************************************************
  
  TH2D* h2_Jpsi_2D = (TH2D*)w->pdf("Jpsi_2D")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
	h2_Jpsi_2D->Scale(nEvents_Jpsi);
  
  TH2D* h2_Template2D = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
  cout << "2D template integral: " << h2_Template2D->Integral() << std::endl;
  
  TH2D* h2_Template2D_diagonal    = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
  TH2D* h2_Template2D_offDiagonal = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
  
  // Projection of signal region on diagonal (m1+m2)/sqrt(2)
  TH1D* h1_BBbar_diagonal = new TH1D("h1_BBbar_diagonal", "h1_BBbar_diagonal", 2.0*m_bins, m_min, m_max);
  TH1D* h1_Jpsi_diagonal = new TH1D("h1_Jpsi_diagonal", "h1_Jpsi_diagonal", 2.0*m_bins, m_min, m_max);
  // Projection of signal region on horizontal axis m1
  TH1D* h1_BBbar_m1 = new TH1D("h1_BBbar_m1", "h1_BBbar_m1", 2.0*m_bins, m_min, m_max);
  TH1D* h1_Jpsi_m1 = new TH1D("h1_Jpsi_m1", "h1_Jpsi_m1", 2.0*m_bins, m_min, m_max);
  // Projection of signal region on vertical axis m2
  TH1D* h1_BBbar_m2 = new TH1D("h1_BBbar_m2", "h1_BBbar_m2", 2.0*m_bins, m_min, m_max);
  TH1D* h1_Jpsi_m2 = new TH1D("h1_Jpsi_m2", "h1_Jpsi_m2", 2.0*m_bins, m_min, m_max);
  
	for(int i=1;i<=m_fine_bins;i++) {
	  for(int j=1;j<=m_fine_bins;j++) {
	    double m_1 = h2_Template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
	    double m_2 = h2_Template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
	    if ( fabs(m_1 - m_2) < 5.*(0.026 + 0.013*(m_1 + m_2)/2.) && m_1 > m_min && m_1 < m_max && m_2 > m_min && m_2 < m_max ) {
	      h2_Template2D_offDiagonal->SetBinContent(i,j,0.);

	      double weight_Template2D = h2_Template2D_diagonal->GetBinContent(i,j);
        double m_diagonal = (m_1 + m_2) / 2.0;
        h1_BBbar_diagonal->Fill(m_diagonal, weight_Template2D);
        h1_BBbar_m1->Fill(m_1, weight_Template2D);
        h1_BBbar_m2->Fill(m_2, weight_Template2D);
        double weight_Jpsi_2D = h2_Jpsi_2D->GetBinContent(i,j);
        h1_Jpsi_diagonal->Fill(m_diagonal, weight_Jpsi_2D);
        h1_Jpsi_m1->Fill(m_1, weight_Jpsi_2D);
        h1_Jpsi_m2->Fill(m_2, weight_Jpsi_2D);
      } else {
        h2_Template2D_diagonal->SetBinContent(i,j,0.);
      }
    }
	}
	cout << "Template2D_offDiagonal integral: " << h2_Template2D_offDiagonal->Integral() << std::endl;
	cout << "Template2D_diagonal integral:    " << h2_Template2D_diagonal->Integral() << std::endl;

	TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
	cout << "Data dimudimu_control_Iso_offDiagonal integral:" << h2_dimudimu_control_Iso_offDiagonal_2D->Integral() << std::endl;

	h2_Template2D->Scale(1.0/h2_Template2D->Integral()*h2_dimudimu_control_Iso_offDiagonal_2D->Integral()*(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral());
	cout << "Scaled 2D template integral: " << h2_Template2D->Integral() << std::endl;
  
  h1_BBbar_diagonal->Scale( 1.0/h1_BBbar_diagonal->Integral()*h2_dimudimu_control_Iso_offDiagonal_2D->Integral()*h2_Template2D_diagonal->Integral() /h2_Template2D_offDiagonal->Integral() );
  cout << "Scaled 1D BBbar diagonal integral: " << h1_BBbar_diagonal->Integral() << endl;
  h1_BBbar_diagonal->GetXaxis()->SetTitle("(m_{#mu#mu_{1}} + m_{#mu#mu_{2}})/2 [GeV/#it{c}^{2}]");
	h1_BBbar_diagonal->GetXaxis()->SetTitleOffset(0.93);
	h1_BBbar_diagonal->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_BBbar_diagonal->GetYaxis()->SetTitleOffset(1.4);
  h1_BBbar_diagonal->SetLineWidth(3);
  h1_BBbar_diagonal->SetLineColor(kMagenta);
  
  h1_BBbar_m1->Scale( 1.0/h1_BBbar_m1->Integral()*h2_dimudimu_control_Iso_offDiagonal_2D->Integral()*h2_Template2D_diagonal->Integral() /h2_Template2D_offDiagonal->Integral() );
  cout << "Scaled 1D BBbar m1 integral: " << h1_BBbar_m1->Integral() << endl;
  h1_BBbar_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h1_BBbar_m1->GetXaxis()->SetTitleOffset(0.93);
	h1_BBbar_m1->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_BBbar_m1->GetYaxis()->SetTitleOffset(1.4);
  h1_BBbar_m1->SetLineWidth(3);
  h1_BBbar_m1->SetLineColor(kMagenta);
  
  h1_BBbar_m2->Scale( 1.0/h1_BBbar_m2->Integral()*h2_dimudimu_control_Iso_offDiagonal_2D->Integral()*h2_Template2D_diagonal->Integral() /h2_Template2D_offDiagonal->Integral() );
  cout << "Scaled 1D BBbar m2 integral: " << h1_BBbar_m2->Integral() << endl;
  h1_BBbar_m2->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h1_BBbar_m2->GetXaxis()->SetTitleOffset(0.93);
	h1_BBbar_m2->GetYaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2}");
	h1_BBbar_m2->GetYaxis()->SetTitleOffset(1.4);
  h1_BBbar_m2->SetLineWidth(3);
  h1_BBbar_m2->SetLineColor(kMagenta);
  
  h1_Jpsi_diagonal->Scale(1.0/h1_Jpsi_diagonal->Integral()*nEvents_Jpsi);
  cout << "Scaled 1D Jpsi diagonal integral: " << h1_Jpsi_diagonal->Integral() << endl;
  h1_Jpsi_diagonal->GetXaxis()->SetTitle("(m_{#mu#mu_{1}} + m_{#mu#mu_{2}})/2 [GeV/#it{c}^{2}]");
	h1_Jpsi_diagonal->GetXaxis()->SetTitleOffset(0.93);
	h1_Jpsi_diagonal->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_Jpsi_diagonal->GetYaxis()->SetTitleOffset(1.4);
  h1_Jpsi_diagonal->SetLineWidth(3);
  h1_Jpsi_diagonal->SetLineColor(kBlue);
  
  h1_Jpsi_m1->Scale( 1.0/h1_Jpsi_m1->Integral()*nEvents_Jpsi );
  cout << "Scaled 1D Jpsi m1 integral: " << h1_Jpsi_m1->Integral() << endl;
  h1_Jpsi_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h1_Jpsi_m1->GetXaxis()->SetTitleOffset(0.93);
	h1_Jpsi_m1->GetYaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2}");
	h1_Jpsi_m1->GetYaxis()->SetTitleOffset(1.4);
  h1_Jpsi_m1->SetLineWidth(3);
  h1_Jpsi_m1->SetLineColor(kBlue);
  
  h1_Jpsi_m2->Scale( 1.0/h1_Jpsi_m2->Integral()*nEvents_Jpsi );
  cout << "Scaled 1D Jpsi m2 integral: " << h1_Jpsi_m2->Integral() << endl;
  h1_Jpsi_m2->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h1_Jpsi_m2->GetXaxis()->SetTitleOffset(0.93);
	h1_Jpsi_m2->GetYaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2}");
	h1_Jpsi_m2->GetYaxis()->SetTitleOffset(1.4);
  h1_Jpsi_m2->SetLineWidth(3);
  h1_Jpsi_m2->SetLineColor(kBlue);
  
  TH1D * h1_background_diagonal = new TH1D( *h1_BBbar_diagonal );
  h1_background_diagonal->Add(h1_Jpsi_diagonal);
  h1_background_diagonal->GetXaxis()->SetTitle("(m_{#mu#mu_{1}} + m_{#mu#mu_{2}})/2 [GeV/#it{c}^{2}]");
	h1_background_diagonal->GetXaxis()->SetTitleOffset(0.93);
	h1_background_diagonal->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_background_diagonal->GetYaxis()->SetTitleOffset(1.4);
  h1_background_diagonal->SetLineWidth(3);
  h1_background_diagonal->SetLineColor(kRed);
  
  TH1D * h1_background_m1 = new TH1D( *h1_BBbar_m1 );
  h1_background_m1->Add(h1_Jpsi_m1);
  h1_background_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h1_background_m1->GetXaxis()->SetTitleOffset(0.93);
	h1_background_m1->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_background_m1->GetYaxis()->SetTitleOffset(1.4);
  h1_background_m1->SetLineWidth(3);
  h1_background_m1->SetLineColor(kRed);
  
  TH1D * h1_background_m2 = new TH1D( *h1_BBbar_m2 );
  h1_background_m2->Add(h1_Jpsi_m2);
  h1_background_m2->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h1_background_m2->GetXaxis()->SetTitleOffset(0.93);
	h1_background_m2->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_background_m2->GetYaxis()->SetTitleOffset(1.4);
  h1_background_m2->SetLineWidth(3);
  h1_background_m2->SetLineColor(kRed);
  
  TCanvas * c_template1D_backgroundInSignalArea_diagonalProjection = new TCanvas("c_template1D_backgroundInSignalArea_diagonalProjection", "c_template1D_backgroundInSignalArea_diagonalProjection",0,1320,1004,928);
  c_template1D_backgroundInSignalArea_diagonalProjection->SetCanvasSize(1000,900);
	c_template1D_backgroundInSignalArea_diagonalProjection->cd();
//	c_template1D_backgroundInSignalArea_diagonalProjection->SetLogy();

  h1_background_diagonal->Draw("C");
//  h1_BBbar_diagonal->Draw("Csame");
//  h1_Jpsi_diagonal->Draw("Csame");
  
  txtHeader->Draw();
  
	c_template1D_backgroundInSignalArea_diagonalProjection->SaveAs("template1D_backgroundInSignalArea_diagonalProjection.pdf");
  c_template1D_backgroundInSignalArea_diagonalProjection->SaveAs("template1D_backgroundInSignalArea_diagonalProjection.png");
  
  TCanvas * c_template1D_backgroundInSignalArea_m1Projection = new TCanvas("c_template1D_backgroundInSignalArea_m1Projection", "c_template1D_backgroundInSignalArea_m1Projection",0,1320,1004,928);
  c_template1D_backgroundInSignalArea_m1Projection->SetCanvasSize(1000,900);
	c_template1D_backgroundInSignalArea_m1Projection->cd();
//	c_template1D_backgroundInSignalArea_m1Projection->SetLogy();
  
  h1_background_m1->Draw("C");
//  h1_BBbar_m1->Draw("Csame");
//  h1_Jpsi_m1->Draw("Csame");
  
  txtHeader->Draw();
  
	c_template1D_backgroundInSignalArea_m1Projection->SaveAs("template1D_backgroundInSignalArea_m1Projection.pdf");
  c_template1D_backgroundInSignalArea_m1Projection->SaveAs("template1D_backgroundInSignalArea_m1Projection.png");
  
  TCanvas * c_template1D_backgroundInSignalArea_m2Projection = new TCanvas("c_template1D_backgroundInSignalArea_m2Projection", "c_template1D_backgroundInSignalArea_m2Projection",0,1320,1004,928);
  c_template1D_backgroundInSignalArea_m2Projection->SetCanvasSize(1000,900);
	c_template1D_backgroundInSignalArea_m2Projection->cd();
//	c_template1D_backgroundInSignalArea_m2Projection->SetLogy();
  
  h1_background_m2->Draw("C");
//  h1_BBbar_m2->Draw("Csame");
//  h1_Jpsi_m2->Draw("Csame");
  
  txtHeader->Draw();
  
	c_template1D_backgroundInSignalArea_m2Projection->SaveAs("template1D_backgroundInSignalArea_m2Projection.pdf");
  c_template1D_backgroundInSignalArea_m2Projection->SaveAs("template1D_backgroundInSignalArea_m2Projection.png");
  
  
  TH2D * h2_background = new TH2D( *h2_Jpsi_2D );
  h2_background->Add( h2_Template2D );
  
  TCanvas * c_template2D_m1_vs_m2 = new TCanvas("c_template2D_m1_vs_m2", "c_template2D_m1_vs_m2",0,1320,1004,928);
  c_template2D_m1_vs_m2->SetCanvasSize(1000,900);
  c_template2D_m1_vs_m2->SetLeftMargin(0.126);
  c_template2D_m1_vs_m2->SetRightMargin(0.154);
	c_template2D_m1_vs_m2->cd();
	c_template2D_m1_vs_m2->SetLogz();
	
	h2_background->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h2_background->GetXaxis()->SetTitleOffset(0.93);
	h2_background->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h2_background->GetYaxis()->SetTitleOffset(0.95);
  h2_background->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_background->GetZaxis()->CenterTitle(true);
  h2_background->GetZaxis()->SetLabelFont(42);
  h2_background->GetZaxis()->SetLabelOffset(-0.005);
  h2_background->GetZaxis()->SetLabelSize(0.044);
  h2_background->GetZaxis()->SetTitleSize(0.044);
  h2_background->GetZaxis()->SetTitleOffset(1.2);
  h2_background->GetZaxis()->SetTitleFont(42);
	
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
	h2_background->SetContour(NCont);
	
  h2_background->Draw("Cont4 Colz");
  
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
  
  
  TH2D* h2_signal = (TH2D*)w->data("ds_dimudimu_signal_2D")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
	cout << "Data signal integral:" << h2_signal->Integral() << std::endl;
  
  h2_signal->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h2_signal->GetXaxis()->SetTitleOffset(0.93);
	h2_signal->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h2_signal->GetYaxis()->SetTitleOffset(0.95);
  h2_signal->GetZaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2} x 0.025 GeV/#it{c}^{2})");
  h2_signal->GetZaxis()->CenterTitle(true);
  h2_signal->GetZaxis()->SetLabelFont(42);
  h2_signal->GetZaxis()->SetLabelOffset(-0.005);
  h2_signal->GetZaxis()->SetLabelSize(0.044);
  h2_signal->GetZaxis()->SetTitleSize(0.044);
  h2_signal->GetZaxis()->SetTitleOffset(1.2);
  h2_signal->GetZaxis()->SetTitleFont(42);
  h2_signal->SetMarkerColor(kBlack);
  h2_signal->SetMarkerStyle(20);
  h2_signal->SetMarkerSize(3.0);
  h2_signal->Draw("same");
  
  TH2D * h2_signal_tmp = new TH2D( *h2_signal);
  h2_signal_tmp->SetMarkerColor(kWhite);
  h2_signal_tmp->SetMarkerStyle(20);
  h2_signal_tmp->SetMarkerSize(2.0);
//  h2_signal_tmp->Draw("same");
  
  double diagonal_x1 = ( (1.0+kB/2.0)*m_min + kA )/( 1.0 - kB/2.0 );
  double diagonal_x2 = ( (1.0-kB/2.0)*m_max - kA )/( 1.0 + kB/2.0 );
  std::cout << "diagonal_x1 " << diagonal_x1 << std::endl;
  std::cout << "diagonal_x2 " << diagonal_x2 << std::endl;
  
  TLine *line1 = new TLine(m_min, diagonal_x1, diagonal_x2, m_max);
  line1->SetLineColor(0);
  line1->SetLineStyle(9);
  line1->SetLineWidth(2);
  line1->Draw();
  line2 = new TLine(diagonal_x1,m_min,m_max,diagonal_x2);
  line2->SetLineColor(0);
  line2->SetLineStyle(9);
  line2->SetLineWidth(2);
  line2->Draw();
  
  txtHeader->Draw();
  
	c_template2D_m1_vs_m2->SaveAs("template2D_signal_and_background_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("template2D_signal_and_background_m1_vs_m2.png");

}
