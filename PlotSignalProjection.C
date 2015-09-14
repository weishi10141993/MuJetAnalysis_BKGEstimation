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

void PlotSignalProjection() {
  
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
	
	RooWorkspace* w_mA04GeV = new RooWorkspace("w_mA04GeV");
	TFile* file = new TFile("ws_H2A4Mu_mA_0.4000_GeV.root");
	RooWorkspace *w_mA04GeV = (RooWorkspace*) file->Get("w_H2A4Mu");
	
	RooWorkspace* w_mA2GeV = new RooWorkspace("w_mA2GeV");
	TFile* file = new TFile("ws_H2A4Mu_mA_2.0000_GeV.root");
	RooWorkspace *w_mA2GeV = (RooWorkspace*) file->Get("w_H2A4Mu");

  const double       m_min  = 0.2113;
	const double       m_max  = 3.5536;
	const unsigned int m_bins = 66;
	const unsigned int m_fine_bins = 1000; // use 4000 for final plot
  
  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

//  double nEvents_Jpsi = 2.0;
  double nEvents_mA04GeV = 1.0;
  double nEvents_mA2GeV  = 1.0;
  
	//****************************************************************************
  //                         Draw 2D template m1 x m2                           
  //****************************************************************************
  
  TH2D* h2_mA04GeV = (TH2D*)w_mA04GeV->pdf("signal")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
	h2_mA04GeV->Scale(nEvents_mA04GeV);
  
  TH1D* h1_diagonal_dummy = new TH1D("h1_diagonal_dummy", "h1_diagonal_dummy", 2.0*m_bins, m_min,m_max);
  
  // Projection of signal region on diagonal (m1+m2)/sqrt(2)
  TH1D* h1_mA04GeV_diagonal = new TH1D("h1_mA04GeV_diagonal", "h1_mA04GeV_diagonal", 20, m_min,0.7);
  // Projection of signal region on horizontal axis m1
  TH1D* h1_mA04GeV_m1 = new TH1D("h1_mA04GeV_m1", "h1_mA04GeV_m1", 2.0*m_bins, m_min, m_max);
  // Projection of signal region on vertical axis m2
  TH1D* h1_mA04GeV_m2 = new TH1D("h1_mA04GeV_m2", "h1_mA04GeV_m2", 2.0*m_bins, m_min, m_max);
  
  TH2D* h2_mA2GeV = (TH2D*)w_mA2GeV->pdf("signal")->createHistogram("m1,m2",m_fine_bins,m_fine_bins);
	h2_mA2GeV->Scale(nEvents_mA2GeV);
  
  // Projection of signal region on diagonal (m1+m2)/sqrt(2)
  TH1D* h1_mA2GeV_diagonal = new TH1D("h1_mA2GeV_diagonal", "h1_mA2GeV_diagonal", 24, 1.6, 2.2);
  // Projection of signal region on horizontal axis m1
  TH1D* h1_mA2GeV_m1 = new TH1D("h1_mA2GeV_m1", "h1_mA2GeV_m1", 2.0*m_bins, m_min, m_max);
  // Projection of signal region on vertical axis m2
  TH1D* h1_mA2GeV_m2 = new TH1D("h1_mA2GeV_m2", "h1_mA2GeV_m2", 2.0*m_bins, m_min, m_max);
  
	for(int i=1;i<=m_fine_bins;i++) {
	  for(int j=1;j<=m_fine_bins;j++) {
	    double m_1 = h2_mA04GeV->GetXaxis()->GetBinCenter(i);
	    double m_2 = h2_mA04GeV->GetYaxis()->GetBinCenter(j);
	    if ( fabs(m_1 - m_2) < 5.*(0.026 + 0.013*(m_1 + m_2)/2.) && m_1 > m_min && m_1 < m_max && m_2 > m_min && m_2 < m_max ) {
	      double m_diagonal = (m_1 + m_2) / 2.0;
	      
	      double weight_mA04GeV = h2_mA04GeV->GetBinContent(i,j);
        h1_mA04GeV_diagonal->Fill(m_diagonal, weight_mA04GeV);
        h1_mA04GeV_m1->Fill(m_1, weight_mA04GeV);
        h1_mA04GeV_m2->Fill(m_2, weight_mA04GeV);
        
        double weight_mA2GeV = h2_mA2GeV->GetBinContent(i,j);
        h1_mA2GeV_diagonal->Fill(m_diagonal, weight_mA2GeV);
        h1_mA2GeV_m1->Fill(m_1, weight_mA2GeV);
        h1_mA2GeV_m2->Fill(m_2, weight_mA2GeV);
      }
    }
	}
  
  h1_diagonal_dummy->GetXaxis()->SetTitle("(m_{#mu#mu_{1}} + m_{#mu#mu_{2}})/2 [GeV/#it{c}^{2}]");
	h1_diagonal_dummy->GetXaxis()->SetTitleOffset(0.93);
	h1_diagonal_dummy->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_diagonal_dummy->GetYaxis()->SetTitleOffset(1.4);
  
  h1_mA04GeV_diagonal->Scale( 1.0 );
  cout << "Scaled 1D mA04GeV diagonal integral: " << h1_mA04GeV_diagonal->Integral() << endl;
  h1_mA04GeV_diagonal->GetXaxis()->SetTitle("(m_{#mu#mu_{1}} + m_{#mu#mu_{2}})/2 [GeV/#it{c}^{2}]");
	h1_mA04GeV_diagonal->GetXaxis()->SetTitleOffset(0.93);
	h1_mA04GeV_diagonal->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_mA04GeV_diagonal->GetYaxis()->SetTitleOffset(1.4);
  h1_mA04GeV_diagonal->SetLineWidth(3);
  h1_mA04GeV_diagonal->SetLineColor(kBlue);
  
  h1_mA04GeV_m1->Scale( 1.0 );
  cout << "Scaled 1D mA04GeV m1 integral: " << h1_mA04GeV_m1->Integral() << endl;
  h1_mA04GeV_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h1_mA04GeV_m1->GetXaxis()->SetTitleOffset(0.93);
	h1_mA04GeV_m1->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_mA04GeV_m1->GetYaxis()->SetTitleOffset(1.4);
  h1_mA04GeV_m1->SetLineWidth(3);
  h1_mA04GeV_m1->SetLineColor(kMagenta);
  
  h1_mA04GeV_m2->Scale( 1.0 );
  cout << "Scaled 1D mA04GeV m2 integral: " << h1_mA04GeV_m2->Integral() << endl;
  h1_mA04GeV_m2->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h1_mA04GeV_m2->GetXaxis()->SetTitleOffset(0.93);
	h1_mA04GeV_m2->GetYaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2}");
	h1_mA04GeV_m2->GetYaxis()->SetTitleOffset(1.4);
  h1_mA04GeV_m2->SetLineWidth(3);
  h1_mA04GeV_m2->SetLineColor(kMagenta);
  
  h1_mA2GeV_diagonal->Scale( 1.0 );
  cout << "Scaled 1D mA2GeV diagonal integral: " << h1_mA2GeV_diagonal->Integral() << endl;
  h1_mA2GeV_diagonal->GetXaxis()->SetTitle("(m_{#mu#mu_{1}} + m_{#mu#mu_{2}})/2 [GeV/#it{c}^{2}]");
	h1_mA2GeV_diagonal->GetXaxis()->SetTitleOffset(0.93);
	h1_mA2GeV_diagonal->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_mA2GeV_diagonal->GetYaxis()->SetTitleOffset(1.4);
  h1_mA2GeV_diagonal->SetLineWidth(3);
  h1_mA2GeV_diagonal->SetLineStyle(7);
  h1_mA2GeV_diagonal->SetLineColor(kMagenta);
  
  h1_mA2GeV_m1->Scale( 1.0 );
  cout << "Scaled 1D mA2GeV m1 integral: " << h1_mA2GeV_m1->Integral() << endl;
  h1_mA2GeV_m1->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV/#it{c}^{2}]");
	h1_mA2GeV_m1->GetXaxis()->SetTitleOffset(0.93);
	h1_mA2GeV_m1->GetYaxis()->SetTitle("Events / 0.025 GeV/#it{c}^{2}");
	h1_mA2GeV_m1->GetYaxis()->SetTitleOffset(1.4);
  h1_mA2GeV_m1->SetLineWidth(3);
  h1_mA2GeV_m1->SetLineColor(kMagenta);
  
  h1_mA2GeV_m2->Scale( 1.0 );
  cout << "Scaled 1D mA2GeV m2 integral: " << h1_mA2GeV_m2->Integral() << endl;
  h1_mA2GeV_m2->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV/#it{c}^{2}]");
	h1_mA2GeV_m2->GetXaxis()->SetTitleOffset(0.93);
	h1_mA2GeV_m2->GetYaxis()->SetTitle("Events / (0.025 GeV/#it{c}^{2}");
	h1_mA2GeV_m2->GetYaxis()->SetTitleOffset(1.4);
  h1_mA2GeV_m2->SetLineWidth(3);
  h1_mA2GeV_m2->SetLineColor(kMagenta);
  
  TCanvas * c_template1D_signals_InSignalArea_diagonalProjection = new TCanvas("c_template1D_signals_InSignalArea_diagonalProjection", "c_template1D_signals_InSignalArea_diagonalProjection",0,1320,1004,928);
  c_template1D_signals_InSignalArea_diagonalProjection->SetCanvasSize(1000,900);
	c_template1D_signals_InSignalArea_diagonalProjection->cd();
  
  h1_diagonal_dummy->SetMaximum(0.5);
  h1_diagonal_dummy->Draw();
  h1_mA04GeV_diagonal->DrawNormalized("Csame");
  h1_mA2GeV_diagonal->DrawNormalized("Csame");
  
  txtHeader->Draw();
  
  TLegend * l_signal = new TLegend(0.5,0.8,0.8,0.9);
  l_signal->SetFillColor(kWhite);
  l_signal->SetMargin(0.4);
  l_signal->SetBorderSize(0);
  l_signal->SetTextFont(42);
//  l_signal->SetTextSize(0.035);
//  l_signal->SetHeader("Masses");
  l_signal->AddEntry(h1_mA04GeV_diagonal, "m_{a} = 0.4 GeV", "L");
  l_signal->AddEntry(h1_mA2GeV_diagonal,  "m_{a} = 2.0 GeV", "L");
  l_signal->Draw();
  
	c_template1D_signals_InSignalArea_diagonalProjection->SaveAs("template1D_signals_InSignalArea_diagonalProjection.pdf");
  c_template1D_signals_InSignalArea_diagonalProjection->SaveAs("template1D_signals_InSignalArea_diagonalProjection.png");
  
}
