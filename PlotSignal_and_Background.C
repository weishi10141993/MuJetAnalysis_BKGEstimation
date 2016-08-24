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
#include "tdrStyle.C"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void PlotSignal_and_Background() {

  setTDRStyle();

  TCanvas * c_template2D_m1_vs_m2 = new TCanvas("c_template2D_m1_vs_m2", "c_template2D_m1_vs_m2",0,1320,1044,928);
  c_template2D_m1_vs_m2->SetCanvasSize(1040,900);
  c_template2D_m1_vs_m2->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2->SetRightMargin(0.148);
  c_template2D_m1_vs_m2->SetTopMargin(0.05);
  c_template2D_m1_vs_m2->cd();
  c_template2D_m1_vs_m2->SetLogz();

  TLegend *txtHeader = new TLegend(.13,.935,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader("CMS Prelim. 2015D  #sqrt{s} = 8 TeV   L_{int} = 2.6 fb^{-1}");


  TLegend *txtHeader_CMS = new TLegend(.2,0.81,0.7,0.86);
  txtHeader_CMS->SetFillColor(kWhite);
  txtHeader_CMS->SetFillStyle(0);
  txtHeader_CMS->SetBorderSize(0);
  txtHeader_CMS->SetTextFont(61);
  txtHeader_CMS->SetTextSize(0.055);
  txtHeader_CMS->SetTextAlign(12);
  txtHeader_CMS->SetHeader("CMS");

  TLegend *txtHeader_lumi = new TLegend(.5,0.95,0.85,1.0);
  txtHeader_lumi->SetFillColor(kWhite);
  txtHeader_lumi->SetFillStyle(0);
  txtHeader_lumi->SetBorderSize(0);
  txtHeader_lumi->SetTextFont(42);
  txtHeader_lumi->SetTextSize(0.042);
  txtHeader_lumi->SetTextAlign(32);
  txtHeader_lumi->SetHeader("2.6 fb^{-1} (13 TeV)");

  TFile* file = new TFile("ws.root");
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
  const double       m_max  = 3.5536;
  const unsigned int m_bins = 66;

  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

  double nEvents_Jpsi = 0.0; //From alfredo: double Jpsi. Not in bb selection. You add it using the Jpsi template from bb, using normalization from Alfredo.

  //****************************************************************************
  //                         Draw 2D template m1 x m2                           
  //****************************************************************************

  TH2D* h2_Jpsi_2D = (TH2D*)w->pdf("Jpsi_2D")->createHistogram("m1,m2",2.0*m_bins,2.0*m_bins);
  h2_Jpsi_2D->Scale(nEvents_Jpsi);  //Double J/Psi
  cout<<"Double J/Psi: "<<h2_Jpsi_2D->Integral()<<endl;

  //bb 2D template
  TH2D* h2_Template2D = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",2.0*m_bins,2.0*m_bins);
  cout << "2D template integral: " << h2_Template2D->Integral() << std::endl;

  TH2D* h2_Template2D_diagonal    = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);
  TH2D* h2_Template2D_offDiagonal = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);

  for(int i=1;i<=1000;i++) {
    for(int j=1;j<=1000;j++) {
	double m_1 = h2_Template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
	double m_2 = h2_Template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
	if ( fabs(m_1 - m_2) < (kA + kB*(m_1 + m_2)/2.) ) {
	  h2_Template2D_offDiagonal->SetBinContent(i,j,0.);
	} else {
	  h2_Template2D_diagonal->SetBinContent(i,j,0.);
	}
    }
  }
  cout << " -> Template2D_offDiagonal integral: " << h2_Template2D_offDiagonal->Integral() << std::endl;
  cout << " -> Template2D_diagonal integral:    " << h2_Template2D_diagonal->Integral() << std::endl;

  //Signal: ISO +off Diag
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2",1000,1000);
  cout << "Signal ISO + offDiag:" << h2_dimudimu_control_Iso_offDiagonal_2D->Integral() << std::endl;

  cout << "Scaled as: "<<h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<" / "<<h2_Template2D->Integral()<<" * "<<(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral()<<endl;
  //Scale to: DimuDimu_iso_offDiag / Template2D_Area (normalize to the off-diag part of the data) * bb_ALL/bb_offDiag (scale factor to pass from a normalization off-diag. to a normalization to the whole area.)
  h2_Template2D->Scale(h2_dimudimu_control_Iso_offDiagonal_2D->Integral()/h2_Template2D->Integral()*(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral());
  cout << "Scaled bb_2D template integral: " << h2_Template2D->Integral() << std::endl;

  TH2D * h2_background = new TH2D( *h2_Jpsi_2D );
  h2_background->Add( h2_Template2D );
  h2_background->GetXaxis()->SetTitle("m_{1 #mu#mu} [GeV]");
  h2_background->GetXaxis()->CenterTitle(true);
  h2_background->GetXaxis()->SetTitleOffset(0.93);
  h2_background->GetYaxis()->SetTitle("m_{2 #mu#mu} [GeV]");
  h2_background->GetYaxis()->CenterTitle(true);
  h2_background->GetYaxis()->SetTitleOffset(0.95);
  h2_background->GetZaxis()->SetTitle("Events / (0.025 GeV x 0.025 GeV)");
  h2_background->GetZaxis()->CenterTitle(true);
  h2_background->GetZaxis()->SetLabelFont(42);
  h2_background->GetZaxis()->SetLabelOffset(-0.005);
  h2_background->GetZaxis()->SetLabelSize(0.044);
  h2_background->GetZaxis()->SetTitleSize(0.044);
  h2_background->GetZaxis()->SetTitleOffset(1.2);
  h2_background->GetZaxis()->SetTitleFont(42);

  gStyle->SetPalette(52);

  //const Int_t NCont = 99;
  //const Int_t NRGBs = 5;
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //Int_t MyPalette[NCont];
  //for ( int i=0; i<NCont; i++ ) MyPalette[i] = FI+i;
  //gStyle->SetPalette(NCont, MyPalette);
  //h2_background->SetContour(NCont);
  h2_background->Draw("Cont4 Colz");

  //double diagonal_x1a = ( (1.0+kB/2.0)*m_min + kA )/( 1.0 - kB/2.0 );
  //double diagonal_x2a = ( (1.0-kB/2.0)*m_max - kA )/( 1.0 + kB/2.0 );
  //std::cout << "diagonal_x1a " << diagonal_x1a <<" "<<m_min<< std::endl;
  //std::cout << "diagonal_x2a " << diagonal_x2a <<" "<<m_max<< std::endl;
  //TLine *line1a = new TLine(m_min, diagonal_x1a, diagonal_x2a, m_max);
  //line1a->SetLineColor(0);
  //line1a->SetLineStyle(9);
  //line1a->SetLineWidth(2);
  //line1a->Draw(" ");
  //TLine *line2a = new TLine(diagonal_x1a,m_min,m_max,diagonal_x2a);
  //line2a->SetLineColor(0);
  //line2a->SetLineStyle(9);
  //line2a->SetLineWidth(2);
  ////line2a->Draw();
  //txtHeader_CMS->Draw();
  //txtHeader_lumi->Draw();

  c_template2D_m1_vs_m2->SaveAs("h2_background.pdf");// LP added
  c_template2D_m1_vs_m2->SaveAs("h2_background.png");

  // This required to draw scatter plot without LogZ 
  TPad* pad = new TPad("pad", "pad",0,0,1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.121);
  pad->SetRightMargin(0.148);
  pad->SetTopMargin(0.05);

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

////  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D_dummy = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_dummy", "h2_dimudimu_control_Iso_offDiagonal_2D_dummy", 1000, m_min, m_max, 1000, m_min, m_max);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetXaxis()->SetTitle("m_{1 #mu#mu} [GeV]");
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetXaxis()->CenterTitle(true);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetXaxis()->SetTitleOffset(0.93);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetYaxis()->SetTitle("m_{2 #mu#mu} [GeV]");
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetYaxis()->CenterTitle(true);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetYaxis()->SetTitleOffset(0.95);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetTitle("Events / (0.025 GeV x 0.025 GeV)");
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->CenterTitle(true);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetLabelFont(42);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetLabelOffset(-0.005);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetLabelSize(0.044);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetTitleSize(0.044);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetTitleOffset(1.2);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->GetZaxis()->SetTitleFont(42);
////  h2_dimudimu_control_Iso_offDiagonal_2D_dummy->Draw("");
////
////  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerColor(kBlack);
////  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerStyle(20);
////  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerSize(3.0);
////  //  h2_dimudimu_control_Iso_offDiagonal_2D->Draw("same");
////
////  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_tmp = new TH2D( *h2_dimudimu_control_Iso_offDiagonal_2D);
////  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerColor(kWhite);
////  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerStyle(20);
////  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerSize(2.0);
  //  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->Draw("same");

  //  -------------------2012----------------------
  // diagonal area
  //************************************
  //*    Row   *     massC *     massF *
  //************************************
  //*        0 * 0.3292791 * 0.2173234 *
  //************************************

  // off-diagonal area
  // Add new recovered points November 24 2014
  //  ************************************************************
  //  *    Row   * diMuonF_M * diMuonC_M * muJetC_mi * muJetF_mi *
  //  ************************************************************
  //  *    46463 * 1.9292701 *  2.414361 * 0.0020144 * 0.0002593 * x
  //  *    58035 * 2.2670404 * 1.7858715 * 0.0017198 * 0.0002601 * <-- new
  //  *   133195 * 1.4747750 * 2.9838371 * 0.0040867 * 0.0082415 * x
  //  *   154301 * 0.2553945 * 0.7071084 * 0.0006586 * 0.0148904 * x
  //  *   229405 * 0.8524319 * 3.1144735 * 0.0073431 * 0.0008649 * x
  //  *   245649 * 2.1131937 * 1.1794409 * 0.0033571 * 0.0056471 * x
  //  *   252508 * 2.6572125 * 0.6805820 * 0.0001614 * 0.0029167 * <-- new
  //  *   317676 * 2.0427842 * 2.3976006 * 0.0017126 * 0.0560831 * x
  //  *   413525 * 2.7330377 * 0.7881325 * 0.0079877 * 0.0034239 * x
  //  ************************************************************
  //  ************************************
  //  *    Row   *     massC *     massF *
  //  ************************************
  //  *        0 *  2.414361 * 1.9292701 *
  //  *        1 * 2.9838371 * 1.4747750 *
  //  *        2 * 1.1207234 * 1.7996511 * <-- lost
  //  *        3 * 0.7071084 * 0.2553945 *
  //  *        4 * 3.1144735 * 0.8524319 *
  //  *        5 * 1.1794409 * 2.1131937 *
  //  *        6 * 0.7881325 * 2.7330377 *
  //  *        7 * 2.3976006 * 2.0427842 *
  //  ************************************
  //  -------------------2015----------------------
  //------OffDiagonal SCAN------
  //************************************
  //*    Row   *     massC *     massF *
  //************************************
  //*        0 * 1.7474905 * 2.7977094 *    260576 * 388414933 *       188 *
  //*        1 * 1.2658947 * 2.3115363 *    256843 * 1.603e+09 *      1166 *
  //*        2 * 1.1979647 * 0.8495020 *    258714 *  14142382 *        10 *

  //************************************
  //------Signal SCAN------
  //************************************
  //*    Row   *     massC *     massF *
  //************************************

  // 2D histogram to nclude all points in new version of the analysis (November 2014). I can not use work space because I currently don't have it. So, all points are hard coded!
  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_points = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_points","h2_dimudimu_control_Iso_offDiagonal_2D_points", m_bins, m_min, m_max, m_bins, m_min, m_max);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->SetMarkerSize(3.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.2658947, 2.3115363);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.7474905, 2.7977094);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.1979647, 0.8495020);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp","h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp", m_bins, m_min, m_max, m_bins, m_min, m_max);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerSize(2.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.2658947, 2.3115363);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.7474905, 2.7977094);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1979647, 0.8495020);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Draw("same");

  TH2D* h2_signal = new TH2D("h2_signal","h2_signal", m_bins, m_min, m_max, m_bins, m_min, m_max);
  //h2_signal->Fill(0.4039067, 0.5626457);

  h2_signal->GetXaxis()->SetTitle("m_{1 #mu#mu} [GeV]");
  h2_signal->GetXaxis()->CenterTitle(true);
  h2_signal->GetXaxis()->SetTitleOffset(0.93);
  h2_signal->GetYaxis()->SetTitle("m_{2 #mu#mu} [GeV]");
  h2_signal->GetYaxis()->CenterTitle(true);
  h2_signal->GetYaxis()->SetTitleOffset(0.95);
  h2_signal->GetZaxis()->SetTitle("Events / (0.025 GeV x 0.025 GeV)");
  h2_signal->GetZaxis()->CenterTitle(true);
  h2_signal->GetZaxis()->SetLabelFont(42);
  h2_signal->GetZaxis()->SetLabelOffset(-0.005);
  h2_signal->GetZaxis()->SetLabelSize(0.044);
  h2_signal->GetZaxis()->SetTitleSize(0.044);
  h2_signal->GetZaxis()->SetTitleOffset(1.2);
  h2_signal->GetZaxis()->SetTitleFont(42);
  h2_signal->SetMarkerColor(kBlack);
  h2_signal->SetMarkerStyle(22);
  h2_signal->SetMarkerSize(3.0);
  h2_signal->Draw("same");

  TH2D * h2_signal_tmp = new TH2D( *h2_signal);
  Int_t col_tmp = TColor::GetColorTransparent(1, 0.1);
  h2_signal_tmp->SetMarkerColor(col_tmp);
  h2_signal_tmp->SetMarkerColor(kYellow);
  h2_signal_tmp->SetMarkerStyle(22);
  h2_signal_tmp->SetMarkerSize(2.0);
  h2_signal_tmp->Draw("same");

  double diagonal_x1 = ( (1.0+kB/2.0)*m_min + kA )/( 1.0 - kB/2.0 );
  double diagonal_x2 = ( (1.0-kB/2.0)*m_max - kA )/( 1.0 + kB/2.0 );
  std::cout << "diagonal_x1 " << diagonal_x1 <<" "<<m_min<< std::endl;
  std::cout << "diagonal_x2 " << diagonal_x2 <<" "<<m_max<< std::endl;

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

  txtHeader_CMS->Draw();
  txtHeader_lumi->Draw();

  c_template2D_m1_vs_m2->SaveAs("template2D_signal_and_background_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("template2D_signal_and_background_m1_vs_m2.png");
}
