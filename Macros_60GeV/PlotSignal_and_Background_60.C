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

void PlotSignal_and_Background_60() {

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
  txtHeader->SetHeader("CMS Prelim. 2015D  #sqrt{s} = 8 TeV   L_{int} = 2.83 fb^{-1}");


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
  txtHeader_lumi->SetHeader("2.83 fb^{-1} (13 TeV)");

  TFile* file = new TFile("ws_60.root");
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
  const double       m_max  = 60.;
  const unsigned int m_bins = 220;

  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

  double nEvents_Jpsi = 0.064; 

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

  float Area_2Jpsi=0., Area_NO2Jpsi=0., Area_2Jpsi_w=0., Area_NO2Jpsi_w=0.;
  float offd_Area_2Jpsi=0., offd_Area_NO2Jpsi=0., offd_Area_2Jpsi_w=0., offd_Area_NO2Jpsi_w=0.;
  for(int i=1;i<=1000;i++) {
    for(int j=1;j<=1000;j++) {
	double m_1 = h2_Template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
	double m_2 = h2_Template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
	if ( fabs(m_1 - m_2) < (kA + kB*(m_1 + m_2)/2.) ) {
	  if( fabs(m_1-3.1)<0.25 && fabs(m_2-3.1)<0.25 ){ Area_2Jpsi++; Area_2Jpsi_w+=h2_Template2D_offDiagonal->GetBinContent(i,j); }
	  else{ Area_NO2Jpsi++; Area_NO2Jpsi_w+=h2_Template2D_offDiagonal->GetBinContent(i,j); }
	  h2_Template2D_offDiagonal->SetBinContent(i,j,0.);
	}
	else {
	  if( fabs(m_1-3.1)<0.25 || fabs(m_2-3.1)<0.25 ){ offd_Area_2Jpsi++; offd_Area_2Jpsi_w+=h2_Template2D_diagonal->GetBinContent(i,j); }
	  else{ offd_Area_NO2Jpsi++; offd_Area_NO2Jpsi_w+=h2_Template2D_diagonal->GetBinContent(i,j); }
	  h2_Template2D_diagonal->SetBinContent(i,j,0.);
	}
    }
  }
  cout<<"Diagonal J/Psi probablity:"<<endl;
  cout<<"2J/Psi Area is " << Area_2Jpsi/(Area_NO2Jpsi+Area_2Jpsi) << " of the rest of the signal region ("<<Area_2Jpsi<<" "<<Area_NO2Jpsi<<")."<<endl;
  cout<<"2J/Psi Area weighted is " << Area_2Jpsi_w/(Area_NO2Jpsi_w+Area_2Jpsi_w) << " of the rest of the signal region ("<<Area_2Jpsi_w<<" "<<Area_NO2Jpsi_w<<")."<<endl;
  cout<<"NON Diagonal J/Psi probablity:"<<endl;
  cout<<"2J/Psi Area is " << offd_Area_2Jpsi/(offd_Area_NO2Jpsi+offd_Area_2Jpsi) << " of the rest of the signal region ("<<offd_Area_2Jpsi<<" "<<offd_Area_NO2Jpsi<<")."<<endl;
  cout<<"2J/Psi Area weighted is " << offd_Area_2Jpsi_w/(offd_Area_NO2Jpsi_w+offd_Area_2Jpsi_w) << " of the rest of the signal region ("<<offd_Area_2Jpsi_w<<" "<<offd_Area_NO2Jpsi_w<<")."<<endl;

  cout<<" -> Template2D_offDiagonal integral: "<<h2_Template2D_offDiagonal->Integral()<<endl;
  cout<<" -> Template2D_diagonal integral:    "<<h2_Template2D_diagonal->Integral()<<endl;

  //Signal: ISO +off Diag
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2",1000,1000);
  cout<<"Signal ISO + offDiag: " << h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<endl;

  cout<<"Scaled as: "<<h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<" / "<<h2_Template2D->Integral()<<" * "<<(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral()<<endl;
  //Scale to: DimuDimu_iso_offDiag / Template2D_Area (normalize to the off-diag part of the data) * bb_ALL/bb_offDiag (scale factor to pass from a normalization off-diag. to a normalization to the whole area.)
  h2_Template2D->Scale(h2_dimudimu_control_Iso_offDiagonal_2D->Integral()/h2_Template2D->Integral()*(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral());
  cout<<"Scaled bb_2D template integral: " << h2_Template2D->Integral() <<" That means " << h2_Template2D->Integral()-h2_dimudimu_control_Iso_offDiagonal_2D->Integral() << " events in signal region "<< std::endl;

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

  //gStyle->SetPalette(52); //Grey Scale
Double_t Red[2]    = { 1.00, 0.00};
Double_t Green[2]  = { 1.00, 0.00};
Double_t Blue[2]   = { 1.00, 0.00};
Double_t Length[2] = { 0.00, 1.00 };
Int_t nb=50;
TColor::CreateGradientColorTable(2,Length,Red,Green,Blue,nb);
h2_background->SetContour(nb);


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

  c_template2D_m1_vs_m2->SaveAs("figures_60/h2_background.pdf");// LP added
  c_template2D_m1_vs_m2->SaveAs("figures_60/h2_background.png");

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

  // 2D histogram to nclude all points in new version of the analysis (November 2014). I can not use work space because I currently don't have it. So, all points are hard coded!
  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_points = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_points","h2_dimudimu_control_Iso_offDiagonal_2D_points", m_bins, m_min, m_max, m_bins, m_min, m_max);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->SetMarkerSize(3.0);
  //2016
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(6.3394618, 4.8828525 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1058311, 21.046278 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(59.704967, 43.594131 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(59.289646, 45.105682 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(39.743816, 35.352249 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(56.168018, 0.7975030 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0989241, 0.4605314 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(11.209195, 34.068927 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.5686466, 23.308450 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(36.335540, 58.007171 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(22.249488, 36.519290 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(52.443889, 3.1008279 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(56.013874, 22.809450 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(28.313312, 20.380514 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.0721776, 39.336704 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(7.8420276, 32.317306 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1085057, 0.8863269 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.1467933, 16.528598 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.7905887, 25.334722 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3275566, 2.9306011 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(49.360717, 37.078338 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(9.0497751, 1.4616435 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(12.672730, 9.3065624 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(48.283596, 39.471145 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(52.593009, 42.147991 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.7733180, 2.9909498 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(15.016396, 0.4810310 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(15.892018, 19.808290 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(16.222017, 9.7123909 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6725933, 2.2082293 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(37.339859, 27.703325 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.9009132, 52.866878 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.9009132, 52.866878 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(46.883014, 21.346233 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(11.538404, 55.400951 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(28.569528, 26.510278 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0867397, 2.5654206 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(35.219913, 23.862997 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(23.124889, 26.847168 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.3725998, 1.5676055 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(25.234773, 41.194061 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(10.526959, 12.168872 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0961875, 2.6408567 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.6235761, 32.032074 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(52.359947, 46.339817 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(50.359432, 42.596439 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(31.921394, 27.564868 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1111078, 0.9851592 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6089491, 34.723232 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0964052, 1.9947856 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(33.847938, 37.534553 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(56.694118, 40.774990 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.5590274, 2.1472802 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(53.186309, 41.733455 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(42.560524, 18.694713 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(28.384697, 26.420057 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(43.457843, 2.2379193 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(42.045887, 21.514329 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(43.762439, 13.936329 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(23.154470, 47.492431 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.2796572, 1.3030284 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(53.643993, 25.872497 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(48.307292, 2.1617939 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.5262008, 47.590011 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(49.299934, 40.728637 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(20.202663, 43.788093 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(56.542530, 0.8674612 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(6.1615929, 47.600956 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.8130171, 15.474525 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.2652084, 19.713544 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.0735573, 3.2055511 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(23.523258, 34.144142 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1052675, 4.2853865 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2261323, 54.650974 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(35.938682, 22.310661 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(43.632278, 37.541214 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(54.621666, 36.506382 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(21.211456, 18.530464 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.9579395, 0.2526102 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(49.482158, 22.623146 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.5465407, 22.858202 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(38.069335, 29.628252 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(49.169578, 42.292602 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(39.564998, 53.640960 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(5.0760717, 12.942659 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.1172916, 0.3317201 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.6937513, 1.3820128 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(36.834030, 41.692997 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0868136, 2.0693941 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2487981, 28.159318 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(58.749683, 38.319442 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(54.396972, 40.095508 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.8807374, 3.0353784 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(10.761882, 1.7307828 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(10.045470, 3.1789214 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(44.929870, 0.9905142 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(59.720779, 30.977441 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(42.936012, 30.577076 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(46.079158, 28.352716 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(44.003917, 52.821025 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3446511, 18.510894 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(17.986646, 27.483783 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(59.111354, 22.785835 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(55.729534, 4.3312573 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(47.617763, 36.275802 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(44.650676, 38.945598 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.6326274, 25.486883 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2141819, 1.3356759 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.3732197, 3.1166579 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(24.886198, 2.5645513 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(54.308445, 25.413944 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(34.049346, 45.058528 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(53.701225, 27.893337 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.7097760, 1.9626674 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.7924485, 19.932649 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0195519, 38.367816 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(41.939609,  35.63274 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(22.228141, 39.524215 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.8375990, 23.524726 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(32.961212, 36.524200 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(9.4377613, 10.287744 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(30.673507, 23.828552 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(38.124122, 46.157466 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(53.382072, 33.468544 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(49.393001, 25.771896 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.3287160, 3.0902001 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(39.955749, 30.323211 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.7560508, 1.4533053 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(46.148002, 1.2616063 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(41.889728, 35.167430 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(54.938133, 36.032436 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3555525, 3.0687544 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(56.321323, 12.863502 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.4548411, 4.9593176 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(43.591377, 12.754740 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(37.563823, 9.5147047 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(33.541244, 1.4301798 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(59.563152, 17.541984 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(40.744159, 57.347049 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(36.621593, 2.5476176 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(46.718303,  36.70187 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6861963, 3.0976259 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(29.370142, 1.4475262 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(55.962074, 43.879081 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(45.212371, 22.505388 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(26.765331, 34.474536 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(59.777263, 32.886989 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(52.592639, 28.080614 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9678230, 21.956131 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(40.993480, 55.651126 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0214359, 55.590107 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(42.890899, 2.7383191 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(58.100662, 2.6796503 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6386389, 20.093421 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(38.900405, 2.5179841 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(29.298503, 35.910564 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(15.751144, 52.449920 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(28.119302, 0.3321138 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(35.297061, 20.921930 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0805194, 1.8598787 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.4601247, 0.5195755 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(22.912105, 1.0601193 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9651384, 1.4651534 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(30.186372, 3.0838427 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(8.0901327, 3.1087918 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(49.411384, 1.8853472 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.1265718, 4.8931441 );
h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(35.915580, 30.104633 ); 
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp","h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp", m_bins, m_min, m_max, m_bins, m_min, m_max);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerSize(2.0);
  //2016
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(6.3394618, 4.8828525 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1058311, 21.046278 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(59.704967, 43.594131 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(59.289646, 45.105682 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(39.743816, 35.352249 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(56.168018, 0.7975030 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0989241, 0.4605314 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(11.209195, 34.068927 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.5686466, 23.308450 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(36.335540, 58.007171 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(22.249488, 36.519290 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(52.443889, 3.1008279 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(56.013874, 22.809450 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(28.313312, 20.380514 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.0721776, 39.336704 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(7.8420276, 32.317306 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1085057, 0.8863269 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.1467933, 16.528598 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.7905887, 25.334722 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3275566, 2.9306011 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(49.360717, 37.078338 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(9.0497751, 1.4616435 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(12.672730, 9.3065624 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(48.283596, 39.471145 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(52.593009, 42.147991 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.7733180, 2.9909498 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(15.016396, 0.4810310 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(15.892018, 19.808290 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(16.222017, 9.7123909 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6725933, 2.2082293 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(37.339859, 27.703325 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.9009132, 52.866878 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.9009132, 52.866878 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(46.883014, 21.346233 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(11.538404, 55.400951 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(28.569528, 26.510278 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0867397, 2.5654206 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(35.219913, 23.862997 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(23.124889, 26.847168 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.3725998, 1.5676055 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(25.234773, 41.194061 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(10.526959, 12.168872 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0961875, 2.6408567 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.6235761, 32.032074 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(52.359947, 46.339817 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(50.359432, 42.596439 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(31.921394, 27.564868 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1111078, 0.9851592 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6089491, 34.723232 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0964052, 1.9947856 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(33.847938, 37.534553 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(56.694118, 40.774990 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.5590274, 2.1472802 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(53.186309, 41.733455 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(42.560524, 18.694713 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(28.384697, 26.420057 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(43.457843, 2.2379193 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(42.045887, 21.514329 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(43.762439, 13.936329 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(23.154470, 47.492431 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.2796572, 1.3030284 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(53.643993, 25.872497 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(48.307292, 2.1617939 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.5262008, 47.590011 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(49.299934, 40.728637 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(20.202663, 43.788093 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(56.542530, 0.8674612 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(6.1615929, 47.600956 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.8130171, 15.474525 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.2652084, 19.713544 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.0735573, 3.2055511 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(23.523258, 34.144142 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1052675, 4.2853865 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2261323, 54.650974 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(35.938682, 22.310661 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(43.632278, 37.541214 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(54.621666, 36.506382 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(21.211456, 18.530464 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.9579395, 0.2526102 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(49.482158, 22.623146 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.5465407, 22.858202 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(38.069335, 29.628252 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(49.169578, 42.292602 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(39.564998, 53.640960 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(5.0760717, 12.942659 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1172916, 0.3317201 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.6937513, 1.3820128 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(36.834030, 41.692997 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0868136, 2.0693941 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2487981, 28.159318 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(58.749683, 38.319442 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(54.396972, 40.095508 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.8807374, 3.0353784 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(10.761882, 1.7307828 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(10.045470, 3.1789214 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(44.929870, 0.9905142 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(59.720779, 30.977441 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(42.936012, 30.577076 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(46.079158, 28.352716 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(44.003917, 52.821025 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3446511, 18.510894 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(17.986646, 27.483783 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(59.111354, 22.785835 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(55.729534, 4.3312573 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(47.617763, 36.275802 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(44.650676, 38.945598 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.6326274, 25.486883 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2141819, 1.3356759 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.3732197, 3.1166579 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(24.886198, 2.5645513 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(54.308445, 25.413944 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(34.049346, 45.058528 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(53.701225, 27.893337 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.7097760, 1.9626674 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.7924485, 19.932649 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0195519, 38.367816 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(41.939609,  35.63274 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(22.228141, 39.524215 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.8375990, 23.524726 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(32.961212, 36.524200 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(9.4377613, 10.287744 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(30.673507, 23.828552 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(38.124122, 46.157466 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(53.382072, 33.468544 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(49.393001, 25.771896 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.3287160, 3.0902001 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(39.955749, 30.323211 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.7560508, 1.4533053 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(46.148002, 1.2616063 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(41.889728, 35.167430 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(54.938133, 36.032436 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3555525, 3.0687544 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(56.321323, 12.863502 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.4548411, 4.9593176 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(43.591377, 12.754740 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(37.563823, 9.5147047 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(33.541244, 1.4301798 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(59.563152, 17.541984 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(40.744159, 57.347049 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(36.621593, 2.5476176 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(46.718303,  36.70187 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6861963, 3.0976259 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(29.370142, 1.4475262 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(55.962074, 43.879081 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(45.212371, 22.505388 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(26.765331, 34.474536 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(59.777263, 32.886989 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(52.592639, 28.080614 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9678230, 21.956131 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(40.993480, 55.651126 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0214359, 55.590107 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(42.890899, 2.7383191 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(58.100662, 2.6796503 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6386389, 20.093421 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(38.900405, 2.5179841 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(29.298503, 35.910564 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(15.751144, 52.449920 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(28.119302, 0.3321138 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(35.297061, 20.921930 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0805194, 1.8598787 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.4601247, 0.5195755 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(22.912105, 1.0601193 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9651384, 1.4651534 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(30.186372, 3.0838427 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(8.0901327, 3.1087918 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(49.411384, 1.8853472 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1265718, 4.8931441 );
h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(35.915580, 30.104633 );  
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Draw("same");

  TH2D* h2_signal = new TH2D("h2_signal","h2_signal", m_bins, m_min, m_max, m_bins, m_min, m_max);
  //2015
  //h2_signal->Fill(0.404, 0.560);
  //2016
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
  line1->SetLineColor(1);
  line1->SetLineStyle(9);
  line1->SetLineWidth(2);
  line1->Draw();
  TLine *line2 = new TLine(diagonal_x1,m_min,m_max,diagonal_x2);
  line2->SetLineColor(1);
  line2->SetLineStyle(9);
  line2->SetLineWidth(2);
  line2->Draw();

  txtHeader_CMS->Draw();
  txtHeader_lumi->Draw();

  c_template2D_m1_vs_m2->SaveAs("figures_60/template2D_signal_and_background_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures_60/template2D_signal_and_background_m1_vs_m2.png");
}
