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
  txtHeader->SetHeader("CMS Prelim. 2015D  #sqrt{s} = 8 TeV   L_{int} = 35.9 fb^{-1}");


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
  txtHeader_lumi->SetHeader("35.9 fb^{-1} (13 TeV)");

  TFile* file = new TFile("ws_FINAL.root");
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

  double nEvents_Jpsi = 0.12; 

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
  cout<<"Probability one event should be in the J/Psi area (SIGNAL REGION):"<<endl;
  cout<<"2J/Psi Area is " << Area_2Jpsi/(Area_NO2Jpsi+Area_2Jpsi) << " of the rest of the signal region ("<<Area_2Jpsi<<" "<<Area_NO2Jpsi<<")."<<endl;
  cout<<"2J/Psi Area (weighted) is " << Area_2Jpsi_w/(Area_NO2Jpsi_w+Area_2Jpsi_w) << " of the rest of the signal region ("<<Area_2Jpsi_w<<" "<<Area_NO2Jpsi_w<<")."<<endl;
  cout<<"Probability one event should be in the J/Psi area (OFF-DIAGONAL REGION):"<<endl;
  cout<<"2J/Psi Area is " << offd_Area_2Jpsi/(offd_Area_NO2Jpsi+offd_Area_2Jpsi) << " of the rest of the signal region ("<<offd_Area_2Jpsi<<" "<<offd_Area_NO2Jpsi<<")."<<endl;
  cout<<"2J/Psi Area (weighted) is " << offd_Area_2Jpsi_w/(offd_Area_NO2Jpsi_w+offd_Area_2Jpsi_w) << " of the rest of the signal region ("<<offd_Area_2Jpsi_w<<" "<<offd_Area_NO2Jpsi_w<<")."<<endl;

  cout<<" -> Template2D_offDiagonal integral: "<<h2_Template2D_offDiagonal->Integral()<<endl;
  cout<<" -> Template2D_diagonal integral:    "<<h2_Template2D_diagonal->Integral()<<endl;

  //Signal: ISO +off Diag
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2",1000,1000);
  cout<<"#Event ISOLATED but offDiag: " << h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<endl;

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

  c_template2D_m1_vs_m2->SaveAs("figures/h2_background.pdf");// LP added
  c_template2D_m1_vs_m2->SaveAs("figures/h2_background.png");

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

  //  -------------------2015----------------------
  //------OffDiagonal SCAN------
  //************************************
  //*    Row   *     massC *     massF *
  //************************************
  //************************************************************************************************
  //*        0 * 0.3789996 * 1.8526362 *         0 * 1.2943311 *    254790 * 859868610 *       611 *
  //*        1 * 1.7474905 * 2.7977094 *         0 * 0.7504828 *    260576 * 388414933 *       188 *
  //*        2 * 1.1979647 * 0.8495020 * 1.5138732 *         0 *    258714 *  14142382 *        10 *
  //*        3 * 1.2658947 * 2.3115363 *         0 *         0 *    256843 * 1.603e+09 *      1166 *
  //************************************
  //------Signal SCAN------
  //************************************
  //*    Row   *     massC *     massF *
  //        0  0.4049646  0.5604345     256843  348750551        241 *
  //************************************************************************

  //  -------------------2016 BCDEF----------------------
  //*        0 * 1.5225365 * 0.8431518 * 0.8962873 *         0 *    274969 * 837317792 *       456 *
  //*        1 * 0.7942560 * 2.3587374 *         0 * 1.8426080 *    274968 * 1.550e+09 *       815 *
  //*        2 * 3.0640931 * 2.6193141 *         0 * 1.8181604 *    274441 * 341665349 *       207 *
  //*        3 * 0.6967173 * 2.0301356 * 1.4252665 * 1.2235757 *    274441 * 659054303 *       394 *
  //*        4 * 3.0594613 * 2.5633280 *         0 * 1.4294159 *    275836 * -1.98e+09 *      1287 *
  //*        5 * 0.3829744 * 0.6517009 *         0 *         0 *    276282 * 2.038e+09 *      1129 *
  //*        6 * 3.0445525 * 0.4301351 * 0.9365223 * 1.5587855 *    276502 * 594211575 *       390 *
  //*        7 * 0.7256926 * 3.0789272 *         0 * 0.8460569 *    276501 * 1.985e+09 *      1182 *
  //*        8 * 0.3180795 * 2.0960705 * 0.7000579 * 0.6592856 *    276525 * 1.635e+09 *       961 *
  //*        9 * 0.2748118 * 3.1098189 *         0 *         0 *    276525 * 757250117 *       490 *
  //*       10 * 1.0963534 * 1.9935944 * 1.8272670 *         0 *    276437 * -67935356 *      2042 *
  //*       11 * 0.2951213 * 7.9970874 *         0 *         0 *    276811 * 985305934 *       541 *
  //*       12 * 3.1160290 * 0.9953680 *         0 *         0 *    276363 * 1.976e+09 *      1103 *
  //*       13 * 2.0859162 * 3.2248451 *         0 *         0 *    276581 * 537441537 *       355 *
  //*       14 * 6.4438476 * 3.0851912 *         0 *         0 *    277070 * 895118089 *       478 *
  //*       15 * 1.9559303 * 0.2534542 * 0.7358005 * 0.7439422 *    276831 * 1.580e+09 *       882 *
  //*       16 * 2.6904547 * 3.0699591 * 1.0419131 * 1.0449538 *    276831 * 373972636 *       244 *
  //*       17 * 1.0854868 * 2.0424716 *         0 * 1.0831496 *    277194 * 2.054e+09 *      1171 *
  //*       18 * 3.4007272 * 2.9169464 * 1.7980103 *         0 *    277194 * 1.170e+09 *       708 *
  //*       19 * 1.6643008 * 1.1302868 *         0 *         0 *    278366 * 768530975 *       384 *
  //*       20 * 2.2167325 * 1.3391919 *         0 *         0 *    278406 * 1.386e+09 *       815 *
  //*       21 * 0.9294350 * 7.1592717 *         0 *         0 *    278509 * -2.00e+09 *      1351 *
  //*       22 * 2.3811967 * 3.1188454 * 1.0696374 *         0 *    278239 * 1.263e+09 *       719 *
  //*       23 * 0.5922431 * 1.2144465 * 0.5427458 * 0.6428673 *    278315 * 551943324 *       362 *
  //*       24 * 3.1027426 * 1.1035094 *         0 *         0 *    278310 *  32888078 *        23 *
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
  //2015
  //h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3789996, 1.8526362);
  //2016
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(6.3394618, 4.8828525);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1085057, 0.8863269);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0989241, 0.4605314);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3275566, 2.9306011);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0705516, 2.6184160);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6725933, 2.2082293);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.7733180, 2.9909498);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0867397, 2.5654206);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0961875, 2.6408567);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.3725998, 1.5676055);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.1172916, 0.3317201);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0868136, 2.0693941);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.6937513, 1.3820128);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1052675, 4.2853865);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.9579395, 0.2526102);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6968700, 3.0599653);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.8807374, 3.0353784);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.6633046, 1.1430702);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.3732197, 3.1166579);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.5893211, 1.2128995);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9293151, 7.1705985);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.0970430, 2.8785300);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.7097760, 1.9626674);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.2659432, 1.8314888);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2141819, 1.3356759);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.4045932, 1.0778753);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.4335794, 3.0383570);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0805194, 1.8598787);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.4548411, 4.9593176);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(8.0901327, 3.1087918);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.3287160, 3.0902001);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.1265718, 4.8931441);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1077470, 2.5274865);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6861963, 3.0976259);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9651384, 1.4651534);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.5701347, 1.2652914);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3555525, 3.0687544);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.4601247, 0.5195755);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.7560508, 1.4533053);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.0998075, 1.3160567);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2375977, 1.3019703);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.8153872, 2.4485802);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.8686915, 2.6802706);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.9575688, 2.4791264);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3246510, 8.8952770);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.2971377, 0.3618577);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.8503810, 3.1030097);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2727637, 3.0747036);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.8011589, 0.3165149); 
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp","h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp", m_bins, m_min, m_max, m_bins, m_min, m_max);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerSize(2.0);
  //2015
  //h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3789996, 1.8526362);
  //h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.2658947, 2.3115363);
  //h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.7474905, 2.7977094);
  //h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1979647, 0.8495020);
  //2016
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(6.3394618, 4.8828525);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1085057, 0.8863269);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0989241, 0.4605314);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3275566, 2.9306011);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0705516, 2.6184160);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6725933, 2.2082293);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.7733180, 2.9909498);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0867397, 2.5654206);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0961875, 2.6408567);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.3725998, 1.5676055);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1172916, 0.3317201);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0868136, 2.0693941);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.6937513, 1.3820128);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1052675, 4.2853865);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.9579395, 0.2526102);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6968700, 3.0599653);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.8807374, 3.0353784);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.6633046, 1.1430702);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.3732197, 3.1166579);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.5893211, 1.2128995);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9293151, 7.1705985);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.0970430, 2.8785300);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.7097760, 1.9626674);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.2659432, 1.8314888);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2141819, 1.3356759);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.4045932, 1.0778753);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.4335794, 3.0383570);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0805194, 1.8598787);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.4548411, 4.9593176);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(8.0901327, 3.1087918);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.3287160, 3.0902001);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1265718, 4.8931441);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1077470, 2.5274865);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6861963, 3.0976259);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9651384, 1.4651534);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.5701347, 1.2652914);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3555525, 3.0687544);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.4601247, 0.5195755);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.7560508, 1.4533053);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.0998075, 1.3160567);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2375977, 1.3019703);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.8153872, 2.4485802);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.8686915, 2.6802706);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.9575688, 2.4791264);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3246510, 8.8952770);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.2971377, 0.3618577);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.8503810, 3.1030097);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2727637, 3.0747036);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.8011589, 0.3165149); 
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Draw("same");

  TH2D* h2_signal = new TH2D("h2_signal","h2_signal", m_bins, m_min, m_max, m_bins, m_min, m_max);
  //2015
  //h2_signal->Fill(0.404, 0.560);
  //2016 BLINDED

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

  c_template2D_m1_vs_m2->SaveAs("figures/template2D_signal_and_background_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures/template2D_signal_and_background_m1_vs_m2.png");
}
