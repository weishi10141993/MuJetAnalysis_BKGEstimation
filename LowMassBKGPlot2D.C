//*****************************************************************************************************
//* cmsenv                                                                                            *
//* To request more time: sintr -t 480                                                                *
//*                                       Wei Shi @Nov 20, 2019, Rice U.                              *
//*****************************************************************************************************
#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"
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

void LowMassBKGPlot2D() {

  setTDRStyle();

  TCanvas * c_template2D_m1_vs_m2 = new TCanvas("c_template2D_m1_vs_m2", "c_template2D_m1_vs_m2",0,1320,1044,928);
  c_template2D_m1_vs_m2->SetCanvasSize(1040,900);
  c_template2D_m1_vs_m2->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2->SetRightMargin(0.17);
  c_template2D_m1_vs_m2->SetTopMargin(0.05);
  c_template2D_m1_vs_m2->cd();
  c_template2D_m1_vs_m2->SetLogz();

  TLegend *txtHeader = new TLegend(.003,.95,0.97,1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader("#bf{CMS} #it{Preliminary}    36.734 fb^{-1} (2017 13 TeV)");

  TFile* file = new TFile("ws_FINAL.root");
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  //double nEvents_Jpsi = 0.12;//from 2016 analysis

  //****************************************************************************
  //                         Draw 2D template m1 x m2
  //****************************************************************************
  //Create and fill ROOT 2D histogram (2*m_bins) with sampling of 2D pdf
  TH2D* h2_Template2D = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",2.0*m_bins,2.0*m_bins);
  cout << "Template2D integral: " << h2_Template2D->Integral() << std::endl;//normalized to 1
  TH2D* h2_Template2D_diagonal    = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);
  cout << "Template2D_diagonal integral: " << h2_Template2D_diagonal->Integral() << std::endl;
  TH2D* h2_Template2D_offDiagonal = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2",1000,1000);
  cout << "Template2D_offDiagonal integral: " << h2_Template2D_offDiagonal->Integral() << std::endl;

  for(int i=1;i<=1000;i++) {
    for(int j=1;j<=1000;j++) {
      double m_1 = h2_Template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
      double m_2 = h2_Template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
      //2017 mass consistency cut
      if ( fabs(m_1 - m_2) < 3*(0.003044 + 0.007025*(m_1+m_2)/2.0 + 0.000053*(m_1+m_2)*(m_1+m_2)/4.0) ) {
        h2_Template2D_offDiagonal->SetBinContent(i,j,0.);
      }
      else {
        h2_Template2D_diagonal->SetBinContent(i,j,0.);
      }
    }
  }

  cout<<" -> Template2D_offDiagonal integral: "<<h2_Template2D_offDiagonal->Integral()<<endl;
  cout<<" -> Template2D_diagonal integral:    "<<h2_Template2D_diagonal->Integral()<<endl;

  //2-dimu events at CR
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2",1000,1000);
  cout<<"#Event ISOLATED but offDiag: " << h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<endl;
  cout<<"Scaled as: "<<h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<" * (1 + "<<h2_Template2D_diagonal->Integral()<<" / "<<h2_Template2D_offDiagonal->Integral()<<") = "<<h2_dimudimu_control_Iso_offDiagonal_2D->Integral()<<" * "<<1+(h2_Template2D_diagonal->Integral()/h2_Template2D_offDiagonal->Integral())<<" = "<<h2_dimudimu_control_Iso_offDiagonal_2D->Integral()*(1+h2_Template2D_diagonal->Integral()/h2_Template2D_offDiagonal->Integral())<<endl;
  //Scale to: DimuDimu_iso_offDiag / Template2D_Area (normalize to the off-diag part of the data) * bb_ALL/bb_offDiag (scale factor to pass from a normalization off-diag. to a normalization to the whole area.)
  h2_Template2D->Scale(h2_dimudimu_control_Iso_offDiagonal_2D->Integral()/h2_Template2D->Integral()*(h2_Template2D_diagonal->Integral() + h2_Template2D_offDiagonal->Integral())/h2_Template2D_offDiagonal->Integral());
  cout<<"Scaled bb_2D template integral: " << h2_Template2D->Integral() <<" That means " << h2_Template2D->Integral()-h2_dimudimu_control_Iso_offDiagonal_2D->Integral() << " events in signal region "<< std::endl;

  TH2D * h2_background = new TH2D( *h2_Template2D );
  h2_background->GetXaxis()->SetTitle("m_{(#mu#mu)_{1}} [GeV]");
  h2_background->GetXaxis()->SetTitleOffset(0.93);
  h2_background->GetYaxis()->SetTitle("m_{(#mu#mu)_{2}} [GeV]");
  h2_background->GetYaxis()->SetTitleOffset(0.85);
  h2_background->GetZaxis()->SetTitle("Events / (0.025 GeV x 0.025 GeV)");
  h2_background->GetZaxis()->CenterTitle(true);
  h2_background->GetZaxis()->SetLabelFont(42);
  h2_background->GetZaxis()->SetLabelSize(0.04);
  h2_background->GetZaxis()->SetTitleSize(0.04);
  h2_background->GetZaxis()->SetTitleOffset(1.42);
  h2_background->GetZaxis()->SetTitleFont(42);

  //gStyle->SetPalette(52); //Grey Scale
  Double_t Red[2]    = { 1.00, 0.00};
  Double_t Green[2]  = { 1.00, 0.00};
  Double_t Blue[2]   = { 1.00, 0.00};
  Double_t Length[2] = { 0.00, 1.00 };
  Int_t nb=50;
  TColor::CreateGradientColorTable(2,Length,Red,Green,Blue,nb);
  h2_background->SetContour(nb);
  h2_background->Draw("Cont4 Colz");

  c_template2D_m1_vs_m2->SaveAs("figures/h2_background.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures/h2_background.png");

  // This required to draw scatter plot without LogZ
  TPad* pad = new TPad("pad", "pad",0,0,1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.121);
  pad->SetRightMargin(0.17);//48);
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
  /*
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3789996, 1.8526362);
  */
  //2016
  /*
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(6.3394618, 4.8828525);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1085057, 0.8863269);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0989241, 0.4605314);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3275566, 2.9306011);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0705516, 2.6184160);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6725933, 2.2082293);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.7733180, 2.9909498);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.3730475, 1.5679047);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3830212, 0.6524042);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0845816, 2.5646343);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.2796572, 1.3030284);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.2739968, 3.1119864);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0964052, 1.9947856);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.0735573, 3.2055511);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3029867, 8.0022907);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3957569, 2.0971429);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1111078, 0.9851592);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.5590274, 2.1472802);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6958830, 3.0606346);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(4.5877881, 0.9336586);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.9578526, 0.2527949);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.8998218, 3.9293952);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.8804173, 3.0341777);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.0866440, 2.0660750);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9268618, 7.1526098);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0994384, 1.1065233);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.3750891, 3.1192662);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.5902231, 1.2115613);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.6644971, 1.1404482);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2108738, 1.3358302);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(8.1372909, 3.1030330);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.4606950, 0.5183745);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.3274743, 3.0891082);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.1260975, 4.8807039);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.4039331, 1.0785657);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.5836787, 0.6993446);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.3601328, 3.0733373);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.6870803, 3.0973193);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.4318099, 3.0343155);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.1111371, 2.5319032);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.7588224, 1.4551333);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9652634, 1.4736373);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(7.6566572, 1.4617062);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.5125466, 3.4488928);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.0834312, 3.5854518);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2351813, 1.2998992);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0159301, 1.7288222);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.8026461, 0.3167102);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.2728815, 3.0707302);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(2.9555747, 2.4789545);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.3244402, 8.8903112);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(0.9500833, 2.3361797);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(1.8660776, 2.6829679);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.2984790, 0.3616573);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(5.3474988, 3.1133358);
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Fill(3.0882096, 1.4727710);
  */
  h2_dimudimu_control_Iso_offDiagonal_2D_points->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp = new TH2D("h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp","h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp", m_bins, m_min, m_max, m_bins, m_min, m_max);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->SetMarkerSize(2.0);
  //2015
  /*
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3789996, 1.8526362);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.2658947, 2.3115363);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.7474905, 2.7977094);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1979647, 0.8495020);
  */
  //2016
  /*
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(6.3394618, 4.8828525);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1085057, 0.8863269);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0989241, 0.4605314);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3275566, 2.9306011);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0705516, 2.6184160);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6725933, 2.2082293);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.7733180, 2.9909498);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.3730475, 1.5679047);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3830212, 0.6524042);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0845816, 2.5646343);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.2796572, 1.3030284);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.2739968, 3.1119864);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0964052, 1.9947856);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.0735573, 3.2055511);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3029867, 8.0022907);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3957569, 2.0971429);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1111078, 0.9851592);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.5590274, 2.1472802);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6958830, 3.0606346);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(4.5877881, 0.9336586);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.9578526, 0.2527949);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.8998218, 3.9293952);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.8804173, 3.0341777);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.0866440, 2.0660750);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9268618, 7.1526098);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0994384, 1.1065233);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.3750891, 3.1192662);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.5902231, 1.2115613);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.6644971, 1.1404482);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2108738, 1.3358302);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(8.1372909, 3.1030330);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.4606950, 0.5183745);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.3274743, 3.0891082);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.1260975, 4.8807039);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.4039331, 1.0785657);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.5836787, 0.6993446);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.3601328, 3.0733373);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.6870803, 3.0973193);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.4318099, 3.0343155);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.1111371, 2.5319032);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.7588224, 1.4551333);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9652634, 1.4736373);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(7.6566572, 1.4617062);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.5125466, 3.4488928);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.0834312, 3.5854518);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2351813, 1.2998992);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0159301, 1.7288222);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.8026461, 0.3167102);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.2728815, 3.0707302);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(2.9555747, 2.4789545);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.3244402, 8.8903112);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(0.9500833, 2.3361797);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(1.8660776, 2.6829679);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.2984790, 0.3616573);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(5.3474988, 3.1133358);
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Fill(3.0882096, 1.4727710);
  */
  h2_dimudimu_control_Iso_offDiagonal_2D_points_tmp->Draw("same");

  TH2D* h2_signal = new TH2D("h2_signal","h2_signal", m_bins, m_min, m_max, m_bins, m_min, m_max);
  //2015
  /*
  h2_signal->Fill(0.404, 0.560);
  */
  //2016
  /*
  h2_signal->Fill(0.8079733, 0.7267103);
  h2_signal->Fill(2.8599584, 3.0017674);
  h2_signal->Fill(0.4258973, 0.5848349);
  h2_signal->Fill(3.0722196, 3.2662851);
  h2_signal->Fill(3.0728187, 3.0538983);
  h2_signal->Fill(3.0950253, 3.3617882);
  h2_signal->Fill(3.1521356, 2.8546791);
  h2_signal->Fill(2.8254406, 2.6496100);
  h2_signal->Fill(1.2541753, 1.1524148);
  h2_signal->Fill(2.3863873, 2.3582603);
  h2_signal->Fill(3.0641751, 3.0972354);
  h2_signal->Fill(1.9403913, 1.8196427);
  h2_signal->Fill(1.3540757, 1.4834892);
  */
  h2_signal->GetXaxis()->SetTitle("m_{(#mu#mu)_{1}} [GeV]");
  h2_signal->GetXaxis()->CenterTitle(true);
  h2_signal->GetXaxis()->SetTitleOffset(0.93);
  h2_signal->GetYaxis()->SetTitle("m_{(#mu#mu)_{2}} [GeV]");
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

  //Pre-calculated m1 and m2 values for drawing the corridor curves:
  //|m1-m2| - 3*(0.003044 + 0.007025*(m1+m2)/2.0 + 0.000053*(m1+m2)*(m1+m2)/4.0)
  //double m1Input[18]={0.25,    0.4,     0.7,     1.0,     2.0,     5.0,     8.0,     10.0,     15.0,     20.0,     25.0,     30.0,     35.0,     40.0,     45.0,     50.0,     55.0,     60.0};
  //double m2Small[18]={0.23574, 0.38259, 0.67629, 0.96995, 1.94863, 4.88284, 7.81428, 9.76704,  14.64356, 19.51244, 24.37369, 29.22732, 34.07335, 38.91180, 43.74269, 48.56604, 53.38186, 58.19017};
  //double m2Large[18]={0.26456, 0.41777, 0.72422, 1.03069, 2.05248, 5.11984, 8.19015, 10.23867, 15.36576, 20.50111, 25.64475, 30.79670, 35.95697, 41.12560, 46.30259, 51.48797, 56.68177, 61.88399};
  double m1Input[8]={0.25,    0.4,     0.7,     1.0,     2.0,     5.0,     8.0,     10.0};
  double m2Small[8]={0.23574, 0.38259, 0.67629, 0.96995, 1.94863, 4.88284, 7.81428, 9.76704};
  double m2Large[8]={0.26456, 0.41777, 0.72422, 1.03069, 2.05248, 5.11984, 8.19015, 10.23867};
  TGraph* corridorDn = new TGraph(8, m1Input, m2Large);
  TGraph* corridorUp = new TGraph(8, m1Input, m2Small);
  corridorDn->SetLineColor(1); corridorDn->SetLineStyle(9); corridorDn->SetLineWidth(2); corridorDn->Draw("C");
  corridorUp->SetLineColor(1); corridorUp->SetLineStyle(9); corridorUp->SetLineWidth(2); corridorUp->Draw("C");
  txtHeader->Draw();

  c_template2D_m1_vs_m2->SaveAs("figures/template2D_signal_and_background_m1_vs_m2.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures/template2D_signal_and_background_m1_vs_m2.png");
}
