//=========================================================================
//= cmsenv                                                                =
//= Run it as: root -l -b -q LowMassBKGPlot2D18.C                         =
//=          Wei Shi @Nov 20, 2019, Rice U.                               =
//=========================================================================
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
#include "Constants.h"
#include "Config.h"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

void LowMassBKGPlot2D18() {

  //Configure inputs for year
  BKG_cfg::ConfigureInput(year);
  setTDRStyle();

  //for template
  TLegend *txtHeader = new TLegend(.00001, 0.95, 0.94, 1.);
  txtHeader->SetFillColor(kWhite);
  txtHeader->SetFillStyle(0);
  txtHeader->SetBorderSize(0);
  txtHeader->SetTextFont(42);
  txtHeader->SetTextSize(0.045);
  txtHeader->SetTextAlign(22);
  txtHeader->SetHeader(headertpl2);

  //for validation plot
  TLegend *txtHeadervld = new TLegend(0.14, 0.93, 0.98, 0.98);
  txtHeadervld->SetFillColor(kWhite);
  txtHeadervld->SetFillStyle(0);
  txtHeadervld->SetBorderSize(0);
  txtHeadervld->SetTextFont(42);
  txtHeadervld->SetTextSize(0.045);
  txtHeadervld->SetTextAlign(22);
  txtHeadervld->SetHeader(headervld);

  //gStyle->SetPalette(52); //Grey Scale
  Double_t Red[2]    = {1.00, 0.00};
  Double_t Green[2]  = {1.00, 0.00};
  Double_t Blue[2]   = {1.00, 0.00};
  Double_t Length[2] = {0.00, 1.00};
  Int_t nb=50;

  //====================================================================================
  //          Pre-calculated m1 and m2 values for drawing the corridor curves
  //====================================================================================
  //Straigt line interpolation, just use the signal mass points
  double m2Small[11];
  double m2Large[11];
  for (int i = 0; i < 11; i++) {
    m2Small[i] = mean[i] - window[i];
    m2Large[i] = mean[i] + window[i];
  }
  TGraph* corridorDn = new TGraph(11, mean, m2Small);
  TGraph* corridorUp = new TGraph(11, mean, m2Large);
  corridorDn->SetLineColor(1); corridorDn->SetLineStyle(9); corridorDn->SetLineWidth(2);
  corridorUp->SetLineColor(1); corridorUp->SetLineStyle(9); corridorUp->SetLineWidth(2);

  /*
  //Const*CB fit sigma: fabs(m_1 - m_2) < Constant*(0.003681 + 0.007583*(m_1 + m_2)/2.0)
  const double Constant = 2.0;
  const double kA = 0.003681*Constant;
  const double kB = 0.007583*Constant;

  //Diagonal lines below J/psi
  double diagonal_x1 = ( (1.0 + kB/2.0)*m_min + kA )/( 1.0 - kB/2.0 );
  double diagonal_x2 = ( (1.0 - kB/2.0)*m_Jpsi_dn - kA )/( 1.0 + kB/2.0 );
  TLine *line1 = new TLine(m_min, diagonal_x1, diagonal_x2, m_Jpsi_dn); line1->SetLineColor(1); line1->SetLineStyle(9); line1->SetLineWidth(2);
  TLine *line2 = new TLine(diagonal_x1, m_min, m_Jpsi_dn, diagonal_x2); line2->SetLineColor(1); line2->SetLineStyle(9); line2->SetLineWidth(2);

  //Diagonal lines above J/psi and below Upsilon
  double diagonal_x3 = ( (1.0 + kB/2.0)*m_Jpsi_up + kA )/( 1.0 - kB/2.0 );
  double diagonal_x4 = ( (1.0 - kB/2.0)*m_max - kA )/( 1.0 + kB/2.0 );
  TLine *line3 = new TLine(m_Jpsi_up, diagonal_x3, diagonal_x4, m_max); line3->SetLineColor(1); line3->SetLineStyle(9); line3->SetLineWidth(2);
  TLine *line4 = new TLine(diagonal_x3, m_Jpsi_up, m_max, diagonal_x4); line4->SetLineColor(1); line4->SetLineStyle(9); line4->SetLineWidth(2);

  //Diagonal lines above Upsilon
  double diagonal_x5 = ( (1.0 + kB/2.0)*m_Upsilon_up + kA )/( 1.0 - kB/2.0 );
  double diagonal_x6 = ( (1.0 - kB/2.0)*m_highmax - kA )/( 1.0 + kB/2.0 );
  TLine *line5 = new TLine(m_Upsilon_up, diagonal_x5, diagonal_x6, m_highmax); line5->SetLineColor(1); line5->SetLineStyle(9); line5->SetLineWidth(2);
  TLine *line6 = new TLine(diagonal_x5, m_Upsilon_up, m_highmax, diagonal_x6); line6->SetLineColor(1); line6->SetLineStyle(9); line6->SetLineWidth(2);
  //Note: Above only works for poly-1 fit
  */

  //-----------------
  //2018 below J/psi
  //-----------------
  //Poly3: fabs(m_1 - m_2) < 5*(0.00849813 + 0.00475107*(m_1 + m_2)/2.0 - 0.00665393*pow((m_1 + m_2)/2.0, 2) + 0.00337777*pow((m_1 + m_2)/2.0, 3) )
  /*
  double m1BelowJpsiInput[6] = {0.25,   0.40,   0.70,   1.00,   2.00,   2.72};
  double m2BelowJpsiSmall[6] = {0.2036, 0.3524, 0.6514, 0.9503, 1.9120, 2.5381};
  double m2BelowJpsiLarge[6] = {0.2968, 0.4479, 0.7487, 1.0500, 2.0967, 2.9469};
  TGraph* corridorDnBelowJpsi = new TGraph(6, m1BelowJpsiInput, m2BelowJpsiSmall);
  TGraph* corridorUpBelowJpsi = new TGraph(6, m1BelowJpsiInput, m2BelowJpsiLarge);
  corridorDnBelowJpsi->SetLineColor(1); corridorDnBelowJpsi->SetLineStyle(9); corridorDnBelowJpsi->SetLineWidth(2);
  corridorUpBelowJpsi->SetLineColor(1); corridorUpBelowJpsi->SetLineStyle(9); corridorUpBelowJpsi->SetLineWidth(2);
  */

  //-----------------
  //2018 above J/psi
  //-----------------
  //Poly4: fabs(m_1 - m_2) < 5*(0.0472738 - 0.00591865*(m_1 + m_2)/2.0 + 0.00113991*pow((m_1 + m_2)/2.0, 2) - 2.62048e-05*pow((m_1 + m_2)/2.0, 3) + 1.92254e-07*pow((m_1 + m_2)/2.0, 4) )
  /*
  double m1AboveJpsiInput[7] = {3.24,   4.00,   5.00,   6.00,   7.00,   8.00,   9.00};
  double m2AboveJpsiSmall[7] = {3.0443, 3.8000, 4.7868, 5.7660, 6.7383, 7.7044, 8.6648};
  double m2AboveJpsiLarge[7] = {3.4363, 4.2021, 5.2171, 6.2401, 7.2703, 8.3072, 9.3503};
  TGraph* corridorDnAboveJpsi = new TGraph(7, m1AboveJpsiInput, m2AboveJpsiSmall);
  TGraph* corridorUpAboveJpsi = new TGraph(7, m1AboveJpsiInput, m2AboveJpsiLarge);
  corridorDnAboveJpsi->SetLineColor(1); corridorDnAboveJpsi->SetLineStyle(9); corridorDnAboveJpsi->SetLineWidth(2);
  corridorUpAboveJpsi->SetLineColor(1); corridorUpAboveJpsi->SetLineStyle(9); corridorUpAboveJpsi->SetLineWidth(2);
  */

  //---------------------------------------
  //2018 above Upsilon (use with CAUTION!)
  //---------------------------------------

  //==================================================================================
  //                                  Use pdf from workspace:
  //           template2D_below_Jpsi, template2D_above_Jpsi, template2D_above_Upsilon
  //==================================================================================

  TFile* file = new TFile(inputFile2);
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  //=================
  //Below Jpsi ONLY
  //=================
  cout<<"****** Part A: Below J/Psi ONLY ******" <<endl;
  TH2D* h2D_template2D_below_Jpsi = (TH2D*)w->pdf("template2D_below_Jpsi")->createHistogram("m1_below_Jpsi,m2_below_Jpsi", m_bins_below_Jpsi, m_bins_below_Jpsi);//normalized
  TH2D *h2D_template2D_below_Jpsi_diagonal = (TH2D*)h2D_template2D_below_Jpsi->Clone();
  TH2D *h2D_template2D_below_Jpsi_offDiagonal = (TH2D*)h2D_template2D_below_Jpsi->Clone();

  for(int i=1;i<=m_bins_below_Jpsi;i++) {
    for(int j=1;j<=m_bins_below_Jpsi;j++) {
      double m_1 = h2D_template2D_below_Jpsi_offDiagonal->GetXaxis()->GetBinCenter(i);
      double m_2 = h2D_template2D_below_Jpsi_offDiagonal->GetYaxis()->GetBinCenter(j);
      //===========================
      //2018 mass window below Jpsi
      //===========================
      //if ( fabs(m_1 - m_2) < 5*(0.00849813 + 0.00475107*(m_1 + m_2)/2.0 - 0.00665393*pow((m_1 + m_2)/2.0, 2) + 0.00337777*pow((m_1 + m_2)/2.0, 3) ) ) {
      //if ( fabs(m_1 - m_2) < (kA + kB*(m_1 + m_2)/2.) ) {
      if ( fabs(m_1 - m_2) < BKG_cfg::My_MassWindow(m_1, m_2) ) {
        h2D_template2D_below_Jpsi_offDiagonal->SetBinContent(i, j, 0.);
      }
      else {
        h2D_template2D_below_Jpsi_diagonal->SetBinContent(i, j, 0.);
      }
    }
  }

  //Fractions of area for diagonal and offdiagonal in 2D template (below Jpsi)
  double Template2D_below_Jpsi_diagonal_integral  = h2D_template2D_below_Jpsi_diagonal->Integral();
  double Template2D_below_Jpsi_offDiagonal_integral  = h2D_template2D_below_Jpsi_offDiagonal->Integral();
  cout<<" -> Template2D (Below J/Psi) diagonal integral:    "<< Template2D_below_Jpsi_diagonal_integral <<endl;
  cout<<" -> Template2D (Below J/Psi) offDiagonal integral: "<< Template2D_below_Jpsi_offDiagonal_integral <<endl;

  //count 2-dimu data events at CR
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi")->createHistogram("m1_below_Jpsi,m2_below_Jpsi", m_bins_below_Jpsi, m_bins_below_Jpsi);
  double Signal_CR_Data_below_Jpsi_integral  = h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->Integral();
  cout<<"2 dimuon events in DATA at CR (Below J/Psi ONLY): " << Signal_CR_Data_below_Jpsi_integral <<endl;
  cout<<"Expected 2 dimuon events in DATA at SR (Below J/Psi ONLY): " << Signal_CR_Data_below_Jpsi_integral*Template2D_below_Jpsi_diagonal_integral/Template2D_below_Jpsi_offDiagonal_integral << std::endl;

  //Scale to actual DATA
  h2D_template2D_below_Jpsi->Scale(Signal_CR_Data_below_Jpsi_integral*(1. + Template2D_below_Jpsi_diagonal_integral/Template2D_below_Jpsi_offDiagonal_integral));
  cout<<" h2D_template2D_below_Jpsi integral (after scale to actual data): "<< h2D_template2D_below_Jpsi->Integral() <<endl;
  TH2D * h2D_background_below_Jpsi = new TH2D( *h2D_template2D_below_Jpsi );
  h2D_background_below_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h2D_background_below_Jpsi->GetXaxis()->SetTitleOffset(0.93);
  h2D_background_below_Jpsi->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h2D_background_below_Jpsi->GetYaxis()->SetTitleOffset(0.85);
  h2D_background_below_Jpsi->GetZaxis()->SetTitle("Events/(0.04 GeV x 0.04 GeV)");
  h2D_background_below_Jpsi->GetZaxis()->CenterTitle(true);
  h2D_background_below_Jpsi->GetZaxis()->SetLabelFont(42);
  h2D_background_below_Jpsi->GetZaxis()->SetLabelSize(0.04);
  h2D_background_below_Jpsi->GetZaxis()->SetTitleSize(0.04);
  h2D_background_below_Jpsi->GetZaxis()->SetTitleOffset(1.42);
  h2D_background_below_Jpsi->GetZaxis()->SetTitleFont(42);

  TCanvas * c_template2D_m1_vs_m2_below_Jpsi = new TCanvas("c_template2D_m1_vs_m2_below_Jpsi", "c_template2D_m1_vs_m2_below_Jpsi", 0, 1320, 1044, 928);
  c_template2D_m1_vs_m2_below_Jpsi->SetCanvasSize(1040, 900);
  c_template2D_m1_vs_m2_below_Jpsi->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2_below_Jpsi->SetRightMargin(0.17);
  c_template2D_m1_vs_m2_below_Jpsi->SetTopMargin(0.05);
  c_template2D_m1_vs_m2_below_Jpsi->cd();
  c_template2D_m1_vs_m2_below_Jpsi->SetLogz();

  TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, nb);
  h2D_background_below_Jpsi->SetContour(nb);
  h2D_background_below_Jpsi->Draw("Cont4 Colz");
  c_template2D_m1_vs_m2_below_Jpsi->SaveAs("figures/Expected_2D_background_below_Jpsi.pdf");
  c_template2D_m1_vs_m2_below_Jpsi->SaveAs("figures/Expected_2D_background_below_Jpsi.png");
  c_template2D_m1_vs_m2_below_Jpsi->SaveAs("figures/Expected_2D_background_below_Jpsi.root");

  //**************************************************************************************
  //        Draw scatter plot at CR from data (SR blinded) Below Jpsi version
  //**************************************************************************************
  //Create pad to draw scatter plot without Logz
  TPad* pad_below_Jpsi = new TPad("pad_below_Jpsi", "pad_below_Jpsi", 0, 0, 1, 1);
  pad_below_Jpsi->Draw();
  pad_below_Jpsi->cd();
  pad_below_Jpsi->SetLeftMargin(0.121);
  pad_below_Jpsi->SetRightMargin(0.17);//48);
  pad_below_Jpsi->SetTopMargin(0.05);
  pad_below_Jpsi->SetFillColor(0);
  pad_below_Jpsi->SetFillStyle(4000);
  pad_below_Jpsi->SetBorderMode(0);
  pad_below_Jpsi->SetBorderSize(2);
  pad_below_Jpsi->SetTickx(1);
  pad_below_Jpsi->SetTicky(1);
  pad_below_Jpsi->SetFrameFillStyle(0);
  pad_below_Jpsi->SetFrameBorderMode(0);
  pad_below_Jpsi->SetFrameFillStyle(0);
  pad_below_Jpsi->SetFrameBorderMode(0);

  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->SetXTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->SetYTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi_tmp = new TH2D( *h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi);
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi_tmp->SetMarkerSize(1.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi_tmp->Draw("same");

  //**************************************************************************************************
  //    !!!BEGIN: Placeholder for unblinding signal (below Jpsi), after green light from pre-approval
  //**************************************************************************************************
  /*
  TH2D* h2D_dimudimu_signal_2D_below_Jpsi = (TH2D*)w->data("ds_dimudimu_signal_2D_below_Jpsi")->createHistogram("m1_below_Jpsi,m2_below_Jpsi", m_bins_below_Jpsi, m_bins_below_Jpsi);
  h2D_dimudimu_signal_2D_below_Jpsi->SetMarkerColor(kBlack);
  h2D_dimudimu_signal_2D_below_Jpsi->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_below_Jpsi->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2D_dimudimu_signal_2D_below_Jpsi->SetXTitle("");
  h2D_dimudimu_signal_2D_below_Jpsi->SetYTitle("");
  h2D_dimudimu_signal_2D_below_Jpsi->Draw("same");

  TH2D * h2D_dimudimu_signal_2D_below_Jpsi_tmp = new TH2D( *h2D_dimudimu_signal_2D_below_Jpsi);
  h2D_dimudimu_signal_2D_below_Jpsi_tmp->SetMarkerColor(kYellow);
  h2D_dimudimu_signal_2D_below_Jpsi_tmp->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_below_Jpsi_tmp->SetMarkerSize(1.0);
  h2D_dimudimu_signal_2D_below_Jpsi_tmp->Draw("same");
  */
  //**************************************************************************************************
  //    !!!END: Placeholder for unblinding signal (below Jpsi), after green light from pre-approval
  //**************************************************************************************************

  //corridorDnBelowJpsi->Draw("L"); corridorUpBelowJpsi->Draw("L"); txtHeader->Draw();
  //line1->Draw(); line2->Draw(); txtHeader->Draw();
  corridorDn->Draw("L"); corridorUp->Draw("L"); txtHeader->Draw();
  c_template2D_m1_vs_m2_below_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_below_Jpsi.pdf");
  c_template2D_m1_vs_m2_below_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_below_Jpsi.png");
  c_template2D_m1_vs_m2_below_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_below_Jpsi.root");

  //Validate for below J/psi ONLY
  //---------------------------------------------------------------------
  cout<<"                                                       " <<endl;
  cout<<"------ Start: validate for m1 (Below J/Psi ONLY) ------" <<endl;
  //---------------------------------------------------------------------
  TH1D *h1_control_Iso_offDiagonal_massC_data_below_Jpsi = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi")->createHistogram("m1_below_Jpsi", m_bins_below_Jpsi);
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->SetStats(0);
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetYaxis()->SetTitleOffset(0.9);
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massC_template_below_Jpsi = new TH1D( *h2D_template2D_below_Jpsi_offDiagonal->ProjectionX() );
  h1_control_Iso_offDiagonal_massC_template_below_Jpsi->Scale( h1_control_Iso_offDiagonal_massC_data_below_Jpsi->Integral()*1./h1_control_Iso_offDiagonal_massC_template_below_Jpsi->Integral() );
  h1_control_Iso_offDiagonal_massC_template_below_Jpsi->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massC_template_below_Jpsi->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massC_template_below_Jpsi->SetMarkerColor(kRed);
  for (Int_t i = 0; i < m_bins; ++i) h1_control_Iso_offDiagonal_massC_template_below_Jpsi->SetBinError(i, 0);

  TCanvas * c_control_Iso_offDiagonal_massC_below_Jpsi = new TCanvas("c_control_Iso_offDiagonal_massC_below_Jpsi", "c_control_Iso_offDiagonal_massC_below_Jpsi");
  c_control_Iso_offDiagonal_massC_below_Jpsi->cd();
  h1_control_Iso_offDiagonal_massC_data_below_Jpsi->Draw("e1");
  h1_control_Iso_offDiagonal_massC_template_below_Jpsi->Draw("HIST same");
  txtHeadervld->Draw();
  c_control_Iso_offDiagonal_massC_below_Jpsi->SaveAs("figures/Validation_m1_CR_below_Jpsi.pdf");
  c_control_Iso_offDiagonal_massC_below_Jpsi->SaveAs("figures/Validation_m1_CR_below_Jpsi.png");
  c_control_Iso_offDiagonal_massC_below_Jpsi->SaveAs("figures/Validation_m1_CR_below_Jpsi.root");
  //--------------------------
  //     Compatibility test
  //--------------------------
  //K-S test
  double KSprob_5 = h1_control_Iso_offDiagonal_massC_data_below_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_below_Jpsi);
  double KSdist_5 = h1_control_Iso_offDiagonal_massC_data_below_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_below_Jpsi, "M");
  cout<<"K-S test prob.: "<< KSprob_5 << "; dist.: "<< KSdist_5 <<endl;
  //Chisquare test
  auto func5 = [&](double *x, double*) { int ibin = h1_control_Iso_offDiagonal_massC_template_below_Jpsi->FindBin(x[0]); return h1_control_Iso_offDiagonal_massC_template_below_Jpsi->GetBinContent(ibin);};
  auto f5 = new TF1("f5", func5, h1_control_Iso_offDiagonal_massC_template_below_Jpsi->GetXaxis()->GetXmin(), h1_control_Iso_offDiagonal_massC_template_below_Jpsi->GetXaxis()->GetXmax(), 0);
  double BCchi2_5 = h1_control_Iso_offDiagonal_massC_data_below_Jpsi->Chisquare(f5, "L");
  double BCprob_5 = TMath::Prob(BCchi2_5, h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetNbinsX());
  cout<<"Chisquare test prob.: "<< BCprob_5 << "; chi2: "<< BCchi2_5 << "; ndof: " << h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetNbinsX() <<endl;
  //Toy experiments for calibration
  auto h_t5 = (TH1*) h1_control_Iso_offDiagonal_massC_data_below_Jpsi->Clone();//placeholder to fill pseudo data from template, same bins as data
  auto hKS_t5 = new TH1D("hKS_t5", "K-S distance", 100, 0, 1);//K-S distance distibution from pseudo exp.
  auto hBC_t5 = new TH1D("hBC_t5", "Baker-Cousins chi2", 100, 0, 300);//chi2 distibution from pseudo exp.
  int nKS_t5 = 0;
  int nBC_t5 = 0;
  int nentries_t5 = h1_control_Iso_offDiagonal_massC_data_below_Jpsi->Integral(1, h1_control_Iso_offDiagonal_massC_data_below_Jpsi->GetNbinsX());
  for (int i = 0; i < ntoys; ++i) {
    h_t5->Reset();
    h_t5->FillRandom("f5", nentries_t5);
    double KSdist_t5 = h_t5->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_below_Jpsi, "M");
    double BCchi2_t5 = h_t5->Chisquare(f5, "L");
    hKS_t5->Fill(KSdist_t5);
    hBC_t5->Fill(BCchi2_t5);
    if (KSdist_t5 > KSdist_5) nKS_t5++;
    if (BCchi2_t5 > BCchi2_5) nBC_t5++;
  }
  std::cout << "Corrected prob. for K-S  test: " << nKS_t5/double(ntoys) << std::endl;
  std::cout << "Corrected prob. for chi2 test: " << nBC_t5/double(ntoys) << std::endl;
  auto c5 = new TCanvas(); c5->Divide(1, 2);
  c5->cd(1); hKS_t5->Draw();
  c5->cd(2); hBC_t5->Draw();
  c5->SaveAs("figures/toys_m1_CR_below_Jpsi.root");
  cout<<"------ End: validate for m1 (Below J/Psi ONLY) ------" <<endl;

  //---------------------------------------------------------------------
  cout<<"                                                       " <<endl;
  cout<<"------ Start: validate for m2 (Below J/Psi ONLY) ------" <<endl;
  //---------------------------------------------------------------------
  TH1D *h1_control_Iso_offDiagonal_massF_data_below_Jpsi = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_below_Jpsi")->createHistogram("m2_below_Jpsi", m_bins_below_Jpsi);
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->SetStats(0);
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetYaxis()->SetTitleOffset(0.9);
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massF_template_below_Jpsi = new TH1D( *h2D_template2D_below_Jpsi_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_massF_template_below_Jpsi->Scale( h1_control_Iso_offDiagonal_massF_data_below_Jpsi->Integral()*1./h1_control_Iso_offDiagonal_massF_template_below_Jpsi->Integral() );
  h1_control_Iso_offDiagonal_massF_template_below_Jpsi->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massF_template_below_Jpsi->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massF_template_below_Jpsi->SetMarkerColor(kRed);
  for (Int_t i = 0; i < m_bins; ++i) h1_control_Iso_offDiagonal_massF_template_below_Jpsi->SetBinError(i, 0);

  TCanvas * c_control_Iso_offDiagonal_massF_below_Jpsi = new TCanvas("c_control_Iso_offDiagonal_massF_below_Jpsi", "c_control_Iso_offDiagonal_massF_below_Jpsi");
  c_control_Iso_offDiagonal_massF_below_Jpsi->cd();
  h1_control_Iso_offDiagonal_massF_data_below_Jpsi->Draw("e1");
  h1_control_Iso_offDiagonal_massF_template_below_Jpsi->Draw("HIST same");
  txtHeadervld->Draw();
  c_control_Iso_offDiagonal_massF_below_Jpsi->SaveAs("figures/Validation_m2_CR_below_Jpsi.pdf");
  c_control_Iso_offDiagonal_massF_below_Jpsi->SaveAs("figures/Validation_m2_CR_below_Jpsi.png");
  c_control_Iso_offDiagonal_massF_below_Jpsi->SaveAs("figures/Validation_m2_CR_below_Jpsi.root");
  //--------------------------
  //     Compatibility test
  //--------------------------
  //K-S test
  double KSprob_6 = h1_control_Iso_offDiagonal_massF_data_below_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_below_Jpsi);
  double KSdist_6 = h1_control_Iso_offDiagonal_massF_data_below_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_below_Jpsi, "M");
  cout<<"K-S test prob.: "<< KSprob_6 << "; dist.: "<< KSdist_6 <<endl;
  //Chisquare test
  auto func6 = [&](double *x, double*) { int ibin = h1_control_Iso_offDiagonal_massF_template_below_Jpsi->FindBin(x[0]); return h1_control_Iso_offDiagonal_massF_template_below_Jpsi->GetBinContent(ibin);};
  auto f6 = new TF1("f6", func6, h1_control_Iso_offDiagonal_massF_template_below_Jpsi->GetXaxis()->GetXmin(), h1_control_Iso_offDiagonal_massF_template_below_Jpsi->GetXaxis()->GetXmax(), 0);
  double BCchi2_6 = h1_control_Iso_offDiagonal_massF_data_below_Jpsi->Chisquare(f6, "L");
  double BCprob_6 = TMath::Prob(BCchi2_6, h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetNbinsX());
  cout<<"Chisquare test prob.: "<< BCprob_6 << "; chi2: "<< BCchi2_6 << "; ndof: " << h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetNbinsX() <<endl;
  //Toy experiments for calibration
  auto h_t6 = (TH1*) h1_control_Iso_offDiagonal_massF_data_below_Jpsi->Clone();//placeholder to fill pseudo data from template, same bins as data
  auto hKS_t6 = new TH1D("hKS_t6", "K-S distance", 100, 0, 1);//K-S distance distibution from pseudo exp.
  auto hBC_t6 = new TH1D("hBC_t6", "Baker-Cousins chi2", 100, 0, 300);//chi2 distibution from pseudo exp.
  int nKS_t6 = 0;
  int nBC_t6 = 0;
  int nentries_t6 = h1_control_Iso_offDiagonal_massF_data_below_Jpsi->Integral(1, h1_control_Iso_offDiagonal_massF_data_below_Jpsi->GetNbinsX());
  for (int i = 0; i < ntoys; ++i) {
    h_t6->Reset();
    h_t6->FillRandom("f6", nentries_t6);
    double KSdist_t6 = h_t6->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_below_Jpsi, "M");
    double BCchi2_t6 = h_t6->Chisquare(f6, "L");
    hKS_t6->Fill(KSdist_t6);
    hBC_t6->Fill(BCchi2_t6);
    if (KSdist_t6 > KSdist_6) nKS_t6++;
    if (BCchi2_t6 > BCchi2_6) nBC_t6++;
  }
  std::cout << "Corrected prob. for K-S  test: " << nKS_t6/double(ntoys) << std::endl;
  std::cout << "Corrected prob. for chi2 test: " << nBC_t6/double(ntoys) << std::endl;
  auto c6 = new TCanvas(); c6->Divide(1, 2);
  c6->cd(1); hKS_t6->Draw();
  c6->cd(2); hBC_t6->Draw();
  c6->SaveAs("figures/toys_m2_CR_below_Jpsi.root");
  cout<<"------ End: validate for m2 (Below J/Psi ONLY) ------" <<endl;
  cout<<"                                                     " <<endl;

  //=================
  //Above Jpsi ONLY
  //=================
  cout<<"****** Part B: Above J/Psi ONLY ******" <<endl;
  TH2D* h2D_template2D_above_Jpsi = (TH2D*)w->pdf("template2D_above_Jpsi")->createHistogram("m1_above_Jpsi,m2_above_Jpsi", m_bins_above_Jpsi, m_bins_above_Jpsi);//normalized
  TH2D *h2D_template2D_above_Jpsi_diagonal = (TH2D*)h2D_template2D_above_Jpsi->Clone();
  TH2D *h2D_template2D_above_Jpsi_offDiagonal = (TH2D*)h2D_template2D_above_Jpsi->Clone();

  for(int i=1;i<=m_bins_above_Jpsi;i++) {
    for(int j=1;j<=m_bins_above_Jpsi;j++) {
      double m_1 = h2D_template2D_above_Jpsi_offDiagonal->GetXaxis()->GetBinCenter(i);
      double m_2 = h2D_template2D_above_Jpsi_offDiagonal->GetYaxis()->GetBinCenter(j);
      //===========================
      //2018 mass window above Jpsi
      //===========================
      //if ( fabs(m_1 - m_2) < 5*(0.0472738 - 0.00591865*(m_1 + m_2)/2.0 + 0.00113991*pow((m_1 + m_2)/2.0, 2) - 2.62048e-05*pow((m_1 + m_2)/2.0, 3) + 1.92254e-07*pow((m_1 + m_2)/2.0, 4) ) ) {
      //if ( fabs(m_1 - m_2) < (kA + kB*(m_1 + m_2)/2.) ) {
      if ( fabs(m_1 - m_2) < BKG_cfg::My_MassWindow(m_1, m_2) ) {
        h2D_template2D_above_Jpsi_offDiagonal->SetBinContent(i, j, 0.);
      }
      else {
        h2D_template2D_above_Jpsi_diagonal->SetBinContent(i, j, 0.);
      }
    }
  }

  //Fractions of area for diagonal and offdiagonal in 2D template (above Jpsi)
  double Template2D_above_Jpsi_diagonal_integral  = h2D_template2D_above_Jpsi_diagonal->Integral();
  double Template2D_above_Jpsi_offDiagonal_integral  = h2D_template2D_above_Jpsi_offDiagonal->Integral();
  cout<<" -> Template2D (Above J/Psi) diagonal integral:    "<< Template2D_above_Jpsi_diagonal_integral <<endl;
  cout<<" -> Template2D (Above J/Psi) offDiagonal integral: "<< Template2D_above_Jpsi_offDiagonal_integral <<endl;

  //count 2-dimu data events at CR
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi")->createHistogram("m1_above_Jpsi,m2_above_Jpsi", m_bins_above_Jpsi, m_bins_above_Jpsi);
  double Signal_CR_Data_above_Jpsi_integral  = h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->Integral();
  cout<<"2 dimuon events in DATA at CR (Above J/Psi ONLY): " << Signal_CR_Data_above_Jpsi_integral <<endl;
  cout<<"Expected 2 dimuon events in DATA at SR (Above J/Psi ONLY): " << Signal_CR_Data_above_Jpsi_integral*Template2D_above_Jpsi_diagonal_integral/Template2D_above_Jpsi_offDiagonal_integral << std::endl;

  //Scale to actual DATA
  h2D_template2D_above_Jpsi->Scale(Signal_CR_Data_above_Jpsi_integral*(1. + Template2D_above_Jpsi_diagonal_integral/Template2D_above_Jpsi_offDiagonal_integral));
  cout<<" h2D_template2D_above_Jpsi integral (after scale to actual data): "<< h2D_template2D_above_Jpsi->Integral() <<endl;
  TH2D * h2D_background_above_Jpsi = new TH2D( *h2D_template2D_above_Jpsi );
  h2D_background_above_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h2D_background_above_Jpsi->GetXaxis()->SetTitleOffset(0.93);
  h2D_background_above_Jpsi->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h2D_background_above_Jpsi->GetYaxis()->SetTitleOffset(0.85);
  h2D_background_above_Jpsi->GetZaxis()->SetTitle("Events/(0.04 GeV x 0.04 GeV)");
  h2D_background_above_Jpsi->GetZaxis()->CenterTitle(true);
  h2D_background_above_Jpsi->GetZaxis()->SetLabelFont(42);
  h2D_background_above_Jpsi->GetZaxis()->SetLabelSize(0.04);
  h2D_background_above_Jpsi->GetZaxis()->SetTitleSize(0.04);
  h2D_background_above_Jpsi->GetZaxis()->SetTitleOffset(1.42);
  h2D_background_above_Jpsi->GetZaxis()->SetTitleFont(42);

  TCanvas * c_template2D_m1_vs_m2_above_Jpsi = new TCanvas("c_template2D_m1_vs_m2_above_Jpsi", "c_template2D_m1_vs_m2_above_Jpsi", 0, 1320, 1044, 928);
  c_template2D_m1_vs_m2_above_Jpsi->SetCanvasSize(1040, 900);
  c_template2D_m1_vs_m2_above_Jpsi->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2_above_Jpsi->SetRightMargin(0.17);
  c_template2D_m1_vs_m2_above_Jpsi->SetTopMargin(0.05);
  c_template2D_m1_vs_m2_above_Jpsi->cd();
  c_template2D_m1_vs_m2_above_Jpsi->SetLogz();

  TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, nb);
  h2D_background_above_Jpsi->SetContour(nb);
  h2D_background_above_Jpsi->Draw("Cont4 Colz");
  c_template2D_m1_vs_m2_above_Jpsi->SaveAs("figures/Expected_2D_background_above_Jpsi.pdf");
  c_template2D_m1_vs_m2_above_Jpsi->SaveAs("figures/Expected_2D_background_above_Jpsi.png");
  c_template2D_m1_vs_m2_above_Jpsi->SaveAs("figures/Expected_2D_background_above_Jpsi.root");

  //**************************************************************************************
  //        Draw scatter plot at CR from data (SR blinded) Above Jpsi version
  //**************************************************************************************
  //Create pad to draw scatter plot without Logz
  TPad* pad_above_Jpsi = new TPad("pad_above_Jpsi", "pad_above_Jpsi", 0, 0, 1, 1);
  pad_above_Jpsi->Draw();
  pad_above_Jpsi->cd();
  pad_above_Jpsi->SetLeftMargin(0.121);
  pad_above_Jpsi->SetRightMargin(0.17);//48);
  pad_above_Jpsi->SetTopMargin(0.05);
  pad_above_Jpsi->SetFillColor(0);
  pad_above_Jpsi->SetFillStyle(4000);
  pad_above_Jpsi->SetBorderMode(0);
  pad_above_Jpsi->SetBorderSize(2);
  pad_above_Jpsi->SetTickx(1);
  pad_above_Jpsi->SetTicky(1);
  pad_above_Jpsi->SetFrameFillStyle(0);
  pad_above_Jpsi->SetFrameBorderMode(0);
  pad_above_Jpsi->SetFrameFillStyle(0);
  pad_above_Jpsi->SetFrameBorderMode(0);

  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->SetXTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->SetYTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi_tmp = new TH2D( *h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi_tmp->SetMarkerSize(1.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi_tmp->Draw("same");

  //**************************************************************************************************
  //    !!!BEGIN: Placeholder for unblinding signal (above Jpsi), after green light from pre-approval
  //**************************************************************************************************
  /*
  TH2D* h2D_dimudimu_signal_2D_above_Jpsi = (TH2D*)w->data("ds_dimudimu_signal_2D_above_Jpsi")->createHistogram("m1_above_Jpsi,m2_above_Jpsi", m_bins_above_Jpsi, m_bins_above_Jpsi);
  h2D_dimudimu_signal_2D_above_Jpsi->SetMarkerColor(kBlack);
  h2D_dimudimu_signal_2D_above_Jpsi->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_above_Jpsi->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2D_dimudimu_signal_2D_above_Jpsi->SetXTitle("");
  h2D_dimudimu_signal_2D_above_Jpsi->SetYTitle("");
  h2D_dimudimu_signal_2D_above_Jpsi->Draw("same");

  TH2D * h2D_dimudimu_signal_2D_above_Jpsi_tmp = new TH2D( *h2D_dimudimu_signal_2D_above_Jpsi);
  h2D_dimudimu_signal_2D_above_Jpsi_tmp->SetMarkerColor(kYellow);
  h2D_dimudimu_signal_2D_above_Jpsi_tmp->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_above_Jpsi_tmp->SetMarkerSize(1.0);
  h2D_dimudimu_signal_2D_above_Jpsi_tmp->Draw("same");
  */
  //**************************************************************************************************
  //    !!!END: Placeholder for unblinding signal (above Jpsi), after green light from pre-approval
  //**************************************************************************************************

  //corridorDnAboveJpsi->Draw("L"); corridorUpAboveJpsi->Draw("L"); txtHeader->Draw();
  //line3->Draw(); line4->Draw(); txtHeader->Draw();
  corridorDn->Draw("L"); corridorUp->Draw("L"); txtHeader->Draw();
  c_template2D_m1_vs_m2_above_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_above_Jpsi.pdf");
  c_template2D_m1_vs_m2_above_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_above_Jpsi.png");
  c_template2D_m1_vs_m2_above_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_above_Jpsi.root");

  //Validate for above J/psi ONLY
  //---------------------------------------------------------------------
  cout<<"                                                       " <<endl;
  cout<<"------ Start: validate for m1 (Above J/Psi ONLY) ------" <<endl;
  //---------------------------------------------------------------------
  TH1D *h1_control_Iso_offDiagonal_massC_data_above_Jpsi = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi")->createHistogram("m1_above_Jpsi", m_bins_above_Jpsi);
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->SetStats(0);
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetYaxis()->SetTitleOffset(0.9);
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massC_template_above_Jpsi = new TH1D( *h2D_template2D_above_Jpsi_offDiagonal->ProjectionX() );
  h1_control_Iso_offDiagonal_massC_template_above_Jpsi->Scale( h1_control_Iso_offDiagonal_massC_data_above_Jpsi->Integral()*1./h1_control_Iso_offDiagonal_massC_template_above_Jpsi->Integral() );
  h1_control_Iso_offDiagonal_massC_template_above_Jpsi->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massC_template_above_Jpsi->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massC_template_above_Jpsi->SetMarkerColor(kRed);
  for (Int_t i = 0; i < m_bins; ++i) h1_control_Iso_offDiagonal_massC_template_above_Jpsi->SetBinError(i, 0);

  TCanvas * c_control_Iso_offDiagonal_massC_above_Jpsi = new TCanvas("c_control_Iso_offDiagonal_massC_above_Jpsi", "c_control_Iso_offDiagonal_massC_above_Jpsi");
  c_control_Iso_offDiagonal_massC_above_Jpsi->cd();
  h1_control_Iso_offDiagonal_massC_data_above_Jpsi->Draw("e1");
  h1_control_Iso_offDiagonal_massC_template_above_Jpsi->Draw("HIST same");
  txtHeadervld->Draw();
  c_control_Iso_offDiagonal_massC_above_Jpsi->SaveAs("figures/Validation_m1_CR_above_Jpsi.pdf");
  c_control_Iso_offDiagonal_massC_above_Jpsi->SaveAs("figures/Validation_m1_CR_above_Jpsi.png");
  c_control_Iso_offDiagonal_massC_above_Jpsi->SaveAs("figures/Validation_m1_CR_above_Jpsi.root");
  //--------------------------
  //     Compatibility test
  //--------------------------
  //K-S test
  double KSprob_7 = h1_control_Iso_offDiagonal_massC_data_above_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_above_Jpsi);
  double KSdist_7 = h1_control_Iso_offDiagonal_massC_data_above_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_above_Jpsi, "M");
  cout<<"K-S test prob.: "<< KSprob_7 << "; dist.: "<< KSdist_7 <<endl;
  //Chisquare test
  auto func7 = [&](double *x, double*) { int ibin = h1_control_Iso_offDiagonal_massC_template_above_Jpsi->FindBin(x[0]); return h1_control_Iso_offDiagonal_massC_template_above_Jpsi->GetBinContent(ibin);};
  auto f7 = new TF1("f7", func7, h1_control_Iso_offDiagonal_massC_template_above_Jpsi->GetXaxis()->GetXmin(), h1_control_Iso_offDiagonal_massC_template_above_Jpsi->GetXaxis()->GetXmax(), 0);
  double BCchi2_7 = h1_control_Iso_offDiagonal_massC_data_above_Jpsi->Chisquare(f7, "L");
  double BCprob_7 = TMath::Prob(BCchi2_7, h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetNbinsX());
  cout<<"Chisquare test prob.: "<< BCprob_7 << "; chi2: "<< BCchi2_7 << "; ndof: " << h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetNbinsX() <<endl;
  //Toy experiments for calibration
  auto h_t7 = (TH1*) h1_control_Iso_offDiagonal_massC_data_above_Jpsi->Clone();//placeholder to fill pseudo data from template, same bins as data
  auto hKS_t7 = new TH1D("hKS_t7", "K-S distance", 100, 0, 1);//K-S distance distibution from pseudo exp.
  auto hBC_t7 = new TH1D("hBC_t7", "Baker-Cousins chi2", 100, 0, 300);//chi2 distibution from pseudo exp.
  int nKS_t7 = 0;
  int nBC_t7 = 0;
  int nentries_t7 = h1_control_Iso_offDiagonal_massC_data_above_Jpsi->Integral(1, h1_control_Iso_offDiagonal_massC_data_above_Jpsi->GetNbinsX());
  for (int i = 0; i < ntoys; ++i) {
    h_t7->Reset();
    h_t7->FillRandom("f7", nentries_t7);
    double KSdist_t7 = h_t7->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_above_Jpsi, "M");
    double BCchi2_t7 = h_t7->Chisquare(f7, "L");
    hKS_t7->Fill(KSdist_t7);
    hBC_t7->Fill(BCchi2_t7);
    if (KSdist_t7 > KSdist_7) nKS_t7++;
    if (BCchi2_t7 > BCchi2_7) nBC_t7++;
  }
  std::cout << "Corrected prob. for K-S  test: " << nKS_t7/double(ntoys) << std::endl;
  std::cout << "Corrected prob. for chi2 test: " << nBC_t7/double(ntoys) << std::endl;
  auto c7 = new TCanvas(); c7->Divide(1, 2);
  c7->cd(1); hKS_t7->Draw();
  c7->cd(2); hBC_t7->Draw();
  c7->SaveAs("figures/toys_m1_CR_above_Jpsi.root");
  cout<<"------ End: validate for m1 (Above J/Psi ONLY) ------" <<endl;

  //---------------------------------------------------------------------
  cout<<"                                                       " <<endl;
  cout<<"------ Start: validate for m2 (Above J/Psi ONLY) ------" <<endl;
  //---------------------------------------------------------------------
  TH1D *h1_control_Iso_offDiagonal_massF_data_above_Jpsi = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_above_Jpsi")->createHistogram("m2_above_Jpsi", m_bins_above_Jpsi);
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->SetStats(0);
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetYaxis()->SetTitleOffset(0.9);
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massF_template_above_Jpsi = new TH1D( *h2D_template2D_above_Jpsi_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_massF_template_above_Jpsi->Scale( h1_control_Iso_offDiagonal_massF_data_above_Jpsi->Integral()*1./h1_control_Iso_offDiagonal_massF_template_above_Jpsi->Integral() );
  h1_control_Iso_offDiagonal_massF_template_above_Jpsi->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massF_template_above_Jpsi->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massF_template_above_Jpsi->SetMarkerColor(kRed);
  for (Int_t i = 0; i < m_bins; ++i) h1_control_Iso_offDiagonal_massF_template_above_Jpsi->SetBinError(i, 0);

  TCanvas * c_control_Iso_offDiagonal_massF_above_Jpsi = new TCanvas("c_control_Iso_offDiagonal_massF_above_Jpsi", "c_control_Iso_offDiagonal_massF_above_Jpsi");
  c_control_Iso_offDiagonal_massF_above_Jpsi->cd();
  h1_control_Iso_offDiagonal_massF_data_above_Jpsi->Draw("e1");
  h1_control_Iso_offDiagonal_massF_template_above_Jpsi->Draw("HIST same");
  txtHeadervld->Draw();
  c_control_Iso_offDiagonal_massF_above_Jpsi->SaveAs("figures/Validation_m2_CR_above_Jpsi.pdf");
  c_control_Iso_offDiagonal_massF_above_Jpsi->SaveAs("figures/Validation_m2_CR_above_Jpsi.png");
  c_control_Iso_offDiagonal_massF_above_Jpsi->SaveAs("figures/Validation_m2_CR_above_Jpsi.root");
  //--------------------------
  //     Compatibility test
  //--------------------------
  //K-S test
  double KSprob_8 = h1_control_Iso_offDiagonal_massF_data_above_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_above_Jpsi);
  double KSdist_8 = h1_control_Iso_offDiagonal_massF_data_above_Jpsi->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_above_Jpsi, "M");
  cout<<"K-S test prob.: "<< KSprob_8 << "; dist.: "<< KSdist_8 <<endl;
  //Chisquare test
  auto func8 = [&](double *x, double*) { int ibin = h1_control_Iso_offDiagonal_massF_template_above_Jpsi->FindBin(x[0]); return h1_control_Iso_offDiagonal_massF_template_above_Jpsi->GetBinContent(ibin);};
  auto f8 = new TF1("f8", func8, h1_control_Iso_offDiagonal_massF_template_above_Jpsi->GetXaxis()->GetXmin(), h1_control_Iso_offDiagonal_massF_template_above_Jpsi->GetXaxis()->GetXmax(), 0);
  double BCchi2_8 = h1_control_Iso_offDiagonal_massF_data_above_Jpsi->Chisquare(f8, "L");
  double BCprob_8 = TMath::Prob(BCchi2_8, h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetNbinsX());
  cout<<"Chisquare test prob.: "<< BCprob_8 << "; chi2: "<< BCchi2_8 << "; ndof: " << h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetNbinsX() <<endl;
  //Toy experiments for calibration
  auto h_t8 = (TH1*) h1_control_Iso_offDiagonal_massF_data_above_Jpsi->Clone();//placeholder to fill pseudo data from template, same bins as data
  auto hKS_t8 = new TH1D("hKS_t8", "K-S distance", 100, 0, 1);//K-S distance distibution from pseudo exp.
  auto hBC_t8 = new TH1D("hBC_t8", "Baker-Cousins chi2", 100, 0, 300);//chi2 distibution from pseudo exp.
  int nKS_t8 = 0;
  int nBC_t8 = 0;
  int nentries_t8 = h1_control_Iso_offDiagonal_massF_data_above_Jpsi->Integral(1, h1_control_Iso_offDiagonal_massF_data_above_Jpsi->GetNbinsX());
  for (int i = 0; i < ntoys; ++i) {
    h_t8->Reset();
    h_t8->FillRandom("f8", nentries_t8);
    double KSdist_t8 = h_t8->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_above_Jpsi, "M");
    double BCchi2_t8 = h_t8->Chisquare(f8, "L");
    hKS_t8->Fill(KSdist_t8);
    hBC_t8->Fill(BCchi2_t8);
    if (KSdist_t8 > KSdist_8) nKS_t8++;
    if (BCchi2_t8 > BCchi2_8) nBC_t8++;
  }
  std::cout << "Corrected prob. for K-S  test: " << nKS_t8/double(ntoys) << std::endl;
  std::cout << "Corrected prob. for chi2 test: " << nBC_t8/double(ntoys) << std::endl;
  auto c8 = new TCanvas(); c8->Divide(1, 2);
  c8->cd(1); hKS_t8->Draw();
  c8->cd(2); hBC_t8->Draw();
  c8->SaveAs("figures/toys_m2_CR_above_Jpsi.root");
  cout<<"------ End: validate for m2 (Above J/Psi ONLY) ------" <<endl;
  cout<<"                                                     " <<endl;

  //==================
  //Above Upsilon ONLY
  //==================
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!   Note: the constructed data-driven template for this region is not accurate,
  //!!!         you should NOT trust the estimated bkg events at this section, it's intended for test ONLY
  //!!!         ONLY the signal dataset (2-dimu) event scatter plot is useful for reference
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  cout<<"****** Part C: Above Upsilon ONLY !!! Use with CAUTION !!! ******" <<endl;
  TH2D* h2D_template2D_above_Upsilon = (TH2D*)w->pdf("template2D_above_Upsilon")->createHistogram("m1_above_Upsilon,m2_above_Upsilon", m_bins_above_Upsilon, m_bins_above_Upsilon);//normalized
  TH2D *h2D_template2D_above_Upsilon_diagonal = (TH2D*)h2D_template2D_above_Upsilon->Clone();
  TH2D *h2D_template2D_above_Upsilon_offDiagonal = (TH2D*)h2D_template2D_above_Upsilon->Clone();

  for (int i = 1; i <= m_bins_above_Upsilon; i++) {
    for (int j = 1;j <= m_bins_above_Upsilon; j++) {
      double m_1 = h2D_template2D_above_Upsilon_offDiagonal->GetXaxis()->GetBinCenter(i);
      double m_2 = h2D_template2D_above_Upsilon_offDiagonal->GetYaxis()->GetBinCenter(j);
      //==============================
      //2018 mass window above Upsilon
      //==============================
      //if ( fabs(m_1 - m_2) < 5*(0.0472738 - 0.00591865*(m_1 + m_2)/2.0 + 0.00113991*pow((m_1 + m_2)/2.0, 2) - 2.62048e-05*pow((m_1 + m_2)/2.0, 3) + 1.92254e-07*pow((m_1 + m_2)/2.0, 4) ) ) {
      //if ( fabs(m_1 - m_2) < (kA + kB*(m_1 + m_2)/2.) ) {
      if ( fabs(m_1 - m_2) < BKG_cfg::My_MassWindow(m_1, m_2) ) {
        h2D_template2D_above_Upsilon_offDiagonal->SetBinContent(i, j, 0.);
      }
      else {
        h2D_template2D_above_Upsilon_diagonal->SetBinContent(i, j, 0.);
      }
    }
  }

  //Fractions of area for diagonal and offdiagonal in 2D template (above Upsilon)
  double Template2D_above_Upsilon_diagonal_integral  = h2D_template2D_above_Upsilon_diagonal->Integral();
  double Template2D_above_Upsilon_offDiagonal_integral  = h2D_template2D_above_Upsilon_offDiagonal->Integral();
  cout<<" -> Template2D (Above Upsilon ONLY) diagonal integral:    "<< Template2D_above_Upsilon_diagonal_integral <<endl;
  cout<<" -> Template2D (Above Upsilon ONLY) offDiagonal integral: "<< Template2D_above_Upsilon_offDiagonal_integral <<endl;

  //count 2-dimu data events at CR
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon")->createHistogram("m1_above_Upsilon,m2_above_Upsilon", m_bins_above_Upsilon, m_bins_above_Upsilon);
  double Signal_CR_Data_above_Upsilon_integral  = h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->Integral();
  cout<<"2 dimuon events in DATA at CR (Above Upsilon ONLY): " << Signal_CR_Data_above_Upsilon_integral <<endl;
  cout<<"Expected 2 dimuon events in DATA at SR (Above Upsilon ONLY): " << Signal_CR_Data_above_Upsilon_integral*Template2D_above_Upsilon_diagonal_integral/Template2D_above_Upsilon_offDiagonal_integral << std::endl;

  //Scale to actual DATA
  h2D_template2D_above_Upsilon->Scale(Signal_CR_Data_above_Upsilon_integral*(1. + Template2D_above_Upsilon_diagonal_integral/Template2D_above_Upsilon_offDiagonal_integral));
  cout<<" h2D_template2D_above_Upsilon integral (after scale to actual data): "<< h2D_template2D_above_Upsilon->Integral() <<endl;
  TH2D * h2D_background_above_Upsilon = new TH2D( *h2D_template2D_above_Upsilon );
  h2D_background_above_Upsilon->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h2D_background_above_Upsilon->GetXaxis()->SetTitleOffset(0.93);
  h2D_background_above_Upsilon->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h2D_background_above_Upsilon->GetYaxis()->SetTitleOffset(0.85);
  h2D_background_above_Upsilon->GetZaxis()->SetTitle("Events/(0.5 GeV x 0.5 GeV)");
  h2D_background_above_Upsilon->GetZaxis()->CenterTitle(true);
  h2D_background_above_Upsilon->GetZaxis()->SetLabelFont(42);
  h2D_background_above_Upsilon->GetZaxis()->SetLabelSize(0.04);
  h2D_background_above_Upsilon->GetZaxis()->SetTitleSize(0.04);
  h2D_background_above_Upsilon->GetZaxis()->SetTitleOffset(1.42);
  h2D_background_above_Upsilon->GetZaxis()->SetTitleFont(42);

  TCanvas * c_template2D_m1_vs_m2_above_Upsilon = new TCanvas("c_template2D_m1_vs_m2_above_Upsilon", "c_template2D_m1_vs_m2_above_Upsilon", 0, 1320, 1044, 928);
  c_template2D_m1_vs_m2_above_Upsilon->SetCanvasSize(1040, 900);
  c_template2D_m1_vs_m2_above_Upsilon->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2_above_Upsilon->SetRightMargin(0.17);
  c_template2D_m1_vs_m2_above_Upsilon->SetTopMargin(0.05);
  c_template2D_m1_vs_m2_above_Upsilon->cd();
  c_template2D_m1_vs_m2_above_Upsilon->SetLogz();

  TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, nb);
  h2D_background_above_Upsilon->SetContour(nb);
  h2D_background_above_Upsilon->Draw("Cont4 Colz");
  c_template2D_m1_vs_m2_above_Upsilon->SaveAs("figures/Expected_2D_background_above_Upsilon.pdf");
  c_template2D_m1_vs_m2_above_Upsilon->SaveAs("figures/Expected_2D_background_above_Upsilon.png");
  c_template2D_m1_vs_m2_above_Upsilon->SaveAs("figures/Expected_2D_background_above_Upsilon.root");

  //**************************************************************************************
  //        Draw scatter plot at CR from data (SR blinded) Above Upsilon version
  //**************************************************************************************
  //Create pad to draw scatter plot without Logz
  TPad* pad_above_Upsilon = new TPad("pad_above_Upsilon", "pad_above_Upsilon", 0, 0, 1, 1);
  pad_above_Upsilon->Draw();
  pad_above_Upsilon->cd();
  pad_above_Upsilon->SetLeftMargin(0.121);
  pad_above_Upsilon->SetRightMargin(0.17);//48);
  pad_above_Upsilon->SetTopMargin(0.05);
  pad_above_Upsilon->SetFillColor(0);
  pad_above_Upsilon->SetFillStyle(4000);
  pad_above_Upsilon->SetBorderMode(0);
  pad_above_Upsilon->SetBorderSize(2);
  pad_above_Upsilon->SetTickx(1);
  pad_above_Upsilon->SetTicky(1);
  pad_above_Upsilon->SetFrameFillStyle(0);
  pad_above_Upsilon->SetFrameBorderMode(0);
  pad_above_Upsilon->SetFrameFillStyle(0);
  pad_above_Upsilon->SetFrameBorderMode(0);

  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->SetXTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->SetYTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon_tmp = new TH2D( *h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon_tmp->SetMarkerSize(1.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon_tmp->Draw("same");

  //**************************************************************************************************
  //    !!!BEGIN: Placeholder for unblinding signal (above Upsilon), after green light from pre-approval
  //**************************************************************************************************
  /*
  TH2D* h2D_dimudimu_signal_2D_above_Upsilon = (TH2D*)w->data("ds_dimudimu_signal_2D_above_Upsilon")->createHistogram("m1_above_Upsilon,m2_above_Upsilon", m_bins_above_Upsilon, m_bins_above_Upsilon);
  h2D_dimudimu_signal_2D_above_Upsilon->SetMarkerColor(kBlack);
  h2D_dimudimu_signal_2D_above_Upsilon->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_above_Upsilon->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2D_dimudimu_signal_2D_above_Upsilon->SetXTitle("");
  h2D_dimudimu_signal_2D_above_Upsilon->SetYTitle("");
  h2D_dimudimu_signal_2D_above_Upsilon->Draw("same");

  TH2D * h2D_dimudimu_signal_2D_above_Upsilon_tmp = new TH2D( *h2D_dimudimu_signal_2D_above_Upsilon);
  h2D_dimudimu_signal_2D_above_Upsilon_tmp->SetMarkerColor(kYellow);
  h2D_dimudimu_signal_2D_above_Upsilon_tmp->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_above_Upsilon_tmp->SetMarkerSize(1.0);
  h2D_dimudimu_signal_2D_above_Upsilon_tmp->Draw("same");
  */
  //**************************************************************************************************
  //    !!!END: Placeholder for unblinding signal (above Upsilon), after green light from pre-approval
  //**************************************************************************************************

  //corridorDnAboveUpsilon->Draw("L"); corridorUpAboveUpsilon->Draw("L"); txtHeader->Draw();
  //line5->Draw(); line6->Draw(); txtHeader->Draw();
  corridorDn->Draw("L"); corridorUp->Draw("L"); txtHeader->Draw();
  c_template2D_m1_vs_m2_above_Upsilon->SaveAs("figures/DATA_and_Expected_2D_background_above_Upsilon.pdf");
  c_template2D_m1_vs_m2_above_Upsilon->SaveAs("figures/DATA_and_Expected_2D_background_above_Upsilon.png");
  c_template2D_m1_vs_m2_above_Upsilon->SaveAs("figures/DATA_and_Expected_2D_background_above_Upsilon.root");

  //-----------------------------------------------------------------------------------------------
  // Here a bit special:
  // just draw signal dataset at CR, don't draw the bkg template as it's not accurate above Upsilon
  // we don't do this for other regions
  //-----------------------------------------------------------------------------------------------
  TCanvas * c_data_m1_vs_m2_above_Upsilon = new TCanvas("c_data_m1_vs_m2_above_Upsilon", "c_data_m1_vs_m2_above_Upsilon", 0, 1320, 1044, 928);
  c_data_m1_vs_m2_above_Upsilon->SetCanvasSize(1040, 900);
  c_data_m1_vs_m2_above_Upsilon->SetLeftMargin(0.121);
  c_data_m1_vs_m2_above_Upsilon->SetRightMargin(0.17);
  c_data_m1_vs_m2_above_Upsilon->SetTopMargin(0.05);
  c_data_m1_vs_m2_above_Upsilon->cd();

  TPad* pad_above_Upsilon_data = new TPad("pad_above_Upsilon_data", "pad_above_Upsilon_data", 0, 0, 1, 1);
  pad_above_Upsilon_data->Draw();
  pad_above_Upsilon_data->cd();
  pad_above_Upsilon_data->SetLeftMargin(0.121);
  pad_above_Upsilon_data->SetRightMargin(0.17);
  pad_above_Upsilon_data->SetTopMargin(0.05);
  pad_above_Upsilon_data->SetFillColor(0);
  pad_above_Upsilon_data->SetFillStyle(4000);
  pad_above_Upsilon_data->SetBorderMode(0);
  pad_above_Upsilon_data->SetBorderSize(2);
  pad_above_Upsilon_data->SetTickx(1);
  pad_above_Upsilon_data->SetTicky(1);
  pad_above_Upsilon_data->SetFrameFillStyle(0);
  pad_above_Upsilon_data->SetFrameBorderMode(0);
  pad_above_Upsilon_data->SetFrameFillStyle(0);
  pad_above_Upsilon_data->SetFrameBorderMode(0);

  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->GetXaxis()->SetTitleOffset(0.93);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->GetYaxis()->SetTitleOffset(0.85);
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon->Draw();
  h2_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon_tmp->Draw("same");

  //corridorDnAboveUpsilon->Draw("L"); corridorUpAboveUpsilon->Draw("L"); txtHeader->Draw();
  //line5->Draw(); line6->Draw(); txtHeader->Draw();
  corridorDn->Draw("L"); corridorUp->Draw("L"); txtHeader->Draw();
  c_data_m1_vs_m2_above_Upsilon->SaveAs("figures/DATA_above_Upsilon.pdf");
  c_data_m1_vs_m2_above_Upsilon->SaveAs("figures/DATA_above_Upsilon.png");
  c_data_m1_vs_m2_above_Upsilon->SaveAs("figures/DATA_above_Upsilon.root");

  //Validate for above Upsilon ONLY, for test, use with caution!!!
  //---------------------------------------------------------------------------------------------
  cout<<"                                                                                  " <<endl;
  cout<<"------ Start: validate for m1 (Above Upsilon ONLY !!! Use with CAUTION !!!) ------" <<endl;
  //---------------------------------------------------------------------------------------------
  TH1D *h1_control_Iso_offDiagonal_massC_data_above_Upsilon = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon")->createHistogram("m1_above_Upsilon", m_bins_above_Upsilon);
  h1_control_Iso_offDiagonal_massC_data_above_Upsilon->SetStats(0);
  h1_control_Iso_offDiagonal_massC_data_above_Upsilon->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massC_data_above_Upsilon->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_Iso_offDiagonal_massC_data_above_Upsilon->GetYaxis()->SetTitle("Events/0.5 GeV");
  h1_control_Iso_offDiagonal_massC_data_above_Upsilon->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massC_template_above_Upsilon = new TH1D( *h2D_template2D_above_Upsilon_offDiagonal->ProjectionX() );
  h1_control_Iso_offDiagonal_massC_template_above_Upsilon->Scale( h1_control_Iso_offDiagonal_massC_data_above_Upsilon->Integral()*1./h1_control_Iso_offDiagonal_massC_template_above_Upsilon->Integral() );
  h1_control_Iso_offDiagonal_massC_template_above_Upsilon->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massC_template_above_Upsilon->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massC_template_above_Upsilon->SetMarkerColor(kRed);
  for (Int_t i = 0; i < m_bins; ++i) h1_control_Iso_offDiagonal_massC_template_above_Upsilon->SetBinError(i, 0);

  TCanvas * c_control_Iso_offDiagonal_massC_above_Upsilon = new TCanvas("c_control_Iso_offDiagonal_massC_above_Upsilon", "c_control_Iso_offDiagonal_massC_above_Upsilon");
  c_control_Iso_offDiagonal_massC_above_Upsilon->cd();
  h1_control_Iso_offDiagonal_massC_data_above_Upsilon->Draw("e1");
  h1_control_Iso_offDiagonal_massC_template_above_Upsilon->Draw("HIST same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massC_above_Upsilon->SaveAs("figures/Validation_m1_CR_above_Upsilon.pdf");
  c_control_Iso_offDiagonal_massC_above_Upsilon->SaveAs("figures/Validation_m1_CR_above_Upsilon.png");
  c_control_Iso_offDiagonal_massC_above_Upsilon->SaveAs("figures/Validation_m1_CR_above_Upsilon.root");
  //--------------------------
  //     Compatibility test
  //--------------------------
  //K-S test
  double KSprob_9 = h1_control_Iso_offDiagonal_massC_data_above_Upsilon->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_above_Upsilon);
  double KSdist_9 = h1_control_Iso_offDiagonal_massC_data_above_Upsilon->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_above_Upsilon, "M");
  cout<<"K-S test prob.: "<< KSprob_9 << "; dist.: "<< KSdist_9 <<endl;
  //Chisquare test
  auto func9 = [&](double *x, double*) { int ibin = h1_control_Iso_offDiagonal_massC_template_above_Upsilon->FindBin(x[0]); return h1_control_Iso_offDiagonal_massC_template_above_Upsilon->GetBinContent(ibin);};
  auto f9 = new TF1("f9", func9, h1_control_Iso_offDiagonal_massC_template_above_Upsilon->GetXaxis()->GetXmin(), h1_control_Iso_offDiagonal_massC_template_above_Upsilon->GetXaxis()->GetXmax(), 0);
  double BCchi2_9 = h1_control_Iso_offDiagonal_massC_data_above_Upsilon->Chisquare(f9, "L");
  double BCprob_9 = TMath::Prob(BCchi2_9, h1_control_Iso_offDiagonal_massC_data_above_Upsilon->GetNbinsX());
  cout<<"Chisquare test prob.: "<< BCprob_9 << "; chi2: "<< BCchi2_9 << "; ndof: " << h1_control_Iso_offDiagonal_massC_data_above_Upsilon->GetNbinsX() <<endl;
  //Toy experiments for calibration
  auto h_t9 = (TH1*) h1_control_Iso_offDiagonal_massC_data_above_Upsilon->Clone();//placeholder to fill pseudo data from template, same bins as data
  auto hKS_t9 = new TH1D("hKS_t9", "K-S distance", 100, 0, 1);//K-S distance distibution from pseudo exp.
  auto hBC_t9 = new TH1D("hBC_t9", "Baker-Cousins chi2", 100, 0, 300);//chi2 distibution from pseudo exp.
  int nKS_t9 = 0;
  int nBC_t9 = 0;
  int nentries_t9 = h1_control_Iso_offDiagonal_massC_data_above_Upsilon->Integral(1, h1_control_Iso_offDiagonal_massC_data_above_Upsilon->GetNbinsX());
  for (int i = 0; i < ntoys; ++i) {
    h_t9->Reset();
    h_t9->FillRandom("f9", nentries_t9);
    double KSdist_t9 = h_t9->KolmogorovTest(h1_control_Iso_offDiagonal_massC_template_above_Upsilon, "M");
    double BCchi2_t9 = h_t9->Chisquare(f9, "L");
    hKS_t9->Fill(KSdist_t9);
    hBC_t9->Fill(BCchi2_t9);
    if (KSdist_t9 > KSdist_9) nKS_t9++;
    if (BCchi2_t9 > BCchi2_9) nBC_t9++;
  }
  std::cout << "Corrected prob. for K-S  test: " << nKS_t9/double(ntoys) << std::endl;
  std::cout << "Corrected prob. for chi2 test: " << nBC_t9/double(ntoys) << std::endl;
  auto c9 = new TCanvas(); c9->Divide(1, 2);
  c9->cd(1); hKS_t9->Draw();
  c9->cd(2); hBC_t9->Draw();
  c9->SaveAs("figures/toys_m1_CR_above_Upsilon.root");
  cout<<"------ End: validate for m1 (Above Upsilon ONLY !!! Use with CAUTION !!!) ------" <<endl;

  //-----------------------------------------------------------------------------------------------
  cout<<"                                                                                  " <<endl;
  cout<<"------ Start: validate for m2 (Above Upsilon ONLY !!! Use with CAUTION !!!) ------" <<endl;
  //-----------------------------------------------------------------------------------------------
  TH1D *h1_control_Iso_offDiagonal_massF_data_above_Upsilon = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_above_Upsilon")->createHistogram("m2_above_Upsilon", m_bins_above_Upsilon);
  h1_control_Iso_offDiagonal_massF_data_above_Upsilon->SetStats(0);
  h1_control_Iso_offDiagonal_massF_data_above_Upsilon->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massF_data_above_Upsilon->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_Iso_offDiagonal_massF_data_above_Upsilon->GetYaxis()->SetTitle("Events/0.5 GeV");
  h1_control_Iso_offDiagonal_massF_data_above_Upsilon->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massF_template_above_Upsilon = new TH1D( *h2D_template2D_above_Upsilon_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_massF_template_above_Upsilon->Scale( h1_control_Iso_offDiagonal_massF_data_above_Upsilon->Integral()*1./h1_control_Iso_offDiagonal_massF_template_above_Upsilon->Integral() );
  h1_control_Iso_offDiagonal_massF_template_above_Upsilon->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massF_template_above_Upsilon->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massF_template_above_Upsilon->SetMarkerColor(kRed);
  for (Int_t i = 0; i < m_bins; ++i) h1_control_Iso_offDiagonal_massF_template_above_Upsilon->SetBinError(i, 0);

  TCanvas * c_control_Iso_offDiagonal_massF_above_Upsilon = new TCanvas("c_control_Iso_offDiagonal_massF_above_Upsilon", "c_control_Iso_offDiagonal_massF_above_Upsilon");
  c_control_Iso_offDiagonal_massF_above_Upsilon->cd();
  h1_control_Iso_offDiagonal_massF_data_above_Upsilon->Draw("e1");
  h1_control_Iso_offDiagonal_massF_template_above_Upsilon->Draw("HIST same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massF_above_Upsilon->SaveAs("figures/Validation_m2_CR_above_Upsilon.pdf");
  c_control_Iso_offDiagonal_massF_above_Upsilon->SaveAs("figures/Validation_m2_CR_above_Upsilon.png");
  c_control_Iso_offDiagonal_massF_above_Upsilon->SaveAs("figures/Validation_m2_CR_above_Upsilon.root");
  //--------------------------
  //     Compatibility test
  //--------------------------
  //K-S test
  double KSprob_0 = h1_control_Iso_offDiagonal_massF_data_above_Upsilon->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_above_Upsilon);
  double KSdist_0 = h1_control_Iso_offDiagonal_massF_data_above_Upsilon->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_above_Upsilon, "M");
  cout<<"K-S test prob.: "<< KSprob_0 << "; dist.: "<< KSdist_0 <<endl;
  //Chisquare test
  auto func0 = [&](double *x, double*) { int ibin = h1_control_Iso_offDiagonal_massF_template_above_Upsilon->FindBin(x[0]); return h1_control_Iso_offDiagonal_massF_template_above_Upsilon->GetBinContent(ibin);};
  auto f0 = new TF1("f0", func0, h1_control_Iso_offDiagonal_massF_template_above_Upsilon->GetXaxis()->GetXmin(), h1_control_Iso_offDiagonal_massF_template_above_Upsilon->GetXaxis()->GetXmax(), 0);
  double BCchi2_0 = h1_control_Iso_offDiagonal_massF_data_above_Upsilon->Chisquare(f0, "L");
  double BCprob_0 = TMath::Prob(BCchi2_0, h1_control_Iso_offDiagonal_massF_data_above_Upsilon->GetNbinsX());
  cout<<"Chisquare test prob.: "<< BCprob_0 << "; chi2: "<< BCchi2_0 << "; ndof: " << h1_control_Iso_offDiagonal_massF_data_above_Upsilon->GetNbinsX() <<endl;
  //Toy experiments for calibration
  auto h_t0 = (TH1*) h1_control_Iso_offDiagonal_massF_data_above_Upsilon->Clone();//placeholder to fill pseudo data from template, same bins as data
  auto hKS_t0 = new TH1D("hKS_t0", "K-S distance", 100, 0, 1);//K-S distance distibution from pseudo exp.
  auto hBC_t0 = new TH1D("hBC_t0", "Baker-Cousins chi2", 100, 0, 300);//chi2 distibution from pseudo exp.
  int nKS_t0 = 0;
  int nBC_t0 = 0;
  int nentries_t0 = h1_control_Iso_offDiagonal_massF_data_above_Upsilon->Integral(1, h1_control_Iso_offDiagonal_massF_data_above_Upsilon->GetNbinsX());
  for (int i = 0; i < ntoys; ++i) {
    h_t0->Reset();
    h_t0->FillRandom("f0", nentries_t0);
    double KSdist_t0 = h_t0->KolmogorovTest(h1_control_Iso_offDiagonal_massF_template_above_Upsilon, "M");
    double BCchi2_t0 = h_t0->Chisquare(f0, "L");
    hKS_t0->Fill(KSdist_t0);
    hBC_t0->Fill(BCchi2_t0);
    if (KSdist_t0 > KSdist_0) nKS_t0++;
    if (BCchi2_t0 > BCchi2_0) nBC_t0++;
  }
  std::cout << "Corrected prob. for K-S  test: " << nKS_t0/double(ntoys) << std::endl;
  std::cout << "Corrected prob. for chi2 test: " << nBC_t0/double(ntoys) << std::endl;
  auto c0 = new TCanvas(); c0->Divide(1, 2);
  c0->cd(1); hKS_t0->Draw();
  c0->cd(2); hBC_t0->Draw();
  c0->SaveAs("figures/toys_m2_CR_above_Upsilon.root");
  cout<<"------ End: validate for m2 (Above Upsilon ONLY !!! Use with CAUTION !!!) ------" <<endl;
  cout<<"                                                                                " <<endl;
}
