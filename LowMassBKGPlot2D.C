//*****************************************************************************************************
//* cmsenv                                                                                            *
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

  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  //Used for excluding the J/psi region when constructing 1D/2D template
  const double       m_Jpsi_dn = 2.72;
  const double       m_Jpsi_up = 3.24;
  const unsigned int m_bins_below_Jpsi = 63;//bin size is ~0.04GeV, as above
  const unsigned int m_bins_above_Jpsi = 144;

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

  //**************************************************************************************
  //                   Draw 2D background (scale 2D template/pdf to data yield)
  //**************************************************************************************
  //Create and fill ROOT 2D histogram with sampling of 2D pdf, normalized to 1
  TH2D* h2D_template2D = (TH2D*)w->pdf("template2D")->createHistogram("m1,m2", m_bins, m_bins);
  TH2D *h2D_template2D_diagonal = (TH2D*)h2D_template2D->Clone();
  TH2D *h2D_template2D_offDiagonal = (TH2D*)h2D_template2D->Clone();

  for(int i=1;i<=m_bins;i++) {
    for(int j=1;j<=m_bins;j++) {
      double m_1 = h2D_template2D_offDiagonal->GetXaxis()->GetBinCenter(i);
      double m_2 = h2D_template2D_offDiagonal->GetYaxis()->GetBinCenter(j);
      //*************************
      //2017 mass consistency cut
      //*************************
      if ( fabs(m_1 - m_2) < 3*(0.003044 + 0.007025*(m_1+m_2)/2.0 + 0.000053*(m_1+m_2)*(m_1+m_2)/4.0) ) {
        h2D_template2D_offDiagonal->SetBinContent(i, j, 0.);
      }
      else {
        h2D_template2D_diagonal->SetBinContent(i, j, 0.);
      }
    }
  }

  //Fractions of area for diagonal and offdiagonal in 2D template
  double Template2D_diagonal_integral  = h2D_template2D_diagonal->Integral();
  double Template2D_offDiagonal_integral  = h2D_template2D_offDiagonal->Integral();
  cout<<" -> Template2D_diagonal integral:    "<< Template2D_diagonal_integral <<endl;
  cout<<" -> Template2D_offDiagonal integral: "<< Template2D_offDiagonal_integral <<endl;

  //count 2-dimu data events at CR
  TH2D* h2_dimudimu_control_Iso_offDiagonal_2D = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1,m2", m_bins, m_bins);
  double Signal_CR_Data_integral  = h2_dimudimu_control_Iso_offDiagonal_2D->Integral();
  cout<<"2 dimuon events in DATA at CR: " << Signal_CR_Data_integral <<endl;
  cout<<"Expected 2 dimuon events in DATA at SR: " << Signal_CR_Data_integral*Template2D_diagonal_integral/Template2D_offDiagonal_integral << std::endl;

  //Scale 2D template to total EXPECTED # of events in 2D plane
  h2D_template2D->Scale(Signal_CR_Data_integral*(1+Template2D_diagonal_integral/Template2D_offDiagonal_integral));
  TH2D * h2_background = new TH2D( *h2D_template2D );
  h2_background->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h2_background->GetXaxis()->SetTitleOffset(0.93);
  h2_background->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h2_background->GetYaxis()->SetTitleOffset(0.85);
  h2_background->GetZaxis()->SetTitle("Events/(0.04 GeV x 0.04 GeV)");
  h2_background->GetZaxis()->CenterTitle(true);
  h2_background->GetZaxis()->SetLabelFont(42);
  h2_background->GetZaxis()->SetLabelSize(0.04);
  h2_background->GetZaxis()->SetTitleSize(0.04);
  h2_background->GetZaxis()->SetTitleOffset(1.42);
  h2_background->GetZaxis()->SetTitleFont(42);

  //gStyle->SetPalette(52); //Grey Scale
  Double_t Red[2]    = {1.00, 0.00};
  Double_t Green[2]  = {1.00, 0.00};
  Double_t Blue[2]   = {1.00, 0.00};
  Double_t Length[2] = {0.00, 1.00};
  Int_t nb=50;
  TColor::CreateGradientColorTable(2,Length,Red,Green,Blue,nb);
  h2_background->SetContour(nb);
  h2_background->Draw("Cont4 Colz");
  c_template2D_m1_vs_m2->SaveAs("figures/Expected_2D_background.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures/Expected_2D_background.png");
  c_template2D_m1_vs_m2->SaveAs("figures/Expected_2D_background.root");

  //**************************************************************************************
  //        Draw scatter plot at CR from data (SR blinded)
  //**************************************************************************************
  //Create pad to draw scatter plot without Logz because c_template2D_m1_vs_m2 is set to Logz
  TPad* pad = new TPad("pad", "pad", 0, 0, 1, 1);
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

  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerColor(kBlack);
  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2_dimudimu_control_Iso_offDiagonal_2D->SetXTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D->SetYTitle("");
  h2_dimudimu_control_Iso_offDiagonal_2D->Draw("same");

  TH2D * h2_dimudimu_control_Iso_offDiagonal_2D_tmp = new TH2D( *h2_dimudimu_control_Iso_offDiagonal_2D);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerColor(kWhite);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerStyle(20);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->SetMarkerSize(1.0);
  h2_dimudimu_control_Iso_offDiagonal_2D_tmp->Draw("same");

  //**************************************************************************************
  //    !!!BEGIN: Placeholder for unblinding signal, after green light from pre-approval
  //**************************************************************************************
  /*
  TH2D* h2_dimudimu_signal_2D = (TH2D*)w->data("ds_dimudimu_signal_2D")->createHistogram("m1,m2", m_bins, m_bins);
  h2_dimudimu_signal_2D->SetMarkerColor(kBlack);
  h2_dimudimu_signal_2D->SetMarkerStyle(22);
  h2_dimudimu_signal_2D->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2_dimudimu_signal_2D->SetXTitle("");
  h2_dimudimu_signal_2D->SetYTitle("");
  h2_dimudimu_signal_2D->Draw("same");

  TH2D * h2_dimudimu_signal_2D_tmp = new TH2D( *h2_dimudimu_signal_2D);
  h2_dimudimu_signal_2D_tmp->SetMarkerColor(kYellow);
  h2_dimudimu_signal_2D_tmp->SetMarkerStyle(22);
  h2_dimudimu_signal_2D_tmp->SetMarkerSize(1.0);
  h2_dimudimu_signal_2D_tmp->Draw("same");
  */
  //************************************************************************************
  //    !!!END: Placeholder for unblinding signal, after green light from pre-approval
  //************************************************************************************

  //************************************************************************************
  //          Pre-calculated m1 and m2 values for drawing the corridor curves
  //************************************************************************************
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
  c_template2D_m1_vs_m2->SaveAs("figures/DATA_and_Expected_2D_background.pdf");
  c_template2D_m1_vs_m2->SaveAs("figures/DATA_and_Expected_2D_background.png");
  c_template2D_m1_vs_m2->SaveAs("figures/DATA_and_Expected_2D_background.root");

  //************************************************************************************
  //           Validate the method with 2 dimu events at CR (iso & no-iso)
  //************************************************************************************
  //Validate m1 with non-iso data at CR
  TH1D *h1_control_offDiagonal_massC_data = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_2D")->createHistogram("m1",m_bins);
  h1_control_offDiagonal_massC_data->SetStats(0);
  h1_control_offDiagonal_massC_data->SetMarkerStyle(20);
  h1_control_offDiagonal_massC_data->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_offDiagonal_massC_data->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_offDiagonal_massC_data->GetYaxis()->SetRangeUser(0.,350.);

  TH1D *h1_control_offDiagonal_massC_template = new TH1D( *h2D_template2D_offDiagonal->ProjectionX() );
  h1_control_offDiagonal_massC_template->Scale( h1_control_offDiagonal_massC_data->Integral() / h1_control_offDiagonal_massC_template->Integral() );
  h1_control_offDiagonal_massC_template->SetLineColor(kRed);
  h1_control_offDiagonal_massC_template->SetLineWidth(2);
  h1_control_offDiagonal_massC_template->SetMarkerColor(kRed);

  TCanvas * c_control_offDiagonal_massC = new TCanvas("c_control_offDiagonal_massC", "c_control_offDiagonal_massC");
  c_control_offDiagonal_massC->cd();
  h1_control_offDiagonal_massC_data->Draw("e1");
  h1_control_offDiagonal_massC_template->Draw("HIST same");
  txtHeader->Draw();
  c_control_offDiagonal_massC->SaveAs("figures/Validation_m1_no_iso_CR.pdf");
  c_control_offDiagonal_massC->SaveAs("figures/Validation_m1_no_iso_CR.png");
  c_control_offDiagonal_massC->SaveAs("figures/Validation_m1_no_iso_CR.root");

  //Validate m1 with iso data at CR: very low stats, just FYI
  TH1D *h1_control_Iso_offDiagonal_massC_data = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m1",m_bins);
  h1_control_Iso_offDiagonal_massC_data->SetStats(0);
  h1_control_Iso_offDiagonal_massC_data->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massC_data->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_Iso_offDiagonal_massC_data->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massC_data->GetYaxis()->SetRangeUser(0.,35.);

  TH1D *h1_control_Iso_offDiagonal_massC_template = new TH1D( *h2D_template2D_offDiagonal->ProjectionX() );
  h1_control_Iso_offDiagonal_massC_template->Scale( h1_control_Iso_offDiagonal_massC_data->Integral() / h1_control_Iso_offDiagonal_massC_template->Integral() );
  h1_control_Iso_offDiagonal_massC_template->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massC_template->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massC_template->SetMarkerColor(kRed);

  TCanvas * c_control_Iso_offDiagonal_massC = new TCanvas("c_control_Iso_offDiagonal_massC", "c_control_Iso_offDiagonal_massC");
  c_control_Iso_offDiagonal_massC->cd();
  h1_control_Iso_offDiagonal_massC_data->Draw("e1");
  h1_control_Iso_offDiagonal_massC_template->Draw("HIST same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massC->SaveAs("figures/Validation_m1_iso_CR.pdf");
  c_control_Iso_offDiagonal_massC->SaveAs("figures/Validation_m1_iso_CR.png");
  c_control_Iso_offDiagonal_massC->SaveAs("figures/Validation_m1_iso_CR.root");

  //Validate m2 with no-iso data at CR
  TH1D *h1_control_offDiagonal_massF_data = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_2D")->createHistogram("m2",m_bins);
  h1_control_offDiagonal_massF_data->SetStats(0);
  h1_control_offDiagonal_massF_data->SetMarkerStyle(20);
  h1_control_offDiagonal_massF_data->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_offDiagonal_massF_data->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_offDiagonal_massF_data->GetYaxis()->SetRangeUser(0.,250.);

  TH1D *h1_control_offDiagonal_massF_template = new TH1D( *h2D_template2D_offDiagonal->ProjectionY() );
  h1_control_offDiagonal_massF_template->Scale( h1_control_offDiagonal_massF_data->Integral() / h1_control_offDiagonal_massF_template->Integral() );
  h1_control_offDiagonal_massF_template->SetLineColor(kRed);
  h1_control_offDiagonal_massF_template->SetLineWidth(2);
  h1_control_offDiagonal_massF_template->SetMarkerColor(kRed);

  TCanvas * c_control_offDiagonal_massF = new TCanvas("c_control_offDiagonal_massF", "c_control_offDiagonal_massF");
  c_control_offDiagonal_massF->cd();
  h1_control_offDiagonal_massF_data->Draw("e1");
  h1_control_offDiagonal_massF_template->Draw("HIST same");
  txtHeader->Draw();
  c_control_offDiagonal_massF->SaveAs("figures/Validation_m2_no_iso_CR.pdf");
  c_control_offDiagonal_massF->SaveAs("figures/Validation_m2_no_iso_CR.png");
  c_control_offDiagonal_massF->SaveAs("figures/Validation_m2_no_iso_CR.root");

  //Validate m2 with iso data at CR: limited stats, just FYI
  TH1D *h1_control_Iso_offDiagonal_massF_data = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D")->createHistogram("m2",m_bins);
  h1_control_Iso_offDiagonal_massF_data->SetStats(0);
  h1_control_Iso_offDiagonal_massF_data->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massF_data->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_Iso_offDiagonal_massF_data->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massF_data->GetYaxis()->SetRangeUser(0.,10.);

  TH1D *h1_control_Iso_offDiagonal_massF_template = new TH1D( *h2D_template2D_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_massF_template->Scale( h1_control_Iso_offDiagonal_massF_data->Integral() / h1_control_Iso_offDiagonal_massF_template->Integral() );
  h1_control_Iso_offDiagonal_massF_template->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massF_template->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massF_template->SetMarkerColor(kRed);

  TCanvas * c_control_Iso_offDiagonal_massF = new TCanvas("c_control_Iso_offDiagonal_massF", "c_control_Iso_offDiagonal_massF");
  c_control_Iso_offDiagonal_massF->cd();
  h1_control_Iso_offDiagonal_massF_data->Draw("e1");
  h1_control_Iso_offDiagonal_massF_template->Draw("HIST same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massF->SaveAs("figures/Validation_m2_iso_CR.pdf");
  c_control_Iso_offDiagonal_massF->SaveAs("figures/Validation_m2_iso_CR.png");
  c_control_Iso_offDiagonal_massF->SaveAs("figures/Validation_m2_iso_CR.root");

  //**************************************************************************************
  //         Similar version to above except Exclude J/psi when construct 2D template
  //**************************************************************************************
  /*
  //TH2D* h2D_template2D_m1_belowJpsi_m2_belowJpsi = (TH2D*)w->pdf("template2D_m1_belowJpsi_m2_belowJpsi")->createHistogram("m1_below_Jpsi,m2_below_Jpsi", m_bins_below_Jpsi, m_bins_below_Jpsi);
  //TH2D* h2D_template2D_m1_belowJpsi_m2_aboveJpsi = (TH2D*)w->pdf("template2D_m1_belowJpsi_m2_aboveJpsi")->createHistogram("m1_below_Jpsi,m2_above_Jpsi", m_bins_below_Jpsi, m_bins_above_Jpsi);
  //TH2D* h2D_template2D_m1_aboveJpsi_m2_belowJpsi = (TH2D*)w->pdf("template2D_m1_aboveJpsi_m2_belowJpsi")->createHistogram("m1_above_Jpsi,m2_below_Jpsi", m_bins_above_Jpsi, m_bins_below_Jpsi);
  //TH2D* h2D_template2D_m1_aboveJpsi_m2_aboveJpsi = (TH2D*)w->pdf("template2D_m1_aboveJpsi_m2_aboveJpsi")->createHistogram("m1_above_Jpsi,m2_above_Jpsi", m_bins_above_Jpsi, m_bins_above_Jpsi);

  //Add entries from above hist into one hist over all range
  TH2D* h2D_template2D_exclude_Jpsi = new TH2D("h2D_template2D_exclude_Jpsi", "", m_bins, m_min, m_max, m_bins, m_min, m_max);
  for(int i=1;i<=m_bins;i++) {
    for(int j=1;j<=m_bins;j++) {
      if(i<=m_bins_below_Jpsi && j<=m_bins_below_Jpsi){ //i:1-63, j:1-63
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h2D_template2D_m1_belowJpsi_m2_belowJpsi->GetBinContent(i, j) );
      }
      else if(i<=m_bins_below_Jpsi && j>= m_bins-m_bins_above_Jpsi+1){ //i:1-63, j:77-220
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h2D_template2D_m1_belowJpsi_m2_aboveJpsi->GetBinContent(i, j-m_bins+m_bins_above_Jpsi ) );
      }
      else if(i>= m_bins-m_bins_above_Jpsi+1 && j<= m_bins_below_Jpsi){//i:77-220, j:1-63
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h2D_template2D_m1_aboveJpsi_m2_belowJpsi->GetBinContent(i-m_bins+m_bins_above_Jpsi, j ) );
      }
      else if(i>= m_bins-m_bins_above_Jpsi+1 && j>= m_bins-m_bins_above_Jpsi+1){//i:77-220, j:77-220
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h2D_template2D_m1_aboveJpsi_m2_aboveJpsi->GetBinContent(i-m_bins+m_bins_above_Jpsi, j-m_bins+m_bins_above_Jpsi ) );
      }
      else{//excluded Jpsi region
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, 0);
      }
    }//end j
  }//end i
  */

  TH1D* h1D_template1D_m1_below_Jpsi = (TH1D*)w->pdf("template1D_m1_below_Jpsi")->createHistogram("m1_below_Jpsi", m_bins_below_Jpsi);//normalized to 1
  TH1D* h1D_template1D_m1_above_Jpsi = (TH1D*)w->pdf("template1D_m1_above_Jpsi")->createHistogram("m1_above_Jpsi", m_bins_above_Jpsi);
  TH1D* h1D_template1D_m2_below_Jpsi = (TH1D*)w->pdf("template1D_m2_below_Jpsi")->createHistogram("m2_below_Jpsi", m_bins_below_Jpsi);
  TH1D* h1D_template1D_m2_above_Jpsi = (TH1D*)w->pdf("template1D_m2_above_Jpsi")->createHistogram("m2_above_Jpsi", m_bins_above_Jpsi);
  //Scale them according to the actual fitted data entry 2017
  cout<<" ds_dimuorphan_bg_m1_below_Jpsi entries: "<< w->data("ds_dimuorphan_bg_m1_below_Jpsi")->sumEntries() <<endl;
  cout<<" ds_dimuorphan_bg_m1_above_Jpsi entries: "<< w->data("ds_dimuorphan_bg_m1_above_Jpsi")->sumEntries() <<endl;
  cout<<" ds_dimuorphan_bg_m2_below_Jpsi entries: "<< w->data("ds_dimuorphan_bg_m2_below_Jpsi")->sumEntries() <<endl;
  cout<<" ds_dimuorphan_bg_m2_above_Jpsi entries: "<< w->data("ds_dimuorphan_bg_m2_above_Jpsi")->sumEntries() <<endl;
  h1D_template1D_m1_below_Jpsi->Scale(w->data("ds_dimuorphan_bg_m1_below_Jpsi")->sumEntries());
  h1D_template1D_m1_above_Jpsi->Scale(w->data("ds_dimuorphan_bg_m1_above_Jpsi")->sumEntries());
  h1D_template1D_m2_below_Jpsi->Scale(w->data("ds_dimuorphan_bg_m2_below_Jpsi")->sumEntries());
  h1D_template1D_m2_above_Jpsi->Scale(w->data("ds_dimuorphan_bg_m2_above_Jpsi")->sumEntries());

  //Combine into 2D histo use 1D histo
  TH2D* h2D_template2D_exclude_Jpsi = new TH2D("h2D_template2D_exclude_Jpsi", "", m_bins, m_min, m_max, m_bins, m_min, m_max);
  for(int i=1;i<=m_bins;i++) {
    for(int j=1;j<=m_bins;j++) {
      if(i<=m_bins_below_Jpsi && j<=m_bins_below_Jpsi){ //i:1-63, j:1-63
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h1D_template1D_m1_below_Jpsi->GetBinContent(i)*h1D_template1D_m2_below_Jpsi->GetBinContent(j) );
      }
      else if(i<=m_bins_below_Jpsi && j>= m_bins-m_bins_above_Jpsi+1){ //i:1-63, j:77-220
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h1D_template1D_m1_below_Jpsi->GetBinContent(i)*h1D_template1D_m2_above_Jpsi->GetBinContent(j-m_bins+m_bins_above_Jpsi) );
      }
      else if(i>= m_bins-m_bins_above_Jpsi+1 && j<= m_bins_below_Jpsi){//i:77-220, j:1-63
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h1D_template1D_m1_above_Jpsi->GetBinContent(i-m_bins+m_bins_above_Jpsi)*h1D_template1D_m2_below_Jpsi->GetBinContent(j) );
      }
      else if(i>= m_bins-m_bins_above_Jpsi+1 && j>= m_bins-m_bins_above_Jpsi+1){//i:77-220, j:77-220
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, h1D_template1D_m1_above_Jpsi->GetBinContent(i-m_bins+m_bins_above_Jpsi)*h1D_template1D_m2_above_Jpsi->GetBinContent(j-m_bins+m_bins_above_Jpsi) );
      }
      else{//excluded Jpsi region
        h2D_template2D_exclude_Jpsi->SetBinContent(i, j, 0);//try not fill anything
      }
    }//end j
  }//end i

  cout<<" h2D_template2D_exclude_Jpsi integral: "<< h2D_template2D_exclude_Jpsi->Integral() <<endl;
  h2D_template2D_exclude_Jpsi->Scale(1./h2D_template2D_exclude_Jpsi->Integral());
  cout<<" h2D_template2D_exclude_Jpsi integral (after scale, should be 1.0): "<< h2D_template2D_exclude_Jpsi->Integral() <<endl;
  TH2D *h2D_template2D_exclude_Jpsi_diagonal = (TH2D*)h2D_template2D_exclude_Jpsi->Clone();
  TH2D *h2D_template2D_exclude_Jpsi_offDiagonal = (TH2D*)h2D_template2D_exclude_Jpsi->Clone();

  for(int i=1;i<=m_bins;i++) {
    for(int j=1;j<=m_bins;j++) {
      double m_1 = h2D_template2D_exclude_Jpsi_offDiagonal->GetXaxis()->GetBinCenter(i);
      double m_2 = h2D_template2D_exclude_Jpsi_offDiagonal->GetYaxis()->GetBinCenter(j);
      //*************************
      //2017 mass consistency cut
      //*************************
      if ( fabs(m_1 - m_2) < 3*(0.003044 + 0.007025*(m_1+m_2)/2.0 + 0.000053*(m_1+m_2)*(m_1+m_2)/4.0) ) {
        h2D_template2D_exclude_Jpsi_offDiagonal->SetBinContent(i, j, 0.);
      }
      else {
        h2D_template2D_exclude_Jpsi_diagonal->SetBinContent(i, j, 0.);
      }
    }
  }

  double h2D_template2D_exclude_Jpsi_diagonal_integral     = h2D_template2D_exclude_Jpsi_diagonal->Integral();
  double h2D_template2D_exclude_Jpsi_offDiagonal_integral  = h2D_template2D_exclude_Jpsi_offDiagonal->Integral();
  cout<<" -> Template2D_diagonal (exclude J/psi) integral:    "<< h2D_template2D_exclude_Jpsi_diagonal_integral <<endl;
  cout<<" -> Template2D_offDiagonal (exclude J/psi) integral: "<< h2D_template2D_exclude_Jpsi_offDiagonal_integral <<endl;

  //count 2-dimu data events at CR
  TH2D* h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi = (TH2D*)w->data("ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi")->createHistogram("m1,m2", m_bins, m_bins);
  double h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_integral  = h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->Integral();
  cout<<"2 dimuon events in DATA at CR (exclude J/psi): " << h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_integral <<endl;
  cout<<"Expected 2 dimuon events in DATA at SR (exclude J/psi): " << h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_integral*h2D_template2D_exclude_Jpsi_diagonal_integral/h2D_template2D_exclude_Jpsi_offDiagonal_integral << std::endl;

  //Scale to actual DATA
  h2D_template2D_exclude_Jpsi->Scale(h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_integral*(1+h2D_template2D_exclude_Jpsi_diagonal_integral/h2D_template2D_exclude_Jpsi_offDiagonal_integral));
  cout<<" h2D_template2D_exclude_Jpsi integral (after scale to actual data): "<< h2D_template2D_exclude_Jpsi->Integral() <<endl;
  TH2D * h2D_background_exclude_Jpsi = new TH2D( *h2D_template2D_exclude_Jpsi );
  h2D_background_exclude_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h2D_background_exclude_Jpsi->GetXaxis()->SetTitleOffset(0.93);
  h2D_background_exclude_Jpsi->GetYaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h2D_background_exclude_Jpsi->GetYaxis()->SetTitleOffset(0.85);
  h2D_background_exclude_Jpsi->GetZaxis()->SetTitle("Events/(0.04 GeV x 0.04 GeV)");
  h2D_background_exclude_Jpsi->GetZaxis()->CenterTitle(true);
  h2D_background_exclude_Jpsi->GetZaxis()->SetLabelFont(42);
  h2D_background_exclude_Jpsi->GetZaxis()->SetLabelSize(0.04);
  h2D_background_exclude_Jpsi->GetZaxis()->SetTitleSize(0.04);
  h2D_background_exclude_Jpsi->GetZaxis()->SetTitleOffset(1.42);
  h2D_background_exclude_Jpsi->GetZaxis()->SetTitleFont(42);

  TCanvas * c_template2D_m1_vs_m2_exclude_Jpsi = new TCanvas("c_template2D_m1_vs_m2_exclude_Jpsi", "c_template2D_m1_vs_m2_exclude_Jpsi", 0, 1320, 1044, 928);
  c_template2D_m1_vs_m2_exclude_Jpsi->SetCanvasSize(1040, 900);
  c_template2D_m1_vs_m2_exclude_Jpsi->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2_exclude_Jpsi->SetRightMargin(0.17);
  c_template2D_m1_vs_m2_exclude_Jpsi->SetTopMargin(0.05);
  c_template2D_m1_vs_m2_exclude_Jpsi->cd();
  c_template2D_m1_vs_m2_exclude_Jpsi->SetLogz();
  TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, nb);
  h2D_background_exclude_Jpsi->SetContour(nb);
  h2D_background_exclude_Jpsi->Draw("Cont4 Colz");
  c_template2D_m1_vs_m2_exclude_Jpsi->SaveAs("figures/Expected_2D_background_exclude_Jpsi.pdf");
  c_template2D_m1_vs_m2_exclude_Jpsi->SaveAs("figures/Expected_2D_background_exclude_Jpsi.png");
  c_template2D_m1_vs_m2_exclude_Jpsi->SaveAs("figures/Expected_2D_background_exclude_Jpsi.root");

  //**************************************************************************************
  //        Draw scatter plot at CR from data (SR blinded) Exclude Jpsi version
  //**************************************************************************************
  //Create pad to draw scatter plot without Logz
  TPad* pad_exclude_Jpsi = new TPad("pad_exclude_Jpsi", "pad_exclude_Jpsi", 0, 0, 1, 1);
  pad_exclude_Jpsi->Draw();
  pad_exclude_Jpsi->cd();
  pad_exclude_Jpsi->SetLeftMargin(0.121);
  pad_exclude_Jpsi->SetRightMargin(0.17);//48);
  pad_exclude_Jpsi->SetTopMargin(0.05);
  pad_exclude_Jpsi->SetFillColor(0);
  pad_exclude_Jpsi->SetFillStyle(4000);
  pad_exclude_Jpsi->SetBorderMode(0);
  pad_exclude_Jpsi->SetBorderSize(2);
  pad_exclude_Jpsi->SetTickx(1);
  pad_exclude_Jpsi->SetTicky(1);
  pad_exclude_Jpsi->SetFrameFillStyle(0);
  pad_exclude_Jpsi->SetFrameBorderMode(0);
  pad_exclude_Jpsi->SetFrameFillStyle(0);
  pad_exclude_Jpsi->SetFrameBorderMode(0);

  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->SetMarkerColor(kBlack);
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->SetMarkerStyle(20);
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->SetXTitle("");
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->SetYTitle("");
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi->Draw("same");

  TH2D * h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_tmp = new TH2D( *h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi);
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_tmp->SetMarkerColor(kWhite);
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_tmp->SetMarkerStyle(20);
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_tmp->SetMarkerSize(1.0);
  h2D_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi_tmp->Draw("same");

  //**************************************************************************************
  //    !!!BEGIN: Placeholder for unblinding signal (exclude Jpsi version), after green light from pre-approval
  //**************************************************************************************
  /*
  TH2D* h2D_dimudimu_signal_2D_no_Jpsi = (TH2D*)w->data("ds_dimudimu_signal_2D_no_Jpsi")->createHistogram("m1,m2", m_bins, m_bins);
  h2D_dimudimu_signal_2D_no_Jpsi->SetMarkerColor(kBlack);
  h2D_dimudimu_signal_2D_no_Jpsi->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_no_Jpsi->SetMarkerSize(1.5);
  //Don't draw titles inheritted from dataset
  h2D_dimudimu_signal_2D_no_Jpsi->SetXTitle("");
  h2D_dimudimu_signal_2D_no_Jpsi->SetYTitle("");
  h2D_dimudimu_signal_2D_no_Jpsi->Draw("same");

  TH2D * h2D_dimudimu_signal_2D_no_Jpsi_tmp = new TH2D( *h2D_dimudimu_signal_2D_no_Jpsi);
  h2D_dimudimu_signal_2D_no_Jpsi_tmp->SetMarkerColor(kYellow);
  h2D_dimudimu_signal_2D_no_Jpsi_tmp->SetMarkerStyle(22);
  h2D_dimudimu_signal_2D_no_Jpsi_tmp->SetMarkerSize(1.0);
  h2D_dimudimu_signal_2D_no_Jpsi_tmp->Draw("same");
  */
  //************************************************************************************
  //    !!!END: Placeholder for unblinding signal, after green light from pre-approval
  //************************************************************************************

  corridorDn->Draw("C"); corridorUp->Draw("C"); txtHeader->Draw();
  c_template2D_m1_vs_m2_exclude_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_exclude_Jpsi.pdf");
  c_template2D_m1_vs_m2_exclude_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_exclude_Jpsi.png");
  c_template2D_m1_vs_m2_exclude_Jpsi->SaveAs("figures/DATA_and_Expected_2D_background_exclude_Jpsi.root");

  //************************************************************************************
  //         Validate the exclude Jpsi version with 2 dimu events at CR (iso & no-iso)
  //************************************************************************************
  //Validate m1 with non-iso data at CR
  TH1D *h1_control_offDiagonal_massC_data_no_Jpsi = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_2D_no_Jpsi")->createHistogram("m1", m_bins);
  h1_control_offDiagonal_massC_data_no_Jpsi->SetStats(0);
  h1_control_offDiagonal_massC_data_no_Jpsi->SetMarkerStyle(20);
  h1_control_offDiagonal_massC_data_no_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_offDiagonal_massC_data_no_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_offDiagonal_massC_data_no_Jpsi->GetYaxis()->SetRangeUser(0., 100.);

  TH1D *h1_control_offDiagonal_massC_template_no_Jpsi = new TH1D( *h2D_template2D_exclude_Jpsi_offDiagonal->ProjectionX() );
  h1_control_offDiagonal_massC_template_no_Jpsi->Scale( h1_control_offDiagonal_massC_data_no_Jpsi->Integral() / h1_control_offDiagonal_massC_template_no_Jpsi->Integral() );
  h1_control_offDiagonal_massC_template_no_Jpsi->SetLineColor(kRed);
  h1_control_offDiagonal_massC_template_no_Jpsi->SetLineWidth(2);
  h1_control_offDiagonal_massC_template_no_Jpsi->SetMarkerColor(kRed);

  TCanvas * c_control_offDiagonal_massC_no_Jpsi = new TCanvas("c_control_offDiagonal_massC_no_Jpsi", "c_control_offDiagonal_massC_no_Jpsi");
  c_control_offDiagonal_massC_no_Jpsi->cd();
  h1_control_offDiagonal_massC_data_no_Jpsi->Draw("e1");
  h1_control_offDiagonal_massC_template_no_Jpsi->Draw("HIST same");
  txtHeader->Draw();
  c_control_offDiagonal_massC_no_Jpsi->SaveAs("figures/Validation_m1_no_iso_CR_exclude_Jpsi.pdf");
  c_control_offDiagonal_massC_no_Jpsi->SaveAs("figures/Validation_m1_no_iso_CR_exclude_Jpsi.png");
  c_control_offDiagonal_massC_no_Jpsi->SaveAs("figures/Validation_m1_no_iso_CR_exclude_Jpsi.root");

  //Validate m1 with iso data at CR: very low stats, just FYI
  TH1D *h1_control_Iso_offDiagonal_massC_data_no_Jpsi = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi")->createHistogram("m1", m_bins);
  h1_control_Iso_offDiagonal_massC_data_no_Jpsi->SetStats(0);
  h1_control_Iso_offDiagonal_massC_data_no_Jpsi->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massC_data_no_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{1}} [GeV]");
  h1_control_Iso_offDiagonal_massC_data_no_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massC_data_no_Jpsi->GetYaxis()->SetRangeUser(0., 10.);

  TH1D *h1_control_Iso_offDiagonal_massC_template_no_Jpsi = new TH1D( *h2D_template2D_exclude_Jpsi_offDiagonal->ProjectionX() );
  h1_control_Iso_offDiagonal_massC_template_no_Jpsi->Scale( h1_control_Iso_offDiagonal_massC_data_no_Jpsi->Integral() / h1_control_Iso_offDiagonal_massC_template_no_Jpsi->Integral() );
  h1_control_Iso_offDiagonal_massC_template_no_Jpsi->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massC_template_no_Jpsi->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massC_template_no_Jpsi->SetMarkerColor(kRed);

  TCanvas * c_control_Iso_offDiagonal_massC_no_Jpsi = new TCanvas("c_control_Iso_offDiagonal_massC_no_Jpsi", "c_control_Iso_offDiagonal_massC_no_Jpsi");
  c_control_Iso_offDiagonal_massC_no_Jpsi->cd();
  h1_control_Iso_offDiagonal_massC_data_no_Jpsi->Draw("e1");
  h1_control_Iso_offDiagonal_massC_template_no_Jpsi->Draw("HIST same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massC_no_Jpsi->SaveAs("figures/Validation_m1_iso_CR_exclude_Jpsi.pdf");
  c_control_Iso_offDiagonal_massC_no_Jpsi->SaveAs("figures/Validation_m1_iso_CR_exclude_Jpsi.png");
  c_control_Iso_offDiagonal_massC_no_Jpsi->SaveAs("figures/Validation_m1_iso_CR_exclude_Jpsi.root");

  //Validate m2 with no-iso data at CR
  TH1D *h1_control_offDiagonal_massF_data_no_Jpsi = (TH1D*) w->data("ds_dimudimu_control_offDiagonal_2D_no_Jpsi")->createHistogram("m2", m_bins);
  h1_control_offDiagonal_massF_data_no_Jpsi->SetStats(0);
  h1_control_offDiagonal_massF_data_no_Jpsi->SetMarkerStyle(20);
  h1_control_offDiagonal_massF_data_no_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_offDiagonal_massF_data_no_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_offDiagonal_massF_data_no_Jpsi->GetYaxis()->SetRangeUser(0., 100.);

  TH1D *h1_control_offDiagonal_massF_template_no_Jpsi = new TH1D( *h2D_template2D_exclude_Jpsi_offDiagonal->ProjectionY() );
  h1_control_offDiagonal_massF_template_no_Jpsi->Scale( h1_control_offDiagonal_massF_data_no_Jpsi->Integral() / h1_control_offDiagonal_massF_template_no_Jpsi->Integral() );
  h1_control_offDiagonal_massF_template_no_Jpsi->SetLineColor(kRed);
  h1_control_offDiagonal_massF_template_no_Jpsi->SetLineWidth(2);
  h1_control_offDiagonal_massF_template_no_Jpsi->SetMarkerColor(kRed);

  TCanvas * c_control_offDiagonal_massF_no_Jpsi = new TCanvas("c_control_offDiagonal_massF_no_Jpsi", "c_control_offDiagonal_massF_no_Jpsi");
  c_control_offDiagonal_massF_no_Jpsi->cd();
  h1_control_offDiagonal_massF_data_no_Jpsi->Draw("e1");
  h1_control_offDiagonal_massF_template_no_Jpsi->Draw("HIST same");
  txtHeader->Draw();
  c_control_offDiagonal_massF_no_Jpsi->SaveAs("figures/Validation_m2_no_iso_CR_exclude_Jpsi.pdf");
  c_control_offDiagonal_massF_no_Jpsi->SaveAs("figures/Validation_m2_no_iso_CR_exclude_Jpsi.png");
  c_control_offDiagonal_massF_no_Jpsi->SaveAs("figures/Validation_m2_no_iso_CR_exclude_Jpsi.root");

  //Validate m2 with iso data at CR: limited stats, just FYI
  TH1D *h1_control_Iso_offDiagonal_massF_data_no_Jpsi = (TH1D*) w->data("ds_dimudimu_control_Iso_offDiagonal_2D_no_Jpsi")->createHistogram("m2", m_bins);
  h1_control_Iso_offDiagonal_massF_data_no_Jpsi->SetStats(0);
  h1_control_Iso_offDiagonal_massF_data_no_Jpsi->SetMarkerStyle(20);
  h1_control_Iso_offDiagonal_massF_data_no_Jpsi->GetXaxis()->SetTitle("m_{#mu#mu_{2}} [GeV]");
  h1_control_Iso_offDiagonal_massF_data_no_Jpsi->GetYaxis()->SetTitle("Events/0.04 GeV");
  h1_control_Iso_offDiagonal_massF_data_no_Jpsi->GetYaxis()->SetRangeUser(0.,10.);

  TH1D *h1_control_Iso_offDiagonal_massF_template_no_Jpsi = new TH1D( *h2D_template2D_exclude_Jpsi_offDiagonal->ProjectionY() );
  h1_control_Iso_offDiagonal_massF_template_no_Jpsi->Scale( h1_control_Iso_offDiagonal_massF_data_no_Jpsi->Integral() / h1_control_Iso_offDiagonal_massF_template_no_Jpsi->Integral() );
  h1_control_Iso_offDiagonal_massF_template_no_Jpsi->SetLineColor(kRed);
  h1_control_Iso_offDiagonal_massF_template_no_Jpsi->SetLineWidth(2);
  h1_control_Iso_offDiagonal_massF_template_no_Jpsi->SetMarkerColor(kRed);

  TCanvas * c_control_Iso_offDiagonal_massF_no_Jpsi = new TCanvas("c_control_Iso_offDiagonal_massF_no_Jpsi", "c_control_Iso_offDiagonal_massF_no_Jpsi");
  c_control_Iso_offDiagonal_massF_no_Jpsi->cd();
  h1_control_Iso_offDiagonal_massF_data_no_Jpsi->Draw("e1");
  h1_control_Iso_offDiagonal_massF_template_no_Jpsi->Draw("HIST same");
  txtHeader->Draw();
  c_control_Iso_offDiagonal_massF_no_Jpsi->SaveAs("figures/Validation_m2_iso_CR_exclude_Jpsi.pdf");
  c_control_Iso_offDiagonal_massF_no_Jpsi->SaveAs("figures/Validation_m2_iso_CR_exclude_Jpsi.png");
  c_control_Iso_offDiagonal_massF_no_Jpsi->SaveAs("figures/Validation_m2_iso_CR_exclude_Jpsi.root");
}
