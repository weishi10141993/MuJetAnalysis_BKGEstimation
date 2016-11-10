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

void PlotSignal_and_Background_ext() {

  setTDRStyle();

  TCanvas * c_template2D_m1_vs_m2 = new TCanvas("c_template2D_m1_vs_m2", "c_template2D_m1_vs_m2",0,1320,1044,928);
  c_template2D_m1_vs_m2->SetCanvasSize(1040,900);
  c_template2D_m1_vs_m2->SetLeftMargin(0.121);
  c_template2D_m1_vs_m2->SetRightMargin(0.148);
  c_template2D_m1_vs_m2->SetTopMargin(0.05);
  c_template2D_m1_vs_m2->cd();
  c_template2D_m1_vs_m2->SetLogz();

  TFile* file = new TFile("ws.root");
  RooWorkspace *w = (RooWorkspace*) file->Get("w");

  const double       m_min  = 0.2113;
  const double       m_max  = 9.;
  const unsigned int m_bins = 220;

  // Diagonal region |m1 - m2| < 5 sigma = kA + kB * (m1 + m2)/2
  const double kA = 0.13;
  const double kB = 0.065;

  double nEvents_Jpsi = 0.064; 

  //****************************************************************************
  //                         Draw 2D template m1 x m2                           
  //****************************************************************************

  RooGenericPdf myJpsi_2D = w->pdf("Jpsi_2D");
  RooGenericPdf mytemplate2D = w->pdf("template2D");
  RooGenericPdf myJpsiAll("myJpsiAll","myJpsiAll","myJpsi_2D+mytemplate2D",RooArgList(myJpsi_2D,mytemplate2D));
  RooRealVar m1("m1","m_{#mu#mu_{1}}",m_min,m_max,"GeV/#it{c}^{2}");
  RooRealVar m2("m2","m_{#mu#mu_{2}}",m_min,m_max,"GeV/#it{c}^{2}");
  m1.setRange("signal",-2.7,3.5);
  m2.setRange("signal",-2.7,3.5);
  RooAbsReal *myintegral = myJpsiAll->createIntegral(m1,m2,NormSet(m1,m2),Range("signal"));

}
