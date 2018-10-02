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
#include "RooAbsReal.h"
#include "RooRealProxy.h"
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
#include "RooCBShape.h"

#ifndef __CINT__
#include "RooCFunction1Binding.h"
#endif

using namespace RooFit;

class RooUserPdf : public RooAbsPdf {

  public:
    RooUserPdf(const char *name, const char *title, RooAbsReal& _x1, RooAbsReal& _x2);
    RooUserPdf(const RooUserPdf& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const {return new RooUserPdf(*this,newname);}
    inline virtual ~RooUserPdf() { }

  protected:
    RooRealProxy x1;
    RooRealProxy x2;
    Double_t evaluate() const ;

  private:
    ClassDef(RooUserPdf,0)
};

RooUserPdf::RooUserPdf(const char *name, const char *title, RooAbsReal& _x1, RooAbsReal& _x2) :
  RooAbsPdf(name,title),
  x1("x1","Dependent1",this,_x1),
  x2("x2","Dependent2",this,_x2)
{}

RooUserPdf::RooUserPdf(const RooUserPdf& other, const char* name) :
  RooAbsPdf(other,name),
  x1("x1", this, other.x1),
  x2("x2", this, other.x2)
{}

Double_t RooUserPdf::evaluate() const
{
  Double_t res;

  if ( fabs(x1 - x2) < 5.*(0.026 + 0.013*(x1 + x2)/2.) ) {
    res = 1.0;
  } else {
    res = 0.0;
  }
  return res;
}

void CheckCorredor() {

  double JpsiF_mean[3]  = {3.0929,     3.09407,  3.09173}; // Mean value from fit, 1sigma up, 1sigma down
  double JpsiF_sigma[3] = {4.0354e-02, 0.041504, 0.039204};
  double JpsiF_alpha[3] = {1.6635,     1.7815,   1.5455};
  const double m_min  = 2.;
  const double m_max  = 6.;

  auto c1 = new TCanvas();
  //(double x, double alpha, double n, double sigma, double mean)
  auto f1  = new TF1("f1","ROOT::Math::crystalball_function(x, 1.6635, 2, 0.040354, 3.0929)",m_min,m_max);
  f1->SetLineColor(kRed);
  f1->Draw();
  auto f2 = new TF1("f2","ROOT::Math::crystalball_function(x, 1.7815, 2, 0.041504, 3.09407)",m_min,m_max);
  f2->SetLineColor(kGreen); 
  f2->Draw("same");
  auto f3 = new TF1("f3","ROOT::Math::crystalball_function(x, 1.5455, 2, 0.039204, 3.09173)",m_min,m_max);
  f3->SetLineColor(kBlue);
  f3->Draw("same");

  auto legend = new TLegend(0.7,0.6,0.9,1.);
  legend->AddEntry(f1,"STD","L");
  legend->AddEntry(f2,"1sigma UP","L");
  legend->AddEntry(f3,"1sigma DOWN","L");
  legend->Draw();

  double int_min = JpsiF_mean[0] - (0.13+0.065*JpsiF_mean[0]);
  double int_max = JpsiF_mean[0] + (0.13+0.065*JpsiF_mean[0]);
  std::cout<<"Integral of withon 5 sigma of the STD function is: "<<f1->Integral(int_min,int_max)/f3->Integral(-999.,999)<<std::endl;
  std::cout<<"Integral of withon 5 sigma of the  UP function is: "<<f2->Integral(int_min,int_max)/f3->Integral(-999.,999)<<std::endl;
  std::cout<<"Integral of withon 5 sigma of the DOW function is: "<<f3->Integral(int_min,int_max)/f3->Integral(-999.,999)<<std::endl;

}
