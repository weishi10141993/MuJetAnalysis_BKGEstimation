//*****************************************************************************************************
//* For estimating the bkg yield at high mass bins based on an ABCD method and check its stability    *
//* Use: . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh  *
//* To request more time and memory: sintr -t 480 -m 102400                                           *
//*                                       Wei Shi @Sep 25, 2019, Rice U.                              *
//*****************************************************************************************************
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>
#include "stdio.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TChain.h"
#include "TBranch.h"

#include "Constants.h"//local const file
#include "Config.h"

void HighMassBKGABCD18() {

  //Configure inputs for year
  BKG_cfg::ConfigureInput(year);

  TChain chain_data_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream Myfile(inputFile1);
  std::string Line;
  if( !Myfile ) std::cout<<"ERROR opening Myfile."<<std::endl;
  while (std::getline(Myfile, Line)){
    TString Line2(Line);
    if( Line2.Contains("root") ){
      chain_data_dimudimu.Add(Line2.Data());
    }
  }

  double VarX[29] = {1000, 60, 50, 40, 30, 20, 10, 5,  60, 50, 40, 30, 20, 10, 5,  60, 50, 40, 30, 20, 10, 5,  60, 50, 40, 30, 20, 10, 5};//Mass inequality
  double VarY[29] = {1000, 80, 80, 80, 80, 80, 80, 80, 60, 60, 60, 60, 60, 60, 60, 40, 40, 40, 40, 40, 40, 40, 20, 20, 20, 20, 20, 20, 20};//Max(Iso1, Iso2)

  double A[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//Estimated BKG yield in region A ~ B*D/C
  double ErrA[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//Error on A ~ |A|*sqrt(1/B+1/C+1/D), assuming poisson stats in BCD
  double B[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//BKG in B
  double C[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//BKG in C
  double D[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//BKG in D

  //Print Table
  cout<<"********************************************************************"<<endl;
  cout<<"* ( VarX, VarY ) *    B    *    C    *    D    *   <A>   *   ErrA  *"<<endl;
  for (int i = 0; i < 29; i++) {
    ostringstream stream_cut_B;
    stream_cut_B << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && recoRePaired2muleading_m < 76 && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && TMath::Max(diMuonC_IsoTk_FittedVtx,diMuonF_IsoTk_FittedVtx)<"<< VarY[i] <<" && TMath::Max(diMuonC_IsoTk_FittedVtx,diMuonF_IsoTk_FittedVtx)>="<< iso_cut <<" && fabs(diMuonC_FittedVtx_m - diMuonF_FittedVtx_m) - 5*( -0.00591865*(diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0 + 0.00113991*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 2) - 2.62048e-05*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 3) + 1.92254e-07*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 4) )<0.236369 && diMuonC_FittedVtx_m>11. && diMuonC_FittedVtx_m<59. && diMuonF_FittedVtx_m>11. && diMuonF_FittedVtx_m<59.";
    TString cut_B = stream_cut_B.str();

    ostringstream stream_cut_C;
    stream_cut_C << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && recoRePaired2muleading_m < 76 && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && TMath::Max(diMuonC_IsoTk_FittedVtx,diMuonF_IsoTk_FittedVtx)<"<< VarY[i] <<" && TMath::Max(diMuonC_IsoTk_FittedVtx,diMuonF_IsoTk_FittedVtx)>="<< iso_cut <<" && fabs(diMuonC_FittedVtx_m - diMuonF_FittedVtx_m) - 5*( -0.00591865*(diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0 + 0.00113991*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 2) - 2.62048e-05*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 3) + 1.92254e-07*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 4) )<"<< VarX[i] <<" && fabs(diMuonC_FittedVtx_m - diMuonF_FittedVtx_m) - 5*( -0.00591865*(diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0 + 0.00113991*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 2) - 2.62048e-05*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 3) + 1.92254e-07*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 4) )>=0.236369 && diMuonC_FittedVtx_m>11. && diMuonC_FittedVtx_m<59. && diMuonF_FittedVtx_m>11. && diMuonF_FittedVtx_m<59.";
    TString cut_C = stream_cut_C.str();

    ostringstream stream_cut_D;
    stream_cut_D << "is1SelMuHighPt && is2SelMuHighPt && is3SelMuLowPt && is4SelMuLowPt && isVertexOK && is2DiMuons && nSAMu <= 1 && diMuonC_FittedVtx_prob > 0.2*(1 - dimuC_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonC_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonC_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonC_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonC_FittedVtx_dR), 4) )*pow(fabs(diMuonC_FittedVtx_Lxy/10.0), 2.0) ) && diMuonF_FittedVtx_prob > 0.2*(1 - dimuF_nSAMu)*exp( -( 8.53647 - 50.4571*(sqrt(diMuonF_FittedVtx_dR)) + 109.83*pow(sqrt(diMuonF_FittedVtx_dR), 2) - 92.7445*pow(sqrt(diMuonF_FittedVtx_dR), 3) + 36.8351*pow(sqrt(diMuonF_FittedVtx_dR), 4) )*pow(fabs(diMuonF_FittedVtx_Lxy/10.0), 2.0) ) && ( nSAMu == 0 || ( nSAMu == 1 && ( diMuonC_FittedVtx_Lxy > 0.1 || diMuonF_FittedVtx_Lxy > 0.1 ) ) ) && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && recoRePaired2muleading_m < 76 && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && isSignalHLTFired && TMath::Max(diMuonC_IsoTk_FittedVtx,diMuonF_IsoTk_FittedVtx)<"<< iso_cut <<" && fabs(diMuonC_FittedVtx_m - diMuonF_FittedVtx_m) - 5*( -0.00591865*(diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0 + 0.00113991*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 2) - 2.62048e-05*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 3) + 1.92254e-07*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 4) )<"<< VarX[i] <<" && fabs(diMuonC_FittedVtx_m - diMuonF_FittedVtx_m) - 5*( -0.00591865*(diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0 + 0.00113991*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 2) - 2.62048e-05*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 3) + 1.92254e-07*pow((diMuonC_FittedVtx_m + diMuonF_FittedVtx_m)/2.0, 4) )>=0.236369 && diMuonC_FittedVtx_m>11. && diMuonC_FittedVtx_m<59. && diMuonF_FittedVtx_m>11. && diMuonF_FittedVtx_m<59.";
    TString cut_D = stream_cut_D.str();

    TTree* tree_dimudimu_B = chain_data_dimudimu.CopyTree(cut_B);
    TTree* tree_dimudimu_C = chain_data_dimudimu.CopyTree(cut_C);
    TTree* tree_dimudimu_D = chain_data_dimudimu.CopyTree(cut_D);

    //Get Entries
    B[i] = tree_dimudimu_B->GetEntries();
    C[i] = tree_dimudimu_C->GetEntries();
    D[i] = tree_dimudimu_D->GetEntries();
    A[i] = B[i] * D[i] / C[i];
    ErrA[i] = A[i] * sqrt( 1/B[i] + 1/C[i] + 1/D[i]);

    cout<<"* ("<<setw(5)<<setprecision(2)<<VarX[i]<<", "<<setw(5)<<setprecision(2)<<VarY[i]<<") * "<<setw(7)<<setprecision(3)<<B[i]<<" * "<<setw(7)<<setprecision(3)<<C[i]<<" * "<<setw(7)<<setprecision(3)<<D[i]<<" * "<<setw(7)<<setprecision(3)<<A[i]<<" * "<<setw(7)<<setprecision(2)<<ErrA[i]<<" *"<<endl;

  }//end for

  cout<<"********************************************************************"<<endl;

  //Draw <A> ErrA in one plot
  TFile myPlot(outFileHMABCD, "RECREATE");

  TCanvas *C1=new TCanvas("C1", "C1", 700, 500);
  C1->cd();
  TH1F *Yield = new TH1F();
  for(unsigned int iB=1; iB<=29; iB++){
    Yield->SetBinContent(iB, A[iB-1] );
    Yield->SetBinError(iB, ErrA[iB-1]);
    Yield->GetXaxis()->SetBinLabel(iB, Form("(%.0f, %.0f)", VarX[iB-1], VarY[iB-1]) );
  }
  Yield->SetMarkerStyle(20);
  Yield->GetXaxis()->SetTitle("(Max VarX, Max VarY)");
  Yield->GetYaxis()->SetTitle("Yield");
  Yield->SetTitle("Estimated Background Yield in Signal Region A");
  Yield->SetStats(0);
  Yield->Draw();
  C1->Write();

  myPlot.Close();
}
