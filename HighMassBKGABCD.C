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

void HighMassBKGABCD() {

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
    stream_cut_B << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)<"<< VarY[i] <<" && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=1.5 && TMath::Min(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=0 && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0<0.009132 && massC>11. && massC<59. && massF>11. && massF<59.";
    TString cut_B = stream_cut_B.str();

    ostringstream stream_cut_C;
    stream_cut_C << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)<"<< VarY[i] <<" && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=1.5 && TMath::Min(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=0 && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0<"<< VarX[i] <<" && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0>=0.009132 && massC>11. && massC<59. && massF>11. && massF<59.";
    TString cut_C = stream_cut_C.str();

    ostringstream stream_cut_D;
    stream_cut_D << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (recoRePaired2mutrailing_dR >= 0.2 || recoRePaired2mutrailing_m >= 3) && (diMuonC_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonC_m2_FittedVtx_hitpix_Phase1 == 1) && (diMuonF_m1_FittedVtx_hitpix_Phase1 == 1 || diMuonF_m2_FittedVtx_hitpix_Phase1 == 1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)<1.5 && TMath::Min(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=0 && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0<"<< VarX[i] <<" && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0>=0.009132 && massC>11. && massC<59. && massF>11. && massF<59.";
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
  Yield->SetTitle("Background yield in A");
  Yield->SetStats(0);
  Yield->Draw();
  C1->Write();

  myPlot.Close();
}
