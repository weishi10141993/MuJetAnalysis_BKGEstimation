//*****************************************************************************************************
//* For estimating the bkg yield at high mass bins based on an ABCD method and check its stability    *
//* Use: . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh  *
//*                                       Wei Shi @Sep 25, 2019, Rice U.                              *
//*****************************************************************************************************
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>
#include "stdio.h"
#include "math.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TChain.h"
#include "TBranch.h"

void HighMassBKGABCD() {

  double VarX[29] = {1000, 60, 50, 40, 30, 20, 10, 5,  60, 50, 40, 30, 20, 10, 5,  60, 50, 40, 30, 20, 10, 5,  60, 50, 40, 30, 20, 10, 5};//Mass inequality
  double VarY[29] = {1000, 80, 80, 80, 80, 80, 80, 80, 60, 60, 60, 60, 60, 60, 60, 40, 40, 40, 40, 40, 40, 40, 20, 20, 20, 20, 20, 20, 20};//Max(Iso1, Iso2)

  double A[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//Estimated BKG yield in region A ~ B*D/C
  double ErrA[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//Error on A ~ |A|*sqrt(1/B+1/C+1/D), assuming poisson stats in BCD
  double B[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//BKG in B
  double C[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//BKG in C
  double D[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};//BKG in D

  //Input file
  TChain chain_data_dimudimu("cutFlowAnalyzerPXBL4PXFL3/Events");
  std::ifstream Myfile( "Input_2017CDEF.txt" );
  std::string Line;
  if( !Myfile ) std::cout<<"ERROR opening Myfile."<<std::endl;
  while (std::getline(Myfile, Line)){
    TString Line2(Line);
    if( Line2.Contains("root") ){
      chain_data_dimudimu.Add(Line2.Data());
    }
  }

  for (int i = 0; i < 29; i++) {
    ostringstream stream_cut_B;
    stream_cut_B << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (diMuonC_m1_FittedVtx_hitpix_Phase1==1 || diMuonC_m2_FittedVtx_hitpix_Phase1) && (diMuonF_m1_FittedVtx_hitpix_Phase1==1 || diMuonF_m2_FittedVtx_hitpix_Phase1==1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)<"<< VarY[i] <<" && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=1.5 && TMath::Min(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=0 && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0<0.009132 && massC>11. && massC<60. && massF>11. && massF<60.";
    TString cut_B = stream_cut_B.str();

    ostringstream stream_cut_C;
    stream_cut_C << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (diMuonC_m1_FittedVtx_hitpix_Phase1==1 || diMuonC_m2_FittedVtx_hitpix_Phase1) && (diMuonF_m1_FittedVtx_hitpix_Phase1==1 || diMuonF_m2_FittedVtx_hitpix_Phase1==1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)<"<< VarY[i] <<" && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=1.5 && TMath::Min(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=0 && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0<"<< VarX[i] <<" && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0>=0.009132 && massC>11. && massC<60. && massF>11. && massF<60.";
    TString cut_C = stream_cut_C.str();

    ostringstream stream_cut_D;
    stream_cut_D << "is1SelMu17 && is2SelMu8 && is3SelMu8 && is4SelMu8 && isVertexOK && is2DiMuons && (diMuonC_m1_FittedVtx_hitpix_Phase1==1 || diMuonC_m2_FittedVtx_hitpix_Phase1) && (diMuonF_m1_FittedVtx_hitpix_Phase1==1 || diMuonF_m2_FittedVtx_hitpix_Phase1==1) && fabs(diMuons_dz_FittedVtx) < 0.1 && isSignalHLTFired && TMath::Max(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)<1.5 && TMath::Min(diMuonCMu0_IsoTk0p3_FittedVtx,diMuonFMu0_IsoTk0p3_FittedVtx)>=0 && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0<"<< VarX[i] <<" && fabs(massC-massF)-3*0.007025*(massC+massF)/2.0-3*0.000053*(massC+massF)*(massC+massF)/4.0>=0.009132 && massC>11. && massC<60. && massF>11. && massF<60.";
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

  }//end for

  //Print Table
  cout<<"*********************************************************************"<<endl;
  cout<<"* (VarX, VarY) *    B    *    C    *    D    *   <A>   *    ErrA    *"<<endl;
  for (int i = 0; i < 29; i++) {
    cout<<"* ("<<setw(4)<<setprecision(2)<<VarX[i]<<", "<<setw(4)<<setprecision(2)<<VarY[i]<<") * "<<setw(9)<<setprecision(3)<<B[i]<<" * "<<setw(9)<<setprecision(3)<<C[i]<<" * "<<setw(8)<<setprecision(3)<<D[i]<<" * "<<setw(8)<<setprecision(3)<<A[i]<<" * "<<setw(13)<<setprecision(3)<<ErrA[i]<<" *"<<endl;
  }
  cout<<"*********************************************************************"<<endl;

}
