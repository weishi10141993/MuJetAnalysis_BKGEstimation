#include <iostream>
#include <fstream>
#include <algorithm>    // std::max
#include <stdlib.h>
#include <iomanip>
#include "TSystemDirectory.h"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!  USER Configure Below        !
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

const int year = 2018;//Configure which year ntuples to run, options: 2017, 2018

//CB fitted mean mass [GeV] for prompt signals
double mean[11] = {0.2560, 0.4012, 0.7003, 1.0000, 1.9990, 4.9980, 8.4920, 14.990, 24.980, 34.970, 57.930};

//80% signal eff needed mass window size for prompt signals
//double window[11] = {0.030922364, 0.0201698988, 0.02337757474, 0.0275968, 0.04615657165, 0.124742502, 0.2110350916, 0.3849085576, 0.67586519, 0.914118934, 1.860449598};

//85% signal eff needed mass window size for prompt signals
//double window[11] = {0.0359823872, 0.0252123735, 0.02787326219, 0.0326656, 0.056518251, 0.1538490858, 0.2586881768, 0.46940068, 0.849659096, 1.102319891, 2.524895883};

//90% signal eff needed mass window size from above prompt signals
double window[11] = {0.0438535344, 0.0336164980, 0.0359654996, 0.0428032000, 0.0753576680, 0.2037460866, 0.3539943472, 0.7041010200, 1.1972469080, 1.6131510600, 3.9866777100};

//95% signal eff needed mass window size for prompt signals
//double window[11] = {0.0629691776, 0.05176940692, 0.05574652438, 0.0664576, 0.1224562105, 0.3617532558, 0.748834196, 1.40820204, 2.70346076, 2.6885851, 9.523730085};
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!  USER Configure Above        !
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!  Initialize variables for macros  !
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Float_t luminosity;//Luminosity of a year
TString header = "";//header of plots

TString inputFile1;//Run 2 data Ntuples used for HighMassBKGABCD and LowMassBKGFit1D
TString outFileLM = "";//output file of LowMassBKGFit1D, will also be used in limit calculation
TString inputFile2;//input workSpace file to LowMassBKGPlot2D, that's also the output file of LowMassBKGFit1D
Float_t validate_m1_iso_Ymax;//Max y in validation plots
Float_t validate_m2_iso_Ymax;

TString store;//Input directory for HighMassBKGShape
TString outFileHMShape = "";//output file of HighMassBKGShape
Float_t MC_ScaleFactors[7];//MC-Data SF for HighMassBKGShape
TString MC_files[7];//Input MC shape for HighMassBKGShape
TString DATA_files[1];//Input DATA shape for HighMassBKGShape

TString outFileHMABCD = "";//output file of HighMassBKGABCD

namespace BKG_cfg {

  inline void ConfigureInput( const int year ) {

    std::cout << "\nConfiguring inputs for year " << year << std::endl;

    if (year == 2017) {
      luminosity = lumi_2017;
      inputFile1 = "Input_2017CDEF_2SAmu_NoVtxProbCut.txt";
      inputFile2 = "ws_2017_FINAL.root";
      store      = store_2017;
      validate_m1_iso_Ymax    = validate_m1_iso_Ymax_2017;
      validate_m2_iso_Ymax    = validate_m2_iso_Ymax_2017;
      for(int i = 0; i < 7; i++){
        MC_ScaleFactors[i] = MC_ScaleFactors_2017[i];
        MC_files[i]        = MC_files_2017[i];
      }
      for(int j = 0; j < 1; j++){
        DATA_files[j] = DATA_files_2017[j];
      }
    }//end 2017
    else if (year == 2018) {
      luminosity = lumi_2018;
      inputFile1 = "Input_2018ABCD_2SAmu_NoVtxProbCut_4HLT.txt";
      inputFile2 = "ws_2018_FINAL.root";
      store      = store_2018;
      validate_m1_iso_Ymax    = validate_m1_iso_Ymax_2018;
      validate_m2_iso_Ymax    = validate_m2_iso_Ymax_2018;
      for(int i = 0; i < 7; i++){
        MC_ScaleFactors[i] = MC_ScaleFactors_2018[i];
        MC_files[i]        = MC_files_2018[i];
      }
      for(int j = 0; j < 1; j++){
        DATA_files[j] = DATA_files_2018[j];
      }
    }//end 2018
    else{
      std::cout << "*** User input year is unknown! Please check. ***" << std::endl;
    }

    header         = header + "#bf{CMS} #it{Preliminary}    " + Form("%.2f", luminosity) + "fb^{-1} (13 TeV)";
    outFileLM      = outFileLM + "ws_" + Form("%d", year) + "_FINAL.root";
    outFileHMABCD  = outFileHMABCD + "ABCD_" + Form("%d", year) + "_FINAL.root";
    outFileHMShape = outFileHMShape + "HighMassBKGShape_" + Form("%d", year) + "_FINAL.root";

  } // End ConfigureInput function

  double My_MassWindow(double m1, double m2) {
    //return this mass window size given m1 and m2
    double mysize = 0.0;

    //Interpolation by drawing staight line, we use (m1+m2)/2 here
    //mass window size = y1 + (x-x1)*(y2-y1)/(x2-x1), x = (m1+m2)/2
    //Start and end with 0.2113 and 60GeV to match bkg analysis
    if ( (m1+m2)/2 >= 0.2113   && (m1+m2)/2 < mean[1] ) mysize = window[0] + ( (m1+m2)/2 - mean[0] )*( window[1]  - window[0] )/( mean[1]  - mean[0] );
    if ( (m1+m2)/2 >= mean[1]  && (m1+m2)/2 < mean[2] ) mysize = window[1] + ( (m1+m2)/2 - mean[1] )*( window[2]  - window[1] )/( mean[2]  - mean[1] );
    if ( (m1+m2)/2 >= mean[2]  && (m1+m2)/2 < mean[3] ) mysize = window[2] + ( (m1+m2)/2 - mean[2] )*( window[3]  - window[2] )/( mean[3]  - mean[2] );
    if ( (m1+m2)/2 >= mean[3]  && (m1+m2)/2 < mean[4] ) mysize = window[3] + ( (m1+m2)/2 - mean[3] )*( window[4]  - window[3] )/( mean[4]  - mean[3] );
    if ( (m1+m2)/2 >= mean[4]  && (m1+m2)/2 < mean[5] ) mysize = window[4] + ( (m1+m2)/2 - mean[4] )*( window[5]  - window[4] )/( mean[5]  - mean[4] );
    if ( (m1+m2)/2 >= mean[5]  && (m1+m2)/2 < mean[6] ) mysize = window[5] + ( (m1+m2)/2 - mean[5] )*( window[6]  - window[5] )/( mean[6]  - mean[5] );
    if ( (m1+m2)/2 >= mean[6]  && (m1+m2)/2 < mean[7] ) mysize = window[6] + ( (m1+m2)/2 - mean[6] )*( window[7]  - window[6] )/( mean[7]  - mean[6] );
    if ( (m1+m2)/2 >= mean[7]  && (m1+m2)/2 < mean[8] ) mysize = window[7] + ( (m1+m2)/2 - mean[7] )*( window[8]  - window[7] )/( mean[8]  - mean[7] );
    if ( (m1+m2)/2 >= mean[8]  && (m1+m2)/2 < mean[9] ) mysize = window[8] + ( (m1+m2)/2 - mean[8] )*( window[9]  - window[8] )/( mean[9]  - mean[8] );
    if ( (m1+m2)/2 >= mean[9]  && (m1+m2)/2 < 60.000  ) mysize = window[9] + ( (m1+m2)/2 - mean[9] )*( window[10] - window[9] )/( mean[10] - mean[9] );

    return mysize;
  }

} // End namespace
