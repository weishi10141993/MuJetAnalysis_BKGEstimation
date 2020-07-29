//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!  USER Configure Below        !
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

const int year = 2018;//Configure which year ntuples to run, options: 2017, 2018

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
TString DATA_files[4];//Input DATA shape for HighMassBKGShape

TString outFileHMABCD = "";//output file of HighMassBKGABCD

namespace BKG_cfg {

  inline void ConfigureInput( const int year ) {

    std::cout << "\nConfiguring inputs for year " << year << std::endl;

    if(year == 2017){
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
      for(int j = 0; j < 4; j++){
        DATA_files[j] = DATA_files_2017[j];
      }
    }//end 2017
    else if(year == 2018){
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
      for(int j = 0; j < 4; j++){
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

  } // End function

} // End namespace
