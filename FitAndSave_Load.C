{
  gROOT->LoadMacro("tdrStyle.C");
  setTDRStyle();
  
  gROOT->LoadMacro("FitAndSave.C");
  FitAndSave();
}
