{
  gROOT->LoadMacro("tdrStyle.C");
  setTDRStyle();
  
  gROOT->LoadMacro("PlotSignal_and_Background.C");
  PlotSignal_and_Background();
}
