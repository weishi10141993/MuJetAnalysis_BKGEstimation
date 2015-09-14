{
  gROOT->LoadMacro("tdrStyle.C");
  setTDRStyle();
  
  gROOT->LoadMacro("PlotSignalProjection.C");
  PlotSignalProjection();
}
