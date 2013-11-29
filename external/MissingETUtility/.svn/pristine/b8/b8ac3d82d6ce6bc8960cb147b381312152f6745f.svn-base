runExample() {

  gSystem->AddIncludePath("-I/usera/khoo/atlas/MetUtil_dev/testarea/Reconstruction/Jet/JetUncertainties");
  gSystem->AddIncludePath("-I/usera/khoo/atlas/MetUtil_dev/testarea/Reconstruction/Jet/JetResolution");
  gSystem->AddIncludePath("-I/usera/khoo/atlas/MetUtil_dev/testarea/Reconstruction/egamma/egammaAnalysis/egammaAnalysisUtils");
  gSystem->AddIncludePath("-I/usera/khoo/atlas/MetUtil_dev/testarea/PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections");

  gSystem->AddLinkedLibs("-L/usera/khoo/atlas/MetUtil_dev/testarea/InstallArea/x86_64-slc5-gcc43-opt/lib");
  gSystem->AddLinkedLibs("-lJetUncertainties");
  gSystem->AddLinkedLibs("-lJetResolution");
  gSystem->AddLinkedLibs("-lMuonMomentumCorrections");
  gSystem->AddLinkedLibs("-legammaAnalysisUtils");
  gSystem->AddLinkedLibs("-lMissingETUtilityLib");

  gROOT->ProcessLine(".L EventReader.C+");

//   TFile* rootfile = TFile::Open("mc12_8TeV.105200.McAtNloJimmy_CT10_ttbar/NTUP_JETMET.00955143._000216.root.1");
  TFile* rootfile = TFile::Open("/r02/atlas/khoo/sample_D3PD/non_SUSY/mc12_8TeV.147807.PowhegPythia8_AU2CT10_Zmumu.merge.NTUP_JETMET.e1169_s1469_s1470_r3542_r3549_p1286_tid01060683_00/NTUP_JETMET.01060683._000963.root.2");
  //  TFile* rootfile = TFile::Open("/r02/atlas/khoo/sample_D3PD/non_SUSY/user.sresconi.mc12_8TeV.147803.PowhegPythia8_AU2CT10_Wminenu.merge.AOD.e1169_s1469_s1470_r3542_r3549.R172test.JETMET.v6_EXT0.120607151801/user.sresconi.006202._00156.qcd.root");
  TTree* qcd = (TTree*) rootfile->Get("qcd");
  qcd->Process("Example.C+");

//   TFile* rootfile = TFile::Open("SUSY_ttbar.root");
//   TTree* susy = (TTree*) rootfile->Get("susy");
//   susy->Process("Example.C+");
}
