

#include "StandardRootIncludes.h"
#include "standardsetup.C"

#include "Tools.C"
void CheckHistos(TString TaggerName = "SV0", TString OPName = "SV050")
{
  

  TFile* file_rel15 = new TFile("EffMapInput.root", "read");
  TFile *file_rel16 = new TFile("ttbar7.root", "read");

  // retrieve scale factor inputs resulting from ptrel fits

  // TFile* fin = new TFile("pTrel.root");
  // TF1* fbSelectionSF = (TF1*) fin->Get("SFFunction");
  // fbSFstatin = (TF1*) fin->Get("StatErrFunction");
  // fbSFsystin = (TF1*) fin->Get("SystErrFunction");

  // retrieve the histograms
  TH2* hb = 0;
  TH2* hc = 0;
  TH2* hlight = 0;

  // rel15:
  // hb = (TH2*) file_rel15->Get("hMCTaggingEff_b_SV0_5_72");
  hc = (TH2*) file_rel15->Get("hMCTaggingEff_c_SV0_5_72");
  // hlight = (TH2*) file_rel15->Get("hMCTaggingEff_light_SV0_5_72");
 
  // rel.16:
  if (OPName == "SV050") {
    TString dirname = "btag/AntiKt4TopoEM/" + TaggerName + "/";
    TString name = "h_btag_AntiKt4TopoEM_SV0_beff_eta_Et_cut1";
    hb = GetHisto2DFromFile(file_rel16, dirname, name);
    name = "h_btag_AntiKt4TopoEM_SV0_ceff_eta_Et_cut1";
    hc = GetHisto2DFromFile(file_rel16, dirname, name);
    dirname = "btag/AntiKt4TopoEM/" + TaggerName + "/";
    name = "h_btag_AntiKt4TopoEM_SV0_udseff_eta_Et_cut1";
    hlight = GetHisto2DFromFile(file_rel16, dirname, name);
    if (!hb) {
      cerr << "Error getting hb!" << endl;
      return;
    }
    if (!hc) {
      cerr << "Error getting hc!" << endl;
      return;
    }
    if (!hlight) {
      cerr << "Error getting hlight!" << endl;
      return;
    }
  } else if (OPName == "JetProb50") {
    hb = (TH2*) file_rel15->Get("btag/AntiKt4TopoEM/" + TaggerName + "h_btag_AntiKt4TopoEM_JetProb_beff_eta_Et_cut1");
    // hc = (TH2*) file_rel15->Get("");
    hlight = (TH2*) file_rel15->Get("btag/AntiKt4TopoEM/" + TaggerName + "h_btag_AntiKt4TopoEM_JetProb_udseff_eta_Et_cut1");
  } else if (OPName == "JetProb70") {
    hb = (TH2*) file_rel15->Get("btag/AntiKt4TopoEM/" + TaggerName + "h_btag_AntiKt4TopoEM_JetProb_beff_eta_Et_cut2");
    // hc = (TH2*) file_rel15->Get("");
    hlight = (TH2*) file_rel15->Get("btag/AntiKt4TopoEM/" + TaggerName + "h_btag_AntiKt4TopoEM_JetProb_udseff_eta_Et_cut2");
   if (!hb) {
      cerr << "Error getting hb!" << endl;
      return;
    }
    if (!hc) {
      cerr << "Error getting hc!" << endl;
      return;
    }
    if (!hlight) {
      cerr << "Error getting hlight!" << endl;
      return;
    }
  }


  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2.0);

  double bptmin = 20.;
  double bptmax = 200.;

  double lightptmin = 20.;
  double lightptmax = 200.;

  double absetamin = 0.;
  double absetamax = 2.5;

  double fitptmin = 15.;
  double fitptmax = 215.;
  double etamin = -2.5;
  double etamax = 2.5;

  double lightSFptmin = 20.;
  double lightSFptmax = 500.;
  int NpFine = 500;

  TCanvas *can = new TCanvas();
  can -> Divide(3,1);

  TString opt = "hist e2";

  can -> cd(1);
  hb -> Draw(opt);
  can -> cd(2);
  hc -> Draw(opt);
  can -> cd(3);
  hlight -> Draw(opt);

}
