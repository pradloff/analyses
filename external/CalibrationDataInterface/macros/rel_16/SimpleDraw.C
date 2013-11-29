

#include "StandardRootIncludes.h"
#include "standardsetup.C"

#include "TF1.h"
#include "TF2.h"

#include "MyTypeDef.h"


void SimpleDraw(TString flavour = "B", bool SF = true)
{

  setup();

  TFile *_file0 = TFile::Open("TopCalibrations_rel16_prelim.root");

  TString ObjNameJP = "JetProb_SF";
  TString ObjNameSV = "SV0_SF";
  TString SFtag = "SF";
  if (flavour == "Light") {
    ObjNameJP = "combined_SF";
    ObjNameSV = "combined_SF";
  }

  if (!SF) {
    ObjNameJP = "top_Eff";
    ObjNameSV = "top_Eff";
    SFtag = "Eff";
  }

  MyType *sf1 = (MyType*) _file0 -> Get("JetProb/AntiKt4Topo/3_25/" + flavour + "/" + ObjNameJP);
  MyType *sf2 = (MyType*) _file0 -> Get("JetProb/AntiKt4Topo/2_05/" + flavour + "/" + ObjNameJP);
  MyType *sf3 = (MyType*) _file0 -> Get("SV0/AntiKt4Topo/5_85/" + flavour + "/" + ObjNameSV);

  MyType *sf1syst = (MyType*) _file0 -> Get("JetProb/AntiKt4Topo/3_25/" + flavour + "/" + ObjNameJP + "_syst");
  MyType *sf2syst = (MyType*) _file0 -> Get("JetProb/AntiKt4Topo/2_05/" + flavour + "/" + ObjNameJP + "_syst");
  MyType *sf3syst = (MyType*) _file0 -> Get("SV0/AntiKt4Topo/5_85/" + flavour + "/" + ObjNameSV + "_syst");
  
  cout << " sf1=" << sf1
       << " sf2=" << sf2
       << " sf3=" << sf3
       << " sf1syst=" << sf1syst
       << " sf2syst=" << sf2syst
       << " sf3syst=" << sf3syst
       << endl;

  TCanvas *can = new TCanvas("check_" + flavour + "_" + SFtag);
  can -> Divide(3,2);

  TString opt = "";
  if (flavour == "Light" || !SF)
    opt = "colz";

  can -> cd(1);
  sf1 -> Draw(opt);
  can -> cd(2);
  sf2 -> Draw(opt);
  can -> cd(3);
  sf3 -> Draw(opt);

  can -> cd(4);
  sf1syst -> Draw(opt);
  can -> cd(5);
  sf2syst -> Draw(opt);
  can -> cd(6);
  sf3syst -> Draw(opt);

  can -> Print(TString(can->GetName()) + ".eps");

}
