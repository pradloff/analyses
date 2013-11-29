//
// File     : fitTopHists.C
// Author   : Frank Filthaut
// Modified : Jiri Kvita 
//            includes, ACLIC-compilable; void->double in SF functions definitions, cloning to 'defaults'; 10.2.2011

// TF1* fbSFstatin = 0;
// TF1* fbSFsystin = 0;


#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMatrixDSym.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
//#include "T.h"

#include "Tools.C"

#include "SF.hpp"
#include "SF_Light.hpp"

#include "standardsetup.C"

const TString JetDirName = "AntiKt4Topo";
const TString FitOpt = "N"; // default" "IN" ; consider also "V" for verbose
// I: use function integral in each bin instead of central value
// N: don't store graphical result uf the fit


// double CombinedUncertainty(Double_t* xx, Double_t* par);
/* see the SF.hpp
  double LightFlavourSF(Double_t* xx, Double_t* par);
  double LightFlavourSFSyst(Double_t* xx, Double_t* par);
  double HeavyFlavourSF(Double_t* xx, Double_t* par);
  double HeavyFlavourSFSyst(Double_t* xx, Double_t* par);
*/


// ______________________________________________________
// ______________________________________________________
// ______________________________________________________


void fitTopHists(TString TaggerName = "SV0",
		 TString OPName = "5_85",
		 TString OutFileName = "TopCalibrations_EPS2011.root",
		 TString OutWriteMode = "recreate",
		 bool cloneToDefault = true, 
		 bool doWrite = true, 
		 bool doDisplay = true,
		 bool doFit = true,
		 bool doZoomZ = true, 
		 bool doOwnRatios = true,
		 bool prelim = false)
{

  //  ROOT::Math::Minimizer::SetDefaultMaxFunctionCalls(10000);

  setup();
  cout << endl;
  cout << "=== Running fitTopHists_rel16 ===" << endl;
  cout << endl;
  cout << " TaggerName: " << TaggerName.Data() << endl 
       << " OPName: " << OPName.Data() << endl 
       << " OutFileName: " << OutFileName.Data() << endl 
       << " OutWriteMode: " << OutWriteMode.Data() << endl 
       << endl;
  cout << endl;
  cout << " cloneToDefault : " << cloneToDefault << endl 
       << " doWrite        : " << doWrite << endl 
       << " doDisplay      : " << doDisplay << endl 
       << " doFit          : " << doFit << endl 
       << " doZoomZ        : " << doZoomZ << endl 
       << " doOwnRatios    : " << doOwnRatios << endl 
       << endl;
  cout << endl;

  TFile *efffile = 0;
  if (OPName == "2_40")
    efffile = new TFile("mergeEt15_2.root", "read");
  else 
    efffile = new TFile("mc10b_merge.root", "read");
  
  // retrieve the histograms
  TH2* hb = 0;
  TH2* hc = 0;
  TH2* hlight = 0;

  // rel.16:
  // produced by Martin zur Nedden et. al.
  TString cutname = "";
  // Legend:

  // new:
  if (!prelim) {
    if (OPName == "5_85") // SV0 50
      cutname = "_cut50";

    if (OPName == "3_25") // JetProb 50
      cutname = "_cut50";
    if (OPName == "2_05") // JetProb 70
      cutname = "_cut70";
    if (OPName == "1_40") // JetProb 80
      cutname = "_cut80";
    if (OPName == "0_60") // JetProb 90
      cutname = "_cut90";

    if (OPName == "7_60") // SV1IP3D 50
      cutname = "_cut50";
    if (OPName == "4_50") // SV1IP3D 60
      cutname = "_cut60";
    if (OPName == "1_55") // SV1IP3D 70
      cutname = "_cut70";
    if (OPName == "-0_85") // SV1IP3D 80
      cutname = "_cut80";

    if (OPName == "3_00") // JetFitterCOMBNN 50
      cutname = "_cut50";
    if (OPName == "2_00") // JetFitterCOMBNN 60
      cutname = "_cut60";
    if (OPName == "0_35") // JetFitterCOMBNN 70
      cutname = "_cut70";
    if (OPName == "-1_25") // JetFitterCOMBNN 80
      cutname = "_cut80";
    if (OPName == "2_40") // JetFitterCOMBNN 80
      cutname = "_cut57";

  }

  TString FileTaggerName = TaggerName;
  if (TaggerName == "SV1IP3D")
    FileTaggerName = "Baseline";

  TString dirname = "btag/AntiKt4TopoEM/" + FileTaggerName + "/";
  if (doOwnRatios) {
    TH2 *hb_num = GetSpecialHisto(efffile, dirname, FileTaggerName, "b", cutname); 
    TH2 *hb_denom = GetSpecialHisto(efffile, dirname, FileTaggerName, "b", "_all"); 
    hb = MakeCustomRatio(hb_num, hb_denom);
    TH2 *hc_num = GetSpecialHisto(efffile, dirname, FileTaggerName, "c", cutname); 
    TH2 *hc_denom = GetSpecialHisto(efffile, dirname, FileTaggerName, "c", "_all"); 
    hc = MakeCustomRatio(hc_num, hc_denom);
    TH2 *hlight_num = GetSpecialHisto(efffile, dirname, FileTaggerName, "uds", cutname); 
    TH2 *hlight_denom = GetSpecialHisto(efffile, dirname, FileTaggerName, "uds", "_all"); 
    hlight = MakeCustomRatio(hlight_num, hlight_denom);
  } else {
    hb = GetSpecialHisto(efffile, dirname, FileTaggerName, "b", cutname, "eff");
    hc = GetSpecialHisto(efffile, dirname, FileTaggerName, "c", cutname, "eff");
    hlight = GetSpecialHisto(efffile, dirname, FileTaggerName, "uds", cutname, "eff");
  }
  
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
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2.0);

  double bptmin = 20.;
  double bptmax = 250.;

  double bptminSF = 20.;
  double bptmaxSF = 200.;

  double lightptmin = 20.;
  double lightptmax = 200.;

  double absetamin = 0.;
  double absetamax = 2.5;

  double fitptmin = 15.;
  double fitptmax = 250.;
  double etamin = -2.5;
  double etamax = 2.5;

  // note that this tries to respect the binning and do 1GeV granularity
  int NpFineHeavy = 180;

  double lightSFptmin = 20.;
  double lightSFptmax = 500.;
  int NpFineLightEta = 250;
  int NpFineLightPt = 480/5; // 480



  // Efficiencies Parametrisations //
  // x[0] ... |eta|
  // x[1] ... pt
  // crashes:
  // TString bFitFormula = "[0]*( 1. + [1]*log(abs(x[0]))*exp(-[2]*x[1]) )*([3]*log(abs(x[0])) + [4]*pow(log(abs(x[0])),2.0)) ";
  //  TString bFitFormula = "[0]*( 1. + [1]*exp(-[2]*x[1]) )*(log(abs(x[0])) + [3]*pow(log(abs(x[0])),2.0))";
    // ORIG:
  TString bFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0])";

  TString cFitFormula = bFitFormula;

  bool MoreDetailedFit1 = OPName == "5_85";
  bool MoreDetailedFit2 = OPName == "2_00" || OPName == "-1_25" || OPName == "0_35";// || OPName == "2_40";
  bool Special1 = OPName == "3_25"; // || OPName == "" || OPName == "" || OPName == "";
  bool Special2 = OPName == "1_40";
  bool Special3 = OPName == "0_60"; // || OPName == "1_40" || OPName == "3_25";
  bool Special4 = OPName == "2_40";

  if (MoreDetailedFit1) {
    bFitFormula += " + [6]*pow(abs(x[0]),3.0) + [7]*pow(log(x[1]),3.0)"; //  +  [8]*pow(log(x[1])*abs(x[0]), 2.0)";
  } else if (MoreDetailedFit2) {
    bFitFormula += " + [6]*pow(abs(x[0]),3.0) + [7]*pow(log(x[1]),3.0)";// + [8]*pow(log(x[1])*abs(x[0]), 2.0)";
  } else if (Special1) {
    bFitFormula =  	"([0] + [1]*exp(-[2]*x[1]) +[3]*abs(x[0]) ) * exp( -(abs(x[0]))^2/(2*([4]*(1+[5]*x[1])^2)) )";
    // SP2: ext
    //    bFitFormula = "[0]*( 1. - [1]*exp(-[2]*x[1]) )*exp(-[4]*abs(x[0])*x[1]) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   ) + [3]*pow(log(abs(x[0])),1.0) + [7]*pow(log(abs(x[0])),2.0)"
    // SM2 ext1: "[0]*( 1. - [1]*exp(-[2]*x[1]) )*(1. - [3]*exp(-[4]*abs(x[0])) ) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   ) + [7]*pow(abs(x[0]),1.0) + [8]*pow(abs(x[0]),2.0)";
    // SP3:
    //    bFitFormula = "[0]*( 1. - [1]*exp(-[2]*x[1]) )*(1. - [3]*exp(-[4]*abs(x[0])) ) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   ) + [7]*log(x[1])*abs(x[0]) + [8]*pow(log(x[1])*abs(x[0]), 2.0)";
    // orig:    bFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*x[1] + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*pow(log(x[1])*abs(x[0]), 2)";
    // bFitFormula = "[0]*( 1. - [1]*exp(-[2]*x[1]) )*(1. - [3]*exp(-[4]*abs(x[0])) ) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   )";
  } else if (Special2) {
    bFitFormula = "[0]*( 1. - [1]*exp(-[2]*x[1]) )*(1. - [3]*exp(-[4]*abs(x[0])) ) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   )"; // + [7]*pow(abs(x[0]),1.0) 
    //    bFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*x[1] + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*pow(log(x[1])*abs(x[0]), 2.0)";
  } else if (Special3) {
    bFitFormula = "[0]*( 1. - [1]*exp(-[2]*x[1]) )*(1. - [3]*exp(-[4]*abs(x[0])) ) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   ) + [7]*log(x[1])*abs(x[0]) + [8]*pow(log(x[1])*abs(x[0]), 2.0)";
    // + ( [5]*log(x[1])*abs(x[0]) + [6]*pow(log(x[1])*abs(x[0]), 2.0)"; //  + [7]*pow(abs(x[0]),2.0))";
  } else if (Special4) {
    // nuclear potential formula to account for efficiency drop around eta=1.5 (barrel/endcap)
    bFitFormula += " + [6]*x[1]/(1+exp((abs(x[0])-[7])/[8]))";
  }


  //  TString cFitFormula = bFitFormula; // "[0] + [1]*log(x[0]) +[2]*abs(x[1]) + [3]*pow(log(x[0]),2.0) + [4]*pow(abs(x[1]),2.0) + [5]*log(x[0])*abs(x[1])";
  // hack
  //if (0 && TaggerName.Contains("COMBNN")) {
  //    bFitFormula = "[0] + [1]*log(x[1]) +[2]*abs(x[0]) + [3]*pow(abs(x[0]),2.0) + [4]*log(x[1])*abs(x[0])";
  //  }

  // first attempt: polynomials in (ln(pt), eta)
  TF2* fbEff = new TF2("fbEff",
		       bFitFormula,
		       GetxIsAbsEta(hb) ? absetamin : etamin, etamax,
		       fitptmin, fitptmax);
  int nbpars = fbEff -> GetNpar();
  double *parb = new double[nbpars];
  // ORIG:
  //  if (0 && TaggerName.Contains("COMBNN")) {
  //    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0; parb[3] = -0.04;  parb[4] = 0;
  //  } else {
  if (MoreDetailedFit1) {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    parb[6] = 0.;    parb[7] = 0.; //   parb[8] = 0.;
    cout << "OK, setting MoreDetailedFit1 params..." << endl;
  } else if (MoreDetailedFit2) {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    parb[6] = 0.;    parb[7] = 0.; // parb[8] = 0.;
    cout << "OK, setting MoreDetailedFit2 params..." << endl;
  } else if (Special1) {
    cout << "OK, setting Special1 params..." << endl;
    //  parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0; parb[6] = 0;
    // like SP2: extended
    /*
    parb[0] = 0.52;  parb[1] = 0.2;  parb[2] = 0.057;  parb[3] = 0.03;  parb[4] = -0.9;  parb[5] = 0.2;
    parb[6] = 0.02;  parb[7] = -0.001; parb[8] = 0.001;
    fbEff -> SetParLimits(1, 0., 1.);
    fbEff -> SetParLimits(3, 0., 1.);
    fbEff -> SetParLimits(5, 0., 1.);
    */
   
    parb[0] = 0.55;  parb[1] = -0.1;  parb[2] = 0.01;  parb[3] = 0.1;  parb[4] = 1.;  parb[5] = 0.01;
    /*
    // like SP3:
    parb[0] = 1.04;  parb[1] = 0.83;  parb[2] = 0.07;  parb[3] = 0.1;  parb[4] = 1.4;  
    parb[5] = 0.0001;    parb[6] = 0.53; parb[7] = -0.02;  parb[8] = 0.0014;
    fbEff -> SetParLimits(1, 0., 1.);
    fbEff -> SetParLimits(3, 0., 1.);
    */

  } else if (Special2) {
    cout << "OK, setting Special2 params..." << endl;
    parb[0] = 0.54;  parb[1] = 0.3;  parb[2] = 0.02;  parb[3] = 0.1;  parb[4] = 0.01;  parb[5] = 0.1;
    parb[6] = 0.1;  // parb[7] = 0.; // parb[8] = 0.;
    fbEff -> SetParLimits(1, 0., 1.);
    fbEff -> SetParLimits(3, 0., 1.);
    fbEff -> SetParLimits(5, 0., 1.);
  } else if (Special3) {
    cout << "OK, setting Special3 params..." << endl;
    parb[0] = 1.04;  parb[1] = 0.83;  parb[2] = 0.07;  parb[3] = 0.1;  parb[4] = 1.4;  
    parb[5] = 0.0001;    parb[6] = 0.53; parb[7] = -0.02;  parb[8] = 0.0014;
    fbEff -> SetParLimits(1, 0., 1.);
    fbEff -> SetParLimits(3, 0., 1.);
  } else if (Special4) {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    parb[6] = 0.0004 ; parb[7] = 1.42; parb[8] = 0.14; 
    // fbEff -> SetParLimits(6, -0.2, 0.2);
    fbEff -> SetParLimits(7, 1.3, 1.7); 
    fbEff -> SetParLimits(8, 0., 0.3);
    cFitFormula = bFitFormula;
  } else {
    cout << "OK, setting default params..." << endl;
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
  }

  fbEff->SetParameters(parb);
  if (doFit) {
    cout << "*** Fitting hb " << hb -> GetName() << ", npars=" << nbpars << endl;
    // check that the x and y axis are what we think they should be:
    cout << "Fitform: " << bFitFormula.Data() << endl;
    for (int ip = 0; ip < nbpars; ++ip)
      cout << " par" << ip << "=" << fbEff -> GetParameter(ip) << " ";
    cout << endl;
    CheckHistoRanges(hb);
    hb->Fit(fbEff, FitOpt);
    if (Special4) {
      for (int haha = 0; haha < 4; haha++)
	hb->Fit(fbEff, FitOpt);
    }

    /*
      if (Special4) {
      hb->Fit(fbEff, FitOpt);
      hb->Fit(fbEff, FitOpt);
      hb->Fit(fbEff, FitOpt);
    }
    */
    for (int ip = 0; ip < nbpars; ++ip)
      cout << " parb" << ip << "=" << fbEff -> GetParameter(ip) << " ";
    cout << endl;
    PrintChi2(fbEff);
  }
  
  //  fbEff->Draw("lego2");
  //return;

  Double_t* covb = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covbEff = new TMatrixDSym(nbpars, covb);
  covbEff -> Print();

  // SV0
  if (MoreDetailedFit1) {
    //cFitFormula = "([0] + [1]*exp(-[2]*x[1]) +[3]*abs(x[0]) ) * exp( -(abs(x[0]))^2/(2*([4]*(1+[5]*x[1])^2)) )"; // SP1
	// cFitFormula = "[0]*( 1. - [1]*exp(-[2]*x[1]) )*(1. - [3]*exp(-[4]*abs(x[0])) ) * (1. + [5]*exp([6]*log(x[1])*abs(x[0]))   )"; // like Special2 bFitFormula;
  }

  //  cFitFormula =  "([0] + [1]*exp(-[2]*x[1]) +[3]*abs(x[0]) ) * exp( -(abs(x[0]))^2/(2*([4]*(1+[5]*x[1])^2)) )";
  // "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0])";
  TF2* fcEff = new TF2("fcEff",
		       cFitFormula,
		       GetxIsAbsEta(hc) ? absetamin : etamin, etamax,
		       fitptmin,fitptmax
		       ); // absetamin, absetamax
  int ncpars = fcEff -> GetNpar();
  double *parc = new double[ncpars];

  //parc[0] = 0.55;  parc[1] = -0.1;  parc[2] = 0.01;  parc[3] = 0.1;  parc[4] = 1.;  parc[5] = 0.01;

  //  if (MoreDetailedFit1) {
    /*    
	  cout << "OK, setting Special1params..." << endl;
	  parc[0] = 0.54;  parc[1] = 0.3;  parc[2] = 0.02;  parc[3] = 0.1;  parc[4] = 0.01;  parc[5] = 0.1;
	  parc[6] = 0.1;  // parc[7] = 0.; // parc[8] = 0.;
	  fcEff -> SetParLimits(1, 0., 1.);
	  fcEff -> SetParLimits(3, 0., 1.);
	  fcEff -> SetParLimits(5, 0., 1.);
    */
    //cout << "OK, setting Special1params..." << endl;
    //parc[0] = 0.55;  parc[1] = -0.1;  parc[2] = 0.01;  parc[3] = 0.1;  parc[4] = 1.;  parc[5] = 0.01;

    //} else {
    //  parc[0] = 0.1;  parc[1] = -1.;  parc[2] = 0.06;  parc[3] = 0.02;  //  parc[4] = 0.1;  //  parc[5] = 0;
  for (int ip = 0; ip < ncpars; ++ip)
    parc[ip] = parb[ip];
  parc[0] /= 3.;

  if (Special4) {
    parc[0] = 0.4/3.;  parc[1] = 0.03;  parc[2] = 0;  parc[3] = 0;  parc[4] = -0.04;  parc[5] = 0;
    parc[6] = 0.0004 ; parc[7] = 1.42; parc[8] = 0.14; 
    // fbEff -> SetParLimits(6, -0.2, 0.2);
    fcEff -> SetParLimits(7, 1.3, 1.7); 
    fcEff -> SetParLimits(8, 0., 0.3);
  }

  fcEff->SetParameters(parc);
  if (doFit) {
    cout << "*** Fitting hc..." << endl;
    CheckHistoRanges(hc);
    hc->Fit(fcEff, FitOpt);
    if (Special4) {
      for (int haha = 0; haha < 4; haha++)
	hc->Fit(fcEff, FitOpt);
    }
    PrintChi2(fcEff);
  }
  Double_t* covc = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covcEff = new TMatrixDSym(ncpars, covc);
  covcEff -> Print();

  // reasonable:
  
  TString lFitFormula = "[0] + [1]*x[1] + [2]*abs(x[0]) + [3]*pow(x[1],2.0) + [4]*pow(abs(x[0]),2.0) + [5]*x[1]*abs(x[0]) + [6]*pow(log(x[1]),4.0) + [7]*pow(abs(x[0]),4.0) + [8]*pow(x[1]*abs(x[0]),2.0)";
  bool MoreDetailedFit =  OPName == "0_60" || OPName == "1_40" || OPName == "-0_85" || OPName == "2_05" || OPName == "1_55";
  if (MoreDetailedFit)
    // reasonable:
    lFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*x[1] + [7]*pow(abs(x[0]),3.0) + [8]*pow(log(x[1])*abs(x[0]),2.0) + [9]*pow(abs(x[0]),4.0) + [10]*pow(x[1], 2.0) + [11]*x[1]*abs(x[0])";
  //  if (OPName == "2_40")
  // lFitFormula = bFitFormula;

  // tried TString lFitFormula = "[0] + [1]*x[1] + [2]*abs(x[0]) + [3]*pow(x[1],2.0) + [4]*pow(abs(x[0]),2.0) + [5]*x[1]*abs(x[0]) + [6]*pow(log(x[1]),4.0) + [7]*pow(abs(x[0]),4.0) + [8]*pow(x[1]*abs(x[0]),2.0) + [9]*pow(abs(x[0]),3.0)"; // + [10]*pow(abs(x[0]),4.0) + [11]*pow(x[1]*abs(x[0]),2.0)

 
  /*
// devel
  TString lFitFormula = "[0] + [1]*x[1] + [2]*abs(x[0]) + [3]*pow(x[1],2.0) + [4]*pow(abs(x[0]),2.0) + [5]*x[1]*abs(x[0]) + [6]*pow(log(x[1]),4.0) + [7]*pow(abs(x[0]),4.0) + [8]*pow(x[1]*abs(x[0]),2.0)";
  MoreDetailedFit =  OPName == "0_60" || OPName == "1_40" || OPName == "-0_85" || OPName == "2_05";
  if (MoreDetailedFit)
    lFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*pow(abs(x[0]),3.0) + [7]*pow(log(x[1])*abs(x[0]),2.0) + [8]*pow(abs(x[0]),4.0) + [9]*x[1]*abs(x[0])";
  */

    // lFitFormula = "[0] + [1]*pow(log(x[1]), [6]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [7]*pow(abs(x[0]),3.0) + [8]*pow(log(x[1])*abs(x[0]),2.0) + [9]*pow(abs(x[0]),4.0) + [10]*pow(x[1], 2.0) + [11]*x[1]*abs(x[0])";

    //if (OPName == "2_05")
  //   lFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*pow(log(x[1]),4.0) + [7]*pow(abs(x[0]),4.0) + [8]*pow(log(x[1])*abs(x[0]),2.0)";

  //  if (OPName == "1_55")
  //    lFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*pow(log(x[1]),4.0) + [7]*pow(abs(x[0]),4.0) + [8]*pow(log(x[1])*abs(x[0]),2.0)";


  // if (Special4) {
  //  lFitFormula = "[0] + [1]*x[1] + [2]*abs(x[0]) + [3]*pow(x[1],2.0) + [4]*pow(abs(x[0]),2.0) + [5]*x[1]*abs(x[0]) + [6]*pow(x[1]*abs(x[0]),2.0)";
  //}
  if (Special4) {
    lFitFormula = bFitFormula;
  } 

  TF2* flightEff = new TF2("flightEff", lFitFormula,
			   GetxIsAbsEta(hlight) ? absetamin : etamin, etamax,
			   fitptmin, fitptmax);
  int nlpars = flightEff -> GetNpar();
  double *parlight = new double[nlpars];
  parlight[0] = 0.04;
  parlight[1] = 0.0;
  parlight[2] = 0.01;
  parlight[3] = 0.01;
  parlight[4] = 0.001;
  parlight[5] = 0.001;
  parlight[6] = 0.0001;
  // flightEff -> SetParLimits(6, 0., 4.);
  

  if (!Special4) {
    parlight[7] = 0.0001;
    parlight[8] = 0.001;
  }
  if (MoreDetailedFit && !Special4) {
      parlight[9] = 0.000;
      parlight[10] = 0.000;
      parlight[11] = 0.000;
  }
  //  parlight[10] = 0.000;
  //  parlight[11] = 0.000;

  if (Special4) {
    parlight[0] = 0.01;  parlight[1] = 0.03;  parlight[2] = 0;  parlight[3] = 0;  parlight[4] = -0.04;  parlight[5] = 0;
    parlight[6] = 0.0004 ; parlight[7] = 1.42; parlight[8] = 0.14; 
    // fbEff -> SetParLimits(6, -0.2, 0.2);
    flightEff -> SetParLimits(7, 1.3, 1.7); 
    flightEff -> SetParLimits(8, 0., 0.3);
  } 

  flightEff->SetParameters(parlight);
  if (doFit) {
    cout << "*** Fitting hlight..." << endl;
    CheckHistoRanges(hlight);
    hlight->Fit(flightEff, FitOpt); // "V"
    if (Special4) {
      for (int haha = 0; haha < 4; haha++)
	hlight->Fit(flightEff, FitOpt); // "V"
    }
    PrintChi2(flightEff);
  }
  Double_t* covlight = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covlightEff = new TMatrixDSym(nlpars, covlight);
  covlightEff -> Print();

  // this is very important, as the axis meaning in the Calib tool as taken from the following axis title names!
  if (GetyIsAbsEta(hb))
    fbEff->SetTitle("b-jet selection efficiency;abseta;pt");
  else
    fbEff->SetTitle("b-jet selection efficiency;eta;pt");
  if (GetyIsAbsEta(hc))
    fcEff->SetTitle("c-jet selection efficiency;abseta;pt");
  else
    fcEff->SetTitle("c-jet selection efficiency;eta;pt");
  if (GetyIsAbsEta(hlight))
    flightEff->SetTitle("light-jet selection efficiency;abseta;pt");
  else
    flightEff->SetTitle("light-jet selection efficiency;eta;pt");
  // Eff. Syst. Parametrisations //

  // dummy systematic uncertainty parametrisations
  TF1* fbEffSyst = new TF1("fbEffSyst", "[0]", bptmin, bptmax);
  fbEffSyst->SetTitle("dummy b-jet selection efficiency syst uncertainty");
  Double_t parbEffSyst[1];
  parbEffSyst[0] = 0.0;
  fbEffSyst->SetParameters(parbEffSyst);

  TF1* fcEffSyst = new TF1("fcEffSyst", "[0]", bptmin, bptmax);
  fcEffSyst->SetTitle("dummy c-jet selection efficiency syst uncertainty");
  Double_t parcEffSyst[1];
  parcEffSyst[0] = 0.0;
  fcEffSyst->SetParameters(parcEffSyst);

  TF1* flightEffSyst = new TF1("flightEffSyst", "[0]", lightptmin, lightptmax);
  flightEffSyst->SetTitle("dummy light-jet selection efficiency syst uncertainty");
  Double_t parlightEffSyst[1];
  parlightEffSyst[0] = 0.0;
  flightEffSyst->SetParameters(parlightEffSyst);

 
  // Scale Factors Parametrisations and their Syst //

  // The TF1's above were provided in sampled form, which implies that it's no use
  // attempting to use a covariance matrix to compute the statistical uncertainties.
  // So instead, we lump the statistical and systematic uncertainties together, into
  // an overall "systematic" uncertainty
  TF1* fbSelectionSF = new TF1("bSFFunction", HeavyFlavourSF, bptminSF, bptmaxSF, 0);
  fbSelectionSF->SetTitle(TaggerName + "-based b-jet selection efficiency SF;pt");
  fbSelectionSF -> SetNpx(NpFineHeavy);
  Double_t dummy = 0; 
  Double_t* dummyCov = &dummy;
  TMatrixDSym* covbSF = new TMatrixDSym(0, dummyCov);
  TF1* fbSelectionSFSyst = new TF1("fbSelectionSFSyst", HeavyFlavourSFSyst, bptminSF, bptmaxSF, 1);
  fbSelectionSFSyst->SetParameter(0, 1.0); // scale to inflate the syst \oplus stat
  fbSelectionSFSyst->SetTitle(TaggerName + "-based b-jet SF systematic uncertainty;pt");
  fbSelectionSFSyst -> SetNpx(NpFineHeavy);

  // c-jet efficiency scale factor: as for b-jets, only with the uncertainty doubled.
  TF1* fcSelectionSF = (TF1*)fbSelectionSF->Clone("cSFFunction");
  fcSelectionSF->SetTitle(TaggerName + "-based c-jet selection efficiency SF;pt");
  fcSelectionSF -> SetNpx(NpFineHeavy);
  TMatrixDSym* covcSF = new TMatrixDSym(0, dummyCov);
  TF1* fcSelectionSFSyst = new TF1("fcSelectionSFSyst", HeavyFlavourSFSyst, bptminSF, bptmaxSF, 1);
  fcSelectionSFSyst->SetParameter(0, 2.0); // c syst 2x of the b syst;-)
  fcSelectionSFSyst->SetTitle(TaggerName + "-based c-jet SF systematic uncertainty;pt");
  fcSelectionSFSyst -> SetNpx(NpFineHeavy);
 
  // light-jet efficiency scale factor: similar fashion
  TF2* flightSelectionSF = new TF2("flightSelectionSF", LightFlavourSF, lightSFptmin, lightSFptmax, absetamin, absetamax, 0);
  flightSelectionSF->SetTitle(TaggerName + "-based light-jet selection efficiency SF;pt;abseta");
  flightSelectionSF -> SetNpx(NpFineLightPt);
  flightSelectionSF -> SetNpy(NpFineLightEta);
  // Double_t parlightSF[1];
  // parlightSF[0] = 1.15;
  // flightSelectionSF->SetParameters(parlightSF);
  // Double_t err2lightSF = TMath::Power(parlightSF[0] * 0.50, 2);
  TMatrixDSym* covlightSF = new TMatrixDSym(0, dummyCov);
  TF2* flightSelectionSFSyst = new TF2("flightSelectionSFSyst", LightFlavourSFSyst, lightSFptmin, lightSFptmax, absetamin, absetamax, 1);
  flightSelectionSFSyst->SetTitle(TaggerName + "-based light-jet selection efficiency SF syst uncertainty;pt;abseta");
  //  flightSelectionSFSyst->SetParameter(0, 1.0);
  flightSelectionSFSyst -> SetNpx(NpFineLightPt);
  flightSelectionSFSyst -> SetNpy(NpFineLightEta);
  // Double_t parlightSFSyst[1];
  // parlightSFSyst[0] = 0.0;
  // flightSelectionSFSyst->SetParameters(parlightSFSyst);

  // write to output file

  if (doWrite) {
    
    TFile* fout = new TFile(OutFileName.Data(), OutWriteMode.Data());
    TDirectory* tdir = 0;
    if (!fout -> cd(TaggerName.Data())) {
      tdir = fout->mkdir(TaggerName.Data());
      tdir->cd();
    } else tdir = gDirectory;

    TDirectory* jdir = 0;
    if (!tdir -> cd(JetDirName.Data())) {
      jdir = tdir->mkdir(JetDirName.Data());
      jdir->cd();
    } else jdir = gDirectory;

    TDirectory* odir = 0;
    if (!jdir -> cd(OPName.Data())) {
      odir = jdir->mkdir(OPName.Data());
      odir->cd();
    }

    TDirectory* bdir = 0;
    if (!odir -> cd("B")) {
      bdir = odir->mkdir("B");
      bdir->cd();
    }
    cout << "now in directory "  << endl; gDirectory->pwd();
    fbSelectionSF->Write(TaggerName + "_SF");
    fbSelectionSFSyst->Write(TaggerName + "_SF_syst");
    covbSF->Write(TaggerName + "_SF_stat");
    fbEff->Write("top_Eff");
    fbEffSyst->Write("top_Eff_syst");
    covbEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClones(fbSelectionSF, fbSelectionSFSyst, covbSF, fbEff, fbEffSyst, covbEff);
    bdir->ls();

    TDirectory* cdir = 0;
    if (!odir -> cd("C")) {
      cdir  = odir->mkdir("C");
      cdir -> cd();
    }
    cout << "now in directory "  << endl; gDirectory->pwd();
    fcSelectionSF->Write(TaggerName + "_SF");
    fcSelectionSFSyst->Write(TaggerName + "_SF_syst");
    covcSF->Write(TaggerName + "_SF_stat");
    fcEff->Write("top_Eff");
    fcEffSyst->Write("top_Eff_syst");
    covcEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClones(fcSelectionSF, fcSelectionSFSyst, covcSF, fcEff, fcEffSyst, covcEff);
    cdir->ls();

    TDirectory* lightdir = 0;
    if (!odir -> cd("Light")) {
      lightdir = odir->mkdir("Light");
      lightdir->cd();
    }
    cout << "now in directory "  << endl; gDirectory->pwd();
    //    TString helpname = "preliminary";
    //    TString helpname = TaggerName;
    // like rel.15:
    TString helpname = "combined";
    flightSelectionSF->Write(helpname + "_SF");
    flightSelectionSFSyst->Write(helpname + "_SF_syst");
    covlightSF->Write(helpname + "_SF_stat");
    flightEff->Write("top_Eff");
    flightEffSyst->Write("top_Eff_syst");
    covlightEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClones2D(flightSelectionSF, flightSelectionSFSyst, covlightSF, flightEff, flightEffSyst, covlightEff);
    lightdir->ls();

    fout->Close();
    delete fout;
  }

  // fin->Close();
  // delete fin;

  // display results

  if (doDisplay) {

    DrawSFs(fbSelectionSF, fbSelectionSFSyst, "b", TaggerName + "_" + OPName, "l");
    DrawSFs(fcSelectionSF, fcSelectionSFSyst, "c", TaggerName + "_" + OPName, "l");
    DrawSFs(flightSelectionSF, flightSelectionSFSyst, "light", TaggerName + "_" + OPName, "colz");

    TString tag = TaggerName + "_" + OPName;
    double zmin = 0.90;// zmin
    double zmax = 1.1; // 1.05
    TString effopt = "lego e2";
    TString checkopt = "colz";

    TCanvas* cb = new TCanvas("cb" + tag, "MC b-jet efficiency", 600, 400);
    hb->DrawCopy(effopt);
    PrintCanvasAs(cb, ".eps");

    // cross-check: divide the histogram by the fit
    TCanvas* cbCheck = new TCanvas("cbCheck" + tag, "MC b-jet efficiency cross-check", 600, 400);
    TH2* hbCheck = (TH2*)hb->Clone("hbCheck");
    hbCheck->Divide(fbEff);
    TH2* hbCheck_copy = (TH2*) hbCheck->DrawCopy(checkopt);
    PrintCanvasAs(cbCheck, ".eps");
    hbCheck_copy -> Draw("boxtext");
    PrintCanvasAs(cbCheck, "_box.eps");
    if (doZoomZ) {
      hbCheck_copy->SetMinimum(zmin);
      hbCheck_copy->SetMaximum(zmax);
      hbCheck_copy->SetDrawOption(checkopt);
      hbCheck_copy->Draw(checkopt);
      gPad -> Update();
      PrintCanvasAs(cbCheck, "_zoom.eps");
    }
    CheckHistoFitDiff2D(hb, fbEff);

    TCanvas* cc = new TCanvas("cc" + tag, "MC c-jet efficiency", 600, 400);
    hc -> SetMaximum(0.25);
    if (OPName.Contains("2_05"))
      hc -> SetMaximum(0.60);
    if (OPName.Contains("1_40"))
      hc -> SetMaximum(1.);
    if (OPName.Contains("0_60"))
      hc -> SetMaximum(1.);
    if (OPName.Contains("7_60"))
      hc -> SetMaximum(0.25);
    if (OPName.Contains("4_50"))
      hc -> SetMaximum(0.25);
    if (OPName.Contains("-0_85"))
      hc -> SetMaximum(0.60);
    if (OPName.Contains("3_00"))
      hc -> SetMaximum(0.20);
    if (OPName.Contains("2_00") || OPName.Contains("2_40"))
      hc -> SetMaximum(0.20);
    if (OPName.Contains("0_35"))
      hc -> SetMaximum(0.25);
    if (OPName.Contains("1_55"))
      hc -> SetMaximum(0.50);
    if (OPName.Contains("-1_25"))
      hc -> SetMaximum(0.50);

   hc->DrawCopy(effopt);
    PrintCanvasAs(cc, ".eps");

    // cross-check: divide the histogram by the fit
    TCanvas* ccCheck = new TCanvas("ccCheck" + tag, "MC c-jet efficiency cross-check", 600, 400);
    TH2* hcCheck = (TH2*)hc->Clone("hcCheck");
    hcCheck->Divide(fcEff);
    TH2* hcCheck_copy = (TH2*) hcCheck->DrawCopy(checkopt);
    PrintCanvasAs(ccCheck, ".eps");
    hcCheck_copy -> Draw("boxtext");
    PrintCanvasAs(ccCheck, "_box.eps");
    if (doZoomZ) {
      hcCheck_copy->SetMinimum(zmin);
      hcCheck_copy->SetMaximum(zmax);
      hcCheck_copy->SetDrawOption(checkopt);
      hcCheck_copy->Draw(checkopt);
      gPad -> Update();
      PrintCanvasAs(ccCheck, "_zoom.eps");
    }
    CheckHistoFitDiff2D(hc, fcEff);

    TCanvas* clight = new TCanvas("clight" + tag, "MC light-jet efficiency", 600, 400);
    hlight -> SetMaximum(0.12);
    if (OPName.Contains("2_05"))
      hlight -> SetMaximum(0.25);
    if (OPName.Contains("1_40"))
      hlight -> SetMaximum(0.55);
    if (OPName.Contains("0_60"))
      hlight -> SetMaximum(0.95);
    if (OPName.Contains("5_85"))
      hlight -> SetMaximum(0.025);
    if (OPName.Contains("7_60"))
      hlight -> SetMaximum(0.025);
    if (OPName.Contains("4_50"))
      hlight -> SetMaximum(0.08);
    if (OPName.Contains("-0_85"))
      hlight -> SetMaximum(0.25);
    if (OPName.Contains("3_00"))
      hlight -> SetMaximum(0.012);
    if (OPName.Contains("2_00") || OPName.Contains("2_40"))
      hlight -> SetMaximum(0.012);
    if (OPName.Contains("0_35"))
      hlight -> SetMaximum(0.05);


    hlight->DrawCopy(effopt);
    PrintCanvasAs(clight, ".eps");

    // cross-check: divide the histogram by the fit
    TCanvas* clightCheck = new TCanvas("clightCheck" + tag, "MC light-jet efficiency cross-check", 600, 400);
    TH2* hlightCheck = (TH2*)hlight->Clone("hlightCheck");
    hlightCheck->Divide(flightEff);
    hlightCheck->SetStats(0);
    hlightCheck->SetDrawOption(checkopt);
    TH2* hlightCheck_copy = (TH2*) hlightCheck->DrawCopy(checkopt);
    PrintCanvasAs(clightCheck, ".eps");
    hlightCheck_copy -> Draw("boxtext");
    PrintCanvasAs(clightCheck, "_box.eps");
    if (doZoomZ) {
      hlightCheck_copy->SetMinimum(zmin);
      hlightCheck_copy->SetMaximum(zmax);
      hlightCheck_copy->SetDrawOption(checkopt);
      gPad -> Update();
      hlightCheck_copy->Draw(checkopt);
      gPad -> Update();
      PrintCanvasAs(clightCheck, "_zoom.eps");
    }
    CheckHistoFitDiff2D(hlight, flightEff);

    TString xtitle = "";
    TString ytitle = "";
    if (GetyIsPt(hb)) {
      ytitle = "jet p_{T} (GeV)";
      if (GetxIsAbsEta(hb))
	xtitle = "jet |#eta|";
      else
	xtitle = "jet #eta";
    } else {
      xtitle = "jet p_{T} (GeV)";
      if (GetyIsAbsEta(hb))
	ytitle = "jet |#eta|";
      else
	ytitle = "jet #eta"; 
    }
    

    // for PR purposes: plots of the parametrisations themselves
    TCanvas* cbParm = new TCanvas("cbParm" + tag, "MC b-jet efficiency parametrisation", 600, 400);
    fbEff->SetTitle(";" + xtitle + ";" + ytitle + ";b-jet tagging efficiency");
    fbEff->GetZaxis()->SetTitleOffset(1.2);
    fbEff->Draw("SURF1");
    MakeText(fbEff);
    PrintCanvasAs(cbParm, ".eps");
    
    TCanvas* ccParm = new TCanvas("ccParm" + tag, "MC c-jet efficiency parametrisation", 600, 400);
    fcEff->SetTitle(";" + xtitle + ";" + ytitle + ";c-jet tagging efficiency");
    fcEff->GetZaxis()->SetTitleOffset(1.2);
    fcEff->Draw("SURF1");
    MakeText(fcEff);
    PrintCanvasAs(ccParm, ".eps");

    TCanvas* clParm = new TCanvas("clightParm" + tag, "MC light-flavour jet efficiency parametrisation", 600, 400);
    flightEff->SetTitle(";" + xtitle + ";" + ytitle + ";light-flavour jet tagging efficiency");
    flightEff->GetZaxis()->SetTitleOffset(1.2);
    flightEff->Draw("SURF1");
    MakeText(flightEff);
    PrintCanvasAs(clParm, ".eps");

  }

  if (efffile)
    efffile->Close();
  //  delete f;
}

// void CombinedUncertainty(Double_t* xx, Double_t* par) {
//   double stat = fbSFstatin->Eval(xx[0]);
//   double syst = fbSFsystin->Eval(xx[0]);
//   return TMath::Sqrt(stat*stat+syst*syst)*par[0];
// }
