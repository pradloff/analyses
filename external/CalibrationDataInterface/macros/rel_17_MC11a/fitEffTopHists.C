//
// File     : fitTopHists.C
// Author   : Frank Filthaut
// Modified : Jiri Kvita 
//            includes, ACLIC-compilable; void->double in SF functions definitions, cloning to 'defaults'; 10.2.2011
//          : Martin zur Nedden
//            separating from filling of scale factors, only to be used for the efficiency fits

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
const TString FitOpt = "NV"; // default" "IN" ; consider also "V" for verbose
// I: use function integral in each bin instead of central value
// N: don't store graphical result uf the fit

void fitEffTopHists(TString TaggerName = "SV0",
		 TString OPName = "5_65",
		 TString OutFileName = "TopCalibrations_rel17_MC11a_Eff_new.root",
		 TString OutWriteMode = "recreate",
		 bool cloneToDefault = false, 
		 bool doWrite = true, 
		 bool doDisplay = true,
		 bool doFit = true,
		 bool doZoomZ = true, 
		 bool doOwnRatios = true,
		 bool prelim = false)
{

  setup();
  cout << endl;
  cout << "**** Running Efficiency for fitTopHists_rel17 ****" << endl;
  cout << endl;
  cout << " TaggerName   : " << TaggerName.Data() << endl 
       << " OPName       : " << OPName.Data() << endl 
       << " OutFileName  : " << OutFileName.Data() << endl 
       << " OutWriteMode : " << OutWriteMode.Data() << endl 
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
  efffile = new TFile("ttbar_rel17_MC11a_new.root", "read");
  //  efffile = new TFile("ttbar_rel17.root", "read");

  // retrieve the histograms
  TH2* hb = 0;
  TH2* hc = 0;
  TH2* hlight = 0;

  TString cutname = "";

  // define the tagger WP for the histogram names
  if (!prelim) {
    if (OPName == "5_65") // SV0 50
      cutname = "_cut50";
    if (OPName == "2_65") // JetProb 50
      cutname = "_cut50";
    if (OPName == "4_55") // SV1IP3D 60
      cutname = "_cut60";
    if (OPName == "1_70") // SV1IP3D 70
      cutname = "_cut70";
    if (OPName == "-0_80") // SV1IP3D 80
      cutname = "_cut80";
    if (OPName == "2_20") // JetFitterCOMBNN 57
      cutname = "_cut57";
    if (OPName == "1_80") // JetFitterCOMBNN 60
      cutname = "_cut60";
    if (OPName == "0_35") // JetFitterCOMBNN 70
      cutname = "_cut70";
    if (OPName == "-1_25") // JetFitterCOMBNN 80
      cutname = "_cut80";
    if (OPName == "0_910") // MV1 60
      cutname = "_cut60";
    if (OPName == "0_614") // MV1 70
      cutname = "_cut70"; 
    if (OPName == "0_416") // MV1 75
      cutname = "_cut75";
    if (OPName == "0_080") // MV1 85
      cutname = "_cut85";
  }

  TString FileTaggerName = TaggerName;
  if (TaggerName == "SV1IP3D")
    FileTaggerName = "Baseline";

  // fill the efficiency histograms with generic path to histogram
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

  // x[0] ... |eta|
  // x[1] ... pt
  // Fitformula for b-Effi
  TString bFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0])";
  // Fitformula for c-Effi
  TString cFitFormula = bFitFormula;
  // special setups for fitformulas and parameters
  bool SV0Fit = OPName == "5_65";
  bool IP3DSV1Fit = OPName == "-0_80" || OPName == "1_70";
  bool JetFitNNFit = OPName == "1_80" || OPName == "-1_25" || OPName == "0_35" || OPName == "2_20";
  bool JetFitNN57Fit = OPName == "0";

  if (SV0Fit) {
    bFitFormula += " + [6]*pow(abs(x[0]),3.0) + [7]*pow(log(x[1]),3.0)"; //  +  [8]*pow(log(x[1])*abs(x[0]), 2.0)";
    cout << "**** setting SV0Fit fit formula ****" << endl;
  } else if (JetFitNNFit) {
    bFitFormula += " + [6]*pow(abs(x[0]),3.0) + [7]*pow(log(x[1]),3.0)";// + [8]*pow(log(x[1])*abs(x[0]), 2.0)";
    cout << "**** setting JetFitNNFit fit formula ****" << endl;
  } else if (JetFitNN57Fit) {
    // nuclear potential formula to account for efficiency drop around eta=1.5 (barrel/endcap)
    bFitFormula += " + [6]*x[1]/(1+exp((abs(x[0])-[7])/[8]))";
    cout << "**** setting JetFitNN57Fit fit formula ****" << endl;
  }

  // first attempt: polynomials in (ln(pt), eta)
  TF2* fbEff = new TF2("fbEff",
		       bFitFormula,
		       GetxIsAbsEta(hb) ? absetamin : etamin, etamax,
		       fitptmin, fitptmax);
  int nbpars = fbEff -> GetNpar();
  double *parb = new double[nbpars];
  
  if (SV0Fit) {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    parb[6] = 0.;    parb[7] = 0.; //   parb[8] = 0.;
    cout << "**** setting SV0Fit params ****" << endl;
  } else if (JetFitNNFit) {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    parb[6] = 0.;    parb[7] = 0.; // parb[8] = 0.;
    cout << "**** setting JetFitNNFit params ****" << endl;
  } else if (JetFitNN57Fit) {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    parb[6] = 0.0004 ; parb[7] = 1.42; parb[8] = 0.14; 
    // fbEff -> SetParLimits(6, -0.2, 0.2);
    fbEff -> SetParLimits(7, 1.3, 1.7); 
    fbEff -> SetParLimits(8, 0., 0.3);
    cFitFormula = bFitFormula;
    cout << "**** setting JetFitNN57Fit params ****" << endl;
  } else {
    parb[0] = 0.4;  parb[1] = 0.03;  parb[2] = 0;  parb[3] = 0;  parb[4] = -0.04;  parb[5] = 0;
    cout << "**** setting default params ****" << endl;
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
    if (JetFitNN57Fit) {
      for (int haha = 0; haha < 4; haha++)
	hb->Fit(fbEff, FitOpt);
    }
    for (int ip = 0; ip < nbpars; ++ip)
      cout << " parb" << ip << "=" << fbEff -> GetParameter(ip) << " ";
    cout << endl;
    PrintChi2(fbEff);
  }
  
  Double_t* covb = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covbEff = new TMatrixDSym(nbpars, covb);
  covbEff -> Print();

  TF2* fcEff = new TF2("fcEff",
		       cFitFormula,
		       GetxIsAbsEta(hc) ? absetamin : etamin, etamax,
		       fitptmin,fitptmax
		       ); // absetamin, absetamax
  int ncpars = fcEff -> GetNpar();
  double *parc = new double[ncpars];
  
  for (int ip = 0; ip < ncpars; ++ip)
    parc[ip] = parb[ip];
  parc[0] /= 3.;
  
  if (JetFitNN57Fit) {
    parc[0] = 0.4/3.;  parc[1] = 0.03;  parc[2] = 0;  parc[3] = 0;  parc[4] = -0.04;  parc[5] = 0;
    parc[6] = 0.0004 ; parc[7] = 1.42; parc[8] = 0.14; 
    // fbEff -> SetParLimits(6, -0.2, 0.2);
    fcEff -> SetParLimits(7, 1.3, 1.7); 
    fcEff -> SetParLimits(8, 0., 0.3);
  }

  fcEff->SetParameters(parc);
  if (doFit) {
    cout << "**** Fitting hc ****" << endl;
    CheckHistoRanges(hc);
    hc->Fit(fcEff, FitOpt);
    if (JetFitNN57Fit) {
      for (int haha = 0; haha < 4; haha++)
	hc->Fit(fcEff, FitOpt);
    }
    PrintChi2(fcEff);
  }
  Double_t* covc = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covcEff = new TMatrixDSym(ncpars, covc);
  covcEff -> Print();
  
  TString lFitFormula = "[0] + [1]*x[1] + [2]*abs(x[0]) + [3]*pow(x[1],2.0) + [4]*pow(abs(x[0]),2.0) + [5]*x[1]*abs(x[0]) + [6]*pow(log(x[1]),4.0) + [7]*pow(abs(x[0]),4.0) + [8]*pow(x[1]*abs(x[0]),2.0)";

  if (IP3DSV1Fit)
    lFitFormula = "[0] + [1]*log(x[1]) + [2]*abs(x[0]) + [3]*pow(log(x[1]),2.0) + [4]*pow(abs(x[0]),2.0) + [5]*log(x[1])*abs(x[0]) + [6]*x[1] + [7]*pow(abs(x[0]),3.0) + [8]*pow(log(x[1])*abs(x[0]),2.0) + [9]*pow(abs(x[0]),4.0) + [10]*pow(x[1], 2.0) + [11]*x[1]*abs(x[0])";
  
  if (JetFitNN57Fit) {
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
  
  if (!JetFitNN57Fit) {
    parlight[7] = 0.0001;
    parlight[8] = 0.001;
  }
  if (IP3DSV1Fit && !JetFitNN57Fit) {
      parlight[9] = 0.000;
      parlight[10] = 0.000;
      parlight[11] = 0.000;
  }
  if (JetFitNN57Fit) {
    parlight[0] = 0.01;  parlight[1] = 0.03;  parlight[2] = 0;  parlight[3] = 0;  parlight[4] = -0.04;  parlight[5] = 0;
    parlight[6] = 0.0004 ; parlight[7] = 1.42; parlight[8] = 0.14; 
    flightEff -> SetParLimits(7, 1.3, 1.7); 
    flightEff -> SetParLimits(8, 0., 0.3);
  } 

  flightEff->SetParameters(parlight);
  if (doFit) {
    cout << "**** Fitting hlight ****" << endl;
    CheckHistoRanges(hlight);
    hlight->Fit(flightEff, FitOpt); // "V"
    if (JetFitNN57Fit) {
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

    
    // SF B
    TDirectory* bdir = 0;
    if (!odir -> cd("B")) {
      bdir = odir->mkdir("B");
      bdir->cd();
    }
    cout << "**** now in directory "  << endl; gDirectory->pwd();
    fbEff->Write("top_Eff");
    fbEffSyst->Write("top_Eff_syst");
    covbEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClonesEff(fbEff, fbEffSyst, covbEff);
    bdir->ls();
    
    //SF C
    TDirectory* cdir = 0;
    if (!odir -> cd("C")) {
      cdir  = odir->mkdir("C");
      cdir -> cd();
    }
    cout << "**** now in directory "  << endl; gDirectory->pwd();
    fcEff->Write("top_Eff");
    fcEffSyst->Write("top_Eff_syst");
    covcEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClonesEff(fcEff, fcEffSyst, covcEff);
    cdir->ls();

    // SF L
    TDirectory* lightdir = 0;
    if (!odir -> cd("Light")) {
      lightdir = odir->mkdir("Light");
      lightdir->cd();
    }
    cout << "**** now in directory "  << endl; gDirectory->pwd();
    TString helpname = "combined";
    flightEff->Write("top_Eff");
    flightEffSyst->Write("top_Eff_syst");
    covlightEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClonesEff2D(flightEff, flightEffSyst, covlightEff);
    lightdir->ls();
    
    fout->Close();
    delete fout;
  }

  // display results

  if (doDisplay) {
    
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
    if (OPName.Contains("4_55"))
      hc -> SetMaximum(0.25);
    if (OPName.Contains("-0_80"))
      hc -> SetMaximum(0.60);
    if (OPName.Contains("1_80") || OPName.Contains("2_20"))
      hc -> SetMaximum(0.20);
    if (OPName.Contains("0_35"))
      hc -> SetMaximum(0.25);
    if (OPName.Contains("1_70"))
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
    if (OPName.Contains("5_65"))
      hlight -> SetMaximum(0.025);
    if (OPName.Contains("4_55"))
      hlight -> SetMaximum(0.08);
    if (OPName.Contains("-0_80"))
      hlight -> SetMaximum(0.25);
    if (OPName.Contains("1_80") || OPName.Contains("2_20"))
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
}
