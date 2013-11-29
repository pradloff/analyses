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




// double CombinedUncertainty(Double_t* xx, Double_t* par);
double LightFlavourSF(Double_t* xx, Double_t* par);
double LightFlavourSFSyst(Double_t* xx, Double_t* par);
double HeavyFlavourSF(Double_t* xx, Double_t* par);
double HeavyFlavourSFSyst(Double_t* xx, Double_t* par);

// ______________________________________________________

void MakeAndWriteClones(TF1 *SelectionSF, TF1 *SelectionSFSyst, TMatrixDSym *covSF, TF2 *Eff, TF1 *EffSyst, TMatrixDSym *covEff)
{
  SelectionSF->Write("default_SF");
  SelectionSFSyst->Write("default_SF_syst");
  covSF->Write("default_SF_stat");
  Eff->Write("default_Eff");
  EffSyst->Write("default_Eff_syst");
  covEff->Write("default_Eff_stat");
}

// ______________________________________________________
// ______________________________________________________
// ______________________________________________________


void fitTopHists(bool cloneToDefault = true, bool doWrite = true, bool doDisplay = true)
{
  TFile* f = TFile::Open("EffMapInput.root");

  // retrieve scale factor inputs resulting from ptrel fits

  // TFile* fin = new TFile("pTrel.root");
  // TF1* fbSelectionSF = (TF1*) fin->Get("SFFunction");
  // fbSFstatin = (TF1*) fin->Get("StatErrFunction");
  // fbSFsystin = (TF1*) fin->Get("SystErrFunction");

  // retrieve the histograms
  TH2* hb = (TH2*) f->Get("hMCTaggingEff_b_SV0_5_72");
  TH2* hc = (TH2*) f->Get("hMCTaggingEff_c_SV0_5_72");
  TH2* hlight = (TH2*) f->Get("hMCTaggingEff_light_SV0_5_72");

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2.0);

  // first attempt: polynomials in (ln(pt), eta)
  TF2* fbEff = new TF2("fbEff",
		       "[0] + [1]*log(x[0]) +[2]*x[1] + [3]*pow(log(x[0]),2.0) + [4]*pow(x[1],2.0) + [5]*log(x[0])*x[1]",
		       20.0, 200.0, 0.0, 2.5);
  Double_t parb[6];
  parb[0] = 0.4;
  parb[1] = 0.03;
  parb[2] = 0;
  parb[3] = 0;
  parb[4] = -0.04;
  parb[5] = 0;
  fbEff->SetParameters(parb);
  hb->Fit(fbEff, "IN");

  Double_t* covb = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covbEff = new TMatrixDSym(6, covb);
  covbEff -> Print();

  TF2* fcEff = new TF2("fcEff",
		       "[0] + [1]*log(x[0]) +[2]*x[1] + [3]*pow(log(x[0]),2.0) + [4]*pow(x[1],2.0) + [5]*log(x[0])*x[1]",
		       20.0, 200.0, 0.0, 2.5);
  Double_t parc[6];
  parc[0] = 0.1;
  parc[1] = 0.01;
  parc[2] = 0;
  parc[3] = 0;
  parc[4] = -0.01;
  parc[5] = 0;
  fcEff->SetParameters(parc);
  hc->Fit(fcEff, "IN");

  Double_t* covc = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covcEff = new TMatrixDSym(6, covc);
  covcEff -> Print();

  TF2* flightEff = new TF2("flightEff",
			   "[0] + [1]*x[0] + [2]*x[1] + [3]*pow(x[0],2.0) + [4]*pow(x[1],2.0) + [5]*x[0]*x[1] + [6]*pow(x[0],4.0) + [7]*pow(x[1],4.0) + [8]*pow(x[0]*x[1],2.0)",
			   20.0, 200.0, 0.0, 2.5);
  Double_t parlight[9];
  parlight[0] = 0.04;
  parlight[1] = 0;
  parlight[2] = 0;
  parlight[3] = 0;
  parlight[4] = 0;
  parlight[5] = 0;
  parlight[6] = 0;
  parlight[7] = 0;
  parlight[8] = 0;
  flightEff->SetParameters(parlight);
  hlight->Fit(flightEff, "IN");

  Double_t* covlight = TVirtualFitter::GetFitter()->GetCovarianceMatrix();
  TMatrixDSym* covlightEff = new TMatrixDSym(9, covlight);
  covlightEff -> Print();

  fbEff->SetTitle("b-jet selection efficiency;pt;abseta");
  fcEff->SetTitle("c-jet selection efficiency;pt;abseta");
  flightEff->SetTitle("light-jet selection efficiency;pt;abseta");

  // dummy systematic uncertainty parametrisations
  TF1* fbEffSyst = new TF1("fbEffSyst", "[0]", 20.0, 200.0);
  fbEffSyst->SetTitle("dummy b-jet selection efficiency syst uncertainty");
  Double_t parbEffSyst[1];
  parbEffSyst[0] = 0.0;
  fbEffSyst->SetParameters(parbEffSyst);

  TF1* fcEffSyst = new TF1("fcEffSyst", "[0]", 20.0, 200.0);
  fcEffSyst->SetTitle("dummy c-jet selection efficiency syst uncertainty");
  Double_t parcEffSyst[1];
  parcEffSyst[0] = 0.0;
  fcEffSyst->SetParameters(parcEffSyst);

  TF1* flightEffSyst = new TF1("flightEffSyst", "[0]", 20.0, 200.0);
  flightEffSyst->SetTitle("dummy light-jet selection efficiency syst uncertainty");
  Double_t parlightEffSyst[1];
  parlightEffSyst[0] = 0.0;
  flightEffSyst->SetParameters(parlightEffSyst);

  // also provide (approximate) scale factor results

  // The TF1's above were provided in sampled form, which implies that it's no use
  // attempting to use a covariance matrix to compute the statistical uncertainties.
  // So instead, we lump the statistical and systematic uncertainties together, into
  // an overall "systematic" uncertainty
  TF1* fbSelectionSF = new TF1("bSFFunction", HeavyFlavourSF, 25.0, 200.0, 0);
  fbSelectionSF->SetTitle("ptrel-based b-jet selection efficiency SF;pt");
  Double_t dummy = 0; 
  Double_t* dummyCov = &dummy;
  TMatrixDSym* covbSF = new TMatrixDSym(0, dummyCov);
  TF1* fbSelectionSFSyst = new TF1("fbSelectionSFSyst", HeavyFlavourSFSyst, 25.0, 200.0, 1);
  fbSelectionSFSyst->SetParameter(0, 1.0); // scale to inflate the syst \oplus stat
  fbSelectionSFSyst->SetTitle("ptrel-based b-jet SF systematic uncertainty;pt");

  // c-jet efficiency scale factor: as for b-jets, only with the uncertainty doubled.
  TF1* fcSelectionSF = (TF1*)fbSelectionSF->Clone("cSFFunction");
  fcSelectionSF->SetTitle("ptrel-based c-jet selection efficiency SF;pt");
  TMatrixDSym* covcSF = new TMatrixDSym(0, dummyCov);
  TF1* fcSelectionSFSyst = new TF1("fcSelectionSFSyst", HeavyFlavourSFSyst, 25.0, 200.0, 1);
  fcSelectionSFSyst->SetParameter(0, 2.0); // c syst 2x of the b syst;-)
  fcSelectionSFSyst->SetTitle("ptrel-based c-jet SF systematic uncertainty;pt");

  // light-jet efficiency scale factor: similar fashion
  TF1* flightSelectionSF = new TF1("flightSelectionSF", LightFlavourSF, 20.0, 200.0, 0);
  flightSelectionSF->SetTitle("light-jet selection efficiency SF;pt");
  // Double_t parlightSF[1];
  // parlightSF[0] = 1.15;
  // flightSelectionSF->SetParameters(parlightSF);
  // Double_t err2lightSF = TMath::Power(parlightSF[0] * 0.50, 2);
  TMatrixDSym* covlightSF = new TMatrixDSym(0, dummyCov);
  TF1* flightSelectionSFSyst = new TF1("flightSelectionSFSyst", LightFlavourSFSyst, 20.0, 200.0, 0);
  flightSelectionSFSyst->SetTitle("preliminary light-jet selection efficiency SF syst uncertainty;pt");
  // Double_t parlightSFSyst[1];
  // parlightSFSyst[0] = 0.0;
  // flightSelectionSFSyst->SetParameters(parlightSFSyst);

  // write to output file

  if (doWrite) {

    TFile* fout = new TFile("TopCalibrations_new.root", "recreate");
    TDirectory* tdir = fout->mkdir("SV0");
    tdir->cd();
    TDirectory* jdir = tdir->mkdir("AntiKt4Topo");
    jdir->cd();
    TDirectory* odir = jdir->mkdir("5_72");
    odir->cd();


    TDirectory* bdir = odir->mkdir("B");
    bdir->cd();
    std::cout << "now in directory "  << std::endl; gDirectory->pwd();
    fbSelectionSF->Write("ptrel_SF");
    fbSelectionSFSyst->Write("ptrel_SF_syst");
    covbSF->Write("ptrel_SF_stat");
    fbEff->Write("top_Eff");
    fbEffSyst->Write("top_Eff_syst");
    covbEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClones(fbSelectionSF, fbSelectionSFSyst, covbSF, fbEff, fbEffSyst, covbEff);
    bdir->ls();

    TDirectory* cdir = odir->mkdir("C");
    cdir->cd();
    std::cout << "now in directory "  << std::endl; gDirectory->pwd();
    fcSelectionSF->Write("ptrel_SF");
    fcSelectionSFSyst->Write("ptrel_SF_syst");
    covcSF->Write("ptrel_SF_stat");
    fcEff->Write("top_Eff");
    fcEffSyst->Write("top_Eff_syst");
    covcEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClones(fcSelectionSF, fcSelectionSFSyst, covcSF, fcEff, fcEffSyst, covcEff);
    cdir->ls();

    TDirectory* lightdir = odir->mkdir("Light");
    lightdir->cd();
    std::cout << "now in directory "  << std::endl; gDirectory->pwd();
    //    TString helpname = "preliminary";
    TString helpname = "combined";
    flightSelectionSF->Write(helpname + "_SF");
    flightSelectionSFSyst->Write(helpname + "_SF_syst");
    covlightSF->Write(helpname + "_SF_stat");
    flightEff->Write("top_Eff");
    flightEffSyst->Write("top_Eff_syst");
    covlightEff->Write("top_Eff_stat");
    if (cloneToDefault) 
      MakeAndWriteClones(flightSelectionSF, flightSelectionSFSyst, covlightSF, flightEff, flightEffSyst, covlightEff);
    lightdir->ls();

    fout->Close();
    delete fout;
  }

  // fin->Close();
  // delete fin;

  // display results

  if (doDisplay) {
    TCanvas* cb = new TCanvas("cb", "MC b-jet efficiency", 600, 400);
    hb->DrawCopy("LEGO");

    // cross-check: divide the histogram by the fit
    TCanvas* cbCheck = new TCanvas("cbCheck", "MC b-jet efficiency cross-check", 600, 400);
    TH2* hbCheck = (TH2*)hb->Clone("hbCheck");
    hbCheck->Divide(fbEff);
    hbCheck->SetMinimum(0.93);
    hbCheck->SetMaximum(1.05);
    hbCheck->DrawCopy("HIST COLZ");

    TCanvas* cc = new TCanvas("cc", "MC c-jet efficiency", 600, 400);
    hc->DrawCopy("LEGO");

    // cross-check: divide the histogram by the fit
    TCanvas* ccCheck = new TCanvas("ccCheck", "MC c-jet efficiency cross-check", 600, 400);
    TH2* hcCheck = (TH2*)hc->Clone("hcCheck");
    hcCheck->Divide(fcEff);
    // hcCheck->SetMinimum(0.93);
    // hcCheck->SetMaximum(1.05);
    hcCheck->DrawCopy("HIST COLZ");

    TCanvas* clight = new TCanvas("clight", "MC light-jet efficiency", 600, 400);
    hlight->DrawCopy("LEGO");

    // cross-check: divide the histogram by the fit
    TCanvas* clightCheck = new TCanvas("clightCheck", "MC light-jet efficiency cross-check", 600, 400);
    TH2* hlightCheck = (TH2*)hlight->Clone("hlightCheck");
    hlightCheck->Divide(flightEff);
    // hlightCheck->SetMinimum(0.93);
    // hlightCheck->SetMaximum(1.05);
    hlightCheck->DrawCopy("HIST COLZ");

    // for PR purposes: plots of the parametrisations themselves
    TCanvas* cbParm = new TCanvas("cbParm", "MC b-jet efficiency parametrisation", 600, 400);
    fbEff->SetTitle(";jet p_{T} (GeV);jet |#eta|;b-jet tagging efficiency");
    fbEff->GetZaxis()->SetTitleOffset(1.2);
    fbEff->Draw("SURF1");
    TCanvas* ccParm = new TCanvas("ccParm", "MC c-jet efficiency parametrisation", 600, 400);
    fcEff->SetTitle(";jet p_{T} (GeV);jet |#eta|;c-jet tagging efficiency");
    fcEff->GetZaxis()->SetTitleOffset(1.2);
    fcEff->Draw("SURF1");
    TCanvas* clParm = new TCanvas("clParm", "MC light-flavour jet efficiency parametrisation", 600, 400);
    flightEff->SetTitle(";jet p_{T} (GeV);jet |#eta|;light-flavour jet tagging efficiency");
    flightEff->GetZaxis()->SetTitleOffset(1.2);
    flightEff->Draw("SURF1");
  }

  f->Close();
  delete f;
}

// void CombinedUncertainty(Double_t* xx, Double_t* par) {
//   double stat = fbSFstatin->Eval(xx[0]);
//   double syst = fbSFsystin->Eval(xx[0]);
//   return TMath::Sqrt(stat*stat+syst*syst)*par[0];
// }

double LightFlavourSF(Double_t* xx, Double_t* par) {
  double pt = xx[0];
  return (pt < 40.0) ? 1.29 : 1.08;
}

double LightFlavourSFSyst(Double_t* xx, Double_t* par) {
  //  double pt = xx[0];
  return 0.50*LightFlavourSF(xx, par);
  // orig. code:
  // return (pt < 40.0) ? 0.26 : 0.25;
}

double HeavyFlavourSF(Double_t* xx, Double_t* par) {
  double pt = xx[0];
  if (pt < 40.0)
    return 1.00;
  else if (pt < 60.0)
    return 0.88;
  else
    return 1.05;
}

double HeavyFlavourSFSyst(Double_t* xx, Double_t* par) {
  double pt = xx[0];
  double scale = par[0];
  double stat = 0., syst = 0.;

  syst = 0.08*HeavyFlavourSF(xx, par);

  /* orig code:
  if (pt < 40.0) {
    stat = 0.03; syst = 0.12;
  } else if (pt < 60.0) {
    stat = 0.04; syst = 0.09;
  } else if (pt < 85.0) {
    stat = 0.11; syst = 0.10;
  } else {
    stat = 0.11; syst = 0.20;
  }
  */
  return TMath::Sqrt(stat*stat+syst*syst)*scale;
}
