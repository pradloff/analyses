#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "TH2.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
// #include "T.h"


TH2* GetHisto2DFromFile(TFile *file, TString dirname, TString hname);

// ______________________________________________________
TH2* GetSpecialHisto(TFile *file, TString dirname, TString TaggerName, 
		     TString flavour, TString tag, TString truthtag = "truth") // eff/truth
{

  if (!file)
    return 0;
  TString name = "h_btag_AntiKt4TopoEM_" + TaggerName + "_" + flavour + truthtag + "_eta_Et" + tag;
  cout << "Getting histo " << name.Data()
       << " from dir " << dirname.Data()
       << " in file " << file -> GetName()
       << endl;

  TH2 *hb = GetHisto2DFromFile(file, dirname, name);
  if (hb)
    hb -> Sumw2();
  return hb;
}

// ______________________________________________________

void MakeAndWriteClones2D(TF2 *SelectionSF, TF2 *SelectionSFSyst, TMatrixDSym *covSF, TF2 *Eff, TF1 *EffSyst, TMatrixDSym *covEff)
{
  SelectionSF->Write("default_SF");
  SelectionSFSyst->Write("default_SF_syst");
  covSF->Write("default_SF_stat");
  Eff->Write("default_Eff");
  EffSyst->Write("default_Eff_syst");
  covEff->Write("default_Eff_stat");
}

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
void PrintChi2(TF1* fun)
{

  double chi2 = fun -> GetChisquare();
  double ndf = fun -> GetNDF();
  cout << "chi2=" << chi2
       << " ndf=" << ndf
       << " chi2/ndf=" << chi2/ndf
       << endl; 

}

// _______________________________________________________________________________
TH2* MakeCustomRatio(TH2* hnum, TH2* hdenom) 
{

  // code by Ford Garberson with modifications by jiri
  cout << "Making a custom ratio..." << endl;
  
  if (!hnum || !hdenom) {
    cerr << "Error getting num or denom!" << endl;
    return 0;
  }

  hnum->Sumw2();
  hdenom->Sumw2();

  TString efficiency_hist_name = hnum -> GetName();
  efficiency_hist_name += "_eff";

  TH2* heff = (TH2*) hnum->Clone(efficiency_hist_name);
  heff->Divide(hdenom);

  for(int bxi = 1; bxi <= heff->GetNbinsX(); bxi++) {
    for(int byi = 1; byi <= heff->GetNbinsY(); byi++) {
      float eff_val = heff->GetBinContent(bxi, byi);
      if( eff_val >= 1.e-6 && eff_val <= 1. ) {
        ///// everything is ok
      } else {
	heff->SetBinContent(bxi, byi, 0.);
	heff->SetBinError(bxi, byi, 2.);    
	/////// So it won't contribute to the fit
      } // else
    } // for
  } // for

  heff -> SetMaximum(1.2);

  return heff;

}

// _______________________________________________________________________________


void CheckHistoFitDiff2D(TH2* hist, TF2 *fun, double tolerance = 0.10) 
{

  if (!hist) {
    cerr << "CheckHistoFitDiff2D: Error getting hist!" << endl;
    return;  
  }
  if (!fun) {
    cerr << "CheckHistoFitDiff2D: Error getting fun!" << endl;
    return;  
  }

  int nBinsX = hist -> GetNbinsX();
  int nBinsY = hist -> GetNbinsY();
  
  int Nneg = 0;
  int Nbad = 0;
  int Nsig = 0;

  for (int j = 1; j <= nBinsY; ++j) {
    for (int i = 1; i <= nBinsX; ++i) {
      double val = hist -> GetBinContent(i, j);
      double x = hist -> GetXaxis() -> GetBinCenter(i);
      double y = hist -> GetYaxis() -> GetBinCenter(j);
      double fit = fun -> Eval(x, y);
      double err = hist -> GetBinError(i, j);
      
      cout << " i=" << i
	   << " j=" << j
	   << " x=" << x
	   << " y=" << y
	   << " val=" << val
	   << " fit=" << fit
	   << " err=" << err
	   << endl;

      if (val > 0.) {
	double ratio = TMath::Abs((val - fit) / val);
	if (ratio > tolerance)
	  Nbad++;
	if (err > 0) {
	  ratio = TMath::Abs((val - fit) / err);
	  if (ratio > 2.5)
	    Nsig++;
	}
      } else {
	if (val < -1.e-6)
	  Nneg++;
      }
    } // x
  } // y

  cout << "Results on " << hist -> GetName()
       << " " << fun -> GetName()
       << " Nneg=" << Nneg
       << " Nbad=" << Nbad
       << " Nsig=" << Nsig
       << endl;

}

// _______________________________________________________________________________
// _______________________________________________________________________________
// _______________________________________________________________________________

  void PrintCanvasAs(TCanvas *can, std::string suffix, bool SaveIndividual = false, std::string name = "")
  {

    if (name == "")
      can -> Print(Form("%s%s", can -> GetName(), suffix.c_str() ));
    else
      can -> Print(Form("%s%s", name.c_str(), suffix.c_str() ));

    if (SaveIndividual) {

      int k = 0;
      TPad *pad = 0;
      while(pad = (TPad*)can -> GetPad(++k)) {
      
	if (name == "")
	  pad -> Print(Form("%s_%i%s", can -> GetName(), k, suffix.c_str() ));
	else
	  pad -> Print(Form("%s_%i%s", name.c_str(), k, suffix.c_str() ));
      }
    
    } // save individual
  
  }
// _______________________________________________________________________________

TH2* GetHisto2DFromFile(TFile *file, TString dirname, TString hname)
{
  

  if (!file) {
    cout << "ERROR: GetHisto2DFromFile::Pointer to file is zero!" << endl;
    return 0;
  }
  // assert(file);
  
  TH2 *h = 0;
  h = (TH2*) file -> Get((dirname + hname).Data());

  /*
    if (h -> InheritsFrom("TH2F")) {
    TH2F *histo = (TH2F*) h;
    //    histo -> Draw("colz");
    } else {
    TH1F *histo = (TH1F*) h;
    //    histo -> Draw("e1hist");
    }
  */
  
  if (!h) 
    cout << "ERROR: GetHisto2DFromFile::Failed getting histogram " << (dirname+hname).Data() << endl;

  return h;

}
// _______________________________________________________________________________

void DrawSFs(TF1 *SF, TF1* SFSyst, TString flavour, TString tag, TString drawopt)
{
  TCanvas* can = new TCanvas("c" + flavour + "_SF_" + tag, "MC SF", 800, 400);
  can -> Divide(2, 1);
  can -> cd(1);
  SF->DrawCopy(drawopt);
  can -> cd(2);
  SFSyst->DrawCopy(drawopt);
  PrintCanvasAs(can, ".eps", true);
}

// _______________________________________________________________________________
// _______________________________________________________________________________
