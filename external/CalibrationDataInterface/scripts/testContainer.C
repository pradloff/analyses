// Either load the library before running ".x testContainer.C+" (i.e. compiled),
// then you'll need the include. Or run ".x testContainer.C" and call
// gSystem->Load("libCalibrationDataInterface.so"); in here - but then the #include
// is in the way: CINT *replaces* elements from the dictionary when loading
// the #include :-(

#ifndef __CINT__
#  include "CalibrationDataContainer.h"
#endif

#include "TMatrixDSym.h"
#include "TFile.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include <string>
#include <iostream>

using Analysis::CalibrationDataContainer;
using Analysis::CalibrationDataHistogramContainer;
using Analysis::CalibrationDataFunctionContainer;
using Analysis::CalibrationDataVariables;
using Analysis::UncertaintyResult;

void testContainer()
{
  // You seem to prefer histograms to behave like all other objects when it
  // comes to object lifetime, i.e. you want to delete them instead of TFile
  // deleting them:
  TH1::AddDirectory(kFALSE);

#ifdef __CINT__
  gSystem->Load("libCalibrationDataInterface.so");
#endif

  // example: open the PLHC file prepared by Jiri
  TFile* fin = TFile::Open("input.root");

  TF1* fB = 0;
  fin->GetObject("SV0/AntiKt4Topo/5_85/B/SV0_SF", fB);
  TF1* fBSyst = 0;
  fin->GetObject("SV0/AntiKt4Topo/5_85/B/SV0_SF_syst", fBSyst);
  TMatrixDSym* fBCov = 0;
  fin->GetObject("SV0/AntiKt4Topo/5_85/B/SV0_SF_stat", fBCov);

  TF1* fBEff = 0;
  fin->GetObject("SV0/AntiKt4Topo/5_85/B/top_Eff", fBEff);
  TF1* fBEffSyst = 0;
  fin->GetObject("SV0/AntiKt4Topo/5_85/B/top_Eff_syst", fBEffSyst);
  TMatrixDSym* fBEffCov = 0;
  fin->GetObject("SV0/AntiKt4Topo/5_85/B/top_Eff_stat", fBEffCov);
  std::cout << "ptrs: fBEff: " << fBEff << ", fBEffSyst: " << fBEffSyst << ", fBEffCov " << fBEffCov << std::endl;

  // bin boundaries for the negative-tag measurements
  Double_t xMisTag[9] = { 20., 25., 40., 60., 90., 140., 200., 300., 500. };
  Double_t yMisTag[3] = { 0.0, 1.2, 2.5 };
  // calibration result for the SV0 5.85 working point
  Double_t valMisTag[2][8] = { { 1.023, 0.891, 0.935, 0.864, 0.866, 0.757, 0.917, 0.888 },
			       { 1.340, 0.877, 1.200, 0.914, 1.024, 0.984, 0.686, 1.044 } };
  Double_t errMisTag[2][8] = { { 0.206, 0.185, 0.182, 0.196, 0.140, 0.142, 0.126, 0.120 },
			       { 0.450, 0.274, 0.275, 0.287, 0.257, 0.244, 0.114, 0.168 } };
  TH2F* hLight = new TH2F("negtag_SF", "SV0 (w>5.85) mistag rate;pt;abseta", 8, xMisTag, 2, yMisTag);
  for (int ix = 1; ix < 9; ++ix)
    for (int iy = 1; iy < 3; ++iy) {
      hLight->SetBinContent(ix, iy, valMisTag[iy-1][ix-1]);
      hLight->SetBinError  (ix, iy, errMisTag[iy-1][ix-1]);
    }
  TH1F* hLightSys = new TH1F("negtag_SF_syst", "dummy corresponding scale factor;pt", 8, xMisTag);
  //  hLight->Draw("COLZ");

  // not sure this can be done here already...
  delete fin;

  // create the new containers using the above information

  CalibrationDataFunctionContainer* cfSF = new CalibrationDataFunctionContainer("ptrel_SF");
  cfSF->setResult(fB);
  cfSF->setUncertainty("statistics", fBCov);
  cfSF->setUncertainty("systematics", fBSyst);
  cfSF->setComment("dummy comment for scale factor object\n\t...and this is a second line");
  cfSF->setLowerBound(CalibrationDataContainer::kPt, 25.0);
  cfSF->setUpperBound(CalibrationDataContainer::kPt, 200.0);
  cfSF->setLowerBound(CalibrationDataContainer::kAbsEta, 0.0);
  cfSF->setUpperBound(CalibrationDataContainer::kAbsEta, 2.5);

  CalibrationDataFunctionContainer* cfEff = new CalibrationDataFunctionContainer("top_Eff");
  cfEff->setResult(fBEff);
  cfEff->setUncertainty("statistics", fBEffCov);
  cfEff->setUncertainty("systematics", fBEffSyst);
  cfEff->setComment("dummy comment for efficiency object");

  CalibrationDataHistogramContainer* chSF = new CalibrationDataHistogramContainer("negtag_SF");
  // ensure that no rebinning is performed to these histograms later
  hLight->ResetBit(TH1::kCanRebin);
  hLightSys->ResetBit(TH1::kCanRebin);
  chSF->setResult(hLight);
  chSF->setUncertainty("systematics", hLightSys);
  chSF->setComment("test comment for custom mistag rate scale object");
  // The statistical uncertainties are expected to be stored with the central values.
  // However, for the purpose of "continuous" b-tagging, the covariance matrix corresponding to
  // the statistical uncertainties becomes nontrivial (in particular, non-diagonal).
  // In this case we store the matrix _in addition to the bin-by-bin uncertainties_.
  // The example below fills a dummy covariance matrix just to show the mechanics
  // (the container is anyway not for the storage of a continuous calibration object).
  Int_t nbx = hLight->GetNbinsX()+2, nby = hLight->GetNbinsY()+2;
  TMatrixDSym* sCov = new TMatrixDSym(nbx*nby); // make sure to keep track of the under- and overflow bins
  for (Int_t binx = 1; binx < nbx-1; ++binx)
    for (Int_t biny = 1; biny < nby-1; ++biny) {
      Int_t bin = ref->GetBin(binx, biny);
      for (Int_t binx2 = 1; binx2 < nbx-1; ++binx2)
	for (Int_t biny2 = 1; biny2 < nby-1; ++biny2) {
	  Int_t bin2 = ref->GetBin(binx2, biny2, binz2);
	  // here it's assumed that the covariance matrix elements have been computed; they would replace the 0 in the following line
	  sCov(bin,bin2) = 0;
	}
    }
  chSF->SetUncertainty("statistics", sCov); // the keyword "statistics" is mandatory (when the statistical uncertainties' covariance matrix is concerned)
  chSF->restrictToRange(true);

  // test to see what happens if first computation is done here...
  CalibrationDataVariables x1;
  x1.jetAuthor = std::string("AntiKt4Topo");
  x1.jetPt     = 18000.;
  x1.jetEta    = 1.1;
  double resTmp; CalibrationDataContainer::CalibrationStatus statTmp = chSF->getResult(x1, resTmp);
  if (statTmp != CalibrationDataContainer::kError) {
    std::cout << "\t hSF result: " << resTmp << ", status = " << statTmp << std::endl;
  } else std::cout << "\t error retrieving hSF result" << std::endl;

  
  // create a new ROOT file
  TString fname("containerFile.root");
  TFile* f = TFile::Open(fname, "RECREATE");
  if (f->IsZombie()) {
    std::cout << "error: cannot create file " << fname << std::endl;
    delete f;
    return;
  }

  f->WriteTObject(cfSF,"ptrel");
  //  cfSF->dummy();
  f->WriteTObject(cfEff);
  //  cfEff->dummy();
  f->WriteTObject(chSF);
  //  chSF->dummy();
  TH1F* hh = new TH1F("hh", "test 1D histogram; example abscissa", 50, 0.0, 5.0);
  hh->Write();
  f->ls();
  // close the file, and delete all objects
  delete cfSF;
  delete cfEff;
  delete chSF;
  delete hh;
  delete f;

  // now re-open the ROOT file and retrieve the objects by
  f = TFile::Open(fname);
  std::cout << "after re-opening the file: " << std::endl;
  f->ls();
  // Type-safe version without (slow) dynamic_cast:
  CalibrationDataContainer* cSF = 0;
  f->GetObject("ptrel", cSF);
  CalibrationDataContainer* cmSF = 0;
  f->GetObject("negtag_SF", cmSF);
  CalibrationDataContainer* cEff = 0;
  f->GetObject("top_Eff", cEff);
  // if (cSF)  cSF->dummy();  else std::cout << " error: null pointer for object cSF" << std::endl;
  // if (cmSF) cmSF->dummy(); else std::cout << " error: null pointer for object cmSF" << std::endl;
  if (!cEff)  std::cout << " error: null pointer for object cEff" << std::endl;
  cSF->restrictToRange(true);
  cmSF->restrictToRange(true);

  double res1; CalibrationDataContainer::CalibrationStatus stat1 = cSF->getResult(x1, res1);
  if (stat1 != CalibrationDataContainer::kError) {
    std::cout << "\t SF result: " << res1 << ", status = " << stat1 << std::endl;
  } else std::cout << "\t error retrieving SF result" << std::endl;
  double statUnc; CalibrationDataContainer::CalibrationStatus stat2 = cSF->getStatUncertainty(x1,statUnc);
  if (stat2 != CalibrationDataContainer::kError) {
    std::cout << "\t SF stat unc: " << statUnc << ", status = " << stat2 << std::endl;
  } else std::cout << "\t error retrieving SF stat unc" << std::endl;
  UncertaintyResult systUnc; CalibrationDataContainer::CalibrationStatus stat3 = cSF->getSystUncertainty(x1,systUnc);
  if (stat3 != CalibrationDataContainer::kError) {
    std::cout << "\t SF syst unc: (+)" << systUnc.first << " (-)" << systUnc.second << ", status = " << stat3 << std::endl;
  } else std::cout << "\t error retrieving SF syst unc" << std::endl;
  std::cout << "\t SF comment: " << cSF->getComment() << std::endl;
  for (unsigned int i = 0; i < 3; ++i) {
    std::cout << "\tboundaries for variable " << i << ": (-)" << cSF->getLowerBound(i)
	      << " (+)" << cSF->getUpperBound(i) << std::endl;
  }

  double resE1; CalibrationDataContainer::CalibrationStatus statE1 = cEff->getResult(x1, resE1);
  if (statE1 != CalibrationDataContainer::kError) {
    std::cout << "\t Eff result: " << resE1 << ", status = " << statE1 << std::endl;
  } else std::cout << "\t error retrieving Eff result" << std::endl;
  double statEUnc; CalibrationDataContainer::CalibrationStatus statE2 = cEff->getStatUncertainty(x1,statEUnc);
  if (statE2 != CalibrationDataContainer::kError) {
    std::cout << "\t Eff stat unc: " << statEUnc << ", status = " << statE2 << std::endl;
  } else std::cout << "\t error retrieving Eff stat unc" << std::endl;
  UncertaintyResult systEUnc; CalibrationDataContainer::CalibrationStatus statE3 = cEff->getSystUncertainty(x1,systEUnc);
  if (statE3 != CalibrationDataContainer::kError) {
    std::cout << "\t Eff syst unc: (+)" << systEUnc.first << " (-)" << systEUnc.second << ", status = " << statE3 << std::endl;
  } else std::cout << "\t error retrieving Eff syst unc" << std::endl;
  std::cout << "\t SF comment: " << cEff->getComment() << std::endl;
  for (unsigned int i = 0; i < 3; ++i) {
    std::cout << "\tboundaries for variable " << i << ": (-)" << cEff->getLowerBound(i)
	      << " (+)" << cEff->getUpperBound(i) << std::endl;
  }

  for (unsigned int i = 0; i < 3; ++i) {
    std::cout << "\tboundaries for variable " << i << ": (-)" << cmSF->getLowerBound(i)
	      << " (+)" << cmSF->getUpperBound(i) << std::endl;
  }
  double resm1; CalibrationDataContainer::CalibrationStatus statm1 = cmSF->getResult(x1, resm1);
  if (statm1 != CalibrationDataContainer::kError) {
    std::cout << "\t mSF result: " << resm1 << ", status = " << statm1 << std::endl;
  } else std::cout << "\t error retrieving mSF result" << std::endl;
  double statmUnc; CalibrationDataContainer::CalibrationStatus statm2 = cmSF->getStatUncertainty(x1,statmUnc);
  if (statm2 != CalibrationDataContainer::kError) {
    std::cout << "\t mSF stat unc: " << statmUnc << ", status = " << statm2 << std::endl;
  } else std::cout << "\t error retrieving mSF stat unc" << std::endl;
  UncertaintyResult systmUnc; CalibrationDataContainer::CalibrationStatus statm3 = cmSF->getSystUncertainty(x1,systmUnc);
  if (statm3 != CalibrationDataContainer::kError) {
    std::cout << "\t mSF syst unc: (+)" << systmUnc.first << " (-)" << systmUnc.second << ", status = " << statm3 << std::endl;
  } else std::cout << "\t error retrieving mSF syst unc" << std::endl;
  std::cout << "\t mSF comment: " << cmSF->getComment() << std::endl;
  for (unsigned int i = 0; i < 3; ++i) {
    std::cout << "\tboundaries for variable " << i << ": (-)" << cmSF->getLowerBound(i)
	      << " (+)" << cmSF->getUpperBound(i) << std::endl;
  }

  std::cout << "full printout for mSF:" << std::endl;
  std::map<std::string, UncertaintyResult> allm;
  if (cmSF->getUncertainties(x1,allm) != CalibrationDataContainer::kError) {
    for (std::map<std::string, UncertaintyResult>::iterator it = allm.begin(); it != allm.end(); ++it)
      std::cout << "\t" << it->first << ": (+)" << it->second.first << " (-)" << it->second.second << std::endl;
  } else
    std::cout << "\terror retrieving results!" << std::endl;
  std::cout << "full printout for cSF:" << std::endl;
  std::map<std::string, UncertaintyResult> allc;
  if (cSF->getUncertainties(x1,allc) != CalibrationDataContainer::kError) {
    for (std::map<std::string, UncertaintyResult>::iterator it = allc.begin(); it != allc.end(); ++it)
      std::cout << "\t" << it->first << ": (+)" << it->second.first << " (-)" << it->second.second << std::endl;
  } else
    std::cout << "\terror retrieving results!" << std::endl;
  std::cout << "full printout for cEff:" << std::endl;
  std::map<std::string, UncertaintyResult> allE;
  if (cEff->getUncertainties(x1,allE) != CalibrationDataContainer::kError) {
    for (std::map<std::string, UncertaintyResult>::iterator it = allE.begin(); it != allE.end(); ++it)
      std::cout << "\t" << it->first << ": (+)" << it->second.first << " (-)" << it->second.second << std::endl;
  } else
    std::cout << "\terror retrieving results!" << std::endl;

  delete cSF;
  delete cmSF;
  delete cEff;
  delete f;
}
