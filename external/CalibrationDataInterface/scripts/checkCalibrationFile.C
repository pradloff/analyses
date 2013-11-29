//
// File     : checkCalibrationFile.C
// Author   : Frank Filthaut
// Purpose  : inspect the calibration ROOT file (new format using CalibrationDataContainers)

// Either load the library before running ".x checkCalibrationFile.C+" (i.e. compiled),
// then you'll need the include. Or run ".x checkCalibrationFile.C" and call
// gSystem->Load("libCalibrationDataInterface.so"); in here - but then the #include
// is in the way: CINT *replaces* elements from the dictionary when loading
// the #include :-(

// One easy way to run this macro is within the RootCore build environment:
// root -l $ROOTCOREDIR/scripts/load_packages.C+
//   .L checkCalibrationFile.C+
//   checkCalibrationFile()

#include <iostream>
#include <vector>
#include <utility>
#include <string.h>
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TMatrixDSym.h"
#include "TDirectory.h"
#include "TCollection.h"
#include "TKey.h"
#include "TString.h"
#include "TClass.h"

using std::cout;
using std::endl;

#ifndef __CINT__
#  include "CalibrationDataInterface/CalibrationDataContainer.h"
#endif

using Analysis::CalibrationDataContainer;
using Analysis::CalibrationDataFunctionContainer;
using Analysis::CalibrationDataHistogramContainer;
using Analysis::CalibrationDataMappedHistogramContainer;

void checkDirectory(int lvl, TString& search);

void checkCalibrationFile(TString fileName = "TopCalibrations_rel17_MC11a.root", int lvl = 0,
			  TString search = "") {

#ifdef __CINT__
  gSystem->Load("libCalibrationDataInterface.so");
#endif

  TFile* f = TFile::Open(fileName.Data(), "READ");
  f->cd();
  checkDirectory(lvl, search);
  f->Close();
}


void checkDirectory(int lvl, TString& search) {

  // TString path( (char*)strstr( dir->GetPath(), ":" ) );
  // path.Remove( 0, 2 );
  // dir->cd();
  TDirectory *current_dir = gDirectory;
  
  // loop over all keys in this directory

  TIter nextkey( current_dir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    TObject *obj = key->ReadObj();

    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      // it's a subdirectory
      if (lvl > 3) cout << "Found subdirectory " << obj->GetName() << endl;
      current_dir->cd(obj->GetName());
      checkDirectory(lvl, search);
    }
    else if (obj->IsA()->InheritsFrom(CalibrationDataContainer::Class())) {
    // else if (! strcmp(obj->GetName(), "top_Eff")) {

      if (obj->IsA()->InheritsFrom(CalibrationDataFunctionContainer::Class())) {
	if (search == "")
	  cout << "found CalibrationDataFunctionContainer object "
	       << current_dir->GetPath() << "/" << key->GetName() << endl;
	else {
	  TString path(current_dir->GetPath()); path += "/"; path += key->GetName();
	  if (path.Contains(search))
	    cout << "found CalibrationDataFunctionContainer object " << path << endl;
	}
      }
      else if (obj->IsA()->InheritsFrom(CalibrationDataMappedHistogramContainer::Class())) {
	if (search == "")
	  cout << "found CalibrationDataMappedHistogramContainer object "
	       << current_dir->GetPath() << "/" << key->GetName() << endl;
	else {
	  TString path(current_dir->GetPath()); path += "/"; path += key->GetName();
	  if (path.Contains(search))
	    cout << "found CalibrationDataMappedHistogramContainer object " << path << endl;
	}
      }
      else if (obj->IsA()->InheritsFrom(CalibrationDataHistogramContainer::Class())) {
	if (search == "")
	  cout << "found CalibrationDataHistogramContainer object "
	       << current_dir->GetPath() << "/" << key->GetName() << endl;
	else {
	  TString path(current_dir->GetPath()); path += "/"; path += key->GetName();
	  if (path.Contains(search))
	    cout << "found CalibrationDataHistogramContainer object " << path << endl;
	}
      }
      else {
	cout << "\tUNKNOWN container type!" << endl;
	continue;
      }
      CalibrationDataContainer* cnt = dynamic_cast<CalibrationDataContainer*>(obj);
      if (lvl > 0) {
	std::string onoff = cnt->isRangeRestricted() ? "on " : "off";
	cout << "\t(1)range restriction switched " << onoff << endl;
	cout << "\t(1)validity bounds: ";
	std::vector<std::pair<double, double> > bounds = cnt->getBounds();
	for (unsigned int t = 0; t < bounds.size(); ++t)
	  cout << " " << t << "(" << bounds[t].first << "," << bounds[t].second << ")";
	cout << endl;
      }
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )
}
