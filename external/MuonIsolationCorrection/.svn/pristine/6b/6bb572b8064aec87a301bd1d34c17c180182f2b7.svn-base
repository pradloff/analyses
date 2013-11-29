// C/C++  
#include <iostream>
#include <math.h>

// Local
#include "MuonIsolationCorrection/CorrectCaloIso.h"

//////////////////////////////////////////////////
// To get the corrected relative Isolation call:
// CorrectEtConeRel(pt, isoConeRel, nvtx, eta, opts)
//    -- author: Doug Schaefer <schae@cern.ch>
//////////////////////////////////////////////////

using namespace std;

//-----------------------------------------------------------------------------
CorrectCaloIso::CorrectCaloIso() : fDebug(false)
{
  //
  // Configuring the eta dependent Isolation corrections
  //
  fConfigTool.Config();
  std::string conf = "cone30Comb";
  Config(conf);
}

//-----------------------------------------------------------------------------
CorrectCaloIso::~CorrectCaloIso()
{
}

//-----------------------------------------------------------------------------    
void CorrectCaloIso::Config(const std::string opts) 
{
  //
  // Setting the Isolation Correction values
  //
  fConfigTool.SetCorrections(opts);
}

//-----------------------------------------------------------------------------
float CorrectCaloIso::CorrectEtCone(const float &isoCone, const float &nvtx, const float &eta, const std::string opts)
{
  // 
  // Calculate corrected isolation variable 
  //    -- isoCone      is the EtCone
  //    -- nVtx         is number of primary vertices as calculated in the HSG3 recommendation
  //    -- eta          is the combined muon eta
  //         e.g. etconecor = (etcone - correction)
  //
  return isoCone - GetCorrectionEtCone(nvtx,eta,opts);
}

//-----------------------------------------------------------------------------
float CorrectCaloIso::CorrectEtConeRel(const float &pt, const float &isoConeRel, const float &nvtx, const float &eta, const std::string opts)
{
  // 
  // Calculate corrected relative isolation   
  //    -- pt           is in MeV and is the combined muon pt                 
  //    -- isoConeRel   is the EtCone/pt
  //    -- nVtx         is number of primary vertices as calculated in the HSG3 recommendation
  //    -- eta          is the combined muon eta
  //
  if(pt==0.0) return isoConeRel;
  float isoCone = isoConeRel*pt;

  //
  // Calculate the cone corrections  
  //   -- (etcone - correction)/pt
  //
  float newIsoCone = (isoCone - GetCorrectionEtCone(nvtx,eta,opts))/pt; 

  if(fDebug){
    std::cout << "Uncorrected Relative Isolation Cone :  " << isoConeRel
	      << " Pileup Corrected Relative Isolation:  " << newIsoCone 
	      << std::endl;
  }

  return newIsoCone;
}

//-----------------------------------------------------------------------------
float CorrectCaloIso::GetCorrectionEtCone(const unsigned &nvtx,const float &eta, const std::string opts)
{
  //
  // Returns the isolation correction
  //
  return fConfigTool.GetCaloCor(opts, nvtx, eta);
}

