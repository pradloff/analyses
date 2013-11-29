#ifndef ANP_CORRECTCALOISO_H
#define ANP_CORRECTCALOISO_H

///////////////////////////////////////////////////
// To get the corrected relative Isolation call:
//     CorrectEtConeRel(pt, isoConeRel, nvtx, eta)
//   - Author: Doug Schaefer <schae@cern.ch>
////////////////////////////////////////////////////
#include <string>

// Local
#include "MuonIsolationCorrection/IsoConfig.h"

class CorrectCaloIso {

 public:

  CorrectCaloIso();
  ~CorrectCaloIso();

  //
  // Returns the fully corrected isolation
  //
  float CorrectEtConeRel     (const float &pt, const float &isoConeRel, const float &nvtx, const float &eta, const std::string opts);
  float CorrectEtCone(const float &isoCone, const float &nvtx, const float &eta, const std::string opts);

  //
  // These functions return the etcone energy correction to the nominal etcone
  //  etconecor = etcone - GetCorrectionEtCone(nvtx);
  //
  float GetCorrectionEtCone(const unsigned &nvtx,const float &eta, const std::string opts);

 private:

  void Config(const std::string opts="");

 private:
  //
  // Declare Vars
  //
  bool        fDebug;
  
  IsoConfig   fConfigTool;
};

#endif
