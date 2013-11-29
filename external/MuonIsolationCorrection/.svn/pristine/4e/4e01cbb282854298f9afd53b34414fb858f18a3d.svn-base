#ifndef ANP_ISOCONFIG_H
#define ANP_ISOCONFIG_H

///////////////////////////////////////////////////
// Reads in isolation corrections
//   - Author: Doug Schaefer <schae@cern.ch>
////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <map>

class IsoConfig {

 public:

  IsoConfig();
  ~IsoConfig();

  //
  // Function to compute the isolation correction
  //
  void Config();

  float GetCaloCor(const std::string opts, const unsigned &nvtx, const float &eta);
  float CaloCor   (const unsigned &nvtx,const float &slope, const float &constant, const float &quad);

  void SetCorrections(const std::string opts="cone30Comb");

 private:

  struct ConeCorr {

    ConeCorr() {}

    ~ConeCorr() {}

    unsigned GetEtaBin(const float &eta) const;    

    std::vector<float> constant;
    std::vector<float> slope;
    std::vector<float> quad;
  };

 private:

  typedef std::map<std::string, ConeCorr> CorMap;

 private:

  ConeCorr ConfigEtCone20Comb();
  ConeCorr ConfigEtCone30Comb();
  ConeCorr ConfigEtCone40Comb();
  ConeCorr ConfigEtCone20MUID();
  ConeCorr ConfigEtCone30MUID();
  ConeCorr ConfigEtCone40MUID();

  ConeCorr ConfigEtCone20Comb2011();
  ConeCorr ConfigEtCone30Comb2011();
  ConeCorr ConfigEtCone40Comb2011();
  ConeCorr ConfigEtCone20Comb2012();
  ConeCorr ConfigEtCone30Comb2012();
  ConeCorr ConfigEtCone40Comb2012();

 private:
  //
  // Declare Vars
  //
  bool         fDebug;
  
  std::string  fOpts;     // Name of current correction

  ConeCorr    *fCor;      // Pointer to the current correction

  CorMap       fCorMap;   // Map of all corrections available
};

#endif
