/*
  TESUncertaintyProvider Winter 2012
---------------------------------------------
  Providing TES uncertainties mentioned on
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TauSystematicsWinterConf2012

  based on the uncertainty estimation using single particle uncertainties by
    - Matteo Volpi <mvolpi@ifae.es>
    - Nils Ruthmann <ruthmann@cern.ch>
*/

#ifndef _TESUNCERTAINTYPROVIDER_H_
#define _TESUNCERTAINTYPROVIDER_H_

#include <vector>
#include <math.h>
#include "stdio.h"
#include <iostream>
class TESUncertaintyProvider {

 public:
  TESUncertaintyProvider();
  ~TESUncertaintyProvider();

  /*
    Provide pt in GeV!
    Returns -1. in case of failure.
    And the uncertainty in % when everything went fine
   */ 
  double GetTESUncertainty(double pt, double eta, int nProng);
  
 private:
  std::vector<double> m_eta_bins;
  std::vector<double> m_pt_bins;
 
  std::vector< std::vector<double> > m_1p_unc;
  std::vector< std::vector<double> > m_mp_unc;

  int m_nEta;
  int m_nPt;
};

#endif
