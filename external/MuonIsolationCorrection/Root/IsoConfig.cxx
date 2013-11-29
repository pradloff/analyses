// C/C++  
#include <iostream>
#include <math.h>

// Local
#include "MuonIsolationCorrection/IsoConfig.h"

//////////////////////////////////////////////////
//    -- author: Doug Schaefer <schae@cern.ch>
//////////////////////////////////////////////////

using namespace std;

//-----------------------------------------------------------------------------    
unsigned IsoConfig::ConeCorr::GetEtaBin(const float &eta) const
{
  //                           
  // Find Eta bin
  //
  unsigned bin = static_cast<unsigned>(fabs(eta)/0.10);

  if(bin<27)  return bin;
  else{
    if(false) std::cout << "ERROR bin is bigger than known isolation corrections!" << std::endl;
    return 26;
  }
}

//-----------------------------------------------------------------------------
IsoConfig::IsoConfig() : fDebug(false), fOpts("cone30Comb"), fCor(0)
{
}

//-----------------------------------------------------------------------------
IsoConfig::~IsoConfig()
{
}

//-----------------------------------------------------------------------------
void IsoConfig::Config()
{
 //
  // Configuring the eta dependent Isolation corrections
  //
  fCorMap.clear();
  fCorMap.insert(CorMap::value_type(std::string("cone30Comb"), ConfigEtCone30Comb()));
  fCorMap.insert(CorMap::value_type(std::string("cone20Comb"), ConfigEtCone20Comb()));
  fCorMap.insert(CorMap::value_type(std::string("cone40Comb"), ConfigEtCone40Comb()));
  fCorMap.insert(CorMap::value_type(std::string("cone20MUID"), ConfigEtCone20MUID()));
  fCorMap.insert(CorMap::value_type(std::string("cone30MUID"), ConfigEtCone30MUID()));
  fCorMap.insert(CorMap::value_type(std::string("cone40MUID"), ConfigEtCone40MUID()));

  // Quadratic corrections
  fCorMap.insert(CorMap::value_type(std::string("cone20Comb2011"), ConfigEtCone20Comb2011()));
  fCorMap.insert(CorMap::value_type(std::string("cone30Comb2011"), ConfigEtCone30Comb2011()));
  fCorMap.insert(CorMap::value_type(std::string("cone40Comb2011"), ConfigEtCone40Comb2011()));
  fCorMap.insert(CorMap::value_type(std::string("cone20Comb2012"), ConfigEtCone20Comb2012()));
  fCorMap.insert(CorMap::value_type(std::string("cone30Comb2012"), ConfigEtCone30Comb2012()));
  fCorMap.insert(CorMap::value_type(std::string("cone40Comb2012"), ConfigEtCone40Comb2012()));

  //
  // Settingt the default corrections to etcone30
  //
  SetCorrections(fOpts);

  if(!fCor) cout << "IsoConfig::IsoConfig - ERROR could not configure default correction!" << endl;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone20Comb2011()
{ //                   
  // Eta dependent isolation from Doug's Ztag and probe  
  //   -- 03/05/2012                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.constant.push_back(21.667);  corr.slope.push_back(23.471); corr.quad.push_back(0.905); //fEta[0] =0.10; 
  corr.constant.push_back(21.667);  corr.slope.push_back(23.471); corr.quad.push_back(0.905); //fEta[1] =0.20; 
  corr.constant.push_back(21.667);  corr.slope.push_back(23.471); corr.quad.push_back(0.905); //fEta[2] =0.30; 
  corr.constant.push_back(21.667);  corr.slope.push_back(23.471); corr.quad.push_back(0.905); //fEta[3] =0.40; 
  corr.constant.push_back(21.667);  corr.slope.push_back(23.471); corr.quad.push_back(0.905); //fEta[4] =0.50; 
  corr.constant.push_back(38.196);  corr.slope.push_back(23.775); corr.quad.push_back(0.791); //fEta[5] =0.60; 
  corr.constant.push_back(38.196);  corr.slope.push_back(23.775); corr.quad.push_back(0.791); //fEta[6] =0.70; 
  corr.constant.push_back(38.196);  corr.slope.push_back(23.775); corr.quad.push_back(0.791); //fEta[7] =0.80; 
  corr.constant.push_back(38.196);  corr.slope.push_back(23.775); corr.quad.push_back(0.791); //fEta[8] =0.90; 
  corr.constant.push_back(38.196);  corr.slope.push_back(23.775); corr.quad.push_back(0.791); //fEta[9] =1.00; 
  corr.constant.push_back(38.196);  corr.slope.push_back(23.775); corr.quad.push_back(0.791); //fEta[10]=1.10; 
  corr.constant.push_back(34.393);  corr.slope.push_back(31.116); corr.quad.push_back(0.3148); //fEta[11]=1.20; 
  corr.constant.push_back(34.393);  corr.slope.push_back(31.116); corr.quad.push_back(0.3148); //fEta[12]=1.30; 
  corr.constant.push_back(34.393);  corr.slope.push_back(31.116); corr.quad.push_back(0.3148); //fEta[13]=1.40; 
  corr.constant.push_back(34.393);  corr.slope.push_back(31.116); corr.quad.push_back(0.3148); //fEta[14]=1.50; 
  corr.constant.push_back(34.393);  corr.slope.push_back(31.116); corr.quad.push_back(0.3148); //fEta[15]=1.60; 
  corr.constant.push_back(34.393);  corr.slope.push_back(31.116); corr.quad.push_back(0.3148); //fEta[16]=1.70; 
  corr.constant.push_back(-22.39);  corr.slope.push_back(29.084); corr.quad.push_back(0.2927); //fEta[17]=1.80; 
  corr.constant.push_back(-22.39);  corr.slope.push_back(29.084); corr.quad.push_back(0.2927); //fEta[18]=1.90; 
  corr.constant.push_back(-22.39);  corr.slope.push_back(29.084); corr.quad.push_back(0.2927); //fEta[19]=2.00; 
  corr.constant.push_back(-22.39);  corr.slope.push_back(29.084); corr.quad.push_back(0.2927); //fEta[20]=2.10; 
  corr.constant.push_back(-22.39);  corr.slope.push_back(29.084); corr.quad.push_back(0.2927); //fEta[21]=2.20; 
  corr.constant.push_back(-22.39);  corr.slope.push_back(29.084); corr.quad.push_back(0.2927); //fEta[22]=2.30; 
  corr.constant.push_back(28.440);  corr.slope.push_back(9.5149); corr.quad.push_back(1.0129); //fEta[23]=2.40; 
  corr.constant.push_back(28.440);  corr.slope.push_back(9.5149); corr.quad.push_back(1.0129); //fEta[24]=2.50; 
  corr.constant.push_back(28.440);  corr.slope.push_back(9.5149); corr.quad.push_back(1.0129); //fEta[25]=2.60; 
  corr.constant.push_back(28.440);  corr.slope.push_back(9.5149); corr.quad.push_back(1.0129); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone30Comb2011()
{ //                   
  // Eta dependent isolation from Doug's Ztag and probe  
  //   -- 03/05/2012                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(48.714);  corr.constant.push_back(10.646);  corr.quad.push_back(1.012); //fEta[0] =0.10; 
  corr.slope.push_back(48.714);  corr.constant.push_back(10.646);  corr.quad.push_back(1.012); //fEta[1] =0.20; 
  corr.slope.push_back(48.714);  corr.constant.push_back(10.646);  corr.quad.push_back(1.012); //fEta[2] =0.30; 
  corr.slope.push_back(48.714);  corr.constant.push_back(10.646);  corr.quad.push_back(1.012); //fEta[3] =0.40; 
  corr.slope.push_back(48.714);  corr.constant.push_back(10.646);  corr.quad.push_back(1.012); //fEta[4] =0.50; 
  corr.slope.push_back(44.330);  corr.constant.push_back(67.8310); corr.quad.push_back(1.708); //fEta[5] =0.60; 
  corr.slope.push_back(44.330);  corr.constant.push_back(67.8310); corr.quad.push_back(1.708); //fEta[6] =0.70; 
  corr.slope.push_back(44.330);  corr.constant.push_back(67.8310); corr.quad.push_back(1.708); //fEta[7] =0.80; 
  corr.slope.push_back(44.330);  corr.constant.push_back(67.8310); corr.quad.push_back(1.708); //fEta[8] =0.90; 
  corr.slope.push_back(44.330);  corr.constant.push_back(67.8310); corr.quad.push_back(1.708); //fEta[9] =1.00; 
  corr.slope.push_back(44.330);  corr.constant.push_back(67.8310); corr.quad.push_back(1.708); //fEta[10]=1.10; 
  corr.slope.push_back(63.287);  corr.constant.push_back(73.3910); corr.quad.push_back(1.291); //fEta[11]=1.20; 
  corr.slope.push_back(63.287);  corr.constant.push_back(73.3910); corr.quad.push_back(1.291); //fEta[12]=1.30; 
  corr.slope.push_back(63.287);  corr.constant.push_back(73.3910); corr.quad.push_back(1.291); //fEta[13]=1.40; 
  corr.slope.push_back(63.287);  corr.constant.push_back(73.3910); corr.quad.push_back(1.291); //fEta[14]=1.50; 
  corr.slope.push_back(63.287);  corr.constant.push_back(73.3910); corr.quad.push_back(1.291); //fEta[15]=1.60; 
  corr.slope.push_back(63.287);  corr.constant.push_back(73.3910); corr.quad.push_back(1.291); //fEta[16]=1.70; 
  corr.slope.push_back(65.637);  corr.constant.push_back(129.261); corr.quad.push_back(1.617); //fEta[17]=1.80; 
  corr.slope.push_back(65.637);  corr.constant.push_back(129.261); corr.quad.push_back(1.617); //fEta[18]=1.90; 
  corr.slope.push_back(65.637);  corr.constant.push_back(129.261); corr.quad.push_back(1.617); //fEta[19]=2.00; 
  corr.slope.push_back(65.637);  corr.constant.push_back(129.261); corr.quad.push_back(1.617); //fEta[20]=2.10; 
  corr.slope.push_back(65.637);  corr.constant.push_back(129.261); corr.quad.push_back(1.617); //fEta[21]=2.20; 
  corr.slope.push_back(65.637);  corr.constant.push_back(129.261); corr.quad.push_back(1.617); //fEta[22]=2.30; 
  corr.slope.push_back(74.510);  corr.constant.push_back(-14.271); corr.quad.push_back(1.012); //fEta[23]=2.40; 
  corr.slope.push_back(74.510);  corr.constant.push_back(-14.271); corr.quad.push_back(1.012); //fEta[24]=2.50; 
  corr.slope.push_back(74.510);  corr.constant.push_back(-14.271); corr.quad.push_back(1.012); //fEta[25]=2.60; 
  corr.slope.push_back(74.510);  corr.constant.push_back(-14.271); corr.quad.push_back(1.012); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone40Comb2011()
{ //                   
  // Eta dependent isolation from Doug's Ztag and probe  
  //   -- 03/05/2012                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(141.47);  corr.constant.push_back(37.790); corr.quad.push_back(2.2497); //fEta[0] =0.10; 
  corr.slope.push_back(141.47);  corr.constant.push_back(37.790); corr.quad.push_back(2.2497); //fEta[1] =0.20; 
  corr.slope.push_back(141.47);  corr.constant.push_back(37.790); corr.quad.push_back(2.2497); //fEta[2] =0.30; 
  corr.slope.push_back(141.47);  corr.constant.push_back(37.790); corr.quad.push_back(2.2497); //fEta[3] =0.40; 
  corr.slope.push_back(141.47);  corr.constant.push_back(37.790); corr.quad.push_back(2.2497); //fEta[4] =0.50; 
  corr.slope.push_back(135.31);  corr.constant.push_back(66.144); corr.quad.push_back(2.2998); //fEta[5] =0.60; 
  corr.slope.push_back(135.31);  corr.constant.push_back(66.144); corr.quad.push_back(2.2998); //fEta[6] =0.70; 
  corr.slope.push_back(135.31);  corr.constant.push_back(66.144); corr.quad.push_back(2.2998); //fEta[7] =0.80; 
  corr.slope.push_back(135.31);  corr.constant.push_back(66.144); corr.quad.push_back(2.2998); //fEta[8] =0.90; 
  corr.slope.push_back(135.31);  corr.constant.push_back(66.144); corr.quad.push_back(2.2998); //fEta[9] =1.00; 
  corr.slope.push_back(135.31);  corr.constant.push_back(66.144); corr.quad.push_back(2.2998); //fEta[10]=1.10; 
  corr.slope.push_back(124.27);  corr.constant.push_back(90.713); corr.quad.push_back(1.9617); //fEta[11]=1.20; 
  corr.slope.push_back(124.27);  corr.constant.push_back(90.713); corr.quad.push_back(1.9617); //fEta[12]=1.30; 
  corr.slope.push_back(124.27);  corr.constant.push_back(90.713); corr.quad.push_back(1.9617); //fEta[13]=1.40; 
  corr.slope.push_back(124.27);  corr.constant.push_back(90.713); corr.quad.push_back(1.9617); //fEta[14]=1.50; 
  corr.slope.push_back(124.27);  corr.constant.push_back(90.713); corr.quad.push_back(1.9617); //fEta[15]=1.60; 
  corr.slope.push_back(124.27);  corr.constant.push_back(90.713); corr.quad.push_back(1.9617); //fEta[16]=1.70; 
  corr.slope.push_back(108.50);  corr.constant.push_back(17.497); corr.quad.push_back(2.8940); //fEta[17]=1.80; 
  corr.slope.push_back(108.50);  corr.constant.push_back(17.497); corr.quad.push_back(2.8940); //fEta[18]=1.90; 
  corr.slope.push_back(108.50);  corr.constant.push_back(17.497); corr.quad.push_back(2.8940); //fEta[19]=2.00; 
  corr.slope.push_back(108.50);  corr.constant.push_back(17.497); corr.quad.push_back(2.8940); //fEta[20]=2.10; 
  corr.slope.push_back(108.50);  corr.constant.push_back(17.497); corr.quad.push_back(2.8940); //fEta[21]=2.20; 
  corr.slope.push_back(108.50);  corr.constant.push_back(17.497); corr.quad.push_back(2.8940); //fEta[22]=2.30; 
  corr.slope.push_back(82.927);  corr.constant.push_back(-22.463); corr.quad.push_back(2.7417); //fEta[23]=2.40; 
  corr.slope.push_back(82.927);  corr.constant.push_back(-22.463); corr.quad.push_back(2.7417); //fEta[24]=2.50; 
  corr.slope.push_back(82.927);  corr.constant.push_back(-22.463); corr.quad.push_back(2.7417); //fEta[25]=2.60; 
  corr.slope.push_back(82.927);  corr.constant.push_back(-22.463); corr.quad.push_back(2.7417); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone20Comb2012()
{ //                   
  // Eta dependent isolation from Doug's Ztag and probe  
  //   -- 05/30/2012                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[0] =0.10; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[1] =0.20; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[2] =0.30; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[3] =0.40; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[4] =0.50; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[5] =0.60; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[6] =0.70; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[7] =0.80; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[8] =0.90; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[9] =1.00; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[10]=1.10; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[11]=1.20; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[12]=1.30; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[13]=1.40; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[14]=1.50; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[15]=1.60; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[16]=1.70; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[17]=1.80; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[18]=1.90; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[19]=2.00; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[20]=2.10; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[21]=2.20; 
  corr.slope.push_back(22.624);  corr.constant.push_back(-57.258); corr.quad.push_back(0.189); //fEta[22]=2.30; 
  corr.slope.push_back(22.624);  corr.constant.push_back(0.0000);  corr.quad.push_back(0.189); //fEta[23]=2.40; 
  corr.slope.push_back(22.624);  corr.constant.push_back(0.0000);  corr.quad.push_back(0.189); //fEta[24]=2.50; 
  corr.slope.push_back(22.624);  corr.constant.push_back(0.0000);  corr.quad.push_back(0.189); //fEta[25]=2.60; 
  corr.slope.push_back(22.624);  corr.constant.push_back(0.0000);  corr.quad.push_back(0.189); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone30Comb2012()
{ //                   
  // Eta dependent isolation from Doug's Ztag and probe  
  //   -- 03/05/2012                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[0] =0.10; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[1] =0.20; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[2] =0.30; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[3] =0.40; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[4] =0.50; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[5] =0.60; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[6] =0.70; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[7] =0.80; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[8] =0.90; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[9] =1.00; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[10]=1.10; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[11]=1.20; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[12]=1.30; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[13]=1.40; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[14]=1.50; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[15]=1.60; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[16]=1.70; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[17]=1.80; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[18]=1.90; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[19]=2.00; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[20]=2.10; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[21]=2.20; 
  corr.slope.push_back(42.455);  corr.constant.push_back(53.039); corr.quad.push_back(1.1265); //fEta[22]=2.30; 
  corr.slope.push_back(42.455);  corr.constant.push_back(0.0000); corr.quad.push_back(1.1265); //fEta[23]=2.40; 
  corr.slope.push_back(42.455);  corr.constant.push_back(0.0000); corr.quad.push_back(1.1265); //fEta[24]=2.50; 
  corr.slope.push_back(42.455);  corr.constant.push_back(0.0000); corr.quad.push_back(1.1265); //fEta[25]=2.60; 
  corr.slope.push_back(42.455);  corr.constant.push_back(0.0000); corr.quad.push_back(1.1265); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone40Comb2012()
{ //                   
  // Eta dependent isolation from Doug's Ztag and probe  
  //   -- 05/30/2012                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[0] =0.10; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[1] =0.20; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[2] =0.30; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[3] =0.40; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[4] =0.50; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[5] =0.60; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[6] =0.70; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[7] =0.80; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[8] =0.90; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[9] =1.00; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[10]=1.10; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[11]=1.20; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[12]=1.30; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[13]=1.40; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[14]=1.50; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[15]=1.60; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[16]=1.70; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[17]=1.80; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[18]=1.90; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[19]=2.00; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[20]=2.10; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[21]=2.20; 
  corr.slope.push_back(88.146);  corr.constant.push_back(65.831); corr.quad.push_back(1.466); //fEta[22]=2.30; 
  corr.slope.push_back(88.146);  corr.constant.push_back(0.0000); corr.quad.push_back(1.466); //fEta[23]=2.40; 
  corr.slope.push_back(88.146);  corr.constant.push_back(0.0000); corr.quad.push_back(1.466); //fEta[24]=2.50; 
  corr.slope.push_back(88.146);  corr.constant.push_back(0.0000); corr.quad.push_back(1.466); //fEta[25]=2.60; 
  corr.slope.push_back(88.146);  corr.constant.push_back(0.0000); corr.quad.push_back(1.466); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone20Comb()
{ //                   
  // Eta dependent isolation from Maria's Ztag and probe  
  //   -- 02/16/2011                                      
  //   -- units are MeV
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[0] =0.10; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[1] =0.20; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[2] =0.30; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[3] =0.40; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[4] =0.50; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[5] =0.60; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[6] =0.70; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[7] =0.80; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[8] =0.90; 
  corr.slope.push_back(44.6);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[9] =1.00; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[10]=1.10; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[11]=1.20; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[12]=1.30; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[13]=1.40; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[14]=1.50; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[15]=1.60; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[16]=1.70; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[17]=1.80; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[18]=1.90; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[19]=2.00; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[20]=2.10; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[21]=2.20; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[22]=2.30; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[23]=2.40; 
  corr.slope.push_back(33.0);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[24]=2.50; 
  corr.slope.push_back(22.8);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[25]=2.60; 
  corr.slope.push_back(22.8);  corr.constant.push_back(207.4); corr.quad.push_back(0.0); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone30Comb()
{ //                   
  // Eta dependent isolation from Peter's Ztag and probe  
  //   -- 11/7/2011                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(81.4380);  corr.constant.push_back(256.798); corr.quad.push_back(0.0); //fEta[0] =0.10; 
  corr.slope.push_back(88.7883);  corr.constant.push_back(222.295); corr.quad.push_back(0.0); //fEta[1] =0.20; 
  corr.slope.push_back(90.4055);  corr.constant.push_back(246.810); corr.quad.push_back(0.0); //fEta[2] =0.30; 
  corr.slope.push_back(88.5415);  corr.constant.push_back(256.253); corr.quad.push_back(0.0); //fEta[3] =0.40; 
  corr.slope.push_back(89.8105);  corr.constant.push_back(253.627); corr.quad.push_back(0.0); //fEta[4] =0.50; 
  corr.slope.push_back(87.2666);  corr.constant.push_back(269.239); corr.quad.push_back(0.0); //fEta[5] =0.60; 
  corr.slope.push_back(89.8660);  corr.constant.push_back(231.288); corr.quad.push_back(0.0); //fEta[6] =0.70; 
  corr.slope.push_back(88.4182);  corr.constant.push_back(219.546); corr.quad.push_back(0.0); //fEta[7] =0.80; 
  corr.slope.push_back(83.7103);  corr.constant.push_back(238.926); corr.quad.push_back(0.0); //fEta[8] =0.90; 
  corr.slope.push_back(88.3318);  corr.constant.push_back(193.716); corr.quad.push_back(0.0); //fEta[9] =1.00; 
  corr.slope.push_back(85.1312);  corr.constant.push_back(207.372); corr.quad.push_back(0.0); //fEta[10]=1.10; 
  corr.slope.push_back(86.0413);  corr.constant.push_back(306.257); corr.quad.push_back(0.0); //fEta[11]=1.20; 
  corr.slope.push_back(81.6753);  corr.constant.push_back(250.002); corr.quad.push_back(0.0); //fEta[12]=1.30; 
  corr.slope.push_back(78.6961);  corr.constant.push_back(377.195); corr.quad.push_back(0.0); //fEta[13]=1.40; 
  corr.slope.push_back(66.8804);  corr.constant.push_back(259.919); corr.quad.push_back(0.0); //fEta[14]=1.50; 
  corr.slope.push_back(71.9913);  corr.constant.push_back(187.506); corr.quad.push_back(0.0); //fEta[15]=1.60; 
  corr.slope.push_back(73.2569);  corr.constant.push_back(275.101); corr.quad.push_back(0.0); //fEta[16]=1.70; 
  corr.slope.push_back(78.2649);  corr.constant.push_back(165.755); corr.quad.push_back(0.0); //fEta[17]=1.80; 
  corr.slope.push_back(77.9631);  corr.constant.push_back(199.264); corr.quad.push_back(0.0); //fEta[18]=1.90; 
  corr.slope.push_back(79.0986);  corr.constant.push_back(198.490); corr.quad.push_back(0.0); //fEta[19]=2.00; 
  corr.slope.push_back(85.7765);  corr.constant.push_back(153.297); corr.quad.push_back(0.0); //fEta[20]=2.10; 
  corr.slope.push_back(77.4238);  corr.constant.push_back(177.392); corr.quad.push_back(0.0); //fEta[21]=2.20; 
  corr.slope.push_back(75.9160);  corr.constant.push_back(120.396); corr.quad.push_back(0.0); //fEta[22]=2.30; 
  corr.slope.push_back(70.5277);  corr.constant.push_back(123.942); corr.quad.push_back(0.0); //fEta[23]=2.40; 
  corr.slope.push_back(70.5277);  corr.constant.push_back(123.942); corr.quad.push_back(0.0); //fEta[24]=2.50; 
  corr.slope.push_back(56.0000);  corr.constant.push_back(123.942); corr.quad.push_back(0.0); //fEta[25]=2.60; 
  corr.slope.push_back(56.0000);  corr.constant.push_back(123.942); corr.quad.push_back(0.0); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone40Comb()
{ //                   
  // Eta dependent isolation from Maria's Ztag and probe  
  //   -- 02/16/2011                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[0] =0.10; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[1] =0.20; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[2] =0.30; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[3] =0.40; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[4] =0.50; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[5] =0.60; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[6] =0.70; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[7] =0.80; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[8] =0.90; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[9] =1.00; 
  corr.slope.push_back(183.2);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[10]=1.10; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[11]=1.20; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[12]=1.30; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[13]=1.40; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[14]=1.50; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[15]=1.60; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[16]=1.70; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[17]=1.80; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[18]=1.90; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[19]=2.00; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[20]=2.10; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[21]=2.20; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[22]=2.30; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[23]=2.40; 
  corr.slope.push_back(154.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[24]=2.50; 
  corr.slope.push_back(114.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[25]=2.60; 
  corr.slope.push_back(114.3);  corr.constant.push_back(365.6);   corr.quad.push_back(0.0); //fEta[26]=2.70; 

  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone20MUID()
{ //                   
  // Eta dependent isolation from Maria's Ztag and probe  
  //   -- 02/16/2011                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[0] =0.10; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[1] =0.20; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[2] =0.30; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[3] =0.40; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[4] =0.50; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[5] =0.60; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[6] =0.70; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[7] =0.80; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[8] =0.90; 
  corr.slope.push_back(43.8);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[9] =1.00; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[10]=1.10; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[11]=1.20; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[12]=1.30; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[13]=1.40; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[14]=1.50; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[15]=1.60; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[16]=1.70; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[17]=1.80; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[18]=1.90; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[19]=2.00; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[20]=2.10; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[21]=2.20; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[22]=2.30; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[23]=2.40; 
  corr.slope.push_back(33.6);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[24]=2.50; 
  corr.slope.push_back(31.5);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[25]=2.60; 
  corr.slope.push_back(31.5);  corr.constant.push_back(200.0); corr.quad.push_back(0.0); //fEta[26]=2.70; 

  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone30MUID()
{ //                   
  // Eta dependent isolation from Maria's Ztag and probe  
  //   -- 02/16/2011                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[0] =0.10; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[1] =0.20; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[2] =0.30; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[3] =0.40; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[4] =0.50; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[5] =0.60; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[6] =0.70; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[7] =0.80; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[8] =0.90; 
  corr.slope.push_back(101.1);  corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[9] =1.00; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[10]=1.10; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[11]=1.20; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[12]=1.30; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[13]=1.40; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[14]=1.50; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[15]=1.60; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[16]=1.70; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[17]=1.80; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[18]=1.90; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[19]=2.00; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[20]=2.10; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[21]=2.20; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[22]=2.30; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[23]=2.40; 
  corr.slope.push_back(85.4);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[24]=2.50; 
  corr.slope.push_back(67.5);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[25]=2.60; 
  corr.slope.push_back(67.5);   corr.constant.push_back(250.);   corr.quad.push_back(0.0); //fEta[26]=2.70; 
  return corr;
}

//-----------------------------------------------------------------------------
IsoConfig::ConeCorr IsoConfig::ConfigEtCone40MUID()
{ //                   
  // Eta dependent isolation from Maria's Ztag and probe  
  //   -- 02/16/2011                                      
  //   -- units are MeV  
  //
  ConeCorr corr;
  corr.slope   .clear();
  corr.constant.clear();
  corr.quad    .clear();
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[0] =0.10; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[1] =0.20; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[2] =0.30; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[3] =0.40; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[4] =0.50; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[5] =0.60; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[6] =0.70; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[7] =0.80; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[8] =0.90; 
  corr.slope.push_back(181.8);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[9] =1.00; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[10]=1.10; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[11]=1.20; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[12]=1.30; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[13]=1.40; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[14]=1.50; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[15]=1.60; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[16]=1.70; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[17]=1.80; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[18]=1.90; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[19]=2.00; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[20]=2.10; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[21]=2.20; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[22]=2.30; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[23]=2.40; 
  corr.slope.push_back(155.4);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[24]=2.50; 
  corr.slope.push_back(149.5);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[25]=2.60; 
  corr.slope.push_back(149.5);  corr.constant.push_back(365.6); corr.quad.push_back(0.0);  //fEta[26]=2.70; 

  return corr;
}

//-----------------------------------------------------------------------------    
void IsoConfig::SetCorrections(const std::string opts) 
{
  //
  // Setting the Isolation Correction values
  //
  if(fOpts==opts && fCor) return;

  CorMap::iterator iter = fCorMap.find(opts);

  // If you find the corrections, then reset the pointer to the requested corrections
  if(iter!=fCorMap.end()){
    fCor = &(iter->second);
    fOpts=iter->first;
  }
  else{
    cout << "IsoConfig::GetCorrections - ERROR could not find corrections" << endl;
  }

  if(!fCor) cout << "IsoConfig::GetCorrections - ERROR invalid corrections for: " << opts << endl;

}

//-----------------------------------------------------------------------------
float IsoConfig::GetCaloCor(const std::string opts, const unsigned &nvtx, const float &eta)
{
  //
  // Returns the correction to the calorimeter isolation
  //   -- this is the amount to subtract from the corrected cone
  //
  SetCorrections(opts);
  if(!fCor) cout << "IsoConfig::GetCaloCor - ERROR invalid corrections for: " << opts << endl;
  const unsigned bin      = fCor->GetEtaBin(eta);

  if(fDebug){
    if(bin>=fCor->quad.size() || bin>=fCor->slope.size() || bin>=fCor->constant.size()){
      cout << "ERROR - eta too big with bin: " << bin << " quad: " << fCor->quad.size() << " slope: " << fCor->slope.size() 
	   << " constant: " << fCor->constant.size()  << " fOpts: " << fOpts << endl;
      return 0.0;
    }
  }
  const float &quad       = fCor->quad    .at(bin);
  const float &slope      = fCor->slope   .at(bin);
  const float &constant   = fCor->constant.at(bin);

  return CaloCor(nvtx,slope,constant,quad);
}

//-----------------------------------------------------------------------------
float IsoConfig::CaloCor(const unsigned &nvtx, const float &slope, const float &constant, const float &quad)
{
  //
  // Returns the correction to the calorimeter isolation
  //   -- this is the amount to subtract from the corrected cone
  //
  return float(nvtx)*slope + constant +  pow(float(nvtx),2.0)*quad;
}
