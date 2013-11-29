#include "MissingETUtility/TESUncertaintyProvider.h"

TESUncertaintyProvider::~TESUncertaintyProvider(){}
TESUncertaintyProvider::TESUncertaintyProvider(){
  m_nEta=5;
  m_eta_bins.push_back(0.);
  m_eta_bins.push_back(0.3);
  m_eta_bins.push_back(0.8);
  m_eta_bins.push_back(1.3);
  m_eta_bins.push_back(1.6);
  m_eta_bins.push_back(2.5);

  m_nPt=5;
  m_pt_bins.push_back(15.);
  m_pt_bins.push_back(20.);
  m_pt_bins.push_back(30.);
  m_pt_bins.push_back(50.);
  m_pt_bins.push_back(150.);
  m_pt_bins.push_back(500.);

  // 1prong
  std::vector<double> tmp;
  //0-0.3
  tmp.push_back(0.03);
  tmp.push_back(0.025);
  tmp.push_back(0.025);
  tmp.push_back(0.03);
  tmp.push_back(0.03);

  m_1p_unc.push_back(tmp);
  
  //0.3-0.8
  tmp.clear();
  tmp.push_back(0.03);
  tmp.push_back(0.025);
  tmp.push_back(0.025);
  tmp.push_back(0.03);
  tmp.push_back(0.03);

  m_1p_unc.push_back(tmp);

  //0.8-1.3
  tmp.clear();
  tmp.push_back(0.03);
  tmp.push_back(0.025);
  tmp.push_back(0.025);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  m_1p_unc.push_back(tmp);

  //1.3-1.6
  tmp.clear();
  tmp.push_back(0.035);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.035);
  tmp.push_back(0.035);
  m_1p_unc.push_back(tmp);

  //>1.6
  tmp.clear();
  tmp.push_back(0.03);
  tmp.push_back(0.025);
  tmp.push_back(0.025);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  m_1p_unc.push_back(tmp);

 // Multi prong
  tmp.clear();
  //0-0.3
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.035);

  m_mp_unc.push_back(tmp);
  
  //0.3-0.8
  tmp.clear();
  tmp.push_back(0.035);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.03);

  m_mp_unc.push_back(tmp);

  //0.8-1.3
  tmp.clear();
  tmp.push_back(0.05);
  tmp.push_back(0.035);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  m_mp_unc.push_back(tmp);
  
  //1.3-1.6
  tmp.clear();
  tmp.push_back(0.045);
  tmp.push_back(0.035);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  tmp.push_back(0.03);
  m_mp_unc.push_back(tmp);
  
  //>1.6
  tmp.clear();
  tmp.push_back(0.03);
  tmp.push_back(0.025);
  tmp.push_back(0.025);
  tmp.push_back(0.025);
  tmp.push_back(0.03);
  m_mp_unc.push_back(tmp);
}

double TESUncertaintyProvider::GetTESUncertainty(double pt, double eta, int nProng){
  if(pt<15. ||fabs(eta)>2.5 || pt>500.){
    return -1.;
  }
  
  int iE=-1;
  int iP=-1;
  
  for(int i=0;i<m_nEta;i++){
    if(m_eta_bins[i]<=fabs(eta) && fabs(eta)<m_eta_bins[i+1]){
      iE=i; break;
    }
  }
  //std::cout<<"EtaBin: "<<iE<<std::endl;

  for(int i=0;i<m_nPt;i++){
    if(m_pt_bins[i]<=pt && pt<m_pt_bins[i+1]){
      iP=i; break;
    }
  }
  //std::cout<<"PtBin: "<<iP<<std::endl;
  
  if(nProng==1)
    return  m_1p_unc[iE][iP];
  else
    return  m_mp_unc[iE][iP];

}
