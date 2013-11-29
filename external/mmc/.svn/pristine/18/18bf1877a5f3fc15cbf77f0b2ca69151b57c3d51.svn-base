#include <iostream>
#include <string>

#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <time.h>
#include <iostream>

#define nDRWEIGHT 16
#define nDTHETAWEIGHT 18
#define M_tau 1.777
#define MaxCorrMET 40	// Limiting the Correction of MET in GeV (Have to be Int_teger)

using namespace std;
using namespace TMath;

class mmc {

public:
  mmc(Int_t wchoice=3);
  void Clear(){m_mz_peak=-1; m_mz_mean=-1; m_mz_maxprob=-1; m_mz_coll=-1; }
  Double_t Scan6d(const TLorentzVector _lep1, const TLorentzVector _lep2, const Double_t _metx, const Double_t _mety, const Double_t _sumet);
  Double_t Scan6dTheta(const TLorentzVector _lep1, const TLorentzVector _lep2, const Double_t _metx, const Double_t _mety, const Double_t _sumet);
  Double_t Scan4dTheta(const TLorentzVector _lep1, const TLorentzVector _lep2, const Double_t _metx, const Double_t _mety);
  Double_t Scan6dAnal(const TLorentzVector _lep1, const TLorentzVector _lep2, const Double_t _metx, const Double_t _mety, const Double_t sigma, const int mode);
  Double_t Scan4dAnal(const TLorentzVector _lep1, const TLorentzVector _lep2, const Double_t _metx, const Double_t _mety, const int mode);
  Double_t GetCollMass(const TLorentzVector _lep1, const TLorentzVector _lep2, const Double_t _metx, const Double_t _mety);

  void SetAngleWeight(Int_t choice=0);
  void SetMassWeight(Int_t choice=1){ m_mm_weight=choice; };
  Double_t GetMZ_peak(){ return m_mz_peak; };
  Double_t GetMZ_mean(){ return m_mz_mean; };
  Double_t GetMZ_maxprob(){ return m_mz_maxprob; };
  Double_t GetMZ_coll(){ return m_mz_coll; };
  Double_t GetExecTime(){ return m_exectime; };
  Double_t GetCorrMEX(){ return Corr_met_ex; };
  Double_t GetCorrMEY(){ return Corr_met_ey; };
  TH1F * GetMassDist(){ return h_mz_dist; };

  void SetPhi_N(Double_t arg1){ m_phi_n=arg1; }
  void SetTheta_N(Double_t arg1){ m_theta_n=arg1; }
  void SetMet_N(Double_t arg1){ m_met_n=arg1; }
  void SetMMa_N(Double_t arg1){ m_mma_n=arg1; }

  void SetPhi_Lim(Double_t arg1){ m_phi_lim=arg1; }
  void SetTheta_Lim(Double_t arg1){ m_theta_lim=arg1; }
  void SetMet_Lim(Double_t arg1){ m_met_lim=arg1; }
  void SetMMa_LimLow(Double_t arg1){ m_mma_limlow=arg1; }
  void SetMMa_LimHigh(Double_t arg1){ m_mma_limhigh=arg1; }

  void SetMet_S(Double_t arg1){ m_met_s=arg1; }
  void SetMet_LimSigma(Double_t arg1){ m_met_limsigma=arg1; }
  void SetMet_ScanMode(Int_t arg1){ m_met_scanmode=arg1; }

  void Set_Alpha(Double_t arg1){ m_alpha=arg1; }
  void Set_METNoise(Double_t arg1){ m_noise=arg1; }

private:
  
  TF1*  ReadFunc( TFile* _file, TString histName, TString taumu);
  
  void  GetDeltaRDistribution();
  void  GetDeltaThetaDistribution();
  void  GetMDistribution();
  Double_t GetMETValue( Double_t met, Double_t sigma);
  Double_t GetDeltaRValue(Double_t dR1, Double_t pt1);
  Double_t GetDeltaThetaValue(Double_t dTheta1, Double_t pt1);
  Double_t GetMValue(Double_t m1);
  Double_t GetAnaWeight_TauMiss(const TLorentzVector tau, const TLorentzVector miss, Double_t &ntheta3d);
  Double_t GetJacobian(const TLorentzVector vis_pos, const TLorentzVector vis_neg, const TLorentzVector miss_pos, const TLorentzVector miss_neg, const Double_t _Ex_miss, const Double_t _Ey_miss);

  void SetMZ_peak(Double_t mz){ m_mz_peak=mz; };
  void SetMZ_mean(Double_t mz){ m_mz_mean=mz; };
  void SetMZ_maxprob(Double_t mz){ m_mz_maxprob=mz; };
  void SetMZ_coll(Double_t mz){ m_mz_coll=mz; };

  //Double_t GetZmass_Max(vector<Double_t> *v_wkeit, vector<Double_t> *v_zm); 
  //Double_t GetZmass_Weighted(vector<Double_t> *v_wkeit, vector<Double_t> *v_zm);

  Int_t _debug;
  Double_t m_mz_peak;
  Double_t m_mz_mean;
  Double_t m_mz_maxprob;
  Double_t m_mz_coll;
  Double_t Corr_met_ex;
  Double_t Corr_met_ey;

  Int_t m_weight;
  Int_t m_mm_weight;
  TF1* f_dR_weight[nDRWEIGHT];
  TF1* f_dTheta_weight[nDTHETAWEIGHT];
  TF1* f_m_weight;
  TH1F * h_mz_dist;

  Double_t m_phi_n;
  Double_t m_theta_n;
  Double_t m_met_n;
  Double_t m_mma_n;
  Double_t m_costheta;
  Double_t theta_star_max;
  Double_t f_cut;		// define f_cut = lep_ptCUT / tauPT
  Double_t m_mma_n_4D;
  Double_t m_mma_n_6D;
  Double_t m_costheta_4D;
  Double_t m_costheta_6D;
  
  Double_t m_phi_lim;
  Double_t m_theta_lim;
  Double_t m_met_lim;
  Double_t m_mma_limlow;
  Double_t m_mma_limhigh;

  Double_t m_met_scanmode;
  Double_t m_met_limsigma;
  Double_t m_met_s;
  Double_t MET_scan_poInt_ts_Phi;
  Double_t MET_scan_poInt_ts_Z;
  
  Double_t m_exectime;
  Double_t m_alpha;
  Double_t m_noise;

};
