#include "HSG4LepLepTriggerSF/HSG4LepLepTriggerSF.h"
#include "TFile.h"
#include "TROOT.h"  //--AM--
#include <iostream>

using namespace std ;

HSG4LepLepTriggerSF::HSG4LepLepTriggerSF( const std::string& path , const Bool_t debug )
    :
  // scale factor histograms

  e12Tvhl1_AB_Hist_SF_eta(NULL),
  e12Tvhl1_AB_Hist_SF_pT(NULL),  
  e12Tvhl1_AB_Hist_SF_eta_errUp(NULL),
  e12Tvhl1_AB_Hist_SF_pT_errUp(NULL),  
  e12Tvhl1_AB_Hist_SF_eta_errDown(NULL),
  e12Tvhl1_AB_Hist_SF_pT_errDown(NULL),
  
  e12Tvhl1_C_Hist_SF_eta(NULL),
  e12Tvhl1_C_Hist_SF_pT(NULL),  
  e12Tvhl1_C_Hist_SF_eta_errUp(NULL),
  e12Tvhl1_C_Hist_SF_pT_errUp(NULL),  
  e12Tvhl1_C_Hist_SF_eta_errDown(NULL),
  e12Tvhl1_C_Hist_SF_pT_errDown(NULL),
  
  e12Tvhl1_D_Hist_SF_eta(NULL),
  e12Tvhl1_D_Hist_SF_pT(NULL),  
  e12Tvhl1_D_Hist_SF_eta_errUp(NULL),
  e12Tvhl1_D_Hist_SF_pT_errUp(NULL),  
  e12Tvhl1_D_Hist_SF_eta_errDown(NULL),
  e12Tvhl1_D_Hist_SF_pT_errDown(NULL),
  
  e12Tvhl1_E_Hist_SF_eta(NULL),
  e12Tvhl1_E_Hist_SF_pT(NULL),  
  e12Tvhl1_E_Hist_SF_eta_errUp(NULL),
  e12Tvhl1_E_Hist_SF_pT_errUp(NULL),  
  e12Tvhl1_E_Hist_SF_eta_errDown(NULL),
  e12Tvhl1_E_Hist_SF_pT_errDown(NULL),
  
  e12Tvhm1_A_Hist_SF_eta(NULL),
  e12Tvhm1_A_Hist_SF_pT(NULL),  
  e12Tvhm1_A_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_A_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_A_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_A_Hist_SF_pT_errDown(NULL),
    
  e12Tvhm1_B_Hist_SF_eta(NULL),
  e12Tvhm1_B_Hist_SF_pT(NULL),  
  e12Tvhm1_B_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_B_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_B_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_B_Hist_SF_pT_errDown(NULL),
    
  e12Tvhm1_C_Hist_SF_eta(NULL),
  e12Tvhm1_C_Hist_SF_pT(NULL),  
  e12Tvhm1_C_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_C_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_C_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_C_Hist_SF_pT_errDown(NULL),
    
  e12Tvhm1_D_Hist_SF_eta(NULL),
  e12Tvhm1_D_Hist_SF_pT(NULL),  
  e12Tvhm1_D_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_D_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_D_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_D_Hist_SF_pT_errDown(NULL),
    
  e12Tvhm1_E_Hist_SF_eta(NULL),
  e12Tvhm1_E_Hist_SF_pT(NULL),  
  e12Tvhm1_E_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_E_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_E_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_E_Hist_SF_pT_errDown(NULL),
  
  e12Tvhm1_G_Hist_SF_eta(NULL),
  e12Tvhm1_G_Hist_SF_pT(NULL),  
  e12Tvhm1_G_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_G_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_G_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_G_Hist_SF_pT_errDown(NULL),
  
  e12Tvhm1_H_Hist_SF_eta(NULL),
  e12Tvhm1_H_Hist_SF_pT(NULL),  
  e12Tvhm1_H_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_H_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_H_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_H_Hist_SF_pT_errDown(NULL),
  
  e12Tvhm1_I_Hist_SF_eta(NULL),
  e12Tvhm1_I_Hist_SF_pT(NULL),  
  e12Tvhm1_I_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_I_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_I_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_I_Hist_SF_pT_errDown(NULL),
  
  e12Tvhm1_J_Hist_SF_eta(NULL),
  e12Tvhm1_J_Hist_SF_pT(NULL),  
  e12Tvhm1_J_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_J_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_J_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_J_Hist_SF_pT_errDown(NULL),
  
  e12Tvhm1_L_Hist_SF_eta(NULL),
  e12Tvhm1_L_Hist_SF_pT(NULL),  
  e12Tvhm1_L_Hist_SF_eta_errUp(NULL),
  e12Tvhm1_L_Hist_SF_pT_errUp(NULL),  
  e12Tvhm1_L_Hist_SF_eta_errDown(NULL),
  e12Tvhm1_L_Hist_SF_pT_errDown(NULL),
  
  mu8_AB_Hist_SF_eta_Barrel(NULL),
  mu8_AB_Hist_SF_pT_Barrel(NULL),  
  mu8_AB_Hist_SF_eta_Barrel_errUp(NULL),
  mu8_AB_Hist_SF_pT_Barrel_errUp(NULL),
  mu8_AB_Hist_SF_eta_Barrel_errDown(NULL),
  mu8_AB_Hist_SF_pT_Barrel_errDown(NULL),

  mu8_AB_Hist_SF_eta_Endcap(NULL),
  mu8_AB_Hist_SF_pT_Endcap(NULL),
  mu8_AB_Hist_SF_eta_Endcap_errUp(NULL),
  mu8_AB_Hist_SF_pT_Endcap_errUp(NULL),  
  mu8_AB_Hist_SF_eta_Endcap_errDown(NULL),
  mu8_AB_Hist_SF_pT_Endcap_errDown(NULL),
  
  mu8_C_Hist_SF_eta_Barrel(NULL),
  mu8_C_Hist_SF_pT_Barrel(NULL),  
  mu8_C_Hist_SF_eta_Barrel_errUp(NULL),
  mu8_C_Hist_SF_pT_Barrel_errUp(NULL),
  mu8_C_Hist_SF_eta_Barrel_errDown(NULL),
  mu8_C_Hist_SF_pT_Barrel_errDown(NULL),

  mu8_C_Hist_SF_eta_Endcap(NULL),
  mu8_C_Hist_SF_pT_Endcap(NULL),
  mu8_C_Hist_SF_eta_Endcap_errUp(NULL),
  mu8_C_Hist_SF_pT_Endcap_errUp(NULL),  
  mu8_C_Hist_SF_eta_Endcap_errDown(NULL),
  mu8_C_Hist_SF_pT_Endcap_errDown(NULL),
  
  mu8_D_Hist_SF_eta_Barrel(NULL),
  mu8_D_Hist_SF_pT_Barrel(NULL),  
  mu8_D_Hist_SF_eta_Barrel_errUp(NULL),
  mu8_D_Hist_SF_pT_Barrel_errUp(NULL),
  mu8_D_Hist_SF_eta_Barrel_errDown(NULL),
  mu8_D_Hist_SF_pT_Barrel_errDown(NULL),

  mu8_D_Hist_SF_eta_Endcap(NULL),
  mu8_D_Hist_SF_pT_Endcap(NULL),
  mu8_D_Hist_SF_eta_Endcap_errUp(NULL),
  mu8_D_Hist_SF_pT_Endcap_errUp(NULL),  
  mu8_D_Hist_SF_eta_Endcap_errDown(NULL),
  mu8_D_Hist_SF_pT_Endcap_errDown(NULL),
  
  mu8_E_Hist_SF_eta_Barrel(NULL),
  mu8_E_Hist_SF_pT_Barrel(NULL),  
  mu8_E_Hist_SF_eta_Barrel_errUp(NULL),
  mu8_E_Hist_SF_pT_Barrel_errUp(NULL),
  mu8_E_Hist_SF_eta_Barrel_errDown(NULL),
  mu8_E_Hist_SF_pT_Barrel_errDown(NULL),

  mu8_E_Hist_SF_eta_Endcap(NULL),
  mu8_E_Hist_SF_pT_Endcap(NULL),
  mu8_E_Hist_SF_eta_Endcap_errUp(NULL),
  mu8_E_Hist_SF_pT_Endcap_errUp(NULL),  
  mu8_E_Hist_SF_eta_Endcap_errDown(NULL),
  mu8_E_Hist_SF_pT_Endcap_errDown(NULL),

  // Data efficiency histograms

  e12Tvhl1_AB_Hist_Eff_Data_eta(NULL),
  e12Tvhl1_AB_Hist_Eff_Data_pT(NULL),  
  e12Tvhl1_AB_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhl1_AB_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhl1_AB_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhl1_AB_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhl1_C_Hist_Eff_Data_eta(NULL),
  e12Tvhl1_C_Hist_Eff_Data_pT(NULL),  
  e12Tvhl1_C_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhl1_C_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhl1_C_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhl1_C_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhl1_D_Hist_Eff_Data_eta(NULL),
  e12Tvhl1_D_Hist_Eff_Data_pT(NULL),  
  e12Tvhl1_D_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhl1_D_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhl1_D_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhl1_D_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhl1_E_Hist_Eff_Data_eta(NULL),
  e12Tvhl1_E_Hist_Eff_Data_pT(NULL),  
  e12Tvhl1_E_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhl1_E_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhl1_E_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhl1_E_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_A_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_A_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_A_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_A_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_A_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_A_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_B_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_B_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_B_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_B_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_B_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_B_Hist_Eff_Data_pT_errDown(NULL),

  e12Tvhm1_C_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_C_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_C_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_C_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_C_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_C_Hist_Eff_Data_pT_errDown(NULL),  

  e12Tvhm1_D_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_D_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_D_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_D_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_D_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_D_Hist_Eff_Data_pT_errDown(NULL),  

  e12Tvhm1_E_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_E_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_E_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_E_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_E_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_E_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_G_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_G_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_G_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_G_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_G_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_G_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_H_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_H_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_H_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_H_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_H_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_H_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_I_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_I_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_I_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_I_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_I_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_I_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_J_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_J_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_J_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_J_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_J_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_J_Hist_Eff_Data_pT_errDown(NULL),
  
  e12Tvhm1_L_Hist_Eff_Data_eta(NULL),
  e12Tvhm1_L_Hist_Eff_Data_pT(NULL),  
  e12Tvhm1_L_Hist_Eff_Data_eta_errUp(NULL),
  e12Tvhm1_L_Hist_Eff_Data_pT_errUp(NULL),
  e12Tvhm1_L_Hist_Eff_Data_eta_errDown(NULL),
  e12Tvhm1_L_Hist_Eff_Data_pT_errDown(NULL),
    
  mu8_AB_Hist_Eff_Data_eta_Barrel(NULL),
  mu8_AB_Hist_Eff_Data_pT_Barrel(NULL),  
  mu8_AB_Hist_Eff_Data_eta_Barrel_errUp(NULL),
  mu8_AB_Hist_Eff_Data_pT_Barrel_errUp(NULL),
  mu8_AB_Hist_Eff_Data_eta_Barrel_errDown(NULL),
  mu8_AB_Hist_Eff_Data_pT_Barrel_errDown(NULL),

  mu8_AB_Hist_Eff_Data_eta_Endcap(NULL),
  mu8_AB_Hist_Eff_Data_pT_Endcap(NULL),
  mu8_AB_Hist_Eff_Data_eta_Endcap_errUp(NULL),
  mu8_AB_Hist_Eff_Data_pT_Endcap_errUp(NULL),  
  mu8_AB_Hist_Eff_Data_eta_Endcap_errDown(NULL),
  mu8_AB_Hist_Eff_Data_pT_Endcap_errDown(NULL),
  
  mu8_C_Hist_Eff_Data_eta_Barrel(NULL),
  mu8_C_Hist_Eff_Data_pT_Barrel(NULL),  
  mu8_C_Hist_Eff_Data_eta_Barrel_errUp(NULL),
  mu8_C_Hist_Eff_Data_pT_Barrel_errUp(NULL),
  mu8_C_Hist_Eff_Data_eta_Barrel_errDown(NULL),
  mu8_C_Hist_Eff_Data_pT_Barrel_errDown(NULL),

  mu8_C_Hist_Eff_Data_eta_Endcap(NULL),
  mu8_C_Hist_Eff_Data_pT_Endcap(NULL),
  mu8_C_Hist_Eff_Data_eta_Endcap_errUp(NULL),
  mu8_C_Hist_Eff_Data_pT_Endcap_errUp(NULL),  
  mu8_C_Hist_Eff_Data_eta_Endcap_errDown(NULL),
  mu8_C_Hist_Eff_Data_pT_Endcap_errDown(NULL),
  
  mu8_D_Hist_Eff_Data_eta_Barrel(NULL),
  mu8_D_Hist_Eff_Data_pT_Barrel(NULL),  
  mu8_D_Hist_Eff_Data_eta_Barrel_errUp(NULL),
  mu8_D_Hist_Eff_Data_pT_Barrel_errUp(NULL),
  mu8_D_Hist_Eff_Data_eta_Barrel_errDown(NULL),
  mu8_D_Hist_Eff_Data_pT_Barrel_errDown(NULL),

  mu8_D_Hist_Eff_Data_eta_Endcap(NULL),
  mu8_D_Hist_Eff_Data_pT_Endcap(NULL),
  mu8_D_Hist_Eff_Data_eta_Endcap_errUp(NULL),
  mu8_D_Hist_Eff_Data_pT_Endcap_errUp(NULL),  
  mu8_D_Hist_Eff_Data_eta_Endcap_errDown(NULL),
  mu8_D_Hist_Eff_Data_pT_Endcap_errDown(NULL),
  
  mu8_E_Hist_Eff_Data_eta_Barrel(NULL),
  mu8_E_Hist_Eff_Data_pT_Barrel(NULL),  
  mu8_E_Hist_Eff_Data_eta_Barrel_errUp(NULL),
  mu8_E_Hist_Eff_Data_pT_Barrel_errUp(NULL),
  mu8_E_Hist_Eff_Data_eta_Barrel_errDown(NULL),
  mu8_E_Hist_Eff_Data_pT_Barrel_errDown(NULL),

  mu8_E_Hist_Eff_Data_eta_Endcap(NULL),
  mu8_E_Hist_Eff_Data_pT_Endcap(NULL),
  mu8_E_Hist_Eff_Data_eta_Endcap_errUp(NULL),
  mu8_E_Hist_Eff_Data_pT_Endcap_errUp(NULL),  
  mu8_E_Hist_Eff_Data_eta_Endcap_errDown(NULL),
  mu8_E_Hist_Eff_Data_pT_Endcap_errDown(NULL),

  // MC efficiency histograms

  e12Tvhl1_AB_Hist_Eff_MC_eta(NULL),
  e12Tvhl1_AB_Hist_Eff_MC_pT(NULL),  
  e12Tvhl1_AB_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhl1_AB_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhl1_AB_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhl1_AB_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhl1_C_Hist_Eff_MC_eta(NULL),
  e12Tvhl1_C_Hist_Eff_MC_pT(NULL),  
  e12Tvhl1_C_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhl1_C_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhl1_C_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhl1_C_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhl1_D_Hist_Eff_MC_eta(NULL),
  e12Tvhl1_D_Hist_Eff_MC_pT(NULL),  
  e12Tvhl1_D_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhl1_D_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhl1_D_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhl1_D_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhl1_E_Hist_Eff_MC_eta(NULL),
  e12Tvhl1_E_Hist_Eff_MC_pT(NULL),  
  e12Tvhl1_E_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhl1_E_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhl1_E_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhl1_E_Hist_Eff_MC_pT_errDown(NULL),

  e12Tvhm1_A_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_A_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_A_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_A_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_A_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_A_Hist_Eff_MC_pT_errDown(NULL),
    
  e12Tvhm1_B_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_B_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_B_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_B_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_B_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_B_Hist_Eff_MC_pT_errDown(NULL),
    
  e12Tvhm1_C_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_C_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_C_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_C_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_C_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_C_Hist_Eff_MC_pT_errDown(NULL),
    
  e12Tvhm1_D_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_D_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_D_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_D_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_D_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_D_Hist_Eff_MC_pT_errDown(NULL),
    
  e12Tvhm1_E_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_E_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_E_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_E_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_E_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_E_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhm1_G_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_G_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_G_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_G_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_G_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_G_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhm1_H_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_H_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_H_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_H_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_H_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_H_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhm1_I_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_I_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_I_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_I_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_I_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_I_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhm1_J_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_J_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_J_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_J_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_J_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_J_Hist_Eff_MC_pT_errDown(NULL),
  
  e12Tvhm1_L_Hist_Eff_MC_eta(NULL),
  e12Tvhm1_L_Hist_Eff_MC_pT(NULL),  
  e12Tvhm1_L_Hist_Eff_MC_eta_errUp(NULL),
  e12Tvhm1_L_Hist_Eff_MC_pT_errUp(NULL),  
  e12Tvhm1_L_Hist_Eff_MC_eta_errDown(NULL),
  e12Tvhm1_L_Hist_Eff_MC_pT_errDown(NULL),
  
  mu8_AB_Hist_Eff_MC_eta_Barrel(NULL),
  mu8_AB_Hist_Eff_MC_pT_Barrel(NULL),  
  mu8_AB_Hist_Eff_MC_eta_Barrel_errUp(NULL),
  mu8_AB_Hist_Eff_MC_pT_Barrel_errUp(NULL),
  mu8_AB_Hist_Eff_MC_eta_Barrel_errDown(NULL),
  mu8_AB_Hist_Eff_MC_pT_Barrel_errDown(NULL),

  mu8_AB_Hist_Eff_MC_eta_Endcap(NULL),
  mu8_AB_Hist_Eff_MC_pT_Endcap(NULL),
  mu8_AB_Hist_Eff_MC_eta_Endcap_errUp(NULL),
  mu8_AB_Hist_Eff_MC_pT_Endcap_errUp(NULL),  
  mu8_AB_Hist_Eff_MC_eta_Endcap_errDown(NULL),
  mu8_AB_Hist_Eff_MC_pT_Endcap_errDown(NULL),
  
  mu8_C_Hist_Eff_MC_eta_Barrel(NULL),
  mu8_C_Hist_Eff_MC_pT_Barrel(NULL),  
  mu8_C_Hist_Eff_MC_eta_Barrel_errUp(NULL),
  mu8_C_Hist_Eff_MC_pT_Barrel_errUp(NULL),
  mu8_C_Hist_Eff_MC_eta_Barrel_errDown(NULL),
  mu8_C_Hist_Eff_MC_pT_Barrel_errDown(NULL),

  mu8_C_Hist_Eff_MC_eta_Endcap(NULL),
  mu8_C_Hist_Eff_MC_pT_Endcap(NULL),
  mu8_C_Hist_Eff_MC_eta_Endcap_errUp(NULL),
  mu8_C_Hist_Eff_MC_pT_Endcap_errUp(NULL),  
  mu8_C_Hist_Eff_MC_eta_Endcap_errDown(NULL),
  mu8_C_Hist_Eff_MC_pT_Endcap_errDown(NULL),
  
  mu8_D_Hist_Eff_MC_eta_Barrel(NULL),
  mu8_D_Hist_Eff_MC_pT_Barrel(NULL),  
  mu8_D_Hist_Eff_MC_eta_Barrel_errUp(NULL),
  mu8_D_Hist_Eff_MC_pT_Barrel_errUp(NULL),
  mu8_D_Hist_Eff_MC_eta_Barrel_errDown(NULL),
  mu8_D_Hist_Eff_MC_pT_Barrel_errDown(NULL),

  mu8_D_Hist_Eff_MC_eta_Endcap(NULL),
  mu8_D_Hist_Eff_MC_pT_Endcap(NULL),
  mu8_D_Hist_Eff_MC_eta_Endcap_errUp(NULL),
  mu8_D_Hist_Eff_MC_pT_Endcap_errUp(NULL),  
  mu8_D_Hist_Eff_MC_eta_Endcap_errDown(NULL),
  mu8_D_Hist_Eff_MC_pT_Endcap_errDown(NULL),
  
  mu8_E_Hist_Eff_MC_eta_Barrel(NULL),
  mu8_E_Hist_Eff_MC_pT_Barrel(NULL),  
  mu8_E_Hist_Eff_MC_eta_Barrel_errUp(NULL),
  mu8_E_Hist_Eff_MC_pT_Barrel_errUp(NULL),
  mu8_E_Hist_Eff_MC_eta_Barrel_errDown(NULL),
  mu8_E_Hist_Eff_MC_pT_Barrel_errDown(NULL),

  mu8_E_Hist_Eff_MC_eta_Endcap(NULL),
  mu8_E_Hist_Eff_MC_pT_Endcap(NULL),
  mu8_E_Hist_Eff_MC_eta_Endcap_errUp(NULL),
  mu8_E_Hist_Eff_MC_pT_Endcap_errUp(NULL),  
  mu8_E_Hist_Eff_MC_eta_Endcap_errDown(NULL),
  mu8_E_Hist_Eff_MC_pT_Endcap_errDown(NULL)

{

    doDebug = debug ;
    if(path.size() > 0 ) loadFiles( path ) ;

}

HSG4LepLepTriggerSF::~HSG4LepLepTriggerSF() 
{

  // scale factor histograms

  delete e12Tvhl1_AB_Hist_SF_eta;
  delete e12Tvhl1_AB_Hist_SF_pT;  
  delete e12Tvhl1_AB_Hist_SF_eta_errUp;
  delete e12Tvhl1_AB_Hist_SF_pT_errUp;  
  delete e12Tvhl1_AB_Hist_SF_eta_errDown;
  delete e12Tvhl1_AB_Hist_SF_pT_errDown;
  
  delete e12Tvhl1_C_Hist_SF_eta;
  delete e12Tvhl1_C_Hist_SF_pT;  
  delete e12Tvhl1_C_Hist_SF_eta_errUp;
  delete e12Tvhl1_C_Hist_SF_pT_errUp;  
  delete e12Tvhl1_C_Hist_SF_eta_errDown;
  delete e12Tvhl1_C_Hist_SF_pT_errDown;
  
  delete e12Tvhl1_D_Hist_SF_eta;
  delete e12Tvhl1_D_Hist_SF_pT;  
  delete e12Tvhl1_D_Hist_SF_eta_errUp;
  delete e12Tvhl1_D_Hist_SF_pT_errUp;  
  delete e12Tvhl1_D_Hist_SF_eta_errDown;
  delete e12Tvhl1_D_Hist_SF_pT_errDown;
  
  delete e12Tvhl1_E_Hist_SF_eta;
  delete e12Tvhl1_E_Hist_SF_pT;  
  delete e12Tvhl1_E_Hist_SF_eta_errUp;
  delete e12Tvhl1_E_Hist_SF_pT_errUp;  
  delete e12Tvhl1_E_Hist_SF_eta_errDown;
  delete e12Tvhl1_E_Hist_SF_pT_errDown;

  delete e12Tvhm1_A_Hist_SF_eta;
  delete e12Tvhm1_A_Hist_SF_pT;  
  delete e12Tvhm1_A_Hist_SF_eta_errUp;
  delete e12Tvhm1_A_Hist_SF_pT_errUp;  
  delete e12Tvhm1_A_Hist_SF_eta_errDown;
  delete e12Tvhm1_A_Hist_SF_pT_errDown;
  
  delete e12Tvhm1_B_Hist_SF_eta;
  delete e12Tvhm1_B_Hist_SF_pT;  
  delete e12Tvhm1_B_Hist_SF_eta_errUp;
  delete e12Tvhm1_B_Hist_SF_pT_errUp;  
  delete e12Tvhm1_B_Hist_SF_eta_errDown;
  delete e12Tvhm1_B_Hist_SF_pT_errDown;
  
  delete e12Tvhm1_C_Hist_SF_eta;
  delete e12Tvhm1_C_Hist_SF_pT;  
  delete e12Tvhm1_C_Hist_SF_eta_errUp;
  delete e12Tvhm1_C_Hist_SF_pT_errUp;  
  delete e12Tvhm1_C_Hist_SF_eta_errDown;
  delete e12Tvhm1_C_Hist_SF_pT_errDown;
  
  delete e12Tvhm1_D_Hist_SF_eta;
  delete e12Tvhm1_D_Hist_SF_pT;  
  delete e12Tvhm1_D_Hist_SF_eta_errUp;
  delete e12Tvhm1_D_Hist_SF_pT_errUp;  
  delete e12Tvhm1_D_Hist_SF_eta_errDown;
  delete e12Tvhm1_D_Hist_SF_pT_errDown;  

  delete e12Tvhm1_E_Hist_SF_eta;
  delete e12Tvhm1_E_Hist_SF_pT;  
  delete e12Tvhm1_E_Hist_SF_eta_errUp;
  delete e12Tvhm1_E_Hist_SF_pT_errUp;  
  delete e12Tvhm1_E_Hist_SF_eta_errDown;
  delete e12Tvhm1_E_Hist_SF_pT_errDown;

  delete e12Tvhm1_G_Hist_SF_eta;
  delete e12Tvhm1_G_Hist_SF_pT;  
  delete e12Tvhm1_G_Hist_SF_eta_errUp;
  delete e12Tvhm1_G_Hist_SF_pT_errUp;  
  delete e12Tvhm1_G_Hist_SF_eta_errDown;
  delete e12Tvhm1_G_Hist_SF_pT_errDown;

  delete e12Tvhm1_H_Hist_SF_eta;
  delete e12Tvhm1_H_Hist_SF_pT;  
  delete e12Tvhm1_H_Hist_SF_eta_errUp;
  delete e12Tvhm1_H_Hist_SF_pT_errUp;  
  delete e12Tvhm1_H_Hist_SF_eta_errDown;
  delete e12Tvhm1_H_Hist_SF_pT_errDown;

  delete e12Tvhm1_I_Hist_SF_eta;
  delete e12Tvhm1_I_Hist_SF_pT;  
  delete e12Tvhm1_I_Hist_SF_eta_errUp;
  delete e12Tvhm1_I_Hist_SF_pT_errUp;  
  delete e12Tvhm1_I_Hist_SF_eta_errDown;
  delete e12Tvhm1_I_Hist_SF_pT_errDown;

  delete e12Tvhm1_J_Hist_SF_eta;
  delete e12Tvhm1_J_Hist_SF_pT;  
  delete e12Tvhm1_J_Hist_SF_eta_errUp;
  delete e12Tvhm1_J_Hist_SF_pT_errUp;  
  delete e12Tvhm1_J_Hist_SF_eta_errDown;
  delete e12Tvhm1_J_Hist_SF_pT_errDown;

  delete e12Tvhm1_L_Hist_SF_eta;
  delete e12Tvhm1_L_Hist_SF_pT;  
  delete e12Tvhm1_L_Hist_SF_eta_errUp;
  delete e12Tvhm1_L_Hist_SF_pT_errUp;  
  delete e12Tvhm1_L_Hist_SF_eta_errDown;
  delete e12Tvhm1_L_Hist_SF_pT_errDown;
  
  delete mu8_AB_Hist_SF_eta_Barrel;
  delete mu8_AB_Hist_SF_pT_Barrel;
  delete mu8_AB_Hist_SF_eta_Barrel_errUp;
  delete mu8_AB_Hist_SF_pT_Barrel_errUp;
  delete mu8_AB_Hist_SF_eta_Barrel_errDown;
  delete mu8_AB_Hist_SF_pT_Barrel_errDown;

  delete mu8_AB_Hist_SF_eta_Endcap;
  delete mu8_AB_Hist_SF_pT_Endcap;
  delete mu8_AB_Hist_SF_eta_Endcap_errUp;
  delete mu8_AB_Hist_SF_pT_Endcap_errUp;
  delete mu8_AB_Hist_SF_eta_Endcap_errDown;
  delete mu8_AB_Hist_SF_pT_Endcap_errDown;

  delete mu8_C_Hist_SF_eta_Barrel;
  delete mu8_C_Hist_SF_pT_Barrel;
  delete mu8_C_Hist_SF_eta_Barrel_errUp;
  delete mu8_C_Hist_SF_pT_Barrel_errUp;
  delete mu8_C_Hist_SF_eta_Barrel_errDown;
  delete mu8_C_Hist_SF_pT_Barrel_errDown;

  delete mu8_C_Hist_SF_eta_Endcap;
  delete mu8_C_Hist_SF_pT_Endcap;
  delete mu8_C_Hist_SF_eta_Endcap_errUp;
  delete mu8_C_Hist_SF_pT_Endcap_errUp;
  delete mu8_C_Hist_SF_eta_Endcap_errDown;
  delete mu8_C_Hist_SF_pT_Endcap_errDown;

  delete mu8_D_Hist_SF_eta_Barrel;
  delete mu8_D_Hist_SF_pT_Barrel;
  delete mu8_D_Hist_SF_eta_Barrel_errUp;
  delete mu8_D_Hist_SF_pT_Barrel_errUp;
  delete mu8_D_Hist_SF_eta_Barrel_errDown;
  delete mu8_D_Hist_SF_pT_Barrel_errDown;

  delete mu8_D_Hist_SF_eta_Endcap;
  delete mu8_D_Hist_SF_pT_Endcap;
  delete mu8_D_Hist_SF_eta_Endcap_errUp;
  delete mu8_D_Hist_SF_pT_Endcap_errUp;
  delete mu8_D_Hist_SF_eta_Endcap_errDown;
  delete mu8_D_Hist_SF_pT_Endcap_errDown;

  delete mu8_E_Hist_SF_eta_Barrel;
  delete mu8_E_Hist_SF_pT_Barrel;
  delete mu8_E_Hist_SF_eta_Barrel_errUp;
  delete mu8_E_Hist_SF_pT_Barrel_errUp;
  delete mu8_E_Hist_SF_eta_Barrel_errDown;
  delete mu8_E_Hist_SF_pT_Barrel_errDown;

  delete mu8_E_Hist_SF_eta_Endcap;
  delete mu8_E_Hist_SF_pT_Endcap;
  delete mu8_E_Hist_SF_eta_Endcap_errUp;
  delete mu8_E_Hist_SF_pT_Endcap_errUp;
  delete mu8_E_Hist_SF_eta_Endcap_errDown;
  delete mu8_E_Hist_SF_pT_Endcap_errDown;

  // data efficiency histograms

  delete e12Tvhl1_AB_Hist_Eff_Data_eta;
  delete e12Tvhl1_AB_Hist_Eff_Data_pT;  
  delete e12Tvhl1_AB_Hist_Eff_Data_eta_errUp;
  delete e12Tvhl1_AB_Hist_Eff_Data_pT_errUp;
  delete e12Tvhl1_AB_Hist_Eff_Data_eta_errDown;
  delete e12Tvhl1_AB_Hist_Eff_Data_pT_errDown;
  
  delete e12Tvhl1_C_Hist_Eff_Data_eta;
  delete e12Tvhl1_C_Hist_Eff_Data_pT;  
  delete e12Tvhl1_C_Hist_Eff_Data_eta_errUp;
  delete e12Tvhl1_C_Hist_Eff_Data_pT_errUp;
  delete e12Tvhl1_C_Hist_Eff_Data_eta_errDown;
  delete e12Tvhl1_C_Hist_Eff_Data_pT_errDown;
  
  delete e12Tvhl1_D_Hist_Eff_Data_eta;
  delete e12Tvhl1_D_Hist_Eff_Data_pT;  
  delete e12Tvhl1_D_Hist_Eff_Data_eta_errUp;
  delete e12Tvhl1_D_Hist_Eff_Data_pT_errUp;
  delete e12Tvhl1_D_Hist_Eff_Data_eta_errDown;
  delete e12Tvhl1_D_Hist_Eff_Data_pT_errDown;
  
  delete e12Tvhl1_E_Hist_Eff_Data_eta;
  delete e12Tvhl1_E_Hist_Eff_Data_pT;  
  delete e12Tvhl1_E_Hist_Eff_Data_eta_errUp;
  delete e12Tvhl1_E_Hist_Eff_Data_pT_errUp;
  delete e12Tvhl1_E_Hist_Eff_Data_eta_errDown;
  delete e12Tvhl1_E_Hist_Eff_Data_pT_errDown;

  delete e12Tvhm1_A_Hist_Eff_Data_eta;
  delete e12Tvhm1_A_Hist_Eff_Data_pT;  
  delete e12Tvhm1_A_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_A_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_A_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_A_Hist_Eff_Data_pT_errDown;
  
  delete e12Tvhm1_B_Hist_Eff_Data_eta;
  delete e12Tvhm1_B_Hist_Eff_Data_pT;  
  delete e12Tvhm1_B_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_B_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_B_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_B_Hist_Eff_Data_pT_errDown;
  
  delete e12Tvhm1_C_Hist_Eff_Data_eta;
  delete e12Tvhm1_C_Hist_Eff_Data_pT;  
  delete e12Tvhm1_C_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_C_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_C_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_C_Hist_Eff_Data_pT_errDown;  

  delete e12Tvhm1_D_Hist_Eff_Data_eta;
  delete e12Tvhm1_D_Hist_Eff_Data_pT;  
  delete e12Tvhm1_D_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_D_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_D_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_D_Hist_Eff_Data_pT_errDown;  

  delete e12Tvhm1_E_Hist_Eff_Data_eta;
  delete e12Tvhm1_E_Hist_Eff_Data_pT;  
  delete e12Tvhm1_E_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_E_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_E_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_E_Hist_Eff_Data_pT_errDown;

  delete e12Tvhm1_G_Hist_Eff_Data_eta;
  delete e12Tvhm1_G_Hist_Eff_Data_pT;  
  delete e12Tvhm1_G_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_G_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_G_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_G_Hist_Eff_Data_pT_errDown;

  delete e12Tvhm1_H_Hist_Eff_Data_eta;
  delete e12Tvhm1_H_Hist_Eff_Data_pT;  
  delete e12Tvhm1_H_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_H_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_H_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_H_Hist_Eff_Data_pT_errDown;

  delete e12Tvhm1_I_Hist_Eff_Data_eta;
  delete e12Tvhm1_I_Hist_Eff_Data_pT;  
  delete e12Tvhm1_I_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_I_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_I_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_I_Hist_Eff_Data_pT_errDown;

  delete e12Tvhm1_J_Hist_Eff_Data_eta;
  delete e12Tvhm1_J_Hist_Eff_Data_pT;  
  delete e12Tvhm1_J_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_J_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_J_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_J_Hist_Eff_Data_pT_errDown;

  delete e12Tvhm1_L_Hist_Eff_Data_eta;
  delete e12Tvhm1_L_Hist_Eff_Data_pT;  
  delete e12Tvhm1_L_Hist_Eff_Data_eta_errUp;
  delete e12Tvhm1_L_Hist_Eff_Data_pT_errUp;
  delete e12Tvhm1_L_Hist_Eff_Data_eta_errDown;
  delete e12Tvhm1_L_Hist_Eff_Data_pT_errDown;
  
  delete mu8_AB_Hist_Eff_Data_eta_Barrel;
  delete mu8_AB_Hist_Eff_Data_pT_Barrel;
  delete mu8_AB_Hist_Eff_Data_eta_Barrel_errUp;
  delete mu8_AB_Hist_Eff_Data_pT_Barrel_errUp;
  delete mu8_AB_Hist_Eff_Data_eta_Barrel_errDown;
  delete mu8_AB_Hist_Eff_Data_pT_Barrel_errDown;

  delete mu8_AB_Hist_Eff_Data_eta_Endcap;
  delete mu8_AB_Hist_Eff_Data_pT_Endcap;
  delete mu8_AB_Hist_Eff_Data_eta_Endcap_errUp;
  delete mu8_AB_Hist_Eff_Data_pT_Endcap_errUp;
  delete mu8_AB_Hist_Eff_Data_eta_Endcap_errDown;
  delete mu8_AB_Hist_Eff_Data_pT_Endcap_errDown;
  
  delete mu8_C_Hist_Eff_Data_eta_Barrel;
  delete mu8_C_Hist_Eff_Data_pT_Barrel;
  delete mu8_C_Hist_Eff_Data_eta_Barrel_errUp;
  delete mu8_C_Hist_Eff_Data_pT_Barrel_errUp;
  delete mu8_C_Hist_Eff_Data_eta_Barrel_errDown;
  delete mu8_C_Hist_Eff_Data_pT_Barrel_errDown;

  delete mu8_C_Hist_Eff_Data_eta_Endcap;
  delete mu8_C_Hist_Eff_Data_pT_Endcap;
  delete mu8_C_Hist_Eff_Data_eta_Endcap_errUp;
  delete mu8_C_Hist_Eff_Data_pT_Endcap_errUp;
  delete mu8_C_Hist_Eff_Data_eta_Endcap_errDown;
  delete mu8_C_Hist_Eff_Data_pT_Endcap_errDown;
  
  delete mu8_D_Hist_Eff_Data_eta_Barrel;
  delete mu8_D_Hist_Eff_Data_pT_Barrel;
  delete mu8_D_Hist_Eff_Data_eta_Barrel_errUp;
  delete mu8_D_Hist_Eff_Data_pT_Barrel_errUp;
  delete mu8_D_Hist_Eff_Data_eta_Barrel_errDown;
  delete mu8_D_Hist_Eff_Data_pT_Barrel_errDown;

  delete mu8_D_Hist_Eff_Data_eta_Endcap;
  delete mu8_D_Hist_Eff_Data_pT_Endcap;
  delete mu8_D_Hist_Eff_Data_eta_Endcap_errUp;
  delete mu8_D_Hist_Eff_Data_pT_Endcap_errUp;
  delete mu8_D_Hist_Eff_Data_eta_Endcap_errDown;
  delete mu8_D_Hist_Eff_Data_pT_Endcap_errDown;
  
  delete mu8_E_Hist_Eff_Data_eta_Barrel;
  delete mu8_E_Hist_Eff_Data_pT_Barrel;
  delete mu8_E_Hist_Eff_Data_eta_Barrel_errUp;
  delete mu8_E_Hist_Eff_Data_pT_Barrel_errUp;
  delete mu8_E_Hist_Eff_Data_eta_Barrel_errDown;
  delete mu8_E_Hist_Eff_Data_pT_Barrel_errDown;

  delete mu8_E_Hist_Eff_Data_eta_Endcap;
  delete mu8_E_Hist_Eff_Data_pT_Endcap;
  delete mu8_E_Hist_Eff_Data_eta_Endcap_errUp;
  delete mu8_E_Hist_Eff_Data_pT_Endcap_errUp;
  delete mu8_E_Hist_Eff_Data_eta_Endcap_errDown;
  delete mu8_E_Hist_Eff_Data_pT_Endcap_errDown;

  // mc efficiency histograms
  
  delete e12Tvhl1_AB_Hist_Eff_MC_eta;
  delete e12Tvhl1_AB_Hist_Eff_MC_pT;  
  delete e12Tvhl1_AB_Hist_Eff_MC_eta_errUp;
  delete e12Tvhl1_AB_Hist_Eff_MC_pT_errUp;
  delete e12Tvhl1_AB_Hist_Eff_MC_eta_errDown;
  delete e12Tvhl1_AB_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhl1_C_Hist_Eff_MC_eta;
  delete e12Tvhl1_C_Hist_Eff_MC_pT;  
  delete e12Tvhl1_C_Hist_Eff_MC_eta_errUp;
  delete e12Tvhl1_C_Hist_Eff_MC_pT_errUp;
  delete e12Tvhl1_C_Hist_Eff_MC_eta_errDown;
  delete e12Tvhl1_C_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhl1_D_Hist_Eff_MC_eta;
  delete e12Tvhl1_D_Hist_Eff_MC_pT;  
  delete e12Tvhl1_D_Hist_Eff_MC_eta_errUp;
  delete e12Tvhl1_D_Hist_Eff_MC_pT_errUp;
  delete e12Tvhl1_D_Hist_Eff_MC_eta_errDown;
  delete e12Tvhl1_D_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhl1_E_Hist_Eff_MC_eta;
  delete e12Tvhl1_E_Hist_Eff_MC_pT;  
  delete e12Tvhl1_E_Hist_Eff_MC_eta_errUp;
  delete e12Tvhl1_E_Hist_Eff_MC_pT_errUp;
  delete e12Tvhl1_E_Hist_Eff_MC_eta_errDown;
  delete e12Tvhl1_E_Hist_Eff_MC_pT_errDown;
 
  delete e12Tvhm1_A_Hist_Eff_MC_eta;
  delete e12Tvhm1_A_Hist_Eff_MC_pT;  
  delete e12Tvhm1_A_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_A_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_A_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_A_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_B_Hist_Eff_MC_eta;
  delete e12Tvhm1_B_Hist_Eff_MC_pT;  
  delete e12Tvhm1_B_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_B_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_B_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_B_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_C_Hist_Eff_MC_eta;
  delete e12Tvhm1_C_Hist_Eff_MC_pT;  
  delete e12Tvhm1_C_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_C_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_C_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_C_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_D_Hist_Eff_MC_eta;
  delete e12Tvhm1_D_Hist_Eff_MC_pT;  
  delete e12Tvhm1_D_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_D_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_D_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_D_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_E_Hist_Eff_MC_eta;
  delete e12Tvhm1_E_Hist_Eff_MC_pT;  
  delete e12Tvhm1_E_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_E_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_E_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_E_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_G_Hist_Eff_MC_eta;
  delete e12Tvhm1_G_Hist_Eff_MC_pT;  
  delete e12Tvhm1_G_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_G_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_G_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_G_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_H_Hist_Eff_MC_eta;
  delete e12Tvhm1_H_Hist_Eff_MC_pT;  
  delete e12Tvhm1_H_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_H_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_H_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_H_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_I_Hist_Eff_MC_eta;
  delete e12Tvhm1_I_Hist_Eff_MC_pT;  
  delete e12Tvhm1_I_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_I_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_I_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_I_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_J_Hist_Eff_MC_eta;
  delete e12Tvhm1_J_Hist_Eff_MC_pT;  
  delete e12Tvhm1_J_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_J_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_J_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_J_Hist_Eff_MC_pT_errDown;
  
  delete e12Tvhm1_L_Hist_Eff_MC_eta;
  delete e12Tvhm1_L_Hist_Eff_MC_pT;  
  delete e12Tvhm1_L_Hist_Eff_MC_eta_errUp;
  delete e12Tvhm1_L_Hist_Eff_MC_pT_errUp;
  delete e12Tvhm1_L_Hist_Eff_MC_eta_errDown;
  delete e12Tvhm1_L_Hist_Eff_MC_pT_errDown;  

  delete mu8_AB_Hist_Eff_MC_eta_Barrel;
  delete mu8_AB_Hist_Eff_MC_pT_Barrel;
  delete mu8_AB_Hist_Eff_MC_eta_Barrel_errUp;
  delete mu8_AB_Hist_Eff_MC_pT_Barrel_errUp;
  delete mu8_AB_Hist_Eff_MC_eta_Barrel_errDown;
  delete mu8_AB_Hist_Eff_MC_pT_Barrel_errDown;

  delete mu8_AB_Hist_Eff_MC_eta_Endcap;
  delete mu8_AB_Hist_Eff_MC_pT_Endcap;
  delete mu8_AB_Hist_Eff_MC_eta_Endcap_errUp;
  delete mu8_AB_Hist_Eff_MC_pT_Endcap_errUp;
  delete mu8_AB_Hist_Eff_MC_eta_Endcap_errDown;
  delete mu8_AB_Hist_Eff_MC_pT_Endcap_errDown;
  
  delete mu8_C_Hist_Eff_MC_eta_Barrel;
  delete mu8_C_Hist_Eff_MC_pT_Barrel;
  delete mu8_C_Hist_Eff_MC_eta_Barrel_errUp;
  delete mu8_C_Hist_Eff_MC_pT_Barrel_errUp;
  delete mu8_C_Hist_Eff_MC_eta_Barrel_errDown;
  delete mu8_C_Hist_Eff_MC_pT_Barrel_errDown;

  delete mu8_C_Hist_Eff_MC_eta_Endcap;
  delete mu8_C_Hist_Eff_MC_pT_Endcap;
  delete mu8_C_Hist_Eff_MC_eta_Endcap_errUp;
  delete mu8_C_Hist_Eff_MC_pT_Endcap_errUp;
  delete mu8_C_Hist_Eff_MC_eta_Endcap_errDown;
  delete mu8_C_Hist_Eff_MC_pT_Endcap_errDown;
  
  delete mu8_D_Hist_Eff_MC_eta_Barrel;
  delete mu8_D_Hist_Eff_MC_pT_Barrel;
  delete mu8_D_Hist_Eff_MC_eta_Barrel_errUp;
  delete mu8_D_Hist_Eff_MC_pT_Barrel_errUp;
  delete mu8_D_Hist_Eff_MC_eta_Barrel_errDown;
  delete mu8_D_Hist_Eff_MC_pT_Barrel_errDown;

  delete mu8_D_Hist_Eff_MC_eta_Endcap;
  delete mu8_D_Hist_Eff_MC_pT_Endcap;
  delete mu8_D_Hist_Eff_MC_eta_Endcap_errUp;
  delete mu8_D_Hist_Eff_MC_pT_Endcap_errUp;
  delete mu8_D_Hist_Eff_MC_eta_Endcap_errDown;
  delete mu8_D_Hist_Eff_MC_pT_Endcap_errDown;
  
  delete mu8_E_Hist_Eff_MC_eta_Barrel;
  delete mu8_E_Hist_Eff_MC_pT_Barrel;
  delete mu8_E_Hist_Eff_MC_eta_Barrel_errUp;
  delete mu8_E_Hist_Eff_MC_pT_Barrel_errUp;
  delete mu8_E_Hist_Eff_MC_eta_Barrel_errDown;
  delete mu8_E_Hist_Eff_MC_pT_Barrel_errDown;

  delete mu8_E_Hist_Eff_MC_eta_Endcap;
  delete mu8_E_Hist_Eff_MC_pT_Endcap;
  delete mu8_E_Hist_Eff_MC_eta_Endcap_errUp;
  delete mu8_E_Hist_Eff_MC_pT_Endcap_errUp;
  delete mu8_E_Hist_Eff_MC_eta_Endcap_errDown;
  delete mu8_E_Hist_Eff_MC_pT_Endcap_errDown;
}

void HSG4LepLepTriggerSF::loadFiles( const TString& path ) 
{

  TFile* e12Tvhl1_AB_File = TFile::Open((path+"SFe12Tvhl1.root")); //--AM-- set to more stable TFile::Open
  TFile* e12Tvhl1_C_File = TFile::Open((path+"SFe12Tvhl1_C.root"));
  TFile* e12Tvhl1_D_File = TFile::Open((path+"SFe12Tvhl1_D.root"));
  TFile* e12Tvhl1_E_File = TFile::Open((path+"SFe12Tvhl1_E.root"));

  TFile* e12Tvhm1_A_File = TFile::Open((path+"SFe12Tvhm1A_p1328.root"));
  TFile* e12Tvhm1_B_File = TFile::Open((path+"SFe12Tvhm1B_p1328.root"));
  TFile* e12Tvhm1_C_File = TFile::Open((path+"SFe12Tvhm1C_p1328.root"));
  TFile* e12Tvhm1_D_File = TFile::Open((path+"SFe12Tvhm1D_p1328.root"));
  TFile* e12Tvhm1_E_File = TFile::Open((path+"SFe12Tvhm1E_p1328.root"));
  TFile* e12Tvhm1_G_File = TFile::Open((path+"SFe12Tvhm1G_p1328.root"));
  TFile* e12Tvhm1_H_File = TFile::Open((path+"SFe12Tvhm1H_p1328.root"));
  TFile* e12Tvhm1_I_File = TFile::Open((path+"SFe12Tvhm1I_p1328.root"));
  TFile* e12Tvhm1_J_File = TFile::Open((path+"SFe12Tvhm1J_p1328.root"));
  TFile* e12Tvhm1_L_File = TFile::Open((path+"SFe12Tvhm1L_p1328.root"));

  TFile* mu8_AB_File = TFile::Open((path+"SFmu8AB.root"));
  TFile* mu8_C_File = TFile::Open((path+"SFmu8C.root"));
  TFile* mu8_D_File = TFile::Open((path+"SFmu8D.root"));
  TFile* mu8_E_File = TFile::Open((path+"SFmu8E.root"));
  
  gROOT->cd(); //--AM-- return to home dir 
    
  // scale factors

  readHisto("e12Tvhl1_AB_", (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_eta"), e12Tvhl1_AB_Hist_SF_eta, e12Tvhl1_AB_Hist_SF_eta_errUp, e12Tvhl1_AB_Hist_SF_eta_errDown);
  readHisto("e12Tvhl1_AB_", (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_pt") , e12Tvhl1_AB_Hist_SF_pT,  e12Tvhl1_AB_Hist_SF_pT_errUp,  e12Tvhl1_AB_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhl1_C_", (TGraphAsymmErrors*)e12Tvhl1_C_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_eta"), e12Tvhl1_C_Hist_SF_eta, e12Tvhl1_C_Hist_SF_eta_errUp, e12Tvhl1_C_Hist_SF_eta_errDown);
  readHisto("e12Tvhl1_C_", (TGraphAsymmErrors*)e12Tvhl1_C_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_pt") , e12Tvhl1_C_Hist_SF_pT,  e12Tvhl1_C_Hist_SF_pT_errUp,  e12Tvhl1_C_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhl1_D_", (TGraphAsymmErrors*)e12Tvhl1_D_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_eta"), e12Tvhl1_D_Hist_SF_eta, e12Tvhl1_D_Hist_SF_eta_errUp, e12Tvhl1_D_Hist_SF_eta_errDown);
  readHisto("e12Tvhl1_D_", (TGraphAsymmErrors*)e12Tvhl1_D_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_pt") , e12Tvhl1_D_Hist_SF_pT,  e12Tvhl1_D_Hist_SF_pT_errUp,  e12Tvhl1_D_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhl1_E_", (TGraphAsymmErrors*)e12Tvhl1_E_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_eta"), e12Tvhl1_E_Hist_SF_eta, e12Tvhl1_E_Hist_SF_eta_errUp, e12Tvhl1_E_Hist_SF_eta_errDown);
  readHisto("e12Tvhl1_E_", (TGraphAsymmErrors*)e12Tvhl1_E_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("SF_err_pt") , e12Tvhl1_E_Hist_SF_pT,  e12Tvhl1_E_Hist_SF_pT_errUp,  e12Tvhl1_E_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_A_", (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("SF_err_eta"), e12Tvhm1_A_Hist_SF_eta, e12Tvhm1_A_Hist_SF_eta_errUp, e12Tvhm1_A_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_A_", (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("SF_err_pt") , e12Tvhm1_A_Hist_SF_pT,  e12Tvhm1_A_Hist_SF_pT_errUp,  e12Tvhm1_A_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_B_", (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("SF_err_eta"), e12Tvhm1_B_Hist_SF_eta, e12Tvhm1_B_Hist_SF_eta_errUp, e12Tvhm1_B_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_B_", (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("SF_err_pt") , e12Tvhm1_B_Hist_SF_pT,  e12Tvhm1_B_Hist_SF_pT_errUp,  e12Tvhm1_B_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_C_", (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("SF_err_eta"), e12Tvhm1_C_Hist_SF_eta, e12Tvhm1_C_Hist_SF_eta_errUp, e12Tvhm1_C_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_C_", (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("SF_err_pt") , e12Tvhm1_C_Hist_SF_pT,  e12Tvhm1_C_Hist_SF_pT_errUp,  e12Tvhm1_C_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_D_", (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("SF_err_eta"), e12Tvhm1_D_Hist_SF_eta, e12Tvhm1_D_Hist_SF_eta_errUp, e12Tvhm1_D_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_D_", (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("SF_err_pt") , e12Tvhm1_D_Hist_SF_pT,  e12Tvhm1_D_Hist_SF_pT_errUp,  e12Tvhm1_D_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_E_", (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("SF_err_eta"), e12Tvhm1_E_Hist_SF_eta, e12Tvhm1_E_Hist_SF_eta_errUp, e12Tvhm1_E_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_E_", (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("SF_err_pt") , e12Tvhm1_E_Hist_SF_pT,  e12Tvhm1_E_Hist_SF_pT_errUp,  e12Tvhm1_E_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_G_", (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("SF_err_eta"), e12Tvhm1_G_Hist_SF_eta, e12Tvhm1_G_Hist_SF_eta_errUp, e12Tvhm1_G_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_G_", (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("SF_err_pt") , e12Tvhm1_G_Hist_SF_pT,  e12Tvhm1_G_Hist_SF_pT_errUp,  e12Tvhm1_G_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_H_", (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("SF_err_eta"), e12Tvhm1_H_Hist_SF_eta, e12Tvhm1_H_Hist_SF_eta_errUp, e12Tvhm1_H_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_H_", (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("SF_err_pt") , e12Tvhm1_H_Hist_SF_pT,  e12Tvhm1_H_Hist_SF_pT_errUp,  e12Tvhm1_H_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_I_", (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("SF_err_eta"), e12Tvhm1_I_Hist_SF_eta, e12Tvhm1_I_Hist_SF_eta_errUp, e12Tvhm1_I_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_I_", (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("SF_err_pt") , e12Tvhm1_I_Hist_SF_pT,  e12Tvhm1_I_Hist_SF_pT_errUp,  e12Tvhm1_I_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_J_", (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("SF_err_eta"), e12Tvhm1_J_Hist_SF_eta, e12Tvhm1_J_Hist_SF_eta_errUp, e12Tvhm1_J_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_J_", (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("SF_err_pt") , e12Tvhm1_J_Hist_SF_pT,  e12Tvhm1_J_Hist_SF_pT_errUp,  e12Tvhm1_J_Hist_SF_pT_errDown);
  
  readHisto("e12Tvhm1_L_", (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("SF_eta"), (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("SF_err_eta"), e12Tvhm1_L_Hist_SF_eta, e12Tvhm1_L_Hist_SF_eta_errUp, e12Tvhm1_L_Hist_SF_eta_errDown);
  readHisto("e12Tvhm1_L_", (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("SF_pt") , (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("SF_err_pt") , e12Tvhm1_L_Hist_SF_pT,  e12Tvhm1_L_Hist_SF_pT_errUp,  e12Tvhm1_L_Hist_SF_pT_errDown);
  
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("SF_eta_Barrel"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Barrel"), mu8_AB_Hist_SF_eta_Barrel, mu8_AB_Hist_SF_eta_Barrel_errUp, mu8_AB_Hist_SF_eta_Barrel_errDown);
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("SF_pt_Barrel"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Barrel"),  mu8_AB_Hist_SF_pT_Barrel,  mu8_AB_Hist_SF_pT_Barrel_errUp,  mu8_AB_Hist_SF_pT_Barrel_errDown);    
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("SF_eta_Endcap"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Endcap"), mu8_AB_Hist_SF_eta_Endcap, mu8_AB_Hist_SF_eta_Endcap_errUp, mu8_AB_Hist_SF_eta_Endcap_errDown);
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("SF_pt_Endcap"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Endcap"),  mu8_AB_Hist_SF_pT_Endcap,  mu8_AB_Hist_SF_pT_Endcap_errUp,  mu8_AB_Hist_SF_pT_Endcap_errDown);
  
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("SF_eta_Barrel"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Barrel"), mu8_C_Hist_SF_eta_Barrel, mu8_C_Hist_SF_eta_Barrel_errUp, mu8_C_Hist_SF_eta_Barrel_errDown);
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("SF_pt_Barrel"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Barrel"),  mu8_C_Hist_SF_pT_Barrel,  mu8_C_Hist_SF_pT_Barrel_errUp,  mu8_C_Hist_SF_pT_Barrel_errDown);    
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("SF_eta_Endcap"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Endcap"), mu8_C_Hist_SF_eta_Endcap, mu8_C_Hist_SF_eta_Endcap_errUp, mu8_C_Hist_SF_eta_Endcap_errDown);
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("SF_pt_Endcap"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Endcap"),  mu8_C_Hist_SF_pT_Endcap,  mu8_C_Hist_SF_pT_Endcap_errUp,  mu8_C_Hist_SF_pT_Endcap_errDown);
  
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("SF_eta_Barrel"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Barrel"), mu8_D_Hist_SF_eta_Barrel, mu8_D_Hist_SF_eta_Barrel_errUp, mu8_D_Hist_SF_eta_Barrel_errDown);
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("SF_pt_Barrel"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Barrel"),  mu8_D_Hist_SF_pT_Barrel,  mu8_D_Hist_SF_pT_Barrel_errUp,  mu8_D_Hist_SF_pT_Barrel_errDown);    
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("SF_eta_Endcap"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Endcap"), mu8_D_Hist_SF_eta_Endcap, mu8_D_Hist_SF_eta_Endcap_errUp, mu8_D_Hist_SF_eta_Endcap_errDown);
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("SF_pt_Endcap"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Endcap"),  mu8_D_Hist_SF_pT_Endcap,  mu8_D_Hist_SF_pT_Endcap_errUp,  mu8_D_Hist_SF_pT_Endcap_errDown);
  
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("SF_eta_Barrel"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Barrel"), mu8_E_Hist_SF_eta_Barrel, mu8_E_Hist_SF_eta_Barrel_errUp, mu8_E_Hist_SF_eta_Barrel_errDown);
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("SF_pt_Barrel"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Barrel"),  mu8_E_Hist_SF_pT_Barrel,  mu8_E_Hist_SF_pT_Barrel_errUp,  mu8_E_Hist_SF_pT_Barrel_errDown);    
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("SF_eta_Endcap"), (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_eta_Endcap"), mu8_E_Hist_SF_eta_Endcap, mu8_E_Hist_SF_eta_Endcap_errUp, mu8_E_Hist_SF_eta_Endcap_errDown);
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("SF_pt_Endcap"),  (TGraphAsymmErrors*)mu8_AB_File->Get("SF_err_pt_Endcap"),  mu8_E_Hist_SF_pT_Endcap,  mu8_E_Hist_SF_pT_Endcap_errUp,  mu8_E_Hist_SF_pT_Endcap_errDown);
  
  // data efficiencies
  
  readHisto("e12Tvhl1_AB_", (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("Eff_Data_eta"), e12Tvhl1_AB_Hist_Eff_Data_eta, e12Tvhl1_AB_Hist_Eff_Data_eta_errUp, e12Tvhl1_AB_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhl1_AB_", (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("Eff_Data_pt") , e12Tvhl1_AB_Hist_Eff_Data_pT,  e12Tvhl1_AB_Hist_Eff_Data_pT_errUp,  e12Tvhl1_AB_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhl1_C_", (TGraphAsymmErrors*)e12Tvhl1_C_File->Get("Eff_Data_eta"), e12Tvhl1_C_Hist_Eff_Data_eta, e12Tvhl1_C_Hist_Eff_Data_eta_errUp, e12Tvhl1_C_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhl1_C_", (TGraphAsymmErrors*)e12Tvhl1_C_File->Get("Eff_Data_pt") , e12Tvhl1_C_Hist_Eff_Data_pT,  e12Tvhl1_C_Hist_Eff_Data_pT_errUp,  e12Tvhl1_C_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhl1_D_", (TGraphAsymmErrors*)e12Tvhl1_D_File->Get("Eff_Data_eta"), e12Tvhl1_D_Hist_Eff_Data_eta, e12Tvhl1_D_Hist_Eff_Data_eta_errUp, e12Tvhl1_D_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhl1_D_", (TGraphAsymmErrors*)e12Tvhl1_D_File->Get("Eff_Data_pt") , e12Tvhl1_D_Hist_Eff_Data_pT,  e12Tvhl1_D_Hist_Eff_Data_pT_errUp,  e12Tvhl1_D_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhl1_E_", (TGraphAsymmErrors*)e12Tvhl1_E_File->Get("Eff_Data_eta"), e12Tvhl1_E_Hist_Eff_Data_eta, e12Tvhl1_E_Hist_Eff_Data_eta_errUp, e12Tvhl1_E_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhl1_E_", (TGraphAsymmErrors*)e12Tvhl1_E_File->Get("Eff_Data_pt") , e12Tvhl1_E_Hist_Eff_Data_pT,  e12Tvhl1_E_Hist_Eff_Data_pT_errUp,  e12Tvhl1_E_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_A_", (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("Eff_Data_eta"), e12Tvhm1_A_Hist_Eff_Data_eta, e12Tvhm1_A_Hist_Eff_Data_eta_errUp, e12Tvhm1_A_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_A_", (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("Eff_Data_pt") , e12Tvhm1_A_Hist_Eff_Data_pT,  e12Tvhm1_A_Hist_Eff_Data_pT_errUp,  e12Tvhm1_A_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_B_", (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("Eff_Data_eta"), e12Tvhm1_B_Hist_Eff_Data_eta, e12Tvhm1_B_Hist_Eff_Data_eta_errUp, e12Tvhm1_B_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_B_", (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("Eff_Data_pt") , e12Tvhm1_B_Hist_Eff_Data_pT,  e12Tvhm1_B_Hist_Eff_Data_pT_errUp,  e12Tvhm1_B_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_C_", (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("Eff_Data_eta"), e12Tvhm1_C_Hist_Eff_Data_eta, e12Tvhm1_C_Hist_Eff_Data_eta_errUp, e12Tvhm1_C_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_C_", (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("Eff_Data_pt") , e12Tvhm1_C_Hist_Eff_Data_pT,  e12Tvhm1_C_Hist_Eff_Data_pT_errUp,  e12Tvhm1_C_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_D_", (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("Eff_Data_eta"), e12Tvhm1_D_Hist_Eff_Data_eta, e12Tvhm1_D_Hist_Eff_Data_eta_errUp, e12Tvhm1_D_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_D_", (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("Eff_Data_pt") , e12Tvhm1_D_Hist_Eff_Data_pT,  e12Tvhm1_D_Hist_Eff_Data_pT_errUp,  e12Tvhm1_D_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_E_", (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("Eff_Data_eta"), e12Tvhm1_E_Hist_Eff_Data_eta, e12Tvhm1_E_Hist_Eff_Data_eta_errUp, e12Tvhm1_E_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_E_", (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("Eff_Data_pt") , e12Tvhm1_E_Hist_Eff_Data_pT,  e12Tvhm1_E_Hist_Eff_Data_pT_errUp,  e12Tvhm1_E_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_G_", (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("Eff_Data_eta"), e12Tvhm1_G_Hist_Eff_Data_eta, e12Tvhm1_G_Hist_Eff_Data_eta_errUp, e12Tvhm1_G_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_G_", (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("Eff_Data_pt") , e12Tvhm1_G_Hist_Eff_Data_pT,  e12Tvhm1_G_Hist_Eff_Data_pT_errUp,  e12Tvhm1_G_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_H_", (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("Eff_Data_eta"), e12Tvhm1_H_Hist_Eff_Data_eta, e12Tvhm1_H_Hist_Eff_Data_eta_errUp, e12Tvhm1_H_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_H_", (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("Eff_Data_pt") , e12Tvhm1_H_Hist_Eff_Data_pT,  e12Tvhm1_H_Hist_Eff_Data_pT_errUp,  e12Tvhm1_H_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_I_", (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("Eff_Data_eta"), e12Tvhm1_I_Hist_Eff_Data_eta, e12Tvhm1_I_Hist_Eff_Data_eta_errUp, e12Tvhm1_I_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_I_", (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("Eff_Data_pt") , e12Tvhm1_I_Hist_Eff_Data_pT,  e12Tvhm1_I_Hist_Eff_Data_pT_errUp,  e12Tvhm1_I_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_J_", (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("Eff_Data_eta"), e12Tvhm1_J_Hist_Eff_Data_eta, e12Tvhm1_J_Hist_Eff_Data_eta_errUp, e12Tvhm1_J_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_J_", (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("Eff_Data_pt") , e12Tvhm1_J_Hist_Eff_Data_pT,  e12Tvhm1_J_Hist_Eff_Data_pT_errUp,  e12Tvhm1_J_Hist_Eff_Data_pT_errDown);
  
  readHisto("e12Tvhm1_L_", (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("Eff_Data_eta"), e12Tvhm1_L_Hist_Eff_Data_eta, e12Tvhm1_L_Hist_Eff_Data_eta_errUp, e12Tvhm1_L_Hist_Eff_Data_eta_errDown);
  readHisto("e12Tvhm1_L_", (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("Eff_Data_pt") , e12Tvhm1_L_Hist_Eff_Data_pT,  e12Tvhm1_L_Hist_Eff_Data_pT_errUp,  e12Tvhm1_L_Hist_Eff_Data_pT_errDown);
  
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_Data_eta_Barrel"), mu8_AB_Hist_Eff_Data_eta_Barrel, mu8_AB_Hist_Eff_Data_eta_Barrel_errUp, mu8_AB_Hist_Eff_Data_eta_Barrel_errDown);
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_Data_pt_Barrel"),  mu8_AB_Hist_Eff_Data_pT_Barrel,  mu8_AB_Hist_Eff_Data_pT_Barrel_errUp,  mu8_AB_Hist_Eff_Data_pT_Barrel_errDown);    
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_Data_eta_Endcap"), mu8_AB_Hist_Eff_Data_eta_Endcap, mu8_AB_Hist_Eff_Data_eta_Endcap_errUp, mu8_AB_Hist_Eff_Data_eta_Endcap_errDown);
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_Data_pt_Endcap"),  mu8_AB_Hist_Eff_Data_pT_Endcap,  mu8_AB_Hist_Eff_Data_pT_Endcap_errUp,  mu8_AB_Hist_Eff_Data_pT_Endcap_errDown);
  
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_Data_eta_Barrel"), mu8_C_Hist_Eff_Data_eta_Barrel, mu8_C_Hist_Eff_Data_eta_Barrel_errUp, mu8_C_Hist_Eff_Data_eta_Barrel_errDown);
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_Data_pt_Barrel"),  mu8_C_Hist_Eff_Data_pT_Barrel,  mu8_C_Hist_Eff_Data_pT_Barrel_errUp,  mu8_C_Hist_Eff_Data_pT_Barrel_errDown);    
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_Data_eta_Endcap"), mu8_C_Hist_Eff_Data_eta_Endcap, mu8_C_Hist_Eff_Data_eta_Endcap_errUp, mu8_C_Hist_Eff_Data_eta_Endcap_errDown);
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_Data_pt_Endcap"),  mu8_C_Hist_Eff_Data_pT_Endcap,  mu8_C_Hist_Eff_Data_pT_Endcap_errUp,  mu8_C_Hist_Eff_Data_pT_Endcap_errDown);
  
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_Data_eta_Barrel"), mu8_D_Hist_Eff_Data_eta_Barrel, mu8_D_Hist_Eff_Data_eta_Barrel_errUp, mu8_D_Hist_Eff_Data_eta_Barrel_errDown);
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_Data_pt_Barrel"),  mu8_D_Hist_Eff_Data_pT_Barrel,  mu8_D_Hist_Eff_Data_pT_Barrel_errUp,  mu8_D_Hist_Eff_Data_pT_Barrel_errDown);    
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_Data_eta_Endcap"), mu8_D_Hist_Eff_Data_eta_Endcap, mu8_D_Hist_Eff_Data_eta_Endcap_errUp, mu8_D_Hist_Eff_Data_eta_Endcap_errDown);
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_Data_pt_Endcap"),  mu8_D_Hist_Eff_Data_pT_Endcap,  mu8_D_Hist_Eff_Data_pT_Endcap_errUp,  mu8_D_Hist_Eff_Data_pT_Endcap_errDown);
  
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_Data_eta_Barrel"), mu8_E_Hist_Eff_Data_eta_Barrel, mu8_E_Hist_Eff_Data_eta_Barrel_errUp, mu8_E_Hist_Eff_Data_eta_Barrel_errDown);
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_Data_pt_Barrel"),  mu8_E_Hist_Eff_Data_pT_Barrel,  mu8_E_Hist_Eff_Data_pT_Barrel_errUp,  mu8_E_Hist_Eff_Data_pT_Barrel_errDown);    
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_Data_eta_Endcap"), mu8_E_Hist_Eff_Data_eta_Endcap, mu8_E_Hist_Eff_Data_eta_Endcap_errUp, mu8_E_Hist_Eff_Data_eta_Endcap_errDown);
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_Data_pt_Endcap"),  mu8_E_Hist_Eff_Data_pT_Endcap,  mu8_E_Hist_Eff_Data_pT_Endcap_errUp,  mu8_E_Hist_Eff_Data_pT_Endcap_errDown);
  
  // MC efficiencies
  
  readHisto("e12Tvhl1_AB_", (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("Eff_MC_eta"), e12Tvhl1_AB_Hist_Eff_MC_eta, e12Tvhl1_AB_Hist_Eff_MC_eta_errUp, e12Tvhl1_AB_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhl1_AB_", (TGraphAsymmErrors*)e12Tvhl1_AB_File->Get("Eff_MC_pt") , e12Tvhl1_AB_Hist_Eff_MC_pT,  e12Tvhl1_AB_Hist_Eff_MC_pT_errUp,  e12Tvhl1_AB_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhl1_C_", (TGraphAsymmErrors*)e12Tvhl1_C_File->Get("Eff_MC_eta"), e12Tvhl1_C_Hist_Eff_MC_eta, e12Tvhl1_C_Hist_Eff_MC_eta_errUp, e12Tvhl1_C_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhl1_C_", (TGraphAsymmErrors*)e12Tvhl1_C_File->Get("Eff_MC_pt") , e12Tvhl1_C_Hist_Eff_MC_pT,  e12Tvhl1_C_Hist_Eff_MC_pT_errUp,  e12Tvhl1_C_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhl1_D_", (TGraphAsymmErrors*)e12Tvhl1_D_File->Get("Eff_MC_eta"), e12Tvhl1_D_Hist_Eff_MC_eta, e12Tvhl1_D_Hist_Eff_MC_eta_errUp, e12Tvhl1_D_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhl1_D_", (TGraphAsymmErrors*)e12Tvhl1_D_File->Get("Eff_MC_pt") , e12Tvhl1_D_Hist_Eff_MC_pT,  e12Tvhl1_D_Hist_Eff_MC_pT_errUp,  e12Tvhl1_D_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhl1_E_", (TGraphAsymmErrors*)e12Tvhl1_E_File->Get("Eff_MC_eta"), e12Tvhl1_E_Hist_Eff_MC_eta, e12Tvhl1_E_Hist_Eff_MC_eta_errUp, e12Tvhl1_E_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhl1_E_", (TGraphAsymmErrors*)e12Tvhl1_E_File->Get("Eff_MC_pt") , e12Tvhl1_E_Hist_Eff_MC_pT,  e12Tvhl1_E_Hist_Eff_MC_pT_errUp,  e12Tvhl1_E_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_A_", (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("Eff_MC_eta"), e12Tvhm1_A_Hist_Eff_MC_eta, e12Tvhm1_A_Hist_Eff_MC_eta_errUp, e12Tvhm1_A_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_A_", (TGraphAsymmErrors*)e12Tvhm1_A_File->Get("Eff_MC_pt") , e12Tvhm1_A_Hist_Eff_MC_pT,  e12Tvhm1_A_Hist_Eff_MC_pT_errUp,  e12Tvhm1_A_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_B_", (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("Eff_MC_eta"), e12Tvhm1_B_Hist_Eff_MC_eta, e12Tvhm1_B_Hist_Eff_MC_eta_errUp, e12Tvhm1_B_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_B_", (TGraphAsymmErrors*)e12Tvhm1_B_File->Get("Eff_MC_pt") , e12Tvhm1_B_Hist_Eff_MC_pT,  e12Tvhm1_B_Hist_Eff_MC_pT_errUp,  e12Tvhm1_B_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_C_", (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("Eff_MC_eta"), e12Tvhm1_C_Hist_Eff_MC_eta, e12Tvhm1_C_Hist_Eff_MC_eta_errUp, e12Tvhm1_C_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_C_", (TGraphAsymmErrors*)e12Tvhm1_C_File->Get("Eff_MC_pt") , e12Tvhm1_C_Hist_Eff_MC_pT,  e12Tvhm1_C_Hist_Eff_MC_pT_errUp,  e12Tvhm1_C_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_D_", (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("Eff_MC_eta"), e12Tvhm1_D_Hist_Eff_MC_eta, e12Tvhm1_D_Hist_Eff_MC_eta_errUp, e12Tvhm1_D_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_D_", (TGraphAsymmErrors*)e12Tvhm1_D_File->Get("Eff_MC_pt") , e12Tvhm1_D_Hist_Eff_MC_pT,  e12Tvhm1_D_Hist_Eff_MC_pT_errUp,  e12Tvhm1_D_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_E_", (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("Eff_MC_eta"), e12Tvhm1_E_Hist_Eff_MC_eta, e12Tvhm1_E_Hist_Eff_MC_eta_errUp, e12Tvhm1_E_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_E_", (TGraphAsymmErrors*)e12Tvhm1_E_File->Get("Eff_MC_pt") , e12Tvhm1_E_Hist_Eff_MC_pT,  e12Tvhm1_E_Hist_Eff_MC_pT_errUp,  e12Tvhm1_E_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_G_", (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("Eff_MC_eta"), e12Tvhm1_G_Hist_Eff_MC_eta, e12Tvhm1_G_Hist_Eff_MC_eta_errUp, e12Tvhm1_G_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_G_", (TGraphAsymmErrors*)e12Tvhm1_G_File->Get("Eff_MC_pt") , e12Tvhm1_G_Hist_Eff_MC_pT,  e12Tvhm1_G_Hist_Eff_MC_pT_errUp,  e12Tvhm1_G_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_H_", (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("Eff_MC_eta"), e12Tvhm1_H_Hist_Eff_MC_eta, e12Tvhm1_H_Hist_Eff_MC_eta_errUp, e12Tvhm1_H_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_H_", (TGraphAsymmErrors*)e12Tvhm1_H_File->Get("Eff_MC_pt") , e12Tvhm1_H_Hist_Eff_MC_pT,  e12Tvhm1_H_Hist_Eff_MC_pT_errUp,  e12Tvhm1_H_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_I_", (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("Eff_MC_eta"), e12Tvhm1_I_Hist_Eff_MC_eta, e12Tvhm1_I_Hist_Eff_MC_eta_errUp, e12Tvhm1_I_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_I_", (TGraphAsymmErrors*)e12Tvhm1_I_File->Get("Eff_MC_pt") , e12Tvhm1_I_Hist_Eff_MC_pT,  e12Tvhm1_I_Hist_Eff_MC_pT_errUp,  e12Tvhm1_I_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_J_", (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("Eff_MC_eta"), e12Tvhm1_J_Hist_Eff_MC_eta, e12Tvhm1_J_Hist_Eff_MC_eta_errUp, e12Tvhm1_J_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_J_", (TGraphAsymmErrors*)e12Tvhm1_J_File->Get("Eff_MC_pt") , e12Tvhm1_J_Hist_Eff_MC_pT,  e12Tvhm1_J_Hist_Eff_MC_pT_errUp,  e12Tvhm1_J_Hist_Eff_MC_pT_errDown);
  
  readHisto("e12Tvhm1_L_", (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("Eff_MC_eta"), e12Tvhm1_L_Hist_Eff_MC_eta, e12Tvhm1_L_Hist_Eff_MC_eta_errUp, e12Tvhm1_L_Hist_Eff_MC_eta_errDown);
  readHisto("e12Tvhm1_L_", (TGraphAsymmErrors*)e12Tvhm1_L_File->Get("Eff_MC_pt") , e12Tvhm1_L_Hist_Eff_MC_pT,  e12Tvhm1_L_Hist_Eff_MC_pT_errUp,  e12Tvhm1_L_Hist_Eff_MC_pT_errDown);
  
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_MC_eta_Barrel"), mu8_AB_Hist_Eff_MC_eta_Barrel, mu8_AB_Hist_Eff_MC_eta_Barrel_errUp, mu8_AB_Hist_Eff_MC_eta_Barrel_errDown);
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_MC_pt_Barrel"),  mu8_AB_Hist_Eff_MC_pT_Barrel,  mu8_AB_Hist_Eff_MC_pT_Barrel_errUp,  mu8_AB_Hist_Eff_MC_pT_Barrel_errDown);    
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_MC_eta_Endcap"), mu8_AB_Hist_Eff_MC_eta_Endcap, mu8_AB_Hist_Eff_MC_eta_Endcap_errUp, mu8_AB_Hist_Eff_MC_eta_Endcap_errDown);
  readHisto("mu8_AB_", (TGraphAsymmErrors*)mu8_AB_File->Get("Eff_MC_pt_Endcap"),  mu8_AB_Hist_Eff_MC_pT_Endcap,  mu8_AB_Hist_Eff_MC_pT_Endcap_errUp,  mu8_AB_Hist_Eff_MC_pT_Endcap_errDown);
  
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_MC_eta_Barrel"), mu8_C_Hist_Eff_MC_eta_Barrel, mu8_C_Hist_Eff_MC_eta_Barrel_errUp, mu8_C_Hist_Eff_MC_eta_Barrel_errDown);
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_MC_pt_Barrel"),  mu8_C_Hist_Eff_MC_pT_Barrel,  mu8_C_Hist_Eff_MC_pT_Barrel_errUp,  mu8_C_Hist_Eff_MC_pT_Barrel_errDown);    
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_MC_eta_Endcap"), mu8_C_Hist_Eff_MC_eta_Endcap, mu8_C_Hist_Eff_MC_eta_Endcap_errUp, mu8_C_Hist_Eff_MC_eta_Endcap_errDown);
  readHisto("mu8_C_", (TGraphAsymmErrors*)mu8_C_File->Get("Eff_MC_pt_Endcap"),  mu8_C_Hist_Eff_MC_pT_Endcap,  mu8_C_Hist_Eff_MC_pT_Endcap_errUp,  mu8_C_Hist_Eff_MC_pT_Endcap_errDown);
  
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_MC_eta_Barrel"), mu8_D_Hist_Eff_MC_eta_Barrel, mu8_D_Hist_Eff_MC_eta_Barrel_errUp, mu8_D_Hist_Eff_MC_eta_Barrel_errDown);
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_MC_pt_Barrel"),  mu8_D_Hist_Eff_MC_pT_Barrel,  mu8_D_Hist_Eff_MC_pT_Barrel_errUp,  mu8_D_Hist_Eff_MC_pT_Barrel_errDown);    
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_MC_eta_Endcap"), mu8_D_Hist_Eff_MC_eta_Endcap, mu8_D_Hist_Eff_MC_eta_Endcap_errUp, mu8_D_Hist_Eff_MC_eta_Endcap_errDown);
  readHisto("mu8_D_", (TGraphAsymmErrors*)mu8_D_File->Get("Eff_MC_pt_Endcap"),  mu8_D_Hist_Eff_MC_pT_Endcap,  mu8_D_Hist_Eff_MC_pT_Endcap_errUp,  mu8_D_Hist_Eff_MC_pT_Endcap_errDown);
  
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_MC_eta_Barrel"), mu8_E_Hist_Eff_MC_eta_Barrel, mu8_E_Hist_Eff_MC_eta_Barrel_errUp, mu8_E_Hist_Eff_MC_eta_Barrel_errDown);
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_MC_pt_Barrel"),  mu8_E_Hist_Eff_MC_pT_Barrel,  mu8_E_Hist_Eff_MC_pT_Barrel_errUp,  mu8_E_Hist_Eff_MC_pT_Barrel_errDown);    
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_MC_eta_Endcap"), mu8_E_Hist_Eff_MC_eta_Endcap, mu8_E_Hist_Eff_MC_eta_Endcap_errUp, mu8_E_Hist_Eff_MC_eta_Endcap_errDown);
  readHisto("mu8_E_", (TGraphAsymmErrors*)mu8_E_File->Get("Eff_MC_pt_Endcap"),  mu8_E_Hist_Eff_MC_pT_Endcap,  mu8_E_Hist_Eff_MC_pT_Endcap_errUp,  mu8_E_Hist_Eff_MC_pT_Endcap_errDown);

  // scale factor histograms
  
  normaliseHisto(e12Tvhl1_AB_Hist_SF_eta);
  normaliseHisto(e12Tvhl1_AB_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhl1_AB_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhl1_C_Hist_SF_eta);
  normaliseHisto(e12Tvhl1_C_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhl1_C_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhl1_D_Hist_SF_eta);
  normaliseHisto(e12Tvhl1_D_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhl1_D_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhl1_E_Hist_SF_eta);
  normaliseHisto(e12Tvhl1_E_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhl1_E_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_A_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_A_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_A_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_B_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_B_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_B_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_C_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_C_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_C_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_D_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_D_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_D_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_E_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_E_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_E_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_G_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_G_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_G_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_H_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_H_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_H_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_I_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_I_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_I_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_J_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_J_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_J_Hist_SF_eta_errDown);
  
  normaliseHisto(e12Tvhm1_L_Hist_SF_eta);
  normaliseHisto(e12Tvhm1_L_Hist_SF_eta_errUp);
  normaliseHisto(e12Tvhm1_L_Hist_SF_eta_errDown);
  
  normaliseHisto(mu8_AB_Hist_SF_eta_Barrel);
  normaliseHisto(mu8_AB_Hist_SF_eta_Barrel_errUp);
  normaliseHisto(mu8_AB_Hist_SF_eta_Barrel_errDown);
  normaliseHisto(mu8_AB_Hist_SF_eta_Endcap);
  normaliseHisto(mu8_AB_Hist_SF_eta_Endcap_errUp);    
  normaliseHisto(mu8_AB_Hist_SF_eta_Endcap_errDown);
  
  normaliseHisto(mu8_C_Hist_SF_eta_Barrel);
  normaliseHisto(mu8_C_Hist_SF_eta_Barrel_errUp);
  normaliseHisto(mu8_C_Hist_SF_eta_Barrel_errDown);
  normaliseHisto(mu8_C_Hist_SF_eta_Endcap);
  normaliseHisto(mu8_C_Hist_SF_eta_Endcap_errUp);    
  normaliseHisto(mu8_C_Hist_SF_eta_Endcap_errDown);
  
  normaliseHisto(mu8_D_Hist_SF_eta_Barrel);
  normaliseHisto(mu8_D_Hist_SF_eta_Barrel_errUp);
  normaliseHisto(mu8_D_Hist_SF_eta_Barrel_errDown);
  normaliseHisto(mu8_D_Hist_SF_eta_Endcap);
  normaliseHisto(mu8_D_Hist_SF_eta_Endcap_errUp);    
  normaliseHisto(mu8_D_Hist_SF_eta_Endcap_errDown);
  
  normaliseHisto(mu8_E_Hist_SF_eta_Barrel);
  normaliseHisto(mu8_E_Hist_SF_eta_Barrel_errUp);
  normaliseHisto(mu8_E_Hist_SF_eta_Barrel_errDown);
  normaliseHisto(mu8_E_Hist_SF_eta_Endcap);
  normaliseHisto(mu8_E_Hist_SF_eta_Endcap_errUp);    
  normaliseHisto(mu8_E_Hist_SF_eta_Endcap_errDown);
  
  // data efficiency histograms
  
  normaliseHisto(e12Tvhl1_AB_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhl1_AB_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhl1_AB_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhl1_C_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhl1_C_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhl1_C_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhl1_D_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhl1_D_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhl1_D_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhl1_E_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhl1_E_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhl1_E_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_A_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_A_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_A_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_B_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_B_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_B_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_C_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_C_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_C_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_D_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_D_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_D_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_E_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_E_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_E_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_G_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_G_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_G_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_H_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_H_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_H_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_I_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_I_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_I_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_J_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_J_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_J_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(e12Tvhm1_L_Hist_Eff_Data_eta);
  normaliseHisto(e12Tvhm1_L_Hist_Eff_Data_eta_errUp);
  normaliseHisto(e12Tvhm1_L_Hist_Eff_Data_eta_errDown);
  
  normaliseHisto(mu8_AB_Hist_Eff_Data_eta_Barrel);    
  normaliseHisto(mu8_AB_Hist_Eff_Data_eta_Barrel_errUp);
  normaliseHisto(mu8_AB_Hist_Eff_Data_eta_Barrel_errDown);
  normaliseHisto(mu8_AB_Hist_Eff_Data_eta_Endcap);
  normaliseHisto(mu8_AB_Hist_Eff_Data_eta_Endcap_errUp);    
  normaliseHisto(mu8_AB_Hist_Eff_Data_eta_Endcap_errDown);
  
  normaliseHisto(mu8_C_Hist_Eff_Data_eta_Barrel);    
  normaliseHisto(mu8_C_Hist_Eff_Data_eta_Barrel_errUp);
  normaliseHisto(mu8_C_Hist_Eff_Data_eta_Barrel_errDown);
  normaliseHisto(mu8_C_Hist_Eff_Data_eta_Endcap);
  normaliseHisto(mu8_C_Hist_Eff_Data_eta_Endcap_errUp);    
  normaliseHisto(mu8_C_Hist_Eff_Data_eta_Endcap_errDown);
  
  normaliseHisto(mu8_D_Hist_Eff_Data_eta_Barrel);    
  normaliseHisto(mu8_D_Hist_Eff_Data_eta_Barrel_errUp);
  normaliseHisto(mu8_D_Hist_Eff_Data_eta_Barrel_errDown);
  normaliseHisto(mu8_D_Hist_Eff_Data_eta_Endcap);
  normaliseHisto(mu8_D_Hist_Eff_Data_eta_Endcap_errUp);    
  normaliseHisto(mu8_D_Hist_Eff_Data_eta_Endcap_errDown);
  
  normaliseHisto(mu8_E_Hist_Eff_Data_eta_Barrel);    
  normaliseHisto(mu8_E_Hist_Eff_Data_eta_Barrel_errUp);
  normaliseHisto(mu8_E_Hist_Eff_Data_eta_Barrel_errDown);
  normaliseHisto(mu8_E_Hist_Eff_Data_eta_Endcap);
  normaliseHisto(mu8_E_Hist_Eff_Data_eta_Endcap_errUp);    
  normaliseHisto(mu8_E_Hist_Eff_Data_eta_Endcap_errDown);
  
  // MC efficiency histograms
  
  normaliseHisto(e12Tvhl1_AB_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhl1_AB_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhl1_AB_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhl1_C_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhl1_C_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhl1_C_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhl1_D_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhl1_D_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhl1_D_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhl1_E_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhl1_E_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhl1_E_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_A_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_A_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_A_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_B_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_B_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_B_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_C_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_C_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_C_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_D_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_D_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_D_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_E_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_E_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_E_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_G_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_G_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_G_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_H_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_H_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_H_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_I_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_I_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_I_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_J_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_J_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_J_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(e12Tvhm1_L_Hist_Eff_MC_eta);
  normaliseHisto(e12Tvhm1_L_Hist_Eff_MC_eta_errUp);
  normaliseHisto(e12Tvhm1_L_Hist_Eff_MC_eta_errDown);
  
  normaliseHisto(mu8_AB_Hist_Eff_MC_eta_Barrel);    
  normaliseHisto(mu8_AB_Hist_Eff_MC_eta_Barrel_errUp);
  normaliseHisto(mu8_AB_Hist_Eff_MC_eta_Barrel_errDown);
  normaliseHisto(mu8_AB_Hist_Eff_MC_eta_Endcap);
  normaliseHisto(mu8_AB_Hist_Eff_MC_eta_Endcap_errUp);    
  normaliseHisto(mu8_AB_Hist_Eff_MC_eta_Endcap_errDown);
  
  normaliseHisto(mu8_C_Hist_Eff_MC_eta_Barrel);    
  normaliseHisto(mu8_C_Hist_Eff_MC_eta_Barrel_errUp);
  normaliseHisto(mu8_C_Hist_Eff_MC_eta_Barrel_errDown);
  normaliseHisto(mu8_C_Hist_Eff_MC_eta_Endcap);
  normaliseHisto(mu8_C_Hist_Eff_MC_eta_Endcap_errUp);    
  normaliseHisto(mu8_C_Hist_Eff_MC_eta_Endcap_errDown);
  
  normaliseHisto(mu8_D_Hist_Eff_MC_eta_Barrel);    
  normaliseHisto(mu8_D_Hist_Eff_MC_eta_Barrel_errUp);
  normaliseHisto(mu8_D_Hist_Eff_MC_eta_Barrel_errDown);
  normaliseHisto(mu8_D_Hist_Eff_MC_eta_Endcap);
  normaliseHisto(mu8_D_Hist_Eff_MC_eta_Endcap_errUp);    
  normaliseHisto(mu8_D_Hist_Eff_MC_eta_Endcap_errDown);
  
  normaliseHisto(mu8_E_Hist_Eff_MC_eta_Barrel);    
  normaliseHisto(mu8_E_Hist_Eff_MC_eta_Barrel_errUp);
  normaliseHisto(mu8_E_Hist_Eff_MC_eta_Barrel_errDown);
  normaliseHisto(mu8_E_Hist_Eff_MC_eta_Endcap);
  normaliseHisto(mu8_E_Hist_Eff_MC_eta_Endcap_errUp);    
  normaliseHisto(mu8_E_Hist_Eff_MC_eta_Endcap_errDown);
  
  //--AM- Closing and deleting file from memory
  e12Tvhl1_AB_File->Close();
  e12Tvhl1_C_File->Close();
  e12Tvhl1_D_File->Close();
  e12Tvhl1_E_File->Close();

  e12Tvhm1_A_File->Close();
  e12Tvhm1_B_File->Close();
  e12Tvhm1_C_File->Close();
  e12Tvhm1_D_File->Close();
  e12Tvhm1_E_File->Close();
  e12Tvhm1_G_File->Close();
  e12Tvhm1_H_File->Close();
  e12Tvhm1_I_File->Close();
  e12Tvhm1_J_File->Close();
  e12Tvhm1_L_File->Close();

  mu8_AB_File->Close();
  mu8_C_File->Close();
  mu8_D_File->Close();
  mu8_E_File->Close();
  
  delete e12Tvhl1_AB_File;
  delete e12Tvhl1_C_File;
  delete e12Tvhl1_D_File;
  delete e12Tvhl1_E_File;

  delete e12Tvhm1_A_File;
  delete e12Tvhm1_B_File;
  delete e12Tvhm1_C_File;
  delete e12Tvhm1_D_File;
  delete e12Tvhm1_E_File;
  delete e12Tvhm1_G_File;
  delete e12Tvhm1_H_File;
  delete e12Tvhm1_I_File;
  delete e12Tvhm1_J_File;
  delete e12Tvhm1_L_File;

  delete mu8_AB_File;
  delete mu8_C_File;
  delete mu8_D_File;
  delete mu8_E_File;
}


void HSG4LepLepTriggerSF::readHisto( const TString& setName, TGraphAsymmErrors* const inGraph , TGraphAsymmErrors* const inGraphSys , TH1F*& outHist , TH1F*& outHistErrUp , TH1F*& outHistErrDown )
{

    TString graphName( setName ) ;
    graphName.Append( inGraph->GetName() ) ;

    // Work out bin boundaries
    // Modified by Yang Heng <Yang.Heng@cern.ch> on 03-Jul-2012
    Int_t nPoints = inGraph->GetN() ; inGraph->Sort() ;

    Double_t binEdges[ nPoints + 1 ] ;
 
    Double_t xBin , xBinErrLow , xBinErrHigh , y ;

    for ( Int_t i = 0 ; i < nPoints ; i++ ){
        inGraph->GetPoint( i , xBin , y ) ;
        xBinErrLow = inGraph->GetErrorXlow( i ) ;
        binEdges[ i ] = xBin - xBinErrLow ;  
    }

    inGraph->GetPoint( nPoints - 1 , xBin , y ) ;
    xBinErrHigh = inGraph->GetErrorXhigh( nPoints - 1) ;
    binEdges[ nPoints ] = xBin + xBinErrHigh ; 

    // Then declare histrograms...
    outHist        = new TH1F( graphName               , graphName               , nPoints , binEdges ) ;				
    outHistErrUp   = new TH1F( (graphName + "ErrUp")   , (graphName + "ErrUp")   , nPoints , binEdges ) ;
    outHistErrDown = new TH1F( (graphName + "ErrDown") , (graphName + "ErrDown") , nPoints , binEdges ) ;

    outHist       ->Sumw2() ; outHist       ->StatOverflows( kTRUE ) ;
    outHistErrUp  ->Sumw2() ; outHistErrUp  ->StatOverflows( kTRUE ) ;
    outHistErrDown->Sumw2() ; outHistErrDown->StatOverflows( kTRUE ) ;

    // ... and fill 'em
    // Modified by Yang Heng <Yang.Heng@cern.ch> on 03-Jul-2012
    for ( Int_t i = 0 ; i < nPoints ; i++ ){
        Double_t x , y , errHigh, errLow , errHighSys, errLowSys ;
        inGraph->GetPoint( i , x , y ) ;
        errHigh = inGraph->GetErrorYhigh( i ) ;
        errLow = inGraph->GetErrorYlow( i ) ;
        errHighSys = inGraphSys->GetErrorYhigh( i ) ;
        errLowSys = inGraphSys->GetErrorYlow( i ) ;

        outHist->SetBinContent( i + 1 , y ) ;
        outHist->SetBinError( i + 1 , 0.0 ) ;

        outHistErrUp->SetBinContent( i + 1 , y ) ;
        outHistErrUp->SetBinError( i + 1 , TMath::Sqrt(TMath::Power(errHigh, 2.0) + TMath::Power(errHighSys, 2.0)) ) ;

        outHistErrDown->SetBinContent( i + 1 , y ) ;
        outHistErrDown->SetBinError( i + 1 , TMath::Sqrt(TMath::Power(errLow, 2.0) + TMath::Power(errLowSys, 2.0)) ) ;
    }

    outHist       ->SetBinContent( 0 , 0 ) ; outHist       ->SetBinError( 0 , 0 ) ;
    outHistErrUp  ->SetBinContent( 0 , 0 ) ; outHistErrUp  ->SetBinError( 0 , 0 ) ;
    outHistErrDown->SetBinContent( 0 , 0 ) ; outHistErrDown->SetBinError( 0 , 0 ) ;

    outHist       ->SetBinContent( nPoints + 1 , 0 ) ; outHist       ->SetBinError( nPoints + 1 , 0 ) ;
    outHistErrUp  ->SetBinContent( nPoints + 1 , 0 ) ; outHistErrUp  ->SetBinError( nPoints + 1 , 0 ) ;
    outHistErrDown->SetBinContent( nPoints + 1 , 0 ) ; outHistErrDown->SetBinError( nPoints + 1 , 0 ) ;

    if(doDebug == kTRUE ){
        cout << "finished read in of : " << inGraph->GetName() << endl ; 
        outHist       ->Print( "all" ) ; 
        outHistErrUp  ->Print( "all" ) ;
        outHistErrDown->Print( "all" ) ;
    }

}

void HSG4LepLepTriggerSF::readHisto( const TString& setName, TGraphAsymmErrors* const inGraph , TH1F*& outHist , TH1F*& outHistErrUp , TH1F*& outHistErrDown )
{

    // readHisto for efficiency values (only stat err, no sys)
    
    TString graphName( setName ) ;
    graphName.Append( inGraph->GetName() ) ;

    // Work out bin boundaries
    // Modified by Yang Heng <Yang.Heng@cern.ch> on 03-Jul-2012
    Int_t nPoints = inGraph->GetN() ; inGraph->Sort() ;

    Double_t binEdges[ nPoints + 1 ] ;
 
    Double_t xBin , xBinErrLow , xBinErrHigh , y ;

    for ( Int_t i = 0 ; i < nPoints ; i++ ){
        inGraph->GetPoint( i , xBin , y ) ;
        xBinErrLow = inGraph->GetErrorXlow( i ) ;
        binEdges[ i ] = xBin - xBinErrLow ;  
    }

    inGraph->GetPoint( nPoints - 1 , xBin , y ) ;
    xBinErrHigh = inGraph->GetErrorXhigh( nPoints - 1) ;
    binEdges[ nPoints ] = xBin + xBinErrHigh ; 

    // Then declare histrograms...
    outHist        = new TH1F( graphName               , graphName               , nPoints , binEdges ) ;				
    outHistErrUp   = new TH1F( (graphName + "ErrUp")   , (graphName + "ErrUp")   , nPoints , binEdges ) ;
    outHistErrDown = new TH1F( (graphName + "ErrDown") , (graphName + "ErrDown") , nPoints , binEdges ) ;

    outHist       ->Sumw2() ; outHist       ->StatOverflows( kTRUE ) ;
    outHistErrUp  ->Sumw2() ; outHistErrUp  ->StatOverflows( kTRUE ) ;
    outHistErrDown->Sumw2() ; outHistErrDown->StatOverflows( kTRUE ) ;

    // ... and fill 'em
    // Modified by Yang Heng <Yang.Heng@cern.ch> on 03-Jul-2012
    for ( Int_t i = 0 ; i < nPoints ; i++ ){
        Double_t x , y , errHigh, errLow ;
        inGraph->GetPoint( i , x , y ) ;
        errHigh = inGraph->GetErrorYhigh( i ) ;
        errLow = inGraph->GetErrorYlow( i ) ;

        outHist->SetBinContent( i + 1 , y ) ;
        outHist->SetBinError( i + 1 , 0.0 ) ;

        outHistErrUp->SetBinContent( i + 1 , y ) ;
        outHistErrUp->SetBinError( i + 1 , errHigh ) ;

        outHistErrDown->SetBinContent( i + 1 , y ) ;
        outHistErrDown->SetBinError( i + 1 , errLow ) ;
    }

    outHist       ->SetBinContent( 0 , 0 ) ; outHist       ->SetBinError( 0 , 0 ) ;
    outHistErrUp  ->SetBinContent( 0 , 0 ) ; outHistErrUp  ->SetBinError( 0 , 0 ) ;
    outHistErrDown->SetBinContent( 0 , 0 ) ; outHistErrDown->SetBinError( 0 , 0 ) ;

    outHist       ->SetBinContent( nPoints + 1 , 0 ) ; outHist       ->SetBinError( nPoints + 1 , 0 ) ;
    outHistErrUp  ->SetBinContent( nPoints + 1 , 0 ) ; outHistErrUp  ->SetBinError( nPoints + 1 , 0 ) ;
    outHistErrDown->SetBinContent( nPoints + 1 , 0 ) ; outHistErrDown->SetBinError( nPoints + 1 , 0 ) ;

    if(doDebug == kTRUE ){
        cout << "finished read in of : " << inGraph->GetName() << endl ; 
        outHist       ->Print( "all" ) ; 
        outHistErrUp  ->Print( "all" ) ;
        outHistErrDown->Print( "all" ) ;
    }

}

void HSG4LepLepTriggerSF::normaliseHisto( TH1F* hist ) 
{

    // this is assuming that the eta, phi distributions are flat
    // but it might not be the case for eta ...
    Int_t nNonZero = 0 ;
    for ( Int_t i = 1 ; i <= hist->GetNbinsX() + 1 ; i++ ){
        if(hist->GetBinContent( i ) > 0 ) nNonZero++ ;
    }
    hist->Scale( 1.0 / ( hist->Integral( 0 , hist->GetNbinsX() + 1) / ( Float_t ) nNonZero ) ) ;
}


Double_t HSG4LepLepTriggerSF::getSFElec( const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const
{

  Double_t sf = 1.0 ; 
  Double_t sfErr = 0.0 ;
  
  if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getSFElec pT " << lepVec.Pt() << " RunNumber = " << runNumber << " Trigger = " << trigger << endl ;
        
  if(trigger == "e12Tvhl1"){
  
    if(lepVec.Pt() < 12e3){
      // Scale factors not valid for pt < 12 GeV
      if(doDebug == kTRUE) cout << " HSG4LepLepTriggerSF::getSFElec pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl;
    } 
    else{
      
//      if(TMath::Abs(lepVec.Eta()) < 1.37 || (TMath::Abs(lepVec.Eta()) > 1.52 && TMath::Abs(lepVec.Eta()) <2.47)){
      if(TMath::Abs(lepVec.Eta()) < 1.37 || TMath::Abs(lepVec.Eta()) > 1.52){
        
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
     	  Int_t bin = e12Tvhl1_AB_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
     	  if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_AB_Hist_SF_pT->FindBin(100-0.0001);
     	  sf =  e12Tvhl1_AB_Hist_SF_pT->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_AB_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhl1_AB_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_AB_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhl1_AB_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }

     	  bin = e12Tvhl1_AB_Hist_SF_eta->FindBin(lepVec.Eta());
     	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
     	    bin = lepVec.Eta() > 0 ? e12Tvhl1_AB_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhl1_AB_Hist_SF_eta->FindBin(-2.47+0.001);
     	  sf *= e12Tvhl1_AB_Hist_SF_eta->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_AB_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhl1_AB_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_AB_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhl1_AB_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }
	}
	
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C
     	  Int_t bin = e12Tvhl1_C_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
     	  if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_C_Hist_SF_pT->FindBin(100-0.0001);
     	  sf =  e12Tvhl1_C_Hist_SF_pT->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_C_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhl1_C_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_C_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhl1_C_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }

     	  bin = e12Tvhl1_C_Hist_SF_eta->FindBin(lepVec.Eta());
     	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
     	    bin = lepVec.Eta() > 0 ? e12Tvhl1_C_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhl1_C_Hist_SF_eta->FindBin(-2.47+0.001);
     	  sf *= e12Tvhl1_C_Hist_SF_eta->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_C_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhl1_C_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_C_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhl1_C_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }
	}

	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
     	  Int_t bin = e12Tvhl1_D_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
     	  if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_D_Hist_SF_pT->FindBin(100-0.0001);
     	  sf =  e12Tvhl1_D_Hist_SF_pT->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_D_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhl1_D_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_D_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhl1_D_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }

     	  bin = e12Tvhl1_D_Hist_SF_eta->FindBin(lepVec.Eta());
     	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
     	    bin = lepVec.Eta() > 0 ? e12Tvhl1_D_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhl1_D_Hist_SF_eta->FindBin(-2.47+0.001);
     	  sf *= e12Tvhl1_D_Hist_SF_eta->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_D_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhl1_D_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_D_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhl1_D_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }
	}
	
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
     	  Int_t bin = e12Tvhl1_E_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
     	  if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_E_Hist_SF_pT->FindBin(100-0.0001);
     	  sf =  e12Tvhl1_E_Hist_SF_pT->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_E_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhl1_E_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_E_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhl1_E_Hist_SF_pT->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }

     	  bin = e12Tvhl1_E_Hist_SF_eta->FindBin(lepVec.Eta());
     	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
     	    bin = lepVec.Eta() > 0 ? e12Tvhl1_E_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhl1_E_Hist_SF_eta->FindBin(-2.47+0.001);
     	  sf *= e12Tvhl1_E_Hist_SF_eta->GetBinContent(bin);
     	  if(mode == 1){
     	    Double_t sfErr_tmp = e12Tvhl1_E_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhl1_E_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  } 
     	  else if(mode == -1){
     	    Double_t sfErr_tmp = e12Tvhl1_E_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhl1_E_Hist_SF_eta->GetBinContent(bin) : 0 ;
     	    sfErr += TMath::Power(sfErr_tmp, 2);
     	  }
	}	
	
  	sfErr = TMath::Sqrt(sfErr) / sf;
	
      }      
    }
  } 
  
  else if(trigger == "e12Tvhm1"){
  
    if(lepVec.Pt() < 12e3){
      // Scale factors not valid for pt < 12 GeV
      if(doDebug == kTRUE) cout << " HSG4LepLepTriggerSF::getSFElec pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl;
    } 
    else{
      
//      if(TMath::Abs(lepVec.Eta()) < 1.37 || (TMath::Abs(lepVec.Eta()) > 1.52 && TMath::Abs(lepVec.Eta()) <2.47)){
      if(TMath::Abs(lepVec.Eta()) < 1.37 || TMath::Abs(lepVec.Eta()) > 1.52){
      
	if(runNumber >= 200804 && runNumber <= 201556){ // period A
          Int_t bin = e12Tvhm1_A_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_A_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_A_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_A_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_A_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_A_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_A_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_A_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_A_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_A_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_A_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_A_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_A_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_A_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_A_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}
	
	if(runNumber >= 202660 && runNumber <= 205113){ // period B
          Int_t bin = e12Tvhm1_B_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_B_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_B_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_B_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_B_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_B_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_B_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_B_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_B_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_B_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_B_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_B_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_B_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_B_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_B_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}
	
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C
	  Int_t bin = e12Tvhm1_C_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_C_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_C_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_C_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_C_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_C_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_C_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_C_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_C_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_C_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_C_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_C_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_C_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_C_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_C_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}

	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
	  Int_t bin = e12Tvhm1_D_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_D_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_D_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_D_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_D_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_D_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_D_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_D_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_D_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_D_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_D_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_D_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_D_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_D_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_D_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}
	
	else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	  Int_t bin = e12Tvhm1_E_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_E_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_E_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_E_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_E_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_E_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_E_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_E_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_E_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_E_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_E_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_E_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_E_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_E_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_E_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}	

	else if(runNumber >= 211522 && runNumber <= 212272){ // period G	
	  Int_t bin = e12Tvhm1_G_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_G_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_G_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_G_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_G_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_G_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_G_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_G_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_G_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_G_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_G_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_G_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_G_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_G_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_G_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}	

	else if(runNumber >= 212619 && runNumber <= 213359){ // period H	
	  Int_t bin = e12Tvhm1_H_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_H_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_H_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_H_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_H_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_H_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_H_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_H_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_H_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_H_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_H_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_H_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_H_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_H_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_H_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}	

	else if(runNumber >= 213431 && runNumber <= 213819){ // period I	
	  Int_t bin = e12Tvhm1_I_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_I_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_I_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_I_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_I_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_I_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_I_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_I_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_I_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_I_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_I_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_I_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_I_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_I_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_I_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}	

	else if(runNumber >= 213900 && runNumber <= 215091){ // period J	
	  Int_t bin = e12Tvhm1_J_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_J_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_J_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_J_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_J_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_J_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_J_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_J_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_J_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_J_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_J_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_J_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_J_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_J_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_J_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}	

	else if(runNumber >= 215414 && runNumber <= 215643){ // period L	
	  Int_t bin = e12Tvhm1_L_Hist_SF_pT->FindBin(lepVec.Pt()/1000.);
      	  if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_L_Hist_SF_pT->FindBin(100-0.0001);
      	  sf =  e12Tvhm1_L_Hist_SF_pT->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_L_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_SF_pT_errUp->GetBinError(bin) / e12Tvhm1_L_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_L_Hist_SF_pT->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_SF_pT_errDown->GetBinError(bin) / e12Tvhm1_L_Hist_SF_pT->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }

      	  bin = e12Tvhm1_L_Hist_SF_eta->FindBin(lepVec.Eta());
      	  if(TMath::Abs(lepVec.Eta()) >= 2.47) 
      	    bin = lepVec.Eta() > 0 ? e12Tvhm1_L_Hist_SF_eta->FindBin(2.47-0.001) : e12Tvhm1_L_Hist_SF_eta->FindBin(-2.47+0.001);
      	  sf *= e12Tvhm1_L_Hist_SF_eta->GetBinContent(bin);
      	  if(mode == 1){
      	    Double_t sfErr_tmp = e12Tvhm1_L_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_SF_eta_errUp->GetBinError(bin) / e12Tvhm1_L_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  } 
      	  else if(mode == -1){
      	    Double_t sfErr_tmp = e12Tvhm1_L_Hist_SF_eta->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_SF_eta_errDown->GetBinError(bin) / e12Tvhm1_L_Hist_SF_eta->GetBinContent(bin) : 0 ;
      	    sfErr += TMath::Power(sfErr_tmp, 2);
      	  }
	}	

  	sfErr = TMath::Sqrt(sfErr) / sf;
	
      }      
    }  
  } 
  
  else{
    if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getSFElec Trigger " << trigger << " not supported... Using SF " << sf << endl ;
  }
  
  if(doDebug == kTRUE ) cout << " SF = " << sf << " err = " << sfErr << endl ;
  
  if(mode == 1){
    sf += sfErr ;
  } 
  else if(mode == -1){
    sf -= sfErr ;
  }
  
  return sf ;
  
}
  
Double_t HSG4LepLepTriggerSF::getDataEffElec( const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const
{

  Double_t sf = 1.0 ; 
  Double_t sfErr = 0.0 ;
  
  if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffElec pT " << lepVec.Pt() << " RunNumber = " << runNumber << " Trigger = " << trigger << endl ;
  
  if(trigger == "e12Tvhl1"){
  
    if(lepVec.Pt() < 12e3){      
      // Scale factors not valid for pt < 12 GeV
      if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffElec pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl ;
    } 
    else{ 
    
//      if(TMath::Abs(lepVec.Eta()) < 1.37 || (TMath::Abs(lepVec.Eta()) > 1.52 && TMath::Abs(lepVec.Eta()) <2.47)){
      if(TMath::Abs(lepVec.Eta()) < 1.37 || TMath::Abs(lepVec.Eta()) > 1.52){
      
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = e12Tvhl1_AB_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_AB_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_AB_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhl1_AB_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhl1_AB_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhl1_AB_Hist_Eff_Data_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_AB_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	else if(runNumber >= 206248 && runNumber <= 207397){ // period C
	  Int_t bin = e12Tvhl1_C_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_C_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_C_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhl1_C_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhl1_C_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhl1_C_Hist_Eff_Data_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_C_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
	  Int_t bin = e12Tvhl1_D_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_D_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_D_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhl1_D_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhl1_D_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhl1_D_Hist_Eff_Data_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_D_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
	  Int_t bin = e12Tvhl1_E_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_E_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_E_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhl1_E_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhl1_E_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhl1_E_Hist_Eff_Data_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_E_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
        
	sfErr = TMath::Sqrt( sfErr ) * sf ;
	
      }
    }
  } 

  else if(trigger == "e12Tvhm1"){
  
    if(lepVec.Pt() < 12e3){
      // Scale factors not valid for pt < 12 GeV
      if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffElec pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl ;
    } 
    else{ 
      
//      if(TMath::Abs(lepVec.Eta()) < 1.37 || (TMath::Abs(lepVec.Eta()) > 1.52 && TMath::Abs(lepVec.Eta()) <2.47)){
      if(TMath::Abs(lepVec.Eta()) < 1.37 || TMath::Abs(lepVec.Eta()) > 1.52){
      
	if(runNumber >= 200804 && runNumber <= 201556){ // period A
          Int_t bin = e12Tvhm1_A_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_A_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_A_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_A_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_A_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_A_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_A_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
        }

	if(runNumber >= 202660 && runNumber <= 205113){ // period B
          Int_t bin = e12Tvhm1_B_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_B_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_B_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_B_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_B_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_B_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_B_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
        }

	else if(runNumber >= 206248 && runNumber <= 207397){ // period C
	  Int_t bin = e12Tvhm1_C_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_C_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_C_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_C_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_C_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_C_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_C_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
	  Int_t bin = e12Tvhm1_D_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_D_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_D_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_D_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_D_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_D_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_D_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	
	else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	  Int_t bin = e12Tvhm1_E_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_E_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_E_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_E_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_E_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_E_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_E_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
	
	else if(runNumber >= 211522 && runNumber <= 212272){ // period G	
	  Int_t bin = e12Tvhm1_G_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_G_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_G_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_G_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_G_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_G_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_G_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
	
	else if(runNumber >= 212619 && runNumber <= 213359){ // period H	
	  Int_t bin = e12Tvhm1_H_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_H_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_H_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_H_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_H_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_H_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_H_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
	
	else if(runNumber >= 213431 && runNumber <= 213819){ // period I	
	  Int_t bin = e12Tvhm1_I_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_I_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_I_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_I_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_I_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_I_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_I_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
	
	else if(runNumber >= 213900 && runNumber <= 215091){ // period J	
	  Int_t bin = e12Tvhm1_J_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_J_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_J_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_J_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_J_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_J_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_J_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
	
	else if(runNumber >= 215414 && runNumber <= 215643){ // period L	
	  Int_t bin = e12Tvhm1_L_Hist_Eff_Data_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_L_Hist_Eff_Data_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_L_Hist_Eff_Data_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_Data_pT_errUp->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_Data_pT->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_Data_pT_errDown->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_Data_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_L_Hist_Eff_Data_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47) 
            bin = lepVec.Eta() > 0 ? e12Tvhm1_L_Hist_Eff_Data_eta->FindBin(2.47-0.001) : e12Tvhm1_L_Hist_Eff_Data_eta->FindBin(-2.47+0.001);    
          sf *= e12Tvhm1_L_Hist_Eff_Data_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_Data_eta_errUp->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_Data_eta->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_Data_eta_errDown->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_Data_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
	
        sfErr = TMath::Sqrt( sfErr ) * sf ;
	
      }	
    }
  } 
  
  else {
    if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffElec Trigger " << trigger << " not supported... Using SF " << sf << endl ;  
  }
      
  if(doDebug == kTRUE ) cout << " eff_data = " << sf << " err = " << sfErr << endl ;
  
  if(mode == 1){
    sf += sfErr ;
  } 
  else if(mode == -1){
    sf -= sfErr ;
  }
  
  return sf ;
  
}

Double_t HSG4LepLepTriggerSF::getMCEffElec( const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const 
{

  Double_t sf = 1.0 ; 
  Double_t sfErr = 0.0 ;
  
  if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffElec pT " << lepVec.Pt() << " RunNumber = " << runNumber << " Trigger = " << trigger << endl ;
  
  if(trigger == "e12Tvhl1"){
  
    if(lepVec.Pt() < 12e3){
      // Scale factors not valid for pt < 12 GeV
      if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffElec pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl ;
    } 
    else{
      
//      if(TMath::Abs(lepVec.Eta()) < 1.37 || (TMath::Abs(lepVec.Eta()) > 1.52 && TMath::Abs(lepVec.Eta()) <2.47)){
      if(TMath::Abs(lepVec.Eta()) < 1.37 || TMath::Abs(lepVec.Eta()) > 1.52){
      
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = e12Tvhl1_AB_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_AB_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_AB_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
    
          bin = e12Tvhl1_AB_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhl1_AB_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhl1_AB_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_AB_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_AB_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_AB_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhl1_AB_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	else if(runNumber >= 206248 && runNumber <= 207397){ // period C
	  Int_t bin = e12Tvhl1_C_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_C_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_C_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
    
          bin = e12Tvhl1_C_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhl1_C_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhl1_C_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_C_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_C_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_C_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhl1_C_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
        }

	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
	  Int_t bin = e12Tvhl1_D_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_D_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_D_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
    
          bin = e12Tvhl1_D_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhl1_D_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhl1_D_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_D_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_D_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_D_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhl1_D_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
	  Int_t bin = e12Tvhl1_E_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhl1_E_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhl1_E_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
    
          bin = e12Tvhl1_E_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhl1_E_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhl1_E_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhl1_E_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhl1_E_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhl1_E_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhl1_E_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
    
        sfErr = TMath::Sqrt( sfErr ) * sf ;
        
      }
    }
  } 
  
  else if(trigger == "e12Tvhm1"){

    if(lepVec.Pt() < 12e3){
      // Scale factors not valid for pt < 12 GeV
      if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffElec pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl ;
    } 
    else{
    
//      if(TMath::Abs(lepVec.Eta()) < 1.37 || (TMath::Abs(lepVec.Eta()) > 1.52 && TMath::Abs(lepVec.Eta()) <2.47)){
      if(TMath::Abs(lepVec.Eta()) < 1.37 || TMath::Abs(lepVec.Eta()) > 1.52){
      
	if(runNumber >= 200804 && runNumber <= 201556){ // period A
          Int_t bin = e12Tvhm1_A_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_A_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_A_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_A_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_A_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_A_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_A_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_A_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_A_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_A_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	if(runNumber >= 202660 && runNumber <= 205113){ // period B
          Int_t bin = e12Tvhm1_B_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_B_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_B_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_B_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_B_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_B_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_B_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_B_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_B_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_B_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	else if(runNumber >= 206248 && runNumber <= 207397){ // period C
	  Int_t bin = e12Tvhm1_C_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_C_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_C_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_C_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_C_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_C_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_C_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_C_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_C_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_C_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}

	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
	  Int_t bin = e12Tvhm1_D_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_D_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_D_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_D_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_D_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_D_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_D_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_D_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_D_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_D_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	
	else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	  Int_t bin = e12Tvhm1_E_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_E_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_E_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_E_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_E_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_E_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_E_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_E_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_E_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_E_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
     
	else if(runNumber >= 211522 && runNumber <= 212272){ // period G	
	  Int_t bin = e12Tvhm1_G_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_G_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_G_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_G_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_G_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_G_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_G_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_G_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_G_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_G_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
     
	else if(runNumber >= 212619 && runNumber <= 213359){ // period H	
	  Int_t bin = e12Tvhm1_H_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_H_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_H_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_H_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_H_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_H_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_H_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_H_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_H_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_H_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
     
	else if(runNumber >= 213431 && runNumber <= 213819){ // period I	
	  Int_t bin = e12Tvhm1_I_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_I_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_I_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_I_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_I_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_I_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_I_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_I_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_I_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_I_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
     
	else if(runNumber >= 213900 && runNumber <= 215091){ // period J	
	  Int_t bin = e12Tvhm1_J_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_J_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_J_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_J_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_J_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_J_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_J_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_J_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_J_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_J_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
     
	else if(runNumber >= 215414 && runNumber <= 215643){ // period L	
	  Int_t bin = e12Tvhm1_L_Hist_Eff_MC_pT->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = e12Tvhm1_L_Hist_Eff_MC_pT->FindBin(100-0.0001);
          sf = e12Tvhm1_L_Hist_Eff_MC_pT->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_MC_pT_errUp->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_MC_pT->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_MC_pT_errDown->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_MC_pT->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
     
          bin = e12Tvhm1_L_Hist_Eff_MC_eta->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.47)
            bin = lepVec.Eta() > 0 ? e12Tvhm1_L_Hist_Eff_MC_eta->FindBin(2.47-0.001) : e12Tvhm1_L_Hist_Eff_MC_eta->FindBin(-2.47+0.001);
          sf *= e12Tvhm1_L_Hist_Eff_MC_eta->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_MC_eta_errUp->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = e12Tvhm1_L_Hist_Eff_MC_eta->GetBinContent(bin) > 0 ? e12Tvhm1_L_Hist_Eff_MC_eta_errDown->GetBinError(bin) / e12Tvhm1_L_Hist_Eff_MC_eta->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}	
     
        sfErr = TMath::Sqrt( sfErr ) * sf ;
      
      }
    }
  } 
  
  else {
    if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffElec Trigger " << trigger << " not supported... Using SF " << sf << endl ;  
  }

  if(doDebug == kTRUE ) cout << " eff_mc = " << sf << " err = " << sfErr << endl ;
  
  if(mode == 1){
    sf += sfErr ;
  } 
  else if(mode == -1){
    sf -= sfErr ;
  }
  
  return sf ;
  
}

Double_t HSG4LepLepTriggerSF::getSFMuon(const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const 
{
  
  if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getSFMuon pT " << lepVec.Pt() << " eta " << lepVec.Eta() << " phi " << lepVec.Phi() << " RunNumber = " << runNumber << " mode = " << mode << endl ;
  
  Double_t sf = 1.0 ; 
  Double_t sfErr = 0.0 ;
  
  if(trigger == "mu8"){

    if(lepVec.Pt() < 8e3){
      // Scale factors not valid for pt < 8 GeV
      if(doDebug == kTRUE) cout << " HSG4LepLepTriggerSF::getSFMuon pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl;
    } 
    else{
      
      if(TMath::Abs(lepVec.Eta()) < 1.05 ){
	
	// Barrel SFs
	
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = mu8_AB_Hist_SF_pT_Barrel->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_AB_Hist_SF_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_AB_Hist_SF_pT_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_pT_Barrel_errUp->GetBinError(bin) / mu8_AB_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_pT_Barrel_errDown->GetBinError(bin) / mu8_AB_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_AB_Hist_SF_eta_Barrel->FindBin(lepVec.Eta());
          sf *= mu8_AB_Hist_SF_eta_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_eta_Barrel_errUp->GetBinError(bin) / mu8_AB_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_eta_Barrel_errDown->GetBinError(bin) / mu8_AB_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
        }
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C	
          Int_t bin = mu8_C_Hist_SF_pT_Barrel->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_C_Hist_SF_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_C_Hist_SF_pT_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_pT_Barrel_errUp->GetBinError(bin) / mu8_C_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_pT_Barrel_errDown->GetBinError(bin) / mu8_C_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_C_Hist_SF_eta_Barrel->FindBin(lepVec.Eta());
          sf *= mu8_C_Hist_SF_eta_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_eta_Barrel_errUp->GetBinError(bin) / mu8_C_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_eta_Barrel_errDown->GetBinError(bin) / mu8_C_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
          Int_t bin = mu8_D_Hist_SF_pT_Barrel->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_D_Hist_SF_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_D_Hist_SF_pT_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_pT_Barrel_errUp->GetBinError(bin) / mu8_D_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_pT_Barrel_errDown->GetBinError(bin) / mu8_D_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_D_Hist_SF_eta_Barrel->FindBin(lepVec.Eta());
          sf *= mu8_D_Hist_SF_eta_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_eta_Barrel_errUp->GetBinError(bin) / mu8_D_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_eta_Barrel_errDown->GetBinError(bin) / mu8_D_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
          Int_t bin = mu8_E_Hist_SF_pT_Barrel->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_E_Hist_SF_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_E_Hist_SF_pT_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_pT_Barrel_errUp->GetBinError(bin) / mu8_E_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_pT_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_pT_Barrel_errDown->GetBinError(bin) / mu8_E_Hist_SF_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_E_Hist_SF_eta_Barrel->FindBin(lepVec.Eta());
          sf *= mu8_E_Hist_SF_eta_Barrel->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_eta_Barrel_errUp->GetBinError(bin) / mu8_E_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_eta_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_eta_Barrel_errDown->GetBinError(bin) / mu8_E_Hist_SF_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
      } 
      
//      else if(TMath::Abs(lepVec.Eta()) < 2.3){
      else{
	
	// Endcap SFs
	
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = mu8_AB_Hist_SF_pT_Endcap->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_AB_Hist_SF_pT_Endcap->FindBin(100-0.0001);
          sf =  mu8_AB_Hist_SF_pT_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_pT_Endcap_errUp->GetBinError(bin) / mu8_AB_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_pT_Endcap_errDown->GetBinError(bin) / mu8_AB_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_AB_Hist_SF_eta_Endcap->FindBin(lepVec.Eta());
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_AB_Hist_SF_eta_Endcap->FindBin(2.2-0.001) : mu8_AB_Hist_SF_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_AB_Hist_SF_eta_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_eta_Endcap_errUp->GetBinError(bin) / mu8_AB_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_SF_eta_Endcap_errDown->GetBinError(bin) / mu8_AB_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C	
          Int_t bin = mu8_C_Hist_SF_pT_Endcap->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_C_Hist_SF_pT_Endcap->FindBin(100-0.0001);
          sf =  mu8_C_Hist_SF_pT_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_pT_Endcap_errUp->GetBinError(bin) / mu8_C_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_pT_Endcap_errDown->GetBinError(bin) / mu8_C_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_C_Hist_SF_eta_Endcap->FindBin(lepVec.Eta());
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_C_Hist_SF_eta_Endcap->FindBin(2.2-0.001) : mu8_C_Hist_SF_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_C_Hist_SF_eta_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_eta_Endcap_errUp->GetBinError(bin) / mu8_C_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_SF_eta_Endcap_errDown->GetBinError(bin) / mu8_C_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
          Int_t bin = mu8_D_Hist_SF_pT_Endcap->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_D_Hist_SF_pT_Endcap->FindBin(100-0.0001);
          sf =  mu8_D_Hist_SF_pT_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_pT_Endcap_errUp->GetBinError(bin) / mu8_D_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_pT_Endcap_errDown->GetBinError(bin) / mu8_D_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_D_Hist_SF_eta_Endcap->FindBin(lepVec.Eta());
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_D_Hist_SF_eta_Endcap->FindBin(2.2-0.001) : mu8_D_Hist_SF_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_D_Hist_SF_eta_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_eta_Endcap_errUp->GetBinError(bin) / mu8_D_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_SF_eta_Endcap_errDown->GetBinError(bin) / mu8_D_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
          Int_t bin = mu8_E_Hist_SF_pT_Endcap->FindBin(lepVec.Pt()/1000.);
          if(lepVec.Pt() >= 100e3) bin = mu8_E_Hist_SF_pT_Endcap->FindBin(100-0.0001);
          sf =  mu8_E_Hist_SF_pT_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_pT_Endcap_errUp->GetBinError(bin) / mu8_E_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_pT_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_pT_Endcap_errDown->GetBinError(bin) / mu8_E_Hist_SF_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_E_Hist_SF_eta_Endcap->FindBin(lepVec.Eta());
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_E_Hist_SF_eta_Endcap->FindBin(2.2-0.001) : mu8_E_Hist_SF_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_E_Hist_SF_eta_Endcap->GetBinContent(bin);
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_eta_Endcap_errUp->GetBinError(bin) / mu8_E_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_SF_eta_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_SF_eta_Endcap_errDown->GetBinError(bin) / mu8_E_Hist_SF_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
      }
      
      sfErr = TMath::Sqrt(sfErr) / sf;
      
    }
  } 
  
  else {    
    if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getSFMuon Trigger " << trigger << " not supported... Using SF " << sf << endl ;  
  }
  
  if(doDebug == kTRUE ) cout << " SF = " << sf << " err = " << sfErr << endl ;
  
  if(mode == 1){
    sf += sfErr ;
  } 
  else if(mode == -1){
    sf -= sfErr ;
  }
  
  return sf ;
  
}

Double_t HSG4LepLepTriggerSF::getDataEffMuon(const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const 
{

  if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffMuon pT " << lepVec.Pt() << " eta " << lepVec.Eta() << " phi " << lepVec.Phi() << " RunNumber = " << runNumber << " mode = " << mode << endl ;
  
  Double_t sf = 1.0 ; 
  Double_t sfErr = 0.0 ;
  
  if(trigger == "mu8"){
    
    if(lepVec.Pt() < 8e3){
      // Scale factors not valid for pt < 8 GeV
      if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffMuon pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl ;
    } 
    else{
      
      if(TMath::Abs(lepVec.Eta()) < 1.05 ){
	
	// Barrel SFs
	
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = mu8_AB_Hist_Eff_Data_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_AB_Hist_Eff_Data_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_AB_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_pT_Barrel_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_pT_Barrel_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_AB_Hist_Eff_Data_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_AB_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_eta_Barrel_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_eta_Barrel_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C	
          Int_t bin = mu8_C_Hist_Eff_Data_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_C_Hist_Eff_Data_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_C_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_pT_Barrel_errUp->GetBinError(bin) / mu8_C_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_pT_Barrel_errDown->GetBinError(bin) / mu8_C_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_C_Hist_Eff_Data_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_C_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_eta_Barrel_errUp->GetBinError(bin) / mu8_C_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_eta_Barrel_errDown->GetBinError(bin) / mu8_C_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
          Int_t bin = mu8_D_Hist_Eff_Data_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_D_Hist_Eff_Data_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_D_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_pT_Barrel_errUp->GetBinError(bin) / mu8_D_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_pT_Barrel_errDown->GetBinError(bin) / mu8_D_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_D_Hist_Eff_Data_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_D_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_eta_Barrel_errUp->GetBinError(bin) / mu8_D_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_eta_Barrel_errDown->GetBinError(bin) / mu8_D_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
          Int_t bin = mu8_E_Hist_Eff_Data_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_E_Hist_Eff_Data_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_E_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_pT_Barrel_errUp->GetBinError(bin) / mu8_E_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_pT_Barrel_errDown->GetBinError(bin) / mu8_E_Hist_Eff_Data_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_E_Hist_Eff_Data_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_E_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_eta_Barrel_errUp->GetBinError(bin) / mu8_E_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_eta_Barrel_errDown->GetBinError(bin) / mu8_E_Hist_Eff_Data_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
      } 
      
//      else if(TMath::Abs(lepVec.Eta()) < 2.3 ){
      else{
	
	// Endcap SFs
	
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = mu8_AB_Hist_Eff_Data_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_AB_Hist_Eff_Data_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_AB_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_pT_Endcap_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_pT_Endcap_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_AB_Hist_Eff_Data_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_AB_Hist_Eff_Data_eta_Endcap->FindBin(2.2-0.001) : mu8_AB_Hist_Eff_Data_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_AB_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_eta_Endcap_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_Data_eta_Endcap_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C	
          Int_t bin = mu8_C_Hist_Eff_Data_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_C_Hist_Eff_Data_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_C_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_pT_Endcap_errUp->GetBinError(bin) / mu8_C_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_pT_Endcap_errDown->GetBinError(bin) / mu8_C_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_C_Hist_Eff_Data_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_C_Hist_Eff_Data_eta_Endcap->FindBin(2.2-0.001) : mu8_C_Hist_Eff_Data_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_C_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_eta_Endcap_errUp->GetBinError(bin) / mu8_C_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_Data_eta_Endcap_errDown->GetBinError(bin) / mu8_C_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
          Int_t bin = mu8_D_Hist_Eff_Data_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_D_Hist_Eff_Data_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_D_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_pT_Endcap_errUp->GetBinError(bin) / mu8_D_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_pT_Endcap_errDown->GetBinError(bin) / mu8_D_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_D_Hist_Eff_Data_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_D_Hist_Eff_Data_eta_Endcap->FindBin(2.2-0.001) : mu8_D_Hist_Eff_Data_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_D_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_eta_Endcap_errUp->GetBinError(bin) / mu8_D_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_Data_eta_Endcap_errDown->GetBinError(bin) / mu8_D_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
 	else if(runNumber >= 209074){ // period E	
         Int_t bin = mu8_E_Hist_Eff_Data_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_E_Hist_Eff_Data_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_E_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_pT_Endcap_errUp->GetBinError(bin) / mu8_E_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_pT_Endcap_errDown->GetBinError(bin) / mu8_E_Hist_Eff_Data_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_E_Hist_Eff_Data_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_E_Hist_Eff_Data_eta_Endcap->FindBin(2.2-0.001) : mu8_E_Hist_Eff_Data_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_E_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_eta_Endcap_errUp->GetBinError(bin) / mu8_E_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_Data_eta_Endcap_errDown->GetBinError(bin) / mu8_E_Hist_Eff_Data_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
      }

      sfErr = TMath::Sqrt( sfErr ) * sf ;
    }
  } 
  
  else{
    if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getDataEffMuon Trigger " << trigger << " not supported... Using SF " << sf << endl ;  
  }
  
  if(doDebug == kTRUE ) cout << " eff_data = " << sf << " err = " << sfErr << endl ;
  
  if(mode == 1){
    sf += sfErr ;
  } 
  else if(mode == -1){
    sf -= sfErr ;
  }
  
  return sf ;
  
}

Double_t HSG4LepLepTriggerSF::getMCEffMuon(const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const 
{

  if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffMuon pT " << lepVec.Pt() << " eta " << lepVec.Eta() << " phi " << lepVec.Phi() << " RunNumber = " << runNumber << " mode = " << mode << endl ;
  
  Double_t sf = 1.0 ; 
  Double_t sfErr = 0.0 ;
  
  if(trigger == "mu8"){

    if(lepVec.Pt() < 8e3){
      // Scale factors not valid for pt < 8 GeV
      if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffMuon pT " << lepVec.Pt() << " out of range.. Using SF: " << sf << endl ;
    } 
    else{
      
      if(TMath::Abs(lepVec.Eta()) < 1.05 ){
	
	// Barrel SFs
	
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = mu8_AB_Hist_Eff_MC_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_AB_Hist_Eff_MC_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_AB_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_pT_Barrel_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_pT_Barrel_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_AB_Hist_Eff_MC_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_AB_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_eta_Barrel_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_eta_Barrel_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C	
          Int_t bin = mu8_C_Hist_Eff_MC_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_C_Hist_Eff_MC_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_C_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_pT_Barrel_errUp->GetBinError(bin) / mu8_C_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_pT_Barrel_errDown->GetBinError(bin) / mu8_C_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_C_Hist_Eff_MC_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_C_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_eta_Barrel_errUp->GetBinError(bin) / mu8_C_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_eta_Barrel_errDown->GetBinError(bin) / mu8_C_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
          Int_t bin = mu8_D_Hist_Eff_MC_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_D_Hist_Eff_MC_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_D_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_pT_Barrel_errUp->GetBinError(bin) / mu8_D_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_pT_Barrel_errDown->GetBinError(bin) / mu8_D_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_D_Hist_Eff_MC_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_D_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_eta_Barrel_errUp->GetBinError(bin) / mu8_D_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_eta_Barrel_errDown->GetBinError(bin) / mu8_D_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
          Int_t bin = mu8_E_Hist_Eff_MC_pT_Barrel->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_E_Hist_Eff_MC_pT_Barrel->FindBin(100-0.0001);
          sf =  mu8_E_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_pT_Barrel_errUp->GetBinError(bin) / mu8_E_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_pT_Barrel_errDown->GetBinError(bin) / mu8_E_Hist_Eff_MC_pT_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_E_Hist_Eff_MC_eta_Barrel->FindBin(lepVec.Eta()) ;
          sf *= mu8_E_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_eta_Barrel_errUp->GetBinError(bin) / mu8_E_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_eta_Barrel_errDown->GetBinError(bin) / mu8_E_Hist_Eff_MC_eta_Barrel->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
      } 
      
//      else if(TMath::Abs(lepVec.Eta()) < 2.3 ){
      else{
	
	// Endcap SFs
	
	if(runNumber >= 200804 && runNumber <= 205113){ // period AB
          Int_t bin = mu8_AB_Hist_Eff_MC_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_AB_Hist_Eff_MC_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_AB_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_pT_Endcap_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_pT_Endcap_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_AB_Hist_Eff_MC_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_AB_Hist_Eff_MC_eta_Endcap->FindBin(2.2-0.001) : mu8_AB_Hist_Eff_MC_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_AB_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_eta_Endcap_errUp->GetBinError(bin) / mu8_AB_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_AB_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_AB_Hist_Eff_MC_eta_Endcap_errDown->GetBinError(bin) / mu8_AB_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 206248 && runNumber <= 207397){ // period C	
          Int_t bin = mu8_C_Hist_Eff_MC_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_C_Hist_Eff_MC_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_C_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_pT_Endcap_errUp->GetBinError(bin) / mu8_C_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_pT_Endcap_errDown->GetBinError(bin) / mu8_C_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_C_Hist_Eff_MC_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_C_Hist_Eff_MC_eta_Endcap->FindBin(2.2-0.001) : mu8_C_Hist_Eff_MC_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_C_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_eta_Endcap_errUp->GetBinError(bin) / mu8_C_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_C_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_C_Hist_Eff_MC_eta_Endcap_errDown->GetBinError(bin) / mu8_C_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	else if(runNumber >= 207447 && runNumber <= 209025){ // period D	
          Int_t bin = mu8_D_Hist_Eff_MC_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_D_Hist_Eff_MC_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_D_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_pT_Endcap_errUp->GetBinError(bin) / mu8_D_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_pT_Endcap_errDown->GetBinError(bin) / mu8_D_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_D_Hist_Eff_MC_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_D_Hist_Eff_MC_eta_Endcap->FindBin(2.2-0.001) : mu8_D_Hist_Eff_MC_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_D_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_eta_Endcap_errUp->GetBinError(bin) / mu8_D_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_D_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_D_Hist_Eff_MC_eta_Endcap_errDown->GetBinError(bin) / mu8_D_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
	//else if(runNumber >= 209074 && runNumber <= 210308){ // period E	
	else if(runNumber >= 209074){ // period E	
          Int_t bin = mu8_E_Hist_Eff_MC_pT_Endcap->FindBin( lepVec.Pt() / 1000. ) ;
          if(lepVec.Pt() >= 100e3) bin = mu8_E_Hist_Eff_MC_pT_Endcap->FindBin(100-0.0001);
          sf = mu8_E_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_pT_Endcap_errUp->GetBinError(bin) / mu8_E_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_pT_Endcap_errDown->GetBinError(bin) / mu8_E_Hist_Eff_MC_pT_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }

          bin = mu8_E_Hist_Eff_MC_eta_Endcap->FindBin(lepVec.Eta()) ;
          if(TMath::Abs(lepVec.Eta()) >= 2.27)
            bin = lepVec.Eta() > 0 ? mu8_E_Hist_Eff_MC_eta_Endcap->FindBin(2.2-0.001) : mu8_E_Hist_Eff_MC_eta_Endcap->FindBin(-2.2+0.001);
          sf *= mu8_E_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) ;
          if(mode == 1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_eta_Endcap_errUp->GetBinError(bin) / mu8_E_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0 ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          } 
          else if(mode == -1){
            Double_t sfErr_tmp = mu8_E_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) > 0 ? mu8_E_Hist_Eff_MC_eta_Endcap_errDown->GetBinError(bin) / mu8_E_Hist_Eff_MC_eta_Endcap->GetBinContent(bin) : 0  ;
            sfErr += TMath::Power(sfErr_tmp, 2);
          }
	}
      }
      
      sfErr = TMath::Sqrt( sfErr ) * sf ;
    }
  } 
  
  else{
    if(doDebug == kTRUE ) cout << " HSG4LepLepTriggerSF::getMCEffMuon Trigger " << trigger << " not supported... Using SF " << sf << endl ;  
  }
  
  if(doDebug == kTRUE ) cout << " eff_data = " << sf << " err = " << sfErr << endl ;
  
  if(mode == 1){
    sf += sfErr ;
  } 
  else if(mode == -1){
    sf -= sfErr ;
  }
  
  return sf ;
  
}

