#ifndef H_HSG4LEPLEPTRIGGERSF
#define H_HSG4LEPLEPTRIGGERSF

#include <string>
#include "TLorentzVector.h"
#include "TString.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"

//
// Class to read and apply HSG4Trigger scale factors
// calculated by Juan Ferrer (juan.valls.ferrer@cern.ch)
//
// Code by Matthew Beckingham (matthew.beckingham@cern.ch) 27.06.12
//
// Code adapted for HSG4 leplep analysis by Lidia Dell'Asta 09.08.12

class HSG4LepLepTriggerSF
{

public:

    HSG4LepLepTriggerSF( const std::string& path = "" , const Bool_t debug = kFALSE) ;
    ~HSG4LepLepTriggerSF() ;

    // mode = 0 (central value), mode = 1 (+1 sigma), mode = -1 (-1 sigma)
    Double_t getSFElec(const TLorentzVector& lepVec, const UInt_t runNumber, const TString trigger, const Int_t mode) const ;
    Double_t getSFMuon(const TLorentzVector& lepVec, const UInt_t runNumber, const TString trigger, const Int_t mode) const ;

    Double_t getDataEffElec(const TLorentzVector& lepVec, const UInt_t runNumber, const TString trigger, const Int_t mode) const ;
    Double_t getDataEffMuon(const TLorentzVector& lepVec, const UInt_t runNumber, const TString trigger, const Int_t mode) const ;

    Double_t getMCEffElec(const TLorentzVector& lepVec, const UInt_t runNumber, const TString trigger, const Int_t mode) const ;
    Double_t getMCEffMuon(const TLorentzVector& lepVec, const UInt_t runNumber, const TString trigger, const Int_t mode) const ;

    void loadFiles( const TString& path) ;

private:

    void readHisto( const TString& setName, TGraphAsymmErrors* const inGraph, TH1F*& outHist, TH1F*& outHistErrUp, TH1F*& outHistErrDown) ;
    void readHisto( const TString& setName, TGraphAsymmErrors* const inGraph, TGraphAsymmErrors* const inGraphSys, TH1F*& outHist, TH1F*& outHistErrUp, TH1F*& outHistErrDown) ;
    void normaliseHisto( TH1F* hist) ;


    // Electron trigger SFs
    // 2012 Triggers
    
    TH1F* e12Tvhl1_AB_Hist_SF_eta;
    TH1F* e12Tvhl1_AB_Hist_SF_pT;    
    TH1F* e12Tvhl1_AB_Hist_SF_eta_errUp;
    TH1F* e12Tvhl1_AB_Hist_SF_pT_errUp;
    TH1F* e12Tvhl1_AB_Hist_SF_eta_errDown;
    TH1F* e12Tvhl1_AB_Hist_SF_pT_errDown;
   
    TH1F* e12Tvhl1_C_Hist_SF_eta;
    TH1F* e12Tvhl1_C_Hist_SF_pT;    
    TH1F* e12Tvhl1_C_Hist_SF_eta_errUp;
    TH1F* e12Tvhl1_C_Hist_SF_pT_errUp;
    TH1F* e12Tvhl1_C_Hist_SF_eta_errDown;
    TH1F* e12Tvhl1_C_Hist_SF_pT_errDown;
    
    TH1F* e12Tvhl1_D_Hist_SF_eta;
    TH1F* e12Tvhl1_D_Hist_SF_pT;    
    TH1F* e12Tvhl1_D_Hist_SF_eta_errUp;
    TH1F* e12Tvhl1_D_Hist_SF_pT_errUp;
    TH1F* e12Tvhl1_D_Hist_SF_eta_errDown;
    TH1F* e12Tvhl1_D_Hist_SF_pT_errDown;
    
    TH1F* e12Tvhl1_E_Hist_SF_eta;
    TH1F* e12Tvhl1_E_Hist_SF_pT;    
    TH1F* e12Tvhl1_E_Hist_SF_eta_errUp;
    TH1F* e12Tvhl1_E_Hist_SF_pT_errUp;
    TH1F* e12Tvhl1_E_Hist_SF_eta_errDown;
    TH1F* e12Tvhl1_E_Hist_SF_pT_errDown;

    TH1F* e12Tvhm1_A_Hist_SF_eta;
    TH1F* e12Tvhm1_A_Hist_SF_pT;
    TH1F* e12Tvhm1_A_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_A_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_A_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_A_Hist_SF_pT_errDown;
    
    TH1F* e12Tvhm1_B_Hist_SF_eta;
    TH1F* e12Tvhm1_B_Hist_SF_pT;
    TH1F* e12Tvhm1_B_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_B_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_B_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_B_Hist_SF_pT_errDown;
        
    TH1F* e12Tvhm1_C_Hist_SF_eta;
    TH1F* e12Tvhm1_C_Hist_SF_pT;
    TH1F* e12Tvhm1_C_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_C_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_C_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_C_Hist_SF_pT_errDown;    
    
    TH1F* e12Tvhm1_D_Hist_SF_eta;
    TH1F* e12Tvhm1_D_Hist_SF_pT;
    TH1F* e12Tvhm1_D_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_D_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_D_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_D_Hist_SF_pT_errDown;    
    
    TH1F* e12Tvhm1_E_Hist_SF_eta;
    TH1F* e12Tvhm1_E_Hist_SF_pT;
    TH1F* e12Tvhm1_E_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_E_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_E_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_E_Hist_SF_pT_errDown;

    TH1F* e12Tvhm1_G_Hist_SF_eta;
    TH1F* e12Tvhm1_G_Hist_SF_pT;
    TH1F* e12Tvhm1_G_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_G_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_G_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_G_Hist_SF_pT_errDown;

    TH1F* e12Tvhm1_H_Hist_SF_eta;
    TH1F* e12Tvhm1_H_Hist_SF_pT;
    TH1F* e12Tvhm1_H_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_H_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_H_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_H_Hist_SF_pT_errDown;

    TH1F* e12Tvhm1_I_Hist_SF_eta;
    TH1F* e12Tvhm1_I_Hist_SF_pT;
    TH1F* e12Tvhm1_I_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_I_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_I_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_I_Hist_SF_pT_errDown;

    TH1F* e12Tvhm1_J_Hist_SF_eta;
    TH1F* e12Tvhm1_J_Hist_SF_pT;
    TH1F* e12Tvhm1_J_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_J_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_J_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_J_Hist_SF_pT_errDown;

    TH1F* e12Tvhm1_L_Hist_SF_eta;
    TH1F* e12Tvhm1_L_Hist_SF_pT;
    TH1F* e12Tvhm1_L_Hist_SF_eta_errUp;
    TH1F* e12Tvhm1_L_Hist_SF_pT_errUp;
    TH1F* e12Tvhm1_L_Hist_SF_eta_errDown;
    TH1F* e12Tvhm1_L_Hist_SF_pT_errDown;

    // Electron data efficiency histograms
    // 2012 Triggers
    
    TH1F* e12Tvhl1_AB_Hist_Eff_Data_eta;
    TH1F* e12Tvhl1_AB_Hist_Eff_Data_pT;    
    TH1F* e12Tvhl1_AB_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhl1_AB_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhl1_AB_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhl1_AB_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhl1_C_Hist_Eff_Data_eta;
    TH1F* e12Tvhl1_C_Hist_Eff_Data_pT;	 
    TH1F* e12Tvhl1_C_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhl1_C_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhl1_C_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhl1_C_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhl1_D_Hist_Eff_Data_eta;
    TH1F* e12Tvhl1_D_Hist_Eff_Data_pT;	
    TH1F* e12Tvhl1_D_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhl1_D_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhl1_D_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhl1_D_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhl1_E_Hist_Eff_Data_eta;
    TH1F* e12Tvhl1_E_Hist_Eff_Data_pT;	
    TH1F* e12Tvhl1_E_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhl1_E_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhl1_E_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhl1_E_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhm1_A_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_A_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_A_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_A_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_A_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_A_Hist_Eff_Data_pT_errDown;

    TH1F* e12Tvhm1_B_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_B_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_B_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_B_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_B_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_B_Hist_Eff_Data_pT_errDown;

    TH1F* e12Tvhm1_C_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_C_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_C_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_C_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_C_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_C_Hist_Eff_Data_pT_errDown;

    TH1F* e12Tvhm1_D_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_D_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_D_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_D_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_D_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_D_Hist_Eff_Data_pT_errDown;

    TH1F* e12Tvhm1_E_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_E_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_E_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_E_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_E_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_E_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhm1_G_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_G_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_G_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_G_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_G_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_G_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhm1_H_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_H_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_H_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_H_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_H_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_H_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhm1_I_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_I_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_I_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_I_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_I_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_I_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhm1_J_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_J_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_J_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_J_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_J_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_J_Hist_Eff_Data_pT_errDown;
    
    TH1F* e12Tvhm1_L_Hist_Eff_Data_eta;
    TH1F* e12Tvhm1_L_Hist_Eff_Data_pT;
    TH1F* e12Tvhm1_L_Hist_Eff_Data_eta_errUp;
    TH1F* e12Tvhm1_L_Hist_Eff_Data_pT_errUp;
    TH1F* e12Tvhm1_L_Hist_Eff_Data_eta_errDown;
    TH1F* e12Tvhm1_L_Hist_Eff_Data_pT_errDown;
    
    // Electron MC efficiency histograms
    // 2012 Triggers
    
    TH1F* e12Tvhl1_AB_Hist_Eff_MC_eta;
    TH1F* e12Tvhl1_AB_Hist_Eff_MC_pT;    
    TH1F* e12Tvhl1_AB_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhl1_AB_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhl1_AB_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhl1_AB_Hist_Eff_MC_pT_errDown;
    
    TH1F* e12Tvhl1_C_Hist_Eff_MC_eta;
    TH1F* e12Tvhl1_C_Hist_Eff_MC_pT;    
    TH1F* e12Tvhl1_C_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhl1_C_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhl1_C_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhl1_C_Hist_Eff_MC_pT_errDown;
   
    TH1F* e12Tvhl1_D_Hist_Eff_MC_eta;
    TH1F* e12Tvhl1_D_Hist_Eff_MC_pT;    
    TH1F* e12Tvhl1_D_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhl1_D_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhl1_D_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhl1_D_Hist_Eff_MC_pT_errDown;
    
    TH1F* e12Tvhl1_E_Hist_Eff_MC_eta;
    TH1F* e12Tvhl1_E_Hist_Eff_MC_pT;    
    TH1F* e12Tvhl1_E_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhl1_E_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhl1_E_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhl1_E_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_A_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_A_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_A_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_A_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_A_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_A_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_B_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_B_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_B_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_B_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_B_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_B_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_C_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_C_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_C_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_C_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_C_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_C_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_D_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_D_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_D_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_D_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_D_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_D_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_E_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_E_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_E_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_E_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_E_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_E_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_G_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_G_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_G_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_G_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_G_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_G_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_H_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_H_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_H_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_H_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_H_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_H_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_I_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_I_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_I_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_I_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_I_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_I_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_J_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_J_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_J_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_J_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_J_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_J_Hist_Eff_MC_pT_errDown;

    TH1F* e12Tvhm1_L_Hist_Eff_MC_eta;
    TH1F* e12Tvhm1_L_Hist_Eff_MC_pT;
    TH1F* e12Tvhm1_L_Hist_Eff_MC_eta_errUp;
    TH1F* e12Tvhm1_L_Hist_Eff_MC_pT_errUp;
    TH1F* e12Tvhm1_L_Hist_Eff_MC_eta_errDown;
    TH1F* e12Tvhm1_L_Hist_Eff_MC_pT_errDown;


    // Muons trigger SFs
    // 2012 Triggers
    
    TH1F* mu8_AB_Hist_SF_eta_Barrel;
    TH1F* mu8_AB_Hist_SF_pT_Barrel;    
    TH1F* mu8_AB_Hist_SF_eta_Barrel_errUp;
    TH1F* mu8_AB_Hist_SF_pT_Barrel_errUp;
    TH1F* mu8_AB_Hist_SF_eta_Barrel_errDown;
    TH1F* mu8_AB_Hist_SF_pT_Barrel_errDown;

    TH1F* mu8_AB_Hist_SF_eta_Endcap;
    TH1F* mu8_AB_Hist_SF_pT_Endcap;
    TH1F* mu8_AB_Hist_SF_eta_Endcap_errUp;
    TH1F* mu8_AB_Hist_SF_pT_Endcap_errUp;    
    TH1F* mu8_AB_Hist_SF_eta_Endcap_errDown;
    TH1F* mu8_AB_Hist_SF_pT_Endcap_errDown;
        
    TH1F* mu8_C_Hist_SF_eta_Barrel;
    TH1F* mu8_C_Hist_SF_pT_Barrel;    
    TH1F* mu8_C_Hist_SF_eta_Barrel_errUp;
    TH1F* mu8_C_Hist_SF_pT_Barrel_errUp;
    TH1F* mu8_C_Hist_SF_eta_Barrel_errDown;
    TH1F* mu8_C_Hist_SF_pT_Barrel_errDown;

    TH1F* mu8_C_Hist_SF_eta_Endcap;
    TH1F* mu8_C_Hist_SF_pT_Endcap;
    TH1F* mu8_C_Hist_SF_eta_Endcap_errUp;
    TH1F* mu8_C_Hist_SF_pT_Endcap_errUp;    
    TH1F* mu8_C_Hist_SF_eta_Endcap_errDown;
    TH1F* mu8_C_Hist_SF_pT_Endcap_errDown;
        
    TH1F* mu8_D_Hist_SF_eta_Barrel;
    TH1F* mu8_D_Hist_SF_pT_Barrel;    
    TH1F* mu8_D_Hist_SF_eta_Barrel_errUp;
    TH1F* mu8_D_Hist_SF_pT_Barrel_errUp;
    TH1F* mu8_D_Hist_SF_eta_Barrel_errDown;
    TH1F* mu8_D_Hist_SF_pT_Barrel_errDown;

    TH1F* mu8_D_Hist_SF_eta_Endcap;
    TH1F* mu8_D_Hist_SF_pT_Endcap;
    TH1F* mu8_D_Hist_SF_eta_Endcap_errUp;
    TH1F* mu8_D_Hist_SF_pT_Endcap_errUp;    
    TH1F* mu8_D_Hist_SF_eta_Endcap_errDown;
    TH1F* mu8_D_Hist_SF_pT_Endcap_errDown;
       
    TH1F* mu8_E_Hist_SF_eta_Barrel;
    TH1F* mu8_E_Hist_SF_pT_Barrel;    
    TH1F* mu8_E_Hist_SF_eta_Barrel_errUp;
    TH1F* mu8_E_Hist_SF_pT_Barrel_errUp;
    TH1F* mu8_E_Hist_SF_eta_Barrel_errDown;
    TH1F* mu8_E_Hist_SF_pT_Barrel_errDown;

    TH1F* mu8_E_Hist_SF_eta_Endcap;
    TH1F* mu8_E_Hist_SF_pT_Endcap;
    TH1F* mu8_E_Hist_SF_eta_Endcap_errUp;
    TH1F* mu8_E_Hist_SF_pT_Endcap_errUp;    
    TH1F* mu8_E_Hist_SF_eta_Endcap_errDown;
    TH1F* mu8_E_Hist_SF_pT_Endcap_errDown;
    
    // Muon data efficiency histograms
    // 2012 Triggers
    
    TH1F* mu8_AB_Hist_Eff_Data_eta_Barrel;
    TH1F* mu8_AB_Hist_Eff_Data_pT_Barrel;
    TH1F* mu8_AB_Hist_Eff_Data_eta_Barrel_errUp;
    TH1F* mu8_AB_Hist_Eff_Data_pT_Barrel_errUp;
    TH1F* mu8_AB_Hist_Eff_Data_eta_Barrel_errDown;
    TH1F* mu8_AB_Hist_Eff_Data_pT_Barrel_errDown;

    TH1F* mu8_AB_Hist_Eff_Data_eta_Endcap;
    TH1F* mu8_AB_Hist_Eff_Data_pT_Endcap;
    TH1F* mu8_AB_Hist_Eff_Data_eta_Endcap_errUp;
    TH1F* mu8_AB_Hist_Eff_Data_pT_Endcap_errUp;
    TH1F* mu8_AB_Hist_Eff_Data_eta_Endcap_errDown;
    TH1F* mu8_AB_Hist_Eff_Data_pT_Endcap_errDown;

    TH1F* mu8_C_Hist_Eff_Data_eta_Barrel;
    TH1F* mu8_C_Hist_Eff_Data_pT_Barrel;
    TH1F* mu8_C_Hist_Eff_Data_eta_Barrel_errUp;
    TH1F* mu8_C_Hist_Eff_Data_pT_Barrel_errUp;
    TH1F* mu8_C_Hist_Eff_Data_eta_Barrel_errDown;
    TH1F* mu8_C_Hist_Eff_Data_pT_Barrel_errDown;

    TH1F* mu8_C_Hist_Eff_Data_eta_Endcap;
    TH1F* mu8_C_Hist_Eff_Data_pT_Endcap;
    TH1F* mu8_C_Hist_Eff_Data_eta_Endcap_errUp;
    TH1F* mu8_C_Hist_Eff_Data_pT_Endcap_errUp;
    TH1F* mu8_C_Hist_Eff_Data_eta_Endcap_errDown;
    TH1F* mu8_C_Hist_Eff_Data_pT_Endcap_errDown;

    TH1F* mu8_D_Hist_Eff_Data_eta_Barrel;
    TH1F* mu8_D_Hist_Eff_Data_pT_Barrel;
    TH1F* mu8_D_Hist_Eff_Data_eta_Barrel_errUp;
    TH1F* mu8_D_Hist_Eff_Data_pT_Barrel_errUp;
    TH1F* mu8_D_Hist_Eff_Data_eta_Barrel_errDown;
    TH1F* mu8_D_Hist_Eff_Data_pT_Barrel_errDown;

    TH1F* mu8_D_Hist_Eff_Data_eta_Endcap;
    TH1F* mu8_D_Hist_Eff_Data_pT_Endcap;
    TH1F* mu8_D_Hist_Eff_Data_eta_Endcap_errUp;
    TH1F* mu8_D_Hist_Eff_Data_pT_Endcap_errUp;
    TH1F* mu8_D_Hist_Eff_Data_eta_Endcap_errDown;
    TH1F* mu8_D_Hist_Eff_Data_pT_Endcap_errDown;

    TH1F* mu8_E_Hist_Eff_Data_eta_Barrel;
    TH1F* mu8_E_Hist_Eff_Data_pT_Barrel;
    TH1F* mu8_E_Hist_Eff_Data_eta_Barrel_errUp;
    TH1F* mu8_E_Hist_Eff_Data_pT_Barrel_errUp;
    TH1F* mu8_E_Hist_Eff_Data_eta_Barrel_errDown;
    TH1F* mu8_E_Hist_Eff_Data_pT_Barrel_errDown;

    TH1F* mu8_E_Hist_Eff_Data_eta_Endcap;
    TH1F* mu8_E_Hist_Eff_Data_pT_Endcap;
    TH1F* mu8_E_Hist_Eff_Data_eta_Endcap_errUp;
    TH1F* mu8_E_Hist_Eff_Data_pT_Endcap_errUp;
    TH1F* mu8_E_Hist_Eff_Data_eta_Endcap_errDown;
    TH1F* mu8_E_Hist_Eff_Data_pT_Endcap_errDown;

    // Muon MC efficiency histograms
    // 2012 Triggers
    
    TH1F* mu8_AB_Hist_Eff_MC_eta_Barrel;
    TH1F* mu8_AB_Hist_Eff_MC_pT_Barrel;    
    TH1F* mu8_AB_Hist_Eff_MC_eta_Barrel_errUp;
    TH1F* mu8_AB_Hist_Eff_MC_pT_Barrel_errUp;
    TH1F* mu8_AB_Hist_Eff_MC_eta_Barrel_errDown;
    TH1F* mu8_AB_Hist_Eff_MC_pT_Barrel_errDown;

    TH1F* mu8_AB_Hist_Eff_MC_eta_Endcap;
    TH1F* mu8_AB_Hist_Eff_MC_pT_Endcap;
    TH1F* mu8_AB_Hist_Eff_MC_eta_Endcap_errUp;
    TH1F* mu8_AB_Hist_Eff_MC_pT_Endcap_errUp;    
    TH1F* mu8_AB_Hist_Eff_MC_eta_Endcap_errDown;
    TH1F* mu8_AB_Hist_Eff_MC_pT_Endcap_errDown;
   
    TH1F* mu8_C_Hist_Eff_MC_eta_Barrel;
    TH1F* mu8_C_Hist_Eff_MC_pT_Barrel;	 
    TH1F* mu8_C_Hist_Eff_MC_eta_Barrel_errUp;
    TH1F* mu8_C_Hist_Eff_MC_pT_Barrel_errUp;
    TH1F* mu8_C_Hist_Eff_MC_eta_Barrel_errDown;
    TH1F* mu8_C_Hist_Eff_MC_pT_Barrel_errDown;

    TH1F* mu8_C_Hist_Eff_MC_eta_Endcap;
    TH1F* mu8_C_Hist_Eff_MC_pT_Endcap;
    TH1F* mu8_C_Hist_Eff_MC_eta_Endcap_errUp;
    TH1F* mu8_C_Hist_Eff_MC_pT_Endcap_errUp;    
    TH1F* mu8_C_Hist_Eff_MC_eta_Endcap_errDown;
    TH1F* mu8_C_Hist_Eff_MC_pT_Endcap_errDown;
   
    TH1F* mu8_D_Hist_Eff_MC_eta_Barrel;
    TH1F* mu8_D_Hist_Eff_MC_pT_Barrel;	 
    TH1F* mu8_D_Hist_Eff_MC_eta_Barrel_errUp;
    TH1F* mu8_D_Hist_Eff_MC_pT_Barrel_errUp;
    TH1F* mu8_D_Hist_Eff_MC_eta_Barrel_errDown;
    TH1F* mu8_D_Hist_Eff_MC_pT_Barrel_errDown;

    TH1F* mu8_D_Hist_Eff_MC_eta_Endcap;
    TH1F* mu8_D_Hist_Eff_MC_pT_Endcap;
    TH1F* mu8_D_Hist_Eff_MC_eta_Endcap_errUp;
    TH1F* mu8_D_Hist_Eff_MC_pT_Endcap_errUp;    
    TH1F* mu8_D_Hist_Eff_MC_eta_Endcap_errDown;
    TH1F* mu8_D_Hist_Eff_MC_pT_Endcap_errDown;
    
    TH1F* mu8_E_Hist_Eff_MC_eta_Barrel;
    TH1F* mu8_E_Hist_Eff_MC_pT_Barrel;	 
    TH1F* mu8_E_Hist_Eff_MC_eta_Barrel_errUp;
    TH1F* mu8_E_Hist_Eff_MC_pT_Barrel_errUp;
    TH1F* mu8_E_Hist_Eff_MC_eta_Barrel_errDown;
    TH1F* mu8_E_Hist_Eff_MC_pT_Barrel_errDown;

    TH1F* mu8_E_Hist_Eff_MC_eta_Endcap;
    TH1F* mu8_E_Hist_Eff_MC_pT_Endcap;
    TH1F* mu8_E_Hist_Eff_MC_eta_Endcap_errUp;
    TH1F* mu8_E_Hist_Eff_MC_pT_Endcap_errUp;    
    TH1F* mu8_E_Hist_Eff_MC_eta_Endcap_errDown;
    TH1F* mu8_E_Hist_Eff_MC_pT_Endcap_errDown;


    
    Bool_t doDebug ;

} ;

#endif //H_HSG4LEPLEPTRIGGERSF
