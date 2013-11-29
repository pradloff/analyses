#include "../HSG4LepLepTriggerSF/HSG4LepLepTriggerSF.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

int main( int argc , char** argv ){

    UInt_t runnum = 200804 ;
    TString period = "A";

    if ( argc > 1 ) {
    	//TString strRunnum = argv[ 1 ] ;
    	//runnum = strRunnum.Atoi() ;
	period = argv[ 1 ];
	if(period == "A") runnum = 200804 ; // start of period A
	else if(period == "B") runnum = 202660 ; // start of period B
	else if(period == "C") runnum = 206248 ; // start of period C
	else if(period == "D") runnum = 207447 ; // start of period D
	else if(period == "E") runnum = 209074 ; // start of period E
	else if(period == "G") runnum = 211522 ; // start of period G
 	else if(period == "H") runnum = 212619 ; // start of period H
	else if(period == "I") runnum = 213431 ; // start of period I
	else if(period == "J") runnum = 213900 ; // start of period J
	else if(period == "L") runnum = 215414 ; // start of period L
	else 
	  return 0;
    }
    
    // **************************************************************	
    // Step 1. Instantiate tool
    // **************************************************************	
    HSG4LepLepTriggerSF* SFTool ;

    // **************************************************************	
    // Step 2. Config the tool
    // HSG4LepLepTriggerSF::HSG4LepLepTriggerSF( const std::string& path , const Bool_t debug )
    // -- "path" is where SF root files are
    // -- "debug" is whether to use debug mode or not
    // **************************************************************	
    SFTool = new HSG4LepLepTriggerSF( "../data/" , kFALSE ) ;
    //SFTool = new HSG4LepLepTriggerSF( "../data/" , kTRUE ) ;

    // **************************************************************	
    // Step 3 Retrive Scale Factors and Effeciency
    // **************************************************************	
    // Available methods:
    //
    // -- get Scale Factor
    //
    // Double_t getSFElec( const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const ;
    //
    // -- get Efficiency (in data or in mc)
    //
    // Double_t getDataEffElec( const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const ;
    //
    // Double_t getMCEffElec( const TLorentzVector& lepVec , const UInt_t runNumber , const TString trigger, const Int_t mode ) const ;
    // 
    // Needed:
    // -- "lepVec" is the TLorentzVector of the offline lepton (energy in MeV)
    // -- "runNumber" (not doing anything at the moment, put 0)
    // -- "trigger" is the trigger chain for which the SF is needed (supported at the moment: mu8, e12Tvhm1 and e12Tvhl1)
    // -- "mode" is 0 for SF, +1 and -1 for systematics
    // 
    // **************************************************************	

    Double_t pt_min;
    UInt_t pt_nstep;
    Double_t pt_step; 

    Double_t eta_min; 
    UInt_t eta_nstep;
    Double_t eta_step;
    
    pt_min = 0. ; 
    pt_nstep = 205. ; 
    pt_step = 500. ;
    
    eta_min = -3.0 ; 
    eta_nstep = 600 ; 
    eta_step = 0.01 ;
	
    Double_t phi_min = 0 ; 
    UInt_t phi_nstep = 0 ; 
    Double_t phi_step = 0 ;

    TH2D* plot_eloose = new TH2D("Plot_Eloose", "Plot_Eloose", eta_nstep, eta_min, eta_min+(eta_step*eta_nstep), pt_nstep, pt_min, pt_min+(pt_step*pt_nstep));
    TH2D* plot_emedium = new TH2D("Plot_Emedium", "Plot_Emedium", eta_nstep, eta_min, eta_min+(eta_step*eta_nstep), pt_nstep, pt_min, pt_min+(pt_step*pt_nstep));
    TH2D* plot_mu = new TH2D("Plot_Mu", "Plot_Mu", eta_nstep, eta_min, eta_min+(eta_step*eta_nstep), pt_nstep, pt_min, pt_min+(pt_step*pt_nstep));

    TLorentzVector lep ;

    for ( UInt_t i = 0 ; i <= pt_nstep ; i++ ) {
        for ( UInt_t j = 0 ; j <= eta_nstep ; j++ ) {
	  //            for ( UInt_t k = 0 ; k <= phi_nstep ; k++ ) {

                Double_t pt = pt_min + pt_step * i ;
                Double_t eta = eta_min + eta_step * j ;
		Double_t phi = 0 ;
                //Double_t phi = phi_min + phi_step * k ;
                
                lep.SetPtEtaPhiM( pt , eta , phi , 0 ) ;
                                                
                Double_t ele_loose_fac    = SFTool->getSFElec( lep , runnum , "e12Tvhl1" , 0 ) ;
                Double_t ele_medium_fac    = SFTool->getSFElec( lep , runnum , "e12Tvhm1" , 0 ) ;
		Double_t mu_fac    = SFTool->getSFMuon( lep , runnum , "mu8" , 0 ) ;

		plot_eloose->SetBinContent(j,i,ele_loose_fac);
		plot_emedium->SetBinContent(j,i,ele_medium_fac);
		plot_mu->SetBinContent(j,i,mu_fac);
		//  }
	}
    }
  
    std::cout << std::endl ;  

    TCanvas *canv = new TCanvas();
    plot_eloose->Draw("colz");
    canv->SaveAs("2D_eLoose_SF_"+period+".pdf");
    plot_emedium->Draw("colz");
    canv->SaveAs("2D_eMedium_SF_"+period+".pdf");
    plot_mu->Draw("colz");
    canv->SaveAs("2D_mu_SF_"+period+".pdf");
    return 0;
}

