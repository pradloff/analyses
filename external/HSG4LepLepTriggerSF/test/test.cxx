#include "../HSG4LepLepTriggerSF/HSG4LepLepTriggerSF.h"
#include "TH1.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

int main( int argc , char** argv ){

    UInt_t runnum = 200804 ;
    TString period = "A";
    TString setup = "pt";

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
    if ( argc > 2 ) {
        setup = argv[ 2 ] ;
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

    TH1F* plot_eloose = new TH1F("Plot_Eloose", "Plot_Eloose", 1000, -5, 5);
    TH1F* plot_emedium = new TH1F("Plot_Emedium", "Plot_Emedium", 1000, -5, 5);
    TH1F* plot_mu = new TH1F("Plot_Mu", "Plot_Mu", 1000, -5, 5);

    Double_t pt_min;
    UInt_t pt_nstep;
    Double_t pt_step; 

    Double_t eta_min; 
    UInt_t eta_nstep;
    Double_t eta_step;
    
    if(setup == "pt"){
	pt_min = 0 * 1000 ; 
	pt_nstep = 205 ; 
	pt_step = 1 * 500 ;

	eta_min = 0 ; 
	eta_nstep = 0 ; 
	eta_step = 0 ;
	//Double_t eta_min = 1 ; 
	//UInt_t eta_nstep = 0 ; 
	//Double_t eta_step = 0 ;
    }
    else if(setup == "eta"){    
	pt_min = 40 * 1000 ; 
	pt_nstep = 0 ; 
	pt_step = 0 * 1000 ;
	//Double_t pt_min = 60 * 1000 ; 
	//UInt_t pt_nstep = 0 ; 
	//Double_t pt_step = 0 * 1000 ;

	eta_min = -3.0 ; 
	eta_nstep = 600 ; 
        eta_step = 0.01 ;
    }
    else{
	cout << "Which setup?" << endl;
    }

    Double_t phi_min = 0 ; 
    UInt_t phi_nstep = 0 ; 
    Double_t phi_step = 0 ;

    TLorentzVector lep ;

    std::cout << std::endl ;
    std::cout << Form( "| *Muon (mu8) Trigger Sf & Eff Table ( Run Number = %6d )* | " , runnum ) << std::endl ;
    std::cout << std::endl ;
    std::cout << Form( "|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|" , 
                       "*pT[GeV]*" , "*eta*" , "*phi*" , 
                       "*Eff(m)*" , "*EffUp(m)*" , "*EffDn(m)*" , "*EffMC(m)*" , "*EffUpMC(m)*" , "*EffDnMC(m)*" , "*Sf(m)*" , "*SfUp(m)*" , "*SfDn(m)*" ) << std::endl ;
  
    for ( UInt_t i = 0 ; i <= pt_nstep ; i++ ) {
        for ( UInt_t j = 0 ; j <= eta_nstep ; j++ ) {
            for ( UInt_t k = 0 ; k <= phi_nstep ; k++ ) {

                Double_t pt = pt_min + pt_step * i ;
                Double_t eta = eta_min + eta_step * j ;
                Double_t phi = phi_min + phi_step * k ;
                
                lep.SetPtEtaPhiM( pt , eta , phi , 0 ) ;
                                
                Double_t muo_eff    = SFTool->getDataEffMuon( lep , runnum , "mu8" , 0 ) ;
                Double_t muo_eff_up = SFTool->getDataEffMuon( lep , runnum , "mu8" ,+1 ) ;
                Double_t muo_eff_dn = SFTool->getDataEffMuon( lep , runnum , "mu8" ,-1 ) ;      
                Double_t muo_mc_eff    = SFTool->getMCEffMuon( lep , runnum , "mu8" , 0 ) ;
                Double_t muo_mc_eff_up = SFTool->getMCEffMuon( lep , runnum , "mu8" ,+1 ) ;
                Double_t muo_mc_eff_dn = SFTool->getMCEffMuon( lep , runnum , "mu8" ,-1 ) ;      
                Double_t muo_fac    = SFTool->getSFMuon( lep , runnum , "mu8" , 0 ) ;
                Double_t muo_fac_up = SFTool->getSFMuon( lep , runnum , "mu8" ,+1 ) ;
                Double_t muo_fac_dn = SFTool->getSFMuon( lep , runnum , "mu8" ,-1 ) ;
 
		plot_mu->Fill(muo_fac);
               
                std::cout << Form( "|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|" , 
                           pt / 1000 , eta , phi , 
                                   muo_eff , muo_eff_up , muo_eff_dn , muo_mc_eff , muo_mc_eff_up , muo_mc_eff_dn , muo_fac , muo_fac_up , muo_fac_dn ) << std::endl ;
                
            }
	}
    }
  
    std::cout << std::endl ;
    std::cout << Form( "| *Electron (e12Tvhl1) Trigger Sf & Eff Table ( Run Number = %6d )* | " , runnum ) << std::endl ;
    std::cout << std::endl ;
    std::cout << Form( "|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|" , 
                       "*pT[GeV]*" , "*eta*" , "*phi*" , 
                       "*Eff(e)*" , "*EffUp(e)*" , "*EffDn(e)*" , "*EffMC(e)*" , "*EffUpMC(e)*" , "*EffDnMC(e)*" , "*Sf(e)*" , "*SfUp(e)*" , "*SfDn(e)*" ) << std::endl ;
  
    for ( UInt_t i = 0 ; i <= pt_nstep ; i++ ) {
        for ( UInt_t j = 0 ; j <= eta_nstep ; j++ ) {
            for ( UInt_t k = 0 ; k <= phi_nstep ; k++ ) {

                Double_t pt = pt_min + pt_step * i ;
                Double_t eta = eta_min + eta_step * j ;
                Double_t phi = phi_min + phi_step * k ;
                
                lep.SetPtEtaPhiM( pt , eta , phi , 0 ) ;
                                                
                Double_t ele_eff    = SFTool->getDataEffElec( lep , runnum , "e12Tvhl1" , 0 ) ;
                Double_t ele_eff_up = SFTool->getDataEffElec( lep , runnum , "e12Tvhl1" ,+1 ) ;
                Double_t ele_eff_dn = SFTool->getDataEffElec( lep , runnum , "e12Tvhl1" ,-1 ) ;      
                Double_t ele_mc_eff    = SFTool->getMCEffElec( lep , runnum , "e12Tvhl1" , 0 ) ;
                Double_t ele_mc_eff_up = SFTool->getMCEffElec( lep , runnum , "e12Tvhl1" ,+1 ) ;
                Double_t ele_mc_eff_dn = SFTool->getMCEffElec( lep , runnum , "e12Tvhl1" ,-1 ) ;      
                Double_t ele_fac    = SFTool->getSFElec( lep , runnum , "e12Tvhl1" , 0 ) ;
                Double_t ele_fac_up = SFTool->getSFElec( lep , runnum , "e12Tvhl1" ,+1 ) ;
                Double_t ele_fac_dn = SFTool->getSFElec( lep , runnum , "e12Tvhl1" ,-1 ) ;
 
		plot_emedium->Fill(ele_fac);
               
                std::cout << Form( "|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|" , 
                           pt / 1000 , eta , phi , 
                                   ele_eff , ele_eff_up , ele_eff_dn , ele_mc_eff , ele_mc_eff_up , ele_mc_eff_dn , ele_fac , ele_fac_up , ele_fac_dn ) << std::endl ;
                
            }
	}
    }

    std::cout << std::endl ;
    std::cout << Form( "| *Electron (e12Tvhm1) Trigger Sf & Eff Table ( Run Number = %6d )* | " , runnum ) << std::endl ;
    std::cout << std::endl ;
    std::cout << Form( "|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|%10s|" , 
                       "*pT[GeV]*" , "*eta*" , "*phi*" , 
                       "*Eff(e)*" , "*EffUp(e)*" , "*EffDn(e)*" , "*EffMC(e)*" , "*EffUpMC(e)*" , "*EffDnMC(e)*" , "*Sf(e)*" , "*SfUp(e)*" , "*SfDn(e)*" ) << std::endl ;
  
    for ( UInt_t i = 0 ; i <= pt_nstep ; i++ ) {
        for ( UInt_t j = 0 ; j <= eta_nstep ; j++ ) {
            for ( UInt_t k = 0 ; k <= phi_nstep ; k++ ) {

                Double_t pt = pt_min + pt_step * i ;
                Double_t eta = eta_min + eta_step * j ;
                Double_t phi = phi_min + phi_step * k ;
                
                lep.SetPtEtaPhiM( pt , eta , phi , 0 ) ;
                                                
                Double_t ele_eff    = SFTool->getDataEffElec( lep , runnum , "e12Tvhm1" , 0 ) ;
                Double_t ele_eff_up = SFTool->getDataEffElec( lep , runnum , "e12Tvhm1" ,+1 ) ;
                Double_t ele_eff_dn = SFTool->getDataEffElec( lep , runnum , "e12Tvhm1" ,-1 ) ;      
                Double_t ele_mc_eff    = SFTool->getMCEffElec( lep , runnum , "e12Tvhm1" , 0 ) ;
                Double_t ele_mc_eff_up = SFTool->getMCEffElec( lep , runnum , "e12Tvhm1" ,+1 ) ;
                Double_t ele_mc_eff_dn = SFTool->getMCEffElec( lep , runnum , "e12Tvhm1" ,-1 ) ;      
                Double_t ele_fac    = SFTool->getSFElec( lep , runnum , "e12Tvhm1" , 0 ) ;
                Double_t ele_fac_up = SFTool->getSFElec( lep , runnum , "e12Tvhm1" ,+1 ) ;
                Double_t ele_fac_dn = SFTool->getSFElec( lep , runnum , "e12Tvhm1" ,-1 ) ;

		plot_eloose->Fill(ele_fac);
                
                std::cout << Form( "|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|%10.7f|" , 
                           pt / 1000 , eta , phi , 
                                   ele_eff , ele_eff_up , ele_eff_dn , ele_mc_eff , ele_mc_eff_up , ele_mc_eff_dn , ele_fac , ele_fac_up , ele_fac_dn ) << std::endl ;
                
            }
	}
    }
  
    std::cout << std::endl ;  

    TCanvas *canv = new TCanvas();
    plot_eloose->Draw();
    canv->SaveAs("eLoose_SF_"+setup+"_"+period+".pdf");
    plot_emedium->Draw();
    canv->SaveAs("eMedium_SF_"+setup+"_"+period+".pdf");
    plot_mu->Draw();
    canv->SaveAs("mu_SF_"+setup+"_"+period+".pdf");
    return 0;
}

