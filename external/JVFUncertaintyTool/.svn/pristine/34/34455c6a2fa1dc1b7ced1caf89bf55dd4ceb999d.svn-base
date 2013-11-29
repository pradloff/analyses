#include "JVFUncertaintyTool/JVFUncertaintyTool.h"


JVFUncertaintyTool::JVFUncertaintyTool(TString jetAlgo):_GeV(1000)
{
	
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Three JVF working points: 0.00, 0.25, 0.50
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (jetAlgo=="AntiKt4TopoEM"){ 
    JVFup_HS[0] = 0.10; JVFup_HS[1] = 0.28; JVFup_HS[2] = 0.53;		
    JVFdw_HS[0] = 0.00; JVFdw_HS[1] = 0.21; JVFdw_HS[2] = 0.47;
  } else if (jetAlgo=="AntiKt6TopoEM"){ 
    JVFup_HS[0] = 0.05; JVFup_HS[1] = 0.27; JVFup_HS[2] = 0.52;		
    JVFdw_HS[0] = 0.00; JVFdw_HS[1] = 0.24; JVFdw_HS[2] = 0.48;
  } else if (jetAlgo=="AntiKt4LCTopo"){ 
    JVFup_HS[0] = 0.10; JVFup_HS[1] = 0.28; JVFup_HS[2] = 0.53;		
    JVFdw_HS[0] = 0.00; JVFdw_HS[1] = 0.22; JVFdw_HS[2] = 0.47;
  } else if (jetAlgo=="AntiKt6LCTopo"){ 
    JVFup_HS[0] = 0.05; JVFup_HS[1] = 0.27; JVFup_HS[2] = 0.52;		
    JVFdw_HS[0] = 0.00; JVFdw_HS[1] = 0.24; JVFdw_HS[2] = 0.48;
  } else std::cout<<"Don't know about jet collection "<<jetAlgo<<std::endl;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Truth-reco jet matching DeltaR cut depends on jet algorithm
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (jetAlgo=="AntiKt4TopoEM" || jetAlgo=="AntiKt4LCTopo") drcut = 0.4;
  else if (jetAlgo=="AntiKt6TopoEM" || jetAlgo=="AntiKt6LCTopo") drcut = 0.6;
  else std::cout<<"Don't know about jet collection "<<jetAlgo<<std::endl;

}

JVFUncertaintyTool::~JVFUncertaintyTool()
{
  //Nothing needed here
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Is the jet matched to a truth jet or not?
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//_________________________________________________________________________________________________

bool JVFUncertaintyTool::isPileUpJet(const TLorentzVector jet, const std::vector<TLorentzVector> &truthJets){

  // truthJets must be ordered by decreasing pt
  for (unsigned int t=0; t < truthJets.size(); t++){
    float truthPt = truthJets.at(t).Pt()/_GeV;
    if(truthPt < 10.) break;
    if((jet.Pt()/_GeV)/truthPt > 3.) break;

    float dr = jet.DeltaR(truthJets.at(t));
    if(dr < drcut) return false;    	
  }
	
  ///////////////////////////////////////////////////////////////
  // no truth jet was found with pt > max(10 GeV, jet.Pt()/3.)
  // and within dR < 0.4 --> jet is labeled as pile-up
  ///////////////////////////////////////////////////////////////
  return true;	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Get JVF cut up/ down
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//_________________________________________________________________________________________________

float JVFUncertaintyTool::getJVFcut (const float JVFcutNominal, const bool isPileUp, const float pT_jet, const float eta_det, const bool isUp){

  double JVFcut = -2.;
  
  if (pT_jet/_GeV>50.){
    std::cout << "WARNING: A JVF cut should not be applied to jets with pT>50 GeV. No JVF cut will be applied!" << std::endl;
    return JVFcut;
  }	

  if (fabs(eta_det)>2.4){
    std::cout << "WARNING: A JVF cut should not be applied to jets with |eta|>2.4. No JVF cut will be applied!" << std::endl;
    return JVFcut;
  } 

  int workingPoint = -1;
  if (fabs(JVFcutNominal - 0.) < 0.001){
    if (isUp){ 
      workingPoint = 0;
    }else{
      std::cout << "WARNING: No JVF down cut for the working point JVF=0.0. Use isUp=true and symmetrize. No JVF cut will be applied!" << std::endl;
      return JVFcut;
    } 
  } else if (fabs(JVFcutNominal - 0.25) < 0.001){
    workingPoint = 1;
  } else if (fabs(JVFcutNominal - 0.50) < 0.001){
    workingPoint = 2;
  } else { 
    std::cout << "WARNING: This JVF working point is not supported! No JVF cut will be applied" << std::endl;
    return JVFcut;
  }	

  if (isPileUp){

    if (isUp)
      JVFcut = 1.1;//Removing all pile-up jets
    else	
      JVFcut = 0.1;
			
  } else {

    if (isUp)
      JVFcut = JVFup_HS[workingPoint];
    else
      JVFcut = JVFdw_HS[workingPoint];	

  }
  
  return JVFcut;

}

