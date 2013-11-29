#ifndef JVFUncertaintyTool_h
#define JVFUncertaintyTool_h 1

/************************************************************************
 * JVFUncertaintyTool.h                                              	*
 *                                                                      *
 * Usage description							*
 
The tool provides the JVF cut up/down variation needed to estimate the 
systematic uncertainty associated to the use of a JVF cut in the 
analysis. Three JVF working points are supported: 0, 0.25, 0.50. 
 
The JVF cut, as well as the up/down systematic variations have to be 
applied only to jets with pT<50 GeV and |eta|<2.4. The output depends 
on whether the jet falls in the PU or HS jet category. Two functions 
are defined: the first one isPileUpJet() classifies the jet as a HS 
or PU jet based in a DeltaR matching between the jet and the truth 
jets collection.

bool isPileUpJet(const TLorentzVector jet, const std::vector<TLorentzVector> &truthJets);
Inputs:
- jet = reconstructed jet after LC/EM+JES calibration, pT<50 GeV and |eta_detector|<2.4
- truthJets = truth jets with pT>10 GeV

The second functions getJVFcut() 
provides the JVF up/down cut variation. These variations are calculated 
separately for HS and PU jets and for the different jet algorithms:
"AntiKt4TopoEM","AntiKt6TopoEM","AntiKt4LCTopo" and "AntiKt6LCTopo". The jet 
algorithm is passed as an argument when initializing the tool.

float getJVFcut(const float JVFcutNominal, const bool isPileUp, const float pT_jet, const float eta_det, const bool isUp );
Inputs:
- pT_jet = reconstructed jet pT after LC/EM+JES calibration (only jets with pT<50 GeV have to be considered)
- eta_det = jet eta detector (only jets with |eta_det|<2.4 GeV have to be considered)
- JVFcutNominal = 0.0, 0.25, 0.50
- isPileUp = output of the previous function
- isUp = true to get the JVF up cut variation, false to get the JVF down cut variation

TLorentzVectors are expected to have units of MeV by default. If you want to use GeV instead of GeV please call the function useGeV(true).

####################################
Especial case: JVF nominal = 0
####################################
Only the JVF up cut variation is provided in this case for obvious reasons. The 
difference between using the JVF up cut and the nominal one has to be symmetrized 
in order to get the down variation.

####################################
How to compile:

The code is Standalone/ATHENA/RootCore compatible

To compile it Standalone:
cd Reconstruction/Jet/JetAnalysisTools/JVFUncertaintyTool/cmt
cmt make -f Makefile.Standalone
It is then possible to link the compiled .so to your standalone executable.

####################################

####################################
Example of use: Getting the JVF cut up/down variation, asumming that 
the nominal JVF cut of the analysis is |JVF|>0.5 and use Antikt R=0.4 
jets
####################################

#include "JVFUncertaintyTool/JVFUncertaintyTool.h"

//Create an instance of JVFUncertaintyTool. It takes as argument the jet collection: "AntiKt4TopoEM","AntiKt6TopoEM","AntiKt4LCTopo","AntiKt6LCTopo"

JVFUncertaintyTool *jvfTool = new JVFUncertaintyTool("AntiKt4TopoEM");

float JVFcutNominal = 0.5;

//Create the TLorentzVector of the reconstructed jet after the EM/LC+JES calibration
TLorentzVector jet;
jet.SetPtEtaPhiE(pT, eta, phi, E)

//Verify is the jet is classified as a PU or HS jet
bool isPU = jvf->isPileUpJet(jet, truthJets);

//Get the JVF cut up 
float JVFcutUp = jvf->getJVFcut(JVFcutNominal, isPU, jet.Pt(), eta_det, true);

//Get the JVF cut down
float JVFcutDown = jvf->getJVFcut(JVFcutNominal, isPU, jet.Pt(), eta_det, false);

 
 * 									*
 * History                                                              *
 *         30 January 2013 -- created by R. Camacho/ J. Backus Mayes    *    
 ************************************************************************/

//Headers

//_________________________________________________________________________________________________
// C-C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

//_________________________________________________________________________________________________
// ROOT specific header files
#include <TLorentzVector.h>
#include "TMath.h"
#include "TString.h"
#include "TNamed.h"
#include "TObject.h"

//Declaration of functions

class JVFUncertaintyTool : public TNamed{

public:
	
  JVFUncertaintyTool(TString jetAlgo="AntiKt4TopoEM");
  virtual ~JVFUncertaintyTool();
  bool isPileUpJet(const TLorentzVector jet, const std::vector<TLorentzVector> &truthJets);
  float getJVFcut (const float JVFcutNominal, const bool isPileUp, const float pT_jet, const float eta_det, const bool isUp);

  void UseGeV(bool useGeV=true) { if (useGeV) _GeV=1; else _GeV=1000; }

  float JVFup_HS[3];
  float JVFdw_HS[3];
	
  float drcut;
	
private:
  double _GeV;


ClassDef(JVFUncertaintyTool,1)

};



#endif
