#include "./METUtilityAthTool.h"

#include "JetUncertainties/MultijetJESUncertaintyProvider.h"
#include "JetResolution/JERProviderAthena.h"
#include "MuonMomentumCorrections/SmearingClass.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "MissingETUtility/TESUncertaintyProvider.h"
#include <TRandom2.h>
#include <TMath.h>

#include "AthenaKernel/errorcheck.h"
#include "FourMomUtils/P4Helpers.h"
#include "EventInfo/EventID.h"
#include "EventInfo/EventType.h"
#include "EventInfo/EventInfo.h"
#include "MissingETEvent/MissingETComposition.h"
#include "JetEvent/JetCollection.h"
#include "CaloEvent/CaloClusterContainer.h"
#include "egammaEvent/ElectronContainer.h"
#include "egammaEvent/PhotonContainer.h"
#include "muonEvent/MuonContainer.h"
#include "tauEvent/TauJetContainer.h"
#include "VxVertex/VxContainer.h"
#include "EventKernel/SignalStateHelper.h"

using namespace MissingETUtility;

// Constructor
METUtilityAthTool::METUtilityAthTool(const string& type, const string& name,
				     const IInterface* pParent) :
  AthAlgTool(type, name, pParent),
  METUtility(true, true, true, true, true, true,
	     true, true, true, true, 20000., false),
  m_lbmuTool( "LumiBlockMuTool" ),
  m_isData(false),
  m_nvtxjet(0)
{

  declareInterface< IMETUtilityAthTool >(this);

  declareProperty("DoEveryEvent",      m_doEveryEvent      =  true                         );
  declareProperty("LumiBlockMuTool",   m_lbmuTool,"Tool used for retrieving the mu values" );

  declareProperty("Verbosity",         m_verbosity         =  false                        );
  declareProperty("DoRefEle",          m_dorefele          =  true                         );
  declareProperty("DoRefGamma",        m_dorefgamma        =  true                         );
  declareProperty("DoRefTau",          m_doreftau          =  true                         );
  declareProperty("DoRefJet",          m_dorefjet          =  true                         );
  declareProperty("DoSoftJets",        m_dosoftjets        =  true                         );
  declareProperty("DoRefMuon",         m_dorefmuon         =  true                         );
  declareProperty("DoMuonTotal",       m_domuontotal       =  true                         );
  declareProperty("DoCellOut",         m_docellout         =  false                        );
  declareProperty("DoCellOutEflow",    m_docellouteflow    =  true                         );
  declareProperty("IsMuid",            m_ismuid            =  false                        );
  declareProperty("SoftJetCut",        m_softjetcut        =  20000.0                      ); 

  declareProperty("DoCloseByJES",      m_doCloseByJES      =  true                         );
  declareProperty("PtCloseByJES",      m_ptCloseByJES      =  7e3                          );
  declareProperty("JetTypeJES",        m_jetTypeJES        =  "AntiKt4LCTopoJets"          );
  declareProperty("MCTypeJES",         m_mcTypeJES         =  "MC11c"                      );
  declareProperty("ConfigFileJES",     m_configfileJES     =  "InsituJES2011_AllNuisanceParameters.config" );
  declareProperty("ConfigMultijet",    m_configMultijet    =  "MultijetJES_Preliminary.config"         );

  declareProperty("JetTypeJER",        m_jetTypeJER        =  "AntiKt6LCTopoJES"           );
  declareProperty("MethodJER",         m_methodJER         =  "Truth"                      );
  declareProperty("InputFileJER",      m_inputfileJER      =  "JERProviderPlots.root"      );

  declareProperty("DoSignif",          m_doSignif          =  false                        );

  // Best if these don't get used.
  // Can't declareProperty the enum, so translate the int instead.
  declareProperty("MuonStatusWord",    m_muonStatusWord    =  0x0001                       );
  declareProperty("ClusterStatusWord", m_clusterStatusWord =  0x0001                       );

  declareProperty("UseSmearedEl",      m_useSmearedEl      =  false                        );
  declareProperty("UseSmearedPh",      m_useSmearedPh      =  false                        );
  declareProperty("UseSmearedMu",      m_useSmearedMu      =  false                        );

  declareProperty("Is2012",            m_flagis2012        =  false                        );
  declareProperty("IsSTVF",            m_flagisSTVF        =  false                        );

//   declareProperty("UseMETComp",        m_useMetComp        =  true                         );
  declareProperty("METCompKey",        m_metCompKey        =  "MET_RefComposition"         );

  declareProperty("ClusterKey",        m_clusterKey        =  "METComposition"             );
  /*"CaloCalTopoCluster"*/
  declareProperty("TrackKey",          m_trackKey          =  "METComposition"             );
  /*"TrackParticleCandidate"*/
  declareProperty("ElectronKey",       m_electronKey       =  "METComposition"             );
  /*"ElectronAODCollection"*/
  declareProperty("JetKey",            m_jetKey            =  "METComposition"             );
  /*"AntiKt4LCTopoJets"*/
  declareProperty("MuonKey",           m_muonKey           =  "METComposition"             );
  /*"StacoMuonCollection"*/
  declareProperty("CaloMuonKey",       m_calomuonKey       =  ""             );
  /*"CaloMuonCollection"*/
  declareProperty("PhotonKey",         m_photonKey         =  "METComposition"             );
  /*"PhotonAODCollection"*/
  declareProperty("TaujetKey",         m_taujetKey         =  "METComposition"             );
  /*"TauRecContainer"*/
  declareProperty("VertexKey",         m_vertexKey         =  "VxPrimaryCandidate"         );
  /*"CaloCalTopoCluster"*/

  declareProperty("RefEleKey",         m_RefEleKey         =  ""                           );
  declareProperty("RefGammaKey",       m_RefGammaKey       =  ""                           );
  declareProperty("RefJetKey",         m_RefJetKey         =  ""                           );
  declareProperty("RefTauKey",         m_RefTauKey         =  ""                           );
  declareProperty("RefMuonKey",        m_RefMuonKey        =  "MET_RefMuon"                );
  declareProperty("SoftJetsKey",       m_SoftJetsKey       =  ""                           );
  declareProperty("CellOutKey",        m_CellOutKey        =  ""                           );
  declareProperty("CellOutEflowKey",   m_CellOutEflowKey   =  ""                           );
  declareProperty("TruthNonIntKey",    m_TruthNonIntKey    =  "MET_Truth_NonInt"           );
}

// Destructor
METUtilityAthTool::~METUtilityAthTool()
{

}

StatusCode METUtilityAthTool::queryInterface( const InterfaceID& riid, void** ppvIf )
{
   if ( riid == IMETUtilityAthTool::interfaceID() )  {
      *ppvIf = (IMETUtilityAthTool*)this;
      addRef();
      return StatusCode::SUCCESS;
   }

   return AthAlgTool::queryInterface( riid, ppvIf );
}


StatusCode METUtilityAthTool::initialize() {
  msg(MSG::INFO) << "Initialising MissingETUtility" << endreq;

  m_multiJES = new MultijetJESUncertaintyProvider(m_configMultijet, m_configfileJES, m_jetTypeJES, m_mcTypeJES);

  m_JERTool = new JetResolution::JERProviderAthena(m_jetTypeJER, m_methodJER, m_inputfileJER);
  m_JERTool->init();

  m_egammaTool = new eg2011::EnergyRescaler;
  m_egammaTool->useDefaultCalibConstants("2011");

  m_muonTool = new MuonSmear::SmearingClass("Data11","staco","q_pT","Rel17");
  m_muonTool->UseScale(1);

  m_TESTool = new TESUncertaintyProvider;

  msg(MSG::INFO) << "Setting up MET systematics "
		  << endreq;
  msg(MSG::INFO) << "Using " << (m_ismuid ? "MUID" : "non-MUID") << " muons"
		  << endreq;

  // Set up the MET details
  doSignificance(m_doSignif);
  setVerbosity(m_verbosity);
  setIsMuid(m_ismuid);
  setSoftJetCut(m_softjetcut);
  defineMissingET(m_dorefele, m_dorefgamma, m_doreftau, m_dorefjet, m_dosoftjets,
		  m_dorefmuon, m_domuontotal, m_docellout, m_docellouteflow);
  configMissingET(m_flagis2012,m_flagisSTVF);

  IIncidentSvc* incsvc;
  StatusCode sc = service("IncidentSvc", incsvc);
  int priority = 100;
  if( sc.isSuccess() ) {
     incsvc->addListener( this, "BeginEvent", priority);
  }

  CHECK( m_lbmuTool.retrieve() );

  return StatusCode::SUCCESS;
}

// not an algorithm, doesn't have an execute
// but we may want to run process every event
StatusCode METUtilityAthTool::process() {
  METUtility::reset();

  // This gets the EventInfo object from StoreGate
  const EventInfo* myEventInfo = 0;
  // Key "MyEvent" is optional, usually not specified for EventInfo because
  // there'll be only one. When not specified, just takes the first container.
  CHECK( evtStore()->retrieve(myEventInfo/*,"MyEvent"*/) );

  if(myEventInfo->event_type()->test(EventType::IS_SIMULATION)) m_isData = false;

  CHECK( evtStore()->retrieve(m_metComp, m_metCompKey) );

  // set the MET terms directly
  const MissingET* MET_RefEle = 0;
  if(m_RefEleKey != "") {
    CHECK( evtStore()->retrieve(MET_RefEle, m_RefEleKey) );
    setMETTerm(METUtil::RefEle, MET_RefEle);
  } else {
    msg(MSG::DEBUG) << "Skip RefEle" << endreq;
  }

  const MissingET* MET_RefGamma = 0;
  if(m_RefGammaKey != "") {
    CHECK( evtStore()->retrieve(MET_RefGamma, m_RefGammaKey) );
    setMETTerm(METUtil::RefGamma, MET_RefGamma);
  } else {
    msg(MSG::DEBUG) << "Skip RefGamma" << endreq;
  }

  const MissingET* MET_RefJet = 0;
  if(m_RefJetKey != "") {
    CHECK( evtStore()->retrieve(MET_RefJet, m_RefJetKey) );
    setMETTerm(METUtil::RefJet, MET_RefJet);
  } else {
    msg(MSG::DEBUG) << "Skip RefJet" << endreq;
  }

  const MissingET* MET_RefTau = 0;
  if(m_RefTauKey != "") {
    CHECK( evtStore()->retrieve(MET_RefTau, m_RefTauKey) );
    setMETTerm(METUtil::RefTau, MET_RefTau);
  } else {
    msg(MSG::DEBUG) << "Skip RefTau" << endreq;
  }

  const MissingET* MET_RefMuon = 0;
  if(m_RefMuonKey != "") {
    CHECK( evtStore()->retrieve(MET_RefMuon, m_RefMuonKey) );
    setMETTerm(METUtil::RefMuon, MET_RefMuon);
  } else {
    msg(MSG::DEBUG) << "Skip RefMuon" << endreq;
  }

  const MissingET* MET_SoftJets = 0;
  if(m_SoftJetsKey != "") {
    CHECK( evtStore()->retrieve(MET_SoftJets, m_SoftJetsKey) );
    setMETTerm(METUtil::SoftJets, MET_SoftJets);
  } else {
    msg(MSG::DEBUG) << "Skip SoftJets" << endreq;
  }

  const MissingET* MET_CellOut = 0;
  if(m_CellOutKey != "") {
    CHECK( evtStore()->retrieve(MET_CellOut, m_CellOutKey) );
    setMETTerm(METUtil::CellOut, MET_CellOut);
  } else {
    msg(MSG::DEBUG) << "Skip CellOut" << endreq;
  }

  msg(MSG::DEBUG) << "CellOut_Eflow key set to " << m_CellOutEflowKey << endreq;
  const MissingET* MET_CellOut_Eflow = 0;
  if(m_CellOutEflowKey != "") {
    CHECK( evtStore()->retrieve(MET_CellOut_Eflow, m_CellOutEflowKey) );
    setMETTerm(METUtil::CellOutEflow, MET_CellOut_Eflow);
  } else {
    msg(MSG::DEBUG) << "Skip CellOut_Eflow" << endreq;
  }

  if(m_isData) {
    msg(MSG::DEBUG) << "Truth_NonInt key set to " << m_TruthNonIntKey << endreq;
    const MissingET* MET_Truth_NonInt = 0;
    if(m_TruthNonIntKey != "") {
      CHECK( evtStore()->retrieve(MET_Truth_NonInt, m_TruthNonIntKey) );
      setMETTerm(METUtil::Truth, MET_Truth_NonInt);
    } else {
      msg(MSG::DEBUG) << "Skip Truth_NonInt" << endreq;
    }
  }

  // This is needed for jet pileup uncertainties
  // as well as for SoftTerms_ptHard uncertainties
  if(m_vertexKey != "") {
    CHECK( setNvtx(m_vertexKey) );
  } else {msg(MSG::DEBUG) << "Skip vertices" << endreq;}

  if(m_jetKey != "") {
    CHECK( setupJets(m_jetKey) );
  } else {msg(MSG::DEBUG) << "Skip jets" << endreq;}

  if(m_electronKey != "") {
    CHECK( setupElectrons(m_electronKey) );
  } else {msg(MSG::DEBUG) << "Skip electrons" << endreq;}

  if(m_photonKey != "") {
    CHECK( setupPhotons(m_photonKey) );
  } else {msg(MSG::DEBUG) << "Skip photons" << endreq;}

  if(m_taujetKey != "") {
    CHECK( setupTaus(m_taujetKey) );
  } else {msg(MSG::DEBUG) << "Skip taus" << endreq;}

  if(m_muonKey != "") {
    if(m_calomuonKey != "") {
      CHECK( setupMuons(m_muonKey, m_calomuonKey, m_muonStatusWord) );
    } else {
      CHECK( setupMuons(m_muonKey, m_muonStatusWord) );
    }
  } else {msg(MSG::DEBUG) << "Skip muons" << endreq;}

  if(m_clusterKey != "") {
    CHECK( setupClusters(m_clusterKey, m_clusterStatusWord) );
  } else {msg(MSG::DEBUG) << "Skip clusters" << endreq;}

  if(m_trackKey != "") {
    CHECK( setupTracks(m_trackKey) );
  } else {msg(MSG::DEBUG) << "Skip tracks" << endreq;}

  msg(MSG::DEBUG) << "End of process" << endreq;
  msg(MSG::DEBUG) << "  MET_RefFinal         = " << getMissingET(METUtil::RefFinal).et() << endreq;
  msg(MSG::DEBUG) << "  MET_RefFinal_JESUP   = " << getMissingET(METUtil::RefFinal,METUtil::JESUp).et() << endreq;
  msg(MSG::DEBUG) << "  MET_RefFinal_JESDOWN = " << getMissingET(METUtil::RefFinal,METUtil::JESDown).et() << endreq;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::finalize(){
  msg(MSG::INFO) << "Finalising MissingETUtility." << endreq;
  METUtility::reset();

  delete m_multiJES;
  delete m_JERTool;
  delete m_egammaTool;
  delete m_muonTool;
  delete m_TESTool;

  msg(MSG::DEBUG) << "Nice to have MET you." << endreq;

  return StatusCode::SUCCESS;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//  methods for setting up the physics objects from AOD containers
//

StatusCode METUtilityAthTool::setupJets(const string &collectionName) { 
  msg(MSG::DEBUG) << "Setting up jets" << endreq;

  // We grab the event info for the purposes of seeding random number generators
  const EventInfo* myEventInfo = 0;
  CHECK( evtStore()->retrieve(myEventInfo/*,"MyEvent"*/) );

  bool useMetComp = collectionName == "METComposition";

  const JetCollection *aod_jets = 0;
  vector<vector<float> > jet_wet, jet_wpx, jet_wpy;
  vector<vector<unsigned int> > jet_statusWord;
  int ijet = 0;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up jet collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the jet collection from the MET composition
    JetCollection *temp_jets = new JetCollection;
    MissingETComposition::const_iterator<Jet> fJet(m_metComp->beginIterator<Jet>());
    MissingETComposition::const_iterator<Jet> lJet(m_metComp->endIterator<Jet>());
    for(;fJet!=lJet;fJet++) {
      msg(MSG::VERBOSE) << "Adding jet " << ijet << " to collection" << endreq;
      if(!(*fJet)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      temp_jets->push_back(new Jet(**fJet));

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fJet.getParameter().wet());
      thewpx.push_back(fJet.getParameter().wpx());
      thewpy.push_back(fJet.getParameter().wpy());

      jet_wet.push_back(thewet);
      jet_wpx.push_back(thewpx);
      jet_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fJet.getParameter().statusWord());
      jet_statusWord.push_back(statusvec);
      ijet++;
    }

    aod_jets = temp_jets;
  } else {
    // Grab the jet collection from StoreGate
    // Less advisable, since the weights will not be properly accounted for.
    msg(MSG::DEBUG) << "Retrieving jet collection " << collectionName << " from StoreGate." << endreq;
    CHECK( evtStore()->retrieve( aod_jets, collectionName ) );
  }
  
  msg(MSG::DEBUG) << "Jet collection contains " << aod_jets->size() << " jets." << endreq;

  // Get the relevant information for the METUtility
  JetCollection::const_iterator iJet = aod_jets->begin();
  JetCollection::const_iterator eJet = aod_jets->end();
  vector<float> jet_pT, jet_eta, jet_phi, jet_E;
  vector<float> jet_jesUp, jet_jesDown, jet_resUp, jet_resDown;

  unsigned int EventNumber = myEventInfo->event_ID()->event_number();
  TRandom2 *jetRandom = new TRandom2;
  jetRandom->SetSeed(EventNumber);

  ijet = 0;
  jetRandom->SetSeed(EventNumber);
  for(;iJet!=eJet;iJet++) {
    msg(MSG::VERBOSE) << ijet;
    if(!(*iJet)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    jet_pT.push_back((*iJet)->pt());
    jet_eta.push_back((*iJet)->eta());
    jet_phi.push_back((*iJet)->phi());
    jet_E.push_back((*iJet)->e());

    // Didn't use MET Composition, so we have to make up some weights instead
    // Put in the default statusWord, which is not used for jets anyway
    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      jet_wet.push_back(univec);
      jet_wpx.push_back(univec);
      jet_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(MissingETTags::DEFAULT);
      jet_statusWord.push_back(statusvec);
    }

    // Now compute the uncertainties using the JES and JER providers
    // Don't apply uncertainties to the jets in SoftJets.
    float jesShift = 0.0;
    float jerShift = 0.0;
    msg(MSG::VERBOSE) << " Getting JES uncertainty... ";
    if((*iJet)->pt() > 20e3 && (*iJet)->pt() < 7000e3 && fabs((*iJet)->eta()) < 4.5){
      double mindR = 99.;
      if(m_doCloseByJES) {
	// Compute deltaR to closest 7 GeV jet
	JetCollection::const_iterator iJet2 = aod_jets->begin();
	
	if((*iJet)->pt() > 20e3 && (*iJet)->pt() < 5000e3) {
	  double pt = jet_pT.at(ijet);
	  double eta = jet_eta.at(ijet);

	  double S = m_JERTool->getRelResolutionMC(pt/1e3,eta);
	  double U = m_JERTool->getResolutionUncert(pt/1e3,eta);
	  double smearingFactorSyst = sqrt(pow(S+U,2)-pow(S,2));

	  msg(MSG::VERBOSE) << " Getting JER uncertainty...";
	  jerShift = jetRandom->Gaus(0, smearingFactorSyst);
	}

	for(;iJet2!=eJet;iJet2++) {
	  // Only take other jets with EM pT of at least 7 GeV
	  if(iJet2==iJet) continue;
	  double empt = (*iJet2)->pt(P4SignalState::JETEMSCALE);
	  if(empt<7e3) continue;

	  mindR = TMath::Min(mindR,P4Helpers::deltaR(*iJet,*iJet2));
	}
      }

      double mu = m_lbmuTool->averageInteractionsPerCrossing();
      // Fill the uncertainty vectors
      // The bool is true for +ve uncertainties
      // Can't use getRelPos/NegUncert because those don't take mu and NPV
      // Unfortuntely the definition of Nvtx is different for the jet and soft terms!
      jesShift = m_multiJES->getRelUncert((*iJet)->pt(), (*iJet)->eta(), mindR, true, m_nvtxjet, mu);
      jet_jesUp.push_back(jesShift);
      jesShift = m_multiJES->getRelUncert((*iJet)->pt(), (*iJet)->eta(), mindR, false, m_nvtxjet, mu);
      jet_jesDown.push_back(-1*jesShift);
    } else {
      jet_jesUp.push_back(jesShift);
      jet_jesDown.push_back(jesShift);
    }
    msg(MSG::VERBOSE) << "JESUP=" << jet_jesUp.at(ijet)
		    << ",JESDOWN=" << jet_jesDown.at(ijet);

    jet_resUp.push_back(jerShift);
    jet_resDown.push_back(-1*jerShift);

    msg(MSG::VERBOSE) << "JER=" << jerShift;

    msg(MSG::VERBOSE) << endreq;
    ijet++;
  }

  msg(MSG::DEBUG) << "Setting jet parameters" << endreq;
  setJetParameters(&jet_pT, &jet_eta, &jet_phi, &jet_E,
		   &jet_wet, &jet_wpx, &jet_wpy, &jet_statusWord);
  setOriJetParameters(&jet_pT);
 

  msg(MSG::DEBUG) << "Setting jet uncertainties" << endreq;
  setObjectEnergyUncertainties(METUtil::Jets, jet_jesUp, jet_jesDown);
  setObjectResolutionShift(METUtil::Jets, jet_resUp, jet_resDown);

  msg(MSG::DEBUG) << "Done setting up jets." << endreq;
  if(useMetComp) delete aod_jets;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setupElectrons(const string &collectionName) { 

  msg(MSG::DEBUG) << "Setting up electrons" << endreq;

  bool useMetComp = collectionName == "METComposition";

  int iel = 0;
  const ElectronContainer *aod_electrons = 0;
  vector<vector<float> > el_wet, el_wpx, el_wpy;
  vector<vector<unsigned int> > el_statusWord;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up electron collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the electron collection from the MET composition
    ElectronContainer *temp_electrons = new ElectronContainer;
    msg(MSG::DEBUG) << "Getting electron iterators" << endreq;
    MissingETComposition::const_iterator<Analysis::Electron> fEl(m_metComp->beginIterator<Analysis::Electron>());
    MissingETComposition::const_iterator<Analysis::Electron> lEl(m_metComp->endIterator<Analysis::Electron>());
    for(;fEl!=lEl;fEl++) {
      msg(MSG::VERBOSE) << "Adding electron " << iel << " to collection" << endreq;
      if(!(*fEl)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      if(!(*fEl)->author(egammaParameters::AuthorElectron)) continue;
      temp_electrons->push_back(new Analysis::Electron(**fEl));

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fEl.getParameter().wet());
      thewpx.push_back(fEl.getParameter().wpx());
      thewpy.push_back(fEl.getParameter().wpy());

      el_wet.push_back(thewet);
      el_wpx.push_back(thewpx);
      el_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fEl.getParameter().statusWord());
      el_statusWord.push_back(statusvec);
      iel++;
    }

    aod_electrons = temp_electrons;
  } else {
    // Grab the electron collection from StoreGate
    msg(MSG::DEBUG) << "Retrieving electron collection " << collectionName << " from StoreGate." << endreq;
    CHECK( evtStore()->retrieve( aod_electrons, collectionName ) );
  }
  
  msg(MSG::DEBUG) << "Electron collection contains " << aod_electrons->size() << " electrons." << endreq;

  iel = 0;
  ElectronContainer::const_iterator iEl = aod_electrons->begin();
  ElectronContainer::const_iterator eEl = aod_electrons->end();
  vector<float> el_pT, el_eta, el_phi, el_E;
  vector<float> el_smeared_pT, el_eesUp, el_eesDown, el_resUp, el_resDown;
  for(;iEl!=eEl;iEl++) {
    msg(MSG::VERBOSE) << iel;
    if(!(*iEl)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    el_pT.push_back((*iEl)->pt());
    el_eta.push_back((*iEl)->eta());
    el_phi.push_back((*iEl)->phi());
    el_E.push_back((*iEl)->e());

    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      el_wet.push_back(univec);
      el_wpx.push_back(univec);
      el_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(MissingETTags::DEFAULT);
      el_statusWord.push_back(statusvec);
    }
    msg(MSG::VERBOSE) << " Applying electron smearing... ";
    // Apply smearing/scaling
    m_egammaTool->SetRandomSeed(int(1.e+5*fabs((*iEl)->phi()))); //set the seed however you like, I'm following the SUSYTools example

    float smear = m_egammaTool->getSmearingCorrectionMeV((*iEl)->cluster()->eta(), (*iEl)->cluster()->e(), 0, true);
    float smearUp = m_egammaTool->getSmearingCorrectionMeV((*iEl)->cluster()->eta(), (*iEl)->cluster()->e(), 2, true);
    float smearDown = m_egammaTool->getSmearingCorrectionMeV((*iEl)->cluster()->eta(), (*iEl)->cluster()->e(), 1, true);

    if(m_useSmearedEl) smear = 1;
    el_smeared_pT.push_back(smear*(*iEl)->pt());
    el_resUp.push_back((smearUp - smear)/smear);
    el_resDown.push_back((smearDown - smear)/smear);

    msg(MSG::VERBOSE) << " Getting energy scale uncertainty" << endreq;

    float correction = 1.;
    if(m_isData)
      correction = m_egammaTool->applyEnergyCorrectionMeV((*iEl)->cluster()->eta(),(*iEl)->cluster()->phi(),
						       (*iEl)->e(),(*iEl)->cluster()->pt(),0,"ELECTRON") / (*iEl)->e();
    el_smeared_pT.at(iel) *= correction;
    double energyUp = m_egammaTool->applyEnergyCorrectionMeV((*iEl)->cluster()->eta(),(*iEl)->cluster()->phi(),
							  (*iEl)->e(),(*iEl)->cluster()->pt(),2,"ELECTRON") / (correction*(*iEl)->e()) - 1;
    double energyDown = m_egammaTool->applyEnergyCorrectionMeV((*iEl)->cluster()->eta(),(*iEl)->cluster()->phi(),
							    (*iEl)->e(),(*iEl)->cluster()->pt(),1,"ELECTRON") / (correction*(*iEl)->e()) - 1;

    el_eesUp.push_back(energyUp);
    el_eesDown.push_back(energyDown);

    iel++;
  }
  

  msg(MSG::DEBUG) << "Setting electron momenta...";
  setElectronParameters(&el_smeared_pT, &el_eta, &el_phi,
			&el_wet, &el_wpx, &el_wpy, &el_statusWord);

  msg(MSG::DEBUG) << "Setting energy scale uncertainty" << endreq;
  setObjectEnergyUncertainties(METUtil::Electrons, el_eesUp, el_eesDown);
  msg(MSG::DEBUG) << "Setting resolution uncertainty" << endreq;
  setObjectResolutionShift(METUtil::Electrons, el_resUp, el_resDown);

  if(useMetComp) delete aod_electrons;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setupPhotons(const string &collectionName) {

  msg(MSG::DEBUG) << "Setting up photons" << endreq;

  bool useMetComp = collectionName == "METComposition";

  const PhotonContainer *aod_photons = 0;
  vector<vector<float> > ph_wet, ph_wpx, ph_wpy;
  vector<vector<unsigned int> > ph_statusWord;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up photon collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the photon collection from the MET composition
    PhotonContainer *temp_photons = new PhotonContainer;
    MissingETComposition::const_iterator<Analysis::Photon> fPh(m_metComp->beginIterator<Analysis::Photon>());
    MissingETComposition::const_iterator<Analysis::Photon> lPh(m_metComp->endIterator<Analysis::Photon>());
    for(;fPh!=lPh;fPh++) {
      if(!(*fPh)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      temp_photons->push_back(new Analysis::Photon(**fPh));

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fPh.getParameter().wet());
      thewpx.push_back(fPh.getParameter().wpx());
      thewpy.push_back(fPh.getParameter().wpy());

      ph_wet.push_back(thewet);
      ph_wpx.push_back(thewpx);
      ph_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fPh.getParameter().statusWord());
      ph_statusWord.push_back(statusvec);
    }

    aod_photons = temp_photons;
  } else {
    // Grab the photon collection from StoreGate
    CHECK( evtStore()->retrieve( aod_photons, collectionName ) );
  }

  msg(MSG::DEBUG) << "Photon collection contains " << aod_photons->size() << " photons." << endreq;

  int iph = 0;
  PhotonContainer::const_iterator iPh = aod_photons->begin();
  PhotonContainer::const_iterator ePh = aod_photons->end();
  vector<float> ph_pT, ph_eta, ph_phi, ph_E;
  vector<float> ph_smeared_pT, ph_eesUp, ph_eesDown, ph_resUp, ph_resDown;
  for(;iPh!=ePh;iPh++) {
    if(!(*iPh)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    ph_pT.push_back((*iPh)->pt());
    ph_eta.push_back((*iPh)->eta());
    ph_phi.push_back((*iPh)->phi());
    ph_E.push_back((*iPh)->e());

    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      ph_wet.push_back(univec);
      ph_wpx.push_back(univec);
      ph_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(MissingETTags::DEFAULT);
      ph_statusWord.push_back(statusvec);
    }
    msg(MSG::VERBOSE) << " Applying electron smearing... ";
    // Apply smearing/scaling
    m_egammaTool->SetRandomSeed(int(1.e+5*fabs((*iPh)->phi()))); //set the seed however you like, I'm following the SUSYTools example

    float smear = m_egammaTool->getSmearingCorrectionMeV((*iPh)->cluster()->eta(), (*iPh)->cluster()->e(), 0, true);
    float smearUp = m_egammaTool->getSmearingCorrectionMeV((*iPh)->cluster()->eta(), (*iPh)->cluster()->e(), 2, true);
    float smearDown = m_egammaTool->getSmearingCorrectionMeV((*iPh)->cluster()->eta(), (*iPh)->cluster()->e(), 1, true);

    if(m_useSmearedPh) smear = 1;
    ph_smeared_pT.push_back(smear*(*iPh)->pt());
    ph_resUp.push_back((smearUp - smear)/smear);
    ph_resDown.push_back((smearDown - smear)/smear);

    msg(MSG::VERBOSE) << " Getting energy scale uncertainty" << endreq;

    string photontype = (*iPh)->nConversions()>0 ? "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON";
    float correction = 1.;
    if(m_isData)
      correction = m_egammaTool->applyEnergyCorrectionMeV((*iPh)->cluster()->eta(),(*iPh)->cluster()->phi(),
						       (*iPh)->e(),(*iPh)->cluster()->pt(),0,photontype) / (*iPh)->e();
    ph_smeared_pT.at(iph)*= correction;
    double energyUp = m_egammaTool->applyEnergyCorrectionMeV((*iPh)->cluster()->eta(),(*iPh)->cluster()->phi(),
							  (*iPh)->e(),(*iPh)->cluster()->pt(),2,photontype) / (correction*(*iPh)->e()) - 1;
    double energyDown = m_egammaTool->applyEnergyCorrectionMeV((*iPh)->cluster()->eta(),(*iPh)->cluster()->phi(),
							    (*iPh)->e(),(*iPh)->cluster()->pt(),1,photontype) / (correction*(*iPh)->e()) - 1;

    ph_eesUp.push_back(energyUp);
    ph_eesDown.push_back(energyDown);

    iph++;
  }

  setPhotonParameters(&ph_smeared_pT, &ph_eta, &ph_phi,
		      &ph_wet, &ph_wpx, &ph_wpy, &ph_statusWord);

  msg(MSG::DEBUG) << "Setting energy scale uncertainty" << endreq;
  setObjectEnergyUncertainties(METUtil::Photons, ph_eesUp, ph_eesDown);
  msg(MSG::DEBUG) << "Setting resolution uncertainty" << endreq;
  setObjectResolutionShift(METUtil::Photons, ph_resUp, ph_resDown);

  if(useMetComp) delete aod_photons;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setupTaus(const string &collectionName) {

  msg(MSG::DEBUG) << "Setting up taus" << endreq;

  bool useMetComp = collectionName == "METComposition";

  const Analysis::TauJetContainer *aod_taus = 0;
  vector<vector<float> > tau_wet, tau_wpx, tau_wpy;
  vector<vector<unsigned int> > tau_statusWord;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up tau collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the tau collection from the MET composition
    Analysis::TauJetContainer *temp_taus = new Analysis::TauJetContainer;
    MissingETComposition::const_iterator<Analysis::TauJet> fTau(m_metComp->beginIterator<Analysis::TauJet>());
    MissingETComposition::const_iterator<Analysis::TauJet> lTau(m_metComp->endIterator<Analysis::TauJet>());
    for(;fTau!=lTau;fTau++) {
      if(!(*fTau)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      temp_taus->push_back(new Analysis::TauJet(**fTau));

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fTau.getParameter().wet());
      thewpx.push_back(fTau.getParameter().wpx());
      thewpy.push_back(fTau.getParameter().wpy());

      tau_wet.push_back(thewet);
      tau_wpx.push_back(thewpx);
      tau_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fTau.getParameter().statusWord());
      tau_statusWord.push_back(statusvec);
    }

    aod_taus = temp_taus;
  } else {
    // Grab the tau collection from StoreGate
    CHECK( evtStore()->retrieve( aod_taus, collectionName ) );
  }

  msg(MSG::DEBUG) << "Tau collection contains " << aod_taus->size() << " taus." << endreq;

  Analysis::TauJetContainer::const_iterator iTau = aod_taus->begin();
  Analysis::TauJetContainer::const_iterator eTau = aod_taus->end();
  vector<float> tau_pT, tau_eta, tau_phi, tau_E;
  vector<float> tesUp, tesDown;
  for(;iTau!=eTau;iTau++) {
    if(!(*iTau)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    tau_pT.push_back((*iTau)->pt());
    tau_eta.push_back((*iTau)->eta());
    tau_phi.push_back((*iTau)->phi());
    tau_E.push_back((*iTau)->e());

    double uncert = m_TESTool->GetTESUncertainty((*iTau)->pt()/1e3,
						 (*iTau)->eta(),
						 (*iTau)->nProng());

      if(uncert < 0) uncert = 0;
      tesUp.push_back(uncert);
      tesDown.push_back(-1*uncert);

    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      tau_wet.push_back(univec);
      tau_wpx.push_back(univec);
      tau_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(MissingETTags::DEFAULT);
      tau_statusWord.push_back(statusvec);
    }
  }

  setTauParameters(&tau_pT, &tau_eta, &tau_phi,
		   &tau_wet, &tau_wpx, &tau_wpy, &tau_statusWord);

  setObjectEnergyUncertainties(METUtil::Taus,tesUp,tesDown);

  if(useMetComp) delete aod_taus;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setupMuons(const string &collectionName,
					 const unsigned short userStatusWord) {

  msg(MSG::DEBUG) << "Setting up muons" << endreq;

  bool useMetComp = collectionName == "METComposition";

  const Analysis::MuonContainer *aod_muons = 0;
  vector<vector<float> > mu_wet, mu_wpx, mu_wpy;
  vector<vector<unsigned int> > mu_statusWord;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up muon collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the muon collection from the MET composition
    Analysis::MuonContainer *temp_muons = new Analysis::MuonContainer;
    MissingETComposition::const_iterator<Analysis::Muon> fMu(m_metComp->beginIterator<Analysis::Muon>());
    MissingETComposition::const_iterator<Analysis::Muon> lMu(m_metComp->endIterator<Analysis::Muon>());
    for(;fMu!=lMu;fMu++) {
      if(!(*fMu)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      temp_muons->push_back(new Analysis::Muon(**fMu));

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fMu.getParameter().wet());
      thewpx.push_back(fMu.getParameter().wpx());
      thewpy.push_back(fMu.getParameter().wpy());

      mu_wet.push_back(thewet);
      mu_wpx.push_back(thewpx);
      mu_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fMu.getParameter().statusWord());
      mu_statusWord.push_back(statusvec);
    }

    aod_muons = temp_muons;
  } else {
    // Grab the muon collection from StoreGate
    CHECK( evtStore()->retrieve( aod_muons, collectionName ) );
  }

  msg(MSG::DEBUG) << "Muon collection contains " << aod_muons->size() << " muons." << endreq;

  int imu = 0;
  Analysis::MuonContainer::const_iterator iMu = aod_muons->begin();
  Analysis::MuonContainer::const_iterator eMu = aod_muons->end();
  vector<float> mu_pT, mu_eta, mu_phi, mu_E;
  vector<float> mu_MS_pT, mu_MS_theta, mu_MS_phi;
  vector<float> mu_smeared_pT, mu_smeared_MS_pT;
  vector<float> cb_meridUp, cb_meridDown;
  vector<float> cb_mermsUp, cb_mermsDown;
  vector<float> mermsUp,    mermsDown;
  vector<float> mesUp,      mesDown;
  for(;iMu!=eMu;iMu++) {
    if(!(*iMu)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    msg(MSG::VERBOSE) << imu;

    mu_pT.push_back((*iMu)->pt());
    mu_eta.push_back((*iMu)->eta());
    mu_phi.push_back((*iMu)->phi());
    mu_E.push_back((*iMu)->e());

    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      mu_wet.push_back(univec);
      mu_wpx.push_back(univec);
      mu_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(userStatusWord);
      mu_statusWord.push_back(statusvec);
    }

    msg(MSG::VERBOSE) << " Getting track particle...";
    const Rec::TrackParticle *thisMSmuon = (*iMu)->muonSpectrometerTrackParticle();
    if(thisMSmuon) {
      mu_MS_pT.push_back(thisMSmuon->pt());
      mu_MS_theta.push_back(asin(thisMSmuon->sinTh()));
      mu_MS_phi.push_back(thisMSmuon->phi());
    } else {
      msg(MSG::VERBOSE) << " FAILED! " << endreq;
      mu_MS_pT.push_back(0.);
      mu_MS_theta.push_back(0.);
      mu_MS_phi.push_back(0.);
    }

    double ptcb = (*iMu)->pt();
    double ptid = (*iMu)->hasInDetTrackParticle() ?
      (*iMu)->inDetTrackParticle()->pt() : 0.;
    double ptms = (*iMu)->hasMuonExtrapolatedTrackParticle() ?
      (*iMu)->muonExtrapolatedTrackParticle()->pt() : 0.;
    m_muonTool->SetSeed(int(1.e+5*fabs((*iMu)->phi())));
    double etaMu = (*iMu)->eta();
    m_muonTool->Event(ptms, ptid, ptcb, etaMu);
    msg(MSG::VERBOSE) << " Original muon pT = " << ptcb << " (cb), " << ptms << " (ms)";

    msg(MSG::VERBOSE) << " Doing muon smearing...";
    float smearedCombinedPt =m_muonTool->pTCB();
    if(!(*iMu)->isCombinedMuon())
      smearedCombinedPt = m_muonTool->pTMS() + m_muonTool->pTID();
    float smearedMSPt = m_muonTool->pTMS();

    msg(MSG::VERBOSE) << " smeared pT = " << ptcb << " (cb), " << ptms << " (ms)";
    mu_smeared_pT.push_back(smearedCombinedPt);
    mu_smeared_MS_pT.push_back(smearedMSPt);

    double ptMS_smeared, ptID_smeared, ptCB_smeared;
    float smearedpTMS, smearedpTID, smearedpTCB;
    smearedpTMS = 0.1; smearedpTID = 0.1; smearedpTCB = 0.1;
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSLOW");
    smearedpTMS = ptMS_smeared/smearedMSPt - 1.0;
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    mermsDown.push_back(smearedpTMS);
    cb_mermsDown.push_back(smearedpTCB);
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSUP");
    smearedpTMS = ptMS_smeared/smearedMSPt - 1.0;
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    mermsUp.push_back(smearedpTMS);
    cb_mermsUp.push_back(smearedpTCB);
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDUP");
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    cb_meridUp.push_back(smearedpTCB);
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDLOW");
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    cb_meridDown.push_back(smearedpTCB);

    // get energies uncertainty from smearing class itself, 
    // it's in there but without a function to retrieve them. And only for CB?
    int detRegion = m_muonTool->DetRegion();

    if(detRegion==-1) detRegion = 3;
    double scalesyst = m_muonTool->getScaleSyst_CB().at(detRegion);    
    mesUp.push_back(scalesyst);
    mesDown.push_back(-scalesyst);

    msg(MSG::VERBOSE) << endreq;
    imu++;
  }

  if(m_useSmearedMu) {
    msg(MSG::DEBUG) << "Setting muon parameters with smeared muons" << endreq;
    setMuonParameters(&mu_smeared_pT, &mu_eta, &mu_phi,
		      &mu_wet, &mu_wpx, &mu_wpy, &mu_statusWord);
    msg(MSG::DEBUG) << "Setting extra muon parameters with smeared muons" << endreq;
    setExtraMuonParameters(&mu_smeared_MS_pT, &mu_MS_theta, &mu_MS_phi);
  } else {
    msg(MSG::DEBUG) << "Setting muon parameters with original muons" << endreq;
    setMuonParameters(&mu_pT, &mu_eta, &mu_phi,
		      &mu_wet, &mu_wpx, &mu_wpy, &mu_statusWord);
    msg(MSG::DEBUG) << "Setting extra muon parameters with original muons" << endreq;
    setExtraMuonParameters(&mu_MS_pT, &mu_MS_theta, &mu_MS_phi);
  } 

  msg(MSG::DEBUG) << "Setting muon resolution" << endreq;
  // muons are complicated, and MET makes use of track, spectro, and combined quantites,
  // so we need all off their resolutions. comboms reflects that it is 
  // the combined muon res affected by shifting ms res up and down.
  // comboid is for shifting the id res up and down
  setObjectResolutionShift(METUtil::MuonsComboMS, cb_mermsUp, cb_mermsDown);
  setObjectResolutionShift(METUtil::MuonsComboID, cb_meridUp, cb_meridDown);
  setObjectResolutionShift(METUtil::SpectroMuons, mermsUp, mermsDown);

  msg(MSG::DEBUG) << "Setting muon energy scale uncertainties" << endreq;
  // for now the mes only affects combined
  setObjectEnergyUncertainties(METUtil::Muons, mesUp, mesDown);

  if(useMetComp) delete aod_muons;

  return StatusCode::SUCCESS;
}

// For setting up muons and including the CaloMuonContainer
// No point in doing this if the METComposition is being used

StatusCode METUtilityAthTool::setupMuons(const string &collectionName,
					 const string &calocollectionName,
					 const unsigned short userStatusWord) {

  msg(MSG::DEBUG) << "Setting up muons" << endreq;

  const Analysis::MuonContainer *aod_muons = 0;
  const Analysis::MuonContainer *aod_calomuons = 0;
  vector<vector<float> > mu_wet, mu_wpx, mu_wpy;
  vector<vector<unsigned int> > mu_statusWord;
    // Grab the muon collection from StoreGate
  CHECK( evtStore()->retrieve( aod_muons, collectionName ) );
  CHECK( evtStore()->retrieve( aod_calomuons, calocollectionName ) );

  Analysis::MuonContainer *all_muons = new Analysis::MuonContainer;
  Analysis::MuonContainer::const_iterator fMu = aod_muons->begin();
  Analysis::MuonContainer::const_iterator lMu = aod_muons->end();
  for(; fMu!=lMu; fMu++) {
    all_muons->push_back(new Analysis::Muon(**fMu));
  }
  Analysis::MuonContainer::const_iterator fCMu = aod_calomuons->begin();
  Analysis::MuonContainer::const_iterator lCMu = aod_calomuons->end();
  for(; fCMu!=lCMu; fCMu++) {
    // Avoid double counting muons -- the caloMuons are just to fill in the central crack
    bool hasmatch;
    if((*fCMu)->caloMuonAlsoFoundByMuonReco() > 0 || fabs((*fCMu)->eta())>0.1) continue;
    for (fMu = aod_muons->begin(); fMu!=lMu; fMu++) {
      // select combined tracks with best match or standalone track (in |eta|>2.5)               
      if ((*fMu)->isCombinedMuon() || ((*fMu)->isStandAloneMuon() &&  fabs((*fMu)->eta())>2.5) ) {
	if (fabs((*fMu)->eta()) < 0.12 ) {
	  double dR = P4Helpers::deltaR(*fMu,*fCMu);
	  if (dR < 0.02) {
	    hasmatch=true;
	  }
	}
      }
    }

    if(!hasmatch) all_muons->push_back(new Analysis::Muon(**fCMu));
  }

  msg(MSG::DEBUG) << "Muon collection contains " << aod_muons->size() << " muons." << endreq;

  Analysis::MuonContainer::const_iterator iMu = all_muons->begin();
  Analysis::MuonContainer::const_iterator eMu = all_muons->end();
  vector<float> mu_pT, mu_eta, mu_phi, mu_E;
  vector<float> mu_MS_pT, mu_MS_theta, mu_MS_phi;
  vector<float> mu_smeared_pT, mu_smeared_MS_pT;
  vector<float> cb_meridUp, cb_meridDown;
  vector<float> cb_mermsUp, cb_mermsDown;
  vector<float> mermsUp,    mermsDown;
  vector<float> mesUp,      mesDown;
  for(;iMu!=eMu;iMu++) {
    if(!(*iMu)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    mu_pT.push_back((*iMu)->pt());
    mu_eta.push_back((*iMu)->eta());
    mu_phi.push_back((*iMu)->phi());
    mu_E.push_back((*iMu)->e());

    vector<float> univec;
    univec.push_back(1.);
    mu_wet.push_back(univec);
    mu_wpx.push_back(univec);
    mu_wpy.push_back(univec);

    vector<unsigned int> statusvec;
    statusvec.push_back(userStatusWord);
    mu_statusWord.push_back(statusvec);

    const Rec::TrackParticle *thisMSmuon = (*iMu)->muonSpectrometerTrackParticle();
    if(thisMSmuon) {
      mu_MS_pT.push_back(thisMSmuon->pt());
      mu_MS_theta.push_back(asin(thisMSmuon->sinTh()));
      mu_MS_phi.push_back(thisMSmuon->phi());
    } else {
      mu_MS_pT.push_back(0.);
      mu_MS_theta.push_back(0.);
      mu_MS_phi.push_back(0.);
    }

    double ptcb = (*iMu)->pt();
    double ptid = (*iMu)->hasInDetTrackParticle() ?
      (*iMu)->inDetTrackParticle()->pt() : 0.;
    double ptms = (*iMu)->hasMuonExtrapolatedTrackParticle() ?
      (*iMu)->muonExtrapolatedTrackParticle()->pt() : 0.;
    m_muonTool->SetSeed(int(1.e+5*fabs((*iMu)->phi())));
    double etaMu = (*iMu)->eta();
    m_muonTool->Event(ptms, ptid, ptcb, etaMu);


    float smearedCombinedPt =m_muonTool->pTCB();
    if(!(*iMu)->isCombinedMuon())
      smearedCombinedPt = m_muonTool->pTMS() + m_muonTool->pTID();

    float smearedMSPt = m_muonTool->pTMS();

    mu_smeared_MS_pT.push_back(ptms);
    mu_smeared_pT.push_back(smearedCombinedPt);

    double ptMS_smeared, ptID_smeared, ptCB_smeared;
    float smearedpTMS, smearedpTID, smearedpTCB;
    smearedpTMS = 0.1; smearedpTID = 0.1; smearedpTCB = 0.1;
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSLOW");
    smearedpTMS = ptMS_smeared/smearedMSPt - 1.0;
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    mermsDown.push_back(smearedpTMS);
    cb_mermsDown.push_back(smearedpTCB);
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSUP");
    smearedpTMS = ptMS_smeared/smearedMSPt - 1.0;
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    mermsUp.push_back(smearedpTMS);
    cb_mermsUp.push_back(smearedpTCB);
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDUP");
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    cb_meridUp.push_back(smearedpTCB);
    m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDLOW");
    smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
    cb_meridDown.push_back(smearedpTCB);
	
    // get energies uncertainty from smearing class itself, 
    // it's in there but without a function to retrieve them. And only for CB?
    if(fabs((*iMu)->eta()) < 1.05){
      mesUp.push_back(0.0002);
      mesDown.push_back(-0.0002);
    }
    else if (fabs((*iMu)->eta()) < 1.7){
      mesUp.push_back(0.0006);
      mesDown.push_back(-0.0006);
    }
    else if (fabs((*iMu)->eta()) < 2.0){
      mesUp.push_back(0.0012);
      mesDown.push_back(-0.0012);
    }
    else {
      mesUp.push_back(0.0007);
      mesDown.push_back(-0.0007);
    }

  }

  if(m_useSmearedMu) {
    setMuonParameters(&mu_smeared_pT, &mu_eta, &mu_phi,
		      &mu_wet, &mu_wpx, &mu_wpy, &mu_statusWord);
    setExtraMuonParameters(&mu_smeared_MS_pT, &mu_MS_theta, &mu_MS_phi);
  } else {
    setMuonParameters(&mu_pT, &mu_eta, &mu_phi,
		      &mu_wet, &mu_wpx, &mu_wpy, &mu_statusWord);
    setExtraMuonParameters(&mu_MS_pT, &mu_MS_theta, &mu_MS_phi);
  } 
  // muons are complicated, and MET makes use of track, spectro, and combined quantites,
  // so we need all off their resolutions. comboms reflects that it is 
  // the combined muon res affected by shifting ms res up and down.
  // comboid is for shifting the id res up and down
  setObjectResolutionShift(METUtil::MuonsComboMS, cb_mermsUp, cb_mermsDown);
  setObjectResolutionShift(METUtil::MuonsComboID, cb_meridUp, cb_meridDown);
  setObjectResolutionShift(METUtil::SpectroMuons, mermsUp, mermsDown);

  // for now the mes only affects combined
  setObjectEnergyUncertainties(METUtil::Muons, mesUp, mesDown);

  delete all_muons;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setupClusters(const string &collectionName,
					    const unsigned short userStatusWord) {
  msg(MSG::DEBUG) << "Setting up clusters" << endreq;

  bool useMetComp = collectionName == "METComposition";

  const CaloClusterContainer *aod_clusters = 0;
  vector<vector<float> > cl_wet, cl_wpx, cl_wpy;
  vector<vector<unsigned int> > cl_statusWord;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up cluster collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the cluster collection from the MET composition
    int iclus = 0;
    CaloClusterContainer *temp_clusters = new CaloClusterContainer;
    MissingETComposition::const_iterator<CaloCluster> fCl(m_metComp->beginIterator<CaloCluster>());
    MissingETComposition::const_iterator<CaloCluster> lCl(m_metComp->endIterator<CaloCluster>());
    for(;fCl!=lCl;fCl++) {
      if(!(*fCl)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      temp_clusters->push_back(new CaloCluster(**fCl));

      msg(MSG::VERBOSE) << "Cluster " << iclus << endl;

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fCl.getParameter().wet());
      thewpx.push_back(fCl.getParameter().wpx());
      thewpy.push_back(fCl.getParameter().wpy());

      cl_wet.push_back(thewet);
      cl_wpx.push_back(thewpx);
      cl_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fCl.getParameter().statusWord());
      cl_statusWord.push_back(statusvec);
//       msg(MSG::VERBOSE) << "statusWord = " << fCl.getParameter().statusWord() << endl;
      iclus++;
    }

    aod_clusters = temp_clusters;
  } else {
    // Grab the cluster collection from StoreGate
    CHECK( evtStore()->retrieve( aod_clusters, collectionName ) );
  }

  msg(MSG::DEBUG) << "Cluster collection contains " << aod_clusters->size() << " clusters." << endreq;

  CaloClusterContainer::const_iterator iCl = aod_clusters->begin();
  CaloClusterContainer::const_iterator eCl = aod_clusters->end();
  vector<float> cl_pT, cl_eta, cl_phi, cl_E;
  for(;iCl!=eCl;iCl++) {
    if(!(*iCl)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    cl_pT.push_back((*iCl)->pt());
    cl_eta.push_back((*iCl)->eta());
    cl_phi.push_back((*iCl)->phi());
    cl_E.push_back((*iCl)->e());

    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      cl_wet.push_back(univec);
      cl_wpx.push_back(univec);
      cl_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(userStatusWord);
      cl_statusWord.push_back(statusvec);
    }
  }

  setClusterParameters(&cl_pT, &cl_eta, &cl_phi, &cl_E,
		       &cl_wet, &cl_wpx, &cl_wpy, &cl_statusWord);

  if(useMetComp) delete aod_clusters;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setupTracks(const string &collectionName) {
  msg(MSG::DEBUG) << "Setting up tracks" << endreq;

  bool useMetComp = collectionName == "METComposition";

  const Rec::TrackParticleContainer *aod_tracks = 0;
  vector<vector<float> > tr_wet, tr_wpx, tr_wpy;
  vector<vector<unsigned int> > tr_statusWord;
  if(useMetComp) {
    msg(MSG::DEBUG) << "Setting up track collection from MET composition" << endreq;
    if(!m_metComp) {
      msg(MSG::WARNING) << "Invalid pointer to MET composition object!" << endreq;
      return StatusCode::FAILURE;
    }
    // Fill the track collection from the MET composition
    int itrk = 0;
    Rec::TrackParticleContainer *temp_tracks = new Rec::TrackParticleContainer;
    MissingETComposition::const_iterator<Rec::TrackParticle> fTr(m_metComp->beginIterator<Rec::TrackParticle>());
    MissingETComposition::const_iterator<Rec::TrackParticle> lTr(m_metComp->endIterator<Rec::TrackParticle>());
    for(;fTr!=lTr;fTr++) {
      if(!(*fTr)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
      temp_tracks->push_back(new Rec::TrackParticle(**fTr));

      msg(MSG::VERBOSE) << "Track " << itrk << endl;

      vector<float> thewet, thewpx, thewpy;
      thewet.push_back(fTr.getParameter().wet());
      thewpx.push_back(fTr.getParameter().wpx());
      thewpy.push_back(fTr.getParameter().wpy());

      tr_wet.push_back(thewet);
      tr_wpx.push_back(thewpx);
      tr_wpy.push_back(thewpy);

      vector<unsigned int> statusvec;
      statusvec.push_back(fTr.getParameter().statusWord());
      tr_statusWord.push_back(statusvec);
      msg(MSG::VERBOSE) << "statusWord = " << fTr.getParameter().statusWord() << endl;
      itrk++;
    }

    aod_tracks = temp_tracks;
  } else {
    // Grab the track collection from StoreGate
    CHECK( evtStore()->retrieve( aod_tracks, collectionName ) );
  }

  msg(MSG::DEBUG) << "Track collection contains " << aod_tracks->size() << " tracks." << endreq;

  Rec::TrackParticleContainer::const_iterator iTr = aod_tracks->begin();
  Rec::TrackParticleContainer::const_iterator eTr = aod_tracks->end();
  vector<float> tr_pT, tr_eta, tr_phi, tr_E;
  for(;iTr!=eTr;iTr++) {
    if(!(*iTr)) msg(MSG::FATAL) << "Oops, didn't get the point(er)" << endreq;
    tr_pT.push_back((*iTr)->pt());
    tr_eta.push_back((*iTr)->eta());
    tr_phi.push_back((*iTr)->phi());
    tr_E.push_back((*iTr)->e());

    if(!useMetComp) {
      vector<float> univec;
      univec.push_back(1.);
      tr_wet.push_back(univec);
      tr_wpx.push_back(univec);
      tr_wpy.push_back(univec);

      vector<unsigned int> statusvec;
      statusvec.push_back(MissingETTags::DEFAULT);
      tr_statusWord.push_back(statusvec);
    }
  }

  setTrackParameters(&tr_pT, &tr_eta, &tr_phi, &tr_E,
		     &tr_wet, &tr_wpx, &tr_wpy, &tr_statusWord);

  if(useMetComp) delete aod_tracks;

  return StatusCode::SUCCESS;
}

StatusCode METUtilityAthTool::setNvtx(const std::string &collectionName)
{
    const VxContainer* aod_vertices;
    CHECK( evtStore()->retrieve( aod_vertices, collectionName ) );

    m_nvtx = 0;
    m_nvtxjet = 0;
    VxContainer::const_iterator iVx;
    VxContainer::const_iterator eVx = aod_vertices->end();
    bool goodPV = false;
    for(iVx = aod_vertices->begin(); iVx != eVx; iVx++) {
      if((*iVx)->vertexType() == Trk::PriVtx) {
	if((*iVx)->vxTrackAtVertex()->size()>2 && fabs((*iVx)->recVertex().position().z()<200.)) {
	  goodPV = true;
	}
      }
    }
    if(goodPV) {
      for(iVx = aod_vertices->begin(); iVx != eVx; iVx++) {
	if((*iVx)->vxTrackAtVertex()->size()>2) m_nvtx++;
	if((*iVx)->vxTrackAtVertex()->size()>2) m_nvtxjet++;
      }
    }

    return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//  methods for getting and setting the desired MET terms
//

void METUtilityAthTool::setMETTerm(const METUtil::Terms term, const MissingET* thisMET)
{
  msg(MSG::DEBUG) << "Set MET term " << METUtil::getTermName(term) << endreq;
  METUtility::setMETTerm(term,thisMET->etx(),thisMET->ety(),thisMET->sumet());
}

MissingET METUtilityAthTool::getMETTerm(const METUtil::Terms term, const METUtil::Systematics systematic)
{
  MissingET thisMET;
  METObject myMET = getMissingET(term,systematic);
  thisMET.setEx(myMET.etx());
  thisMET.setEy(myMET.ety());
  thisMET.setEtSum(myMET.sumet());
  return thisMET;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
//  trigger process on hitting a new event
//

void METUtilityAthTool::handle(const Incident& inc) {
  // Get the messaging service, print where you are

  msg(MSG::DEBUG) << "entering handle(), incidence type " << inc.type()
		  << " from " << inc.source() << endreq;

  // Only call fillIOV for EndEvent incident
  if (inc.type() != "BeginEvent") return;
  if(m_doEveryEvent) 
    if(process().isFailure()) msg(MSG::FATAL) << "Setup of METUtility failed" << endreq;

  msg(MSG::DEBUG) << "end event handle" << endreq;

}
