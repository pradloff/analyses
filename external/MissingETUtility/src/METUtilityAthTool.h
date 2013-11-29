/// Dear emacs, this is -*-c++-*-
///***********************************************************************
///
/// METUtilityAthTool
/// Authors: J. Goodson, T.J. Khoo
///
/// METUtility class, for use in ATHENA
///
///#### Usage:
///
///
///***********************************************************************

#ifndef _METUTILITYATHTOOL_
#define _METUTILITYATHTOOL_

#include "MissingETUtility/IMETUtilityAthTool.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/ToolHandle.h"
#include "LumiBlockComps/ILumiBlockMuTool.h"

namespace eg2011 {
  class EnergyRescaler;
}

namespace JetResolution {
  class JERProviderAthena;
}

namespace MuonSmear{
  class SmearingClass;
}

class MultijetJESUncertaintyProvider;
class TESUncertaintyProvider;
class MissingET;
class MissingETComposition;

namespace MissingETUtility {

  class METUtilityAthTool : virtual public AthAlgTool, 
			    virtual public IMETUtilityAthTool,
			    virtual public IIncidentListener {

  public:
  
    // Constructors, destructor
    METUtilityAthTool(const std::string& type, const std::string& name,
		      const IInterface* pParent);

    virtual ~METUtilityAthTool();

    StatusCode queryInterface( const InterfaceID& riid, void** ppvIf );

    virtual StatusCode initialize();
    virtual StatusCode process();
    virtual StatusCode finalize();

    virtual void handle(const Incident& inc);

    virtual StatusCode setupJets(const std::string &collectionName);
    virtual StatusCode setupElectrons(const std::string &collectionName);
    virtual StatusCode setupPhotons(const std::string &collectionName);
    virtual StatusCode setupTaus(const std::string &collectionName);
    virtual StatusCode setupMuons(const std::string &collectionName,
				  const unsigned short userStatusWord);
    virtual StatusCode setupMuons(const std::string &collectionName,
				  const std::string &calocollectionName,
				  const unsigned short userStatusWord);
    virtual StatusCode setupClusters(const std::string &collectionName,
				     const unsigned short userStatusWord);
    virtual StatusCode setupTracks(const std::string &collectionName);

    virtual StatusCode setNvtx(const std::string &collectionName);

    virtual void       setMETTerm(const METUtil::Terms term, const MissingET* thisMET);
    virtual MissingET  getMETTerm(const METUtil::Terms term,
				  const METUtil::Systematics systematic = METUtil::None);

  private:  
    ToolHandle< ILumiBlockMuTool > m_lbmuTool; ///< Tool returning the mu values
    MultijetJESUncertaintyProvider* m_multiJES;
    eg2011::EnergyRescaler* m_egammaTool;
    MuonSmear::SmearingClass* m_muonTool;
    JetResolution::JERProviderAthena* m_JERTool;
    TESUncertaintyProvider* m_TESTool;

    bool   m_isData;
    bool   m_doEveryEvent;

    int    m_nvtxjet;
    bool   m_doCloseByJES;
    double m_ptCloseByJES;
    string m_jetTypeJES;
    string m_mcTypeJES;
    string m_configfileJES;
    string m_configMultijet;

    bool m_verbosity;
    bool m_dorefele;
    bool m_dorefgamma;
    bool m_doreftau;
    bool m_dorefjet;
    bool m_dosoftjets;
    bool m_dorefmuon;
    bool m_domuontotal;
    bool m_docellout;
    bool m_docellouteflow;
    bool m_ismuid;
    double m_softjetcut; 

    string m_jetTypeJER;
    string m_methodJER;
    string m_inputfileJER;

    bool   m_useSmearedEl;
    bool   m_useSmearedPh;
    bool   m_useSmearedMu;
    bool   m_doSignif;

    bool m_flagis2012;
    bool m_flagisSTVF;

    unsigned short m_muonStatusWord;
    unsigned short m_clusterStatusWord;

    string m_clusterKey;
    string m_trackKey;
    string m_electronKey;
    string m_jetKey;
    string m_muonKey;
    string m_calomuonKey;
    string m_photonKey;
    string m_taujetKey;
    string m_vertexKey;

    string m_RefEleKey;
    string m_RefGammaKey;
    string m_RefJetKey;
    string m_RefTauKey;
    string m_RefMuonKey;
    string m_SoftJetsKey;
    string m_CellOutKey;
    string m_CellOutEflowKey;
    string m_TruthNonIntKey;

//     bool m_useMetComp;
    string m_metCompKey;
    const MissingETComposition *m_metComp;
   
  };
}

#endif
