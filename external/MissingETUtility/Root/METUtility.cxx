//see header for instruction
#include<cmath>
#include "MissingETUtility/METUtility.h"
#include "TRandom.h"

using std::fabs;

#ifdef METUTIL_STANDALONE
ClassImp(METObject);
ClassImp(PhysicsObject);
ClassImp(METUtility);
#endif //METUTIL_STANDALONE


///////////////////////////////////////////////////////////////////////////////////////////////
///functions for PhysicsObjects helper class///////////////////////////////////////////////////   
///////////////////////////////////////////////////////////////////////////////////////////////

PhysicsObject::PhysicsObject() : 
  m_momenergy(0.,0.,0.,0.),
  m_secondaryMomenergy(0.,0.,0.,0.),
  m_statusWord(0),
  m_relativeResolution(0.),
  m_index(0),
  m_relativeResolutionShift(0.,0.),
  m_relativeResolutionComboMS(0.,0.),
  m_relativeResolutionSpectro(0.,0.)
{

  m_weights[0] = 0.;
  m_weights[1] = 0.;
  m_weights[2] = 0.;

}

PhysicsObject::~PhysicsObject()
{

  m_relEnergyUncertainty.clear();

} // end 

void PhysicsObject::setMomenergy(float _pT, float _eta, float _phi, float _E, bool isSecondary)
{

  // TLorentzVector::SetPtEtaPhiE() applies a fabs() to the pT
  // This screws up the angle in the case of -ve pT clusters/objects
  double pi = atan(1)*4;
  if(_pT < 0) _phi += pi;
 
  if(isSecondary)
    m_secondaryMomenergy.SetPtEtaPhiE(_pT, _eta, _phi, _E);
  
  else m_momenergy.SetPtEtaPhiE(_pT, _eta, _phi, _E);

}// end of setMomenergy


void PhysicsObject::setMomenergy(const TLorentzVector& _vec, bool isSecondary)
{

  if(isSecondary)
    m_secondaryMomenergy = _vec;
  
  else m_momenergy = _vec;

}// end of setMomenergy


void PhysicsObject::setEnergyUncertainty(float energyUp, float energyDown, EnergyShift::ShiftType iShift)
{

  if ((int) m_relEnergyUncertainty.size()<iShift+1)
  {
    m_relEnergyUncertainty.resize(iShift+1);
  }

  m_relEnergyUncertainty[iShift].first = energyUp; 
  m_relEnergyUncertainty[iShift].second = energyDown;

}// end of setEnergyUncertainty

void PhysicsObject::setResShift(float resUp, float resDown, unsigned short type)
{
  
  if(type == 2) { // Spectro muons
    m_relativeResolutionSpectro.first = resUp; 
    m_relativeResolutionSpectro.second = resDown;    
  }
  else if(type == 1) { // Spectro for combined muons
    m_relativeResolutionComboMS.first = resUp; 
    m_relativeResolutionComboMS.second = resDown;
  }
  else{ // ID for combined muons or anything else
    m_relativeResolutionShift.first = resUp; 
    m_relativeResolutionShift.second = resDown;
  }
  
}// end of setResUncertainty

pair<float, float> PhysicsObject::resShift(unsigned short type) const
{

  if(type == 3) { // Cluster standard
    float  _resUp = 0.0;
    float  _resDown = 0.0;
    pair<float, float> _resolution;
    _resolution.first = _resUp;
    _resolution.second = _resDown;
    return _resolution;
  } else if(type == 2) { // Spectro muons
    return m_relativeResolutionSpectro;
  }
  else if(type == 1) { // Spectro for combined muons
    return m_relativeResolutionComboMS;
  }

  // ID for combined muons, or anything else
  return m_relativeResolutionShift;
}// end of resShift


pair<float, float> PhysicsObject::energyShift(unsigned short type, EnergyShift::ShiftType iShift) const
{
  if(type == 3) { // Cluster standard
    // This is currently unused, but could be if a default
    // cluster resolution uncertainty is determined.
    float a = (fabs(m_momenergy.Eta()) < 3.2 ) ? 0.03 : 0.1;
    float shift =  a * (1200.0/fabs(Pt()*wet()));
    pair<float, float> uncertainty;
    uncertainty.first = shift;
    uncertainty.second = -1.0*shift;
    return uncertainty;
  }

  if (iShift>=(int)m_relEnergyUncertainty.size())
    cerr << "METUTILITY: Requested energy uncertainty " << iShift
	 << " but only " << m_relEnergyUncertainty.size() << " available." << endl;

  return m_relEnergyUncertainty.at(iShift);
}// end of energyShift

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////METUtility classes////////////////////
///////////////////////////////////////////////////////////////////////////////


METUtility::METUtility(bool doRefEle, bool doRefGamma, bool doRefTau, bool doRefJet,
		       bool doSoftJets, bool doRefMuon, bool doMuonTotal, bool doCellOut,
		       bool doCellOutEflow, bool isMuid, double softJetCut, bool verbose) :
  m_useStandardClusterRes(false),
  m_useStandardClusterEnergySigma(false),
  m_useJesEigenvectors(false),
  m_addJesPileupSyst(false),
  m_doSignificance(false),
  m_nvtx(0),
  m_averageIntPerXing(0.),
  m_rndGaus(0.),
  m_jet_pu_statusWord(MissingETTags::DEFAULT)
{

  m_verbose = verbose;
  m_etaCut[0] = 0;
  m_etaCut[1] = 10.;

  m_doTerm[METUtil::RefEle]    = doRefEle;
  m_doTerm[METUtil::RefGamma]  = doRefGamma;
  m_doTerm[METUtil::RefTau]    = doRefTau;
  m_doTerm[METUtil::RefJet]    = doRefJet;
  m_doTerm[METUtil::SoftJets]  = doSoftJets;
  m_doTerm[METUtil::RefMuon]   = doRefMuon;
  m_doTerm[METUtil::MuonTotal] = doMuonTotal;
  m_doTerm[METUtil::CellOut]   = doCellOut;
  m_doTerm[METUtil::CellOutEflow] = doCellOutEflow;

  m_softJetCut = softJetCut;
  
  if(m_doTerm[METUtil::SoftJets] == false) m_softJetCut = 0.0; // so smeared jets don't get lost from the RefJet term
  
  m_isMuid = isMuid;

  // is2012, isSTVF
  configMissingET(true,false);

  reset();
}// end constructor


METUtility::~METUtility() {

  reset();

}// end destructor


void METUtility::reset()
{

  // Reset object vectors
  for(int iObj = METUtil::Electrons; iObj<=METUtil::Tracks; iObj++) {
    m_objects[iObj].clear();
    m_momentaSet[iObj] = false;
    m_scaleShiftsSet[iObj] = false;
    m_resShiftsSet[iObj] = false;
    m_resolutionsSet[iObj] = false;
  }
  // A few extra cases for muons
  for(int iObj = METUtil::SpectroMuons; iObj<=METUtil::MuonsComboID; iObj++) {
    m_resShiftsSet[iObj] = false;
    m_resolutionsSet[iObj] = false;
  }
  // And jets
  m_oriJetsSet = false;

  // Reset stored MET terms
  for(int iTerm = METUtil::RefEle; iTerm<=METUtil::RefFinal; iTerm++) {
    m_terms[iTerm].setBase(0.,0.,0.);
    m_termsScaled[iTerm].setBase(0.,0.,0.);
  }

  m_nvtx = 0;
  m_averageIntPerXing = 0;

}// end of reset

METObject METUtility::deltaMissingET(const METUtil::Terms term, METUtil::SystType systtype)
{

  vector<METUtil::SystType> sources;
  if(systtype == METUtil::All) {
    if(m_doTerm[METUtil::Jets]) {
      if (!m_useJesEigenvectors) {
	sources.push_back(METUtil::JES);
	//        sources.push_back(METUtil::bJES);
      } else {
        sources.push_back(METUtil::JES1);
        sources.push_back(METUtil::JES2);
        sources.push_back(METUtil::JES3);
        sources.push_back(METUtil::JES4);
        sources.push_back(METUtil::JES5);
        sources.push_back(METUtil::JES6);
        sources.push_back(METUtil::JESEtaIntercalibration);
        sources.push_back(METUtil::JESHighPt);
        sources.push_back(METUtil::JESRelNonClosure);
        sources.push_back(METUtil::bJES);
      }
      if (m_addJesPileupSyst)
      {
         sources.push_back(METUtil::JESNPV);
         sources.push_back(METUtil::JESMu);
      }
      sources.push_back(METUtil::JER);
    }
    if(m_doTerm[METUtil::Electrons]) {
      sources.push_back(METUtil::EES);
      sources.push_back(METUtil::EER);
    }
    if(m_doTerm[METUtil::Photons]) {
      sources.push_back(METUtil::PES);
      sources.push_back(METUtil::PER);
    }
    if(m_doTerm[METUtil::Taus]) {
      sources.push_back(METUtil::TES);
      sources.push_back(METUtil::TER);
    }
    if(m_doTerm[METUtil::Muons]) {
      sources.push_back(METUtil::MES);
      sources.push_back(METUtil::MERID);
      sources.push_back(METUtil::MERMS);
    }
    if(m_doTerm[METUtil::CellOut] || m_doTerm[METUtil::CellOutEflow] || m_doTerm[METUtil::SoftJets]) { 
      // if(m_momentaSet[METUtil::Clusters]) {
      // Enable only if case anybody feels like playing with clusters
      // sources.push_back(METUtil::CES);
      // sources.push_back(METUtil::CER);
      // } else {
      sources.push_back(METUtil::ResoSoftTermsSyst);
      sources.push_back(METUtil::ScaleSoftTermsSyst);
      // }
    }
    if(m_doTerm[METUtil::CellOutEflow] && m_momentaSet[METUtil::Tracks]) {
      // Just in case anybody feels like playing with tracks in CellOutEflow
      sources.push_back(METUtil::TrkES);
      sources.push_back(METUtil::TrkER);
    }
  } else if(systtype == METUtil::All_ptHard) {
    if(m_doTerm[METUtil::Jets]) {
      sources.push_back(METUtil::JES);
      sources.push_back(METUtil::JER);
    }
    if(m_doTerm[METUtil::Electrons]) {
      sources.push_back(METUtil::EES);
      sources.push_back(METUtil::EER);
    }
    if(m_doTerm[METUtil::Photons]) {
      sources.push_back(METUtil::PES);
      sources.push_back(METUtil::PER);
    }
    if(m_doTerm[METUtil::Taus]) {
      sources.push_back(METUtil::TES);
      sources.push_back(METUtil::TER);
    }
    if(m_doTerm[METUtil::Muons]) {
      sources.push_back(METUtil::MES);
      sources.push_back(METUtil::MERID);
      sources.push_back(METUtil::MERMS);
    }
    if(m_doTerm[METUtil::CellOut] || m_doTerm[METUtil::CellOutEflow] || m_doTerm[METUtil::SoftJets]) { 
      sources.push_back(METUtil::ResoSoftTermsCo_ptHard);
      sources.push_back(METUtil::ResoSoftTermsAnti_ptHard);
      sources.push_back(METUtil::ScaleSoftTermsSyst_ptHard);
    }
    if(m_doTerm[METUtil::CellOutEflow] && m_momentaSet[METUtil::Tracks]) {
      sources.push_back(METUtil::TrkES);
      sources.push_back(METUtil::TrkER);
    }
  } else if(systtype == METUtil::ResoSoftTermsSyst_ptHard) {
    sources.push_back(METUtil::ResoSoftTermsCo_ptHard);
    sources.push_back(METUtil::ResoSoftTermsAnti_ptHard);
  } else {
    sources.push_back(systtype);
  }

  return deltaMissingET(term, sources);
}

METObject METUtility::deltaMissingET(const METUtil::Terms term, const vector<METUtil::SystType> &sources)
{

  //use commas as delimiters

  Float_t summedSquares_et = 0.0;
  Float_t summedSquares_etx = 0.0;
  Float_t summedSquares_ety = 0.0;
  Float_t summedSquares_sumet = 0.0;
  METObject metTerm;

  for(unsigned int j = 0; j < sources.size(); ++j) {
    METUtil::SystType systo = sources.at(j);
    metTerm = getMissingETDiff(term, systo);
    summedSquares_et += metTerm.et()*metTerm.et();
    summedSquares_etx += metTerm.etx()*metTerm.etx();
    summedSquares_ety += metTerm.ety()*metTerm.ety();
    summedSquares_sumet += metTerm.sumet()*metTerm.sumet();
  }
  summedSquares_et = sqrt(summedSquares_et);
  summedSquares_etx = sqrt(summedSquares_etx);
  summedSquares_ety = sqrt(summedSquares_ety);
  summedSquares_sumet = sqrt(summedSquares_sumet);
  
  METObject deltaMET(summedSquares_etx, summedSquares_ety, summedSquares_sumet);
  deltaMET.setEt(summedSquares_et); //substitutes our et diff for what would be summed in quadrature from etx/ety
  return deltaMET;
}// end of deltaMissingET

METObject METUtility::absDeltaMissingET(METUtil::Terms term, METUtil::SystType systtype)
{
  METObject relMET = deltaMissingET(term, systtype);
  METObject theMET = getMissingET(term);

  Float_t _et = relMET.et()*theMET.et();
  Float_t _etx = relMET.etx()*theMET.etx();
  Float_t _ety = relMET.ety()*theMET.ety();
  Float_t _sumet = relMET.sumet()*theMET.sumet();
  
  METObject absDeltaMET(_etx, _ety, _sumet);
  absDeltaMET.setEt(_et);
  return absDeltaMET;

}// end of absDeltaMissingET

METObject METUtility::getMissingETDiff(METUtil::Terms term, METUtil::SystType systtype)
{
  
  METUtil::Systematics systup = METUtil::None;
  METUtil::Systematics systdown = METUtil::None;
  switch(systtype) {
  case METUtil::EES: systup = METUtil::EESUp; systdown = METUtil::EESDown; break;
  case METUtil::EER: systup = METUtil::EERUp; systdown = METUtil::EERDown; break;
    // Single JES systematic
  case METUtil::JES: systup = METUtil::JESUp; systdown = METUtil::JESDown; break;
    // Uncorrelated JES components
  case METUtil::JES1: systup = METUtil::JES1Up; systdown = METUtil::JES1Down; break;
  case METUtil::JES2: systup = METUtil::JES2Up; systdown = METUtil::JES2Down; break;
  case METUtil::JES3: systup = METUtil::JES3Up; systdown = METUtil::JES3Down; break;
  case METUtil::JES4: systup = METUtil::JES4Up; systdown = METUtil::JES4Down; break;
  case METUtil::JES5: systup = METUtil::JES5Up; systdown = METUtil::JES5Down; break;
  case METUtil::JES6: systup = METUtil::JES6Up; systdown = METUtil::JES6Down; break;
  case METUtil::JESEtaIntercalibration: systup = METUtil::JESEtaIntercalibrationUp; systdown = METUtil::JESEtaIntercalibrationDown; break;
  case METUtil::JESHighPt: systup = METUtil::JESHighPtUp; systdown = METUtil::JESHighPtDown; break;
  case METUtil::JESRelNonClosure: systup = METUtil::JESRelNonClosureUp; systdown = METUtil::JESRelNonClosureDown; break;
  case METUtil::JESNPV: systup = METUtil::JESNPVUp; systdown = METUtil::JESNPVDown; break;
  case METUtil::JESMu: systup = METUtil::JESMuUp; systdown = METUtil::JESMuDown; break;
  case METUtil::bJES: systup = METUtil::bJESUp; systdown = METUtil::bJESDown; break;
    // End JES systematics
  case METUtil::JER: systup = METUtil::JERUp; systdown = METUtil::JERDown; break;
  case METUtil::PES: systup = METUtil::PESUp; systdown = METUtil::PESDown; break;
  case METUtil::PER: systup = METUtil::PERUp; systdown = METUtil::PERDown; break;
  case METUtil::TES: systup = METUtil::TESUp; systdown = METUtil::TESDown; break;
  case METUtil::TER: systup = METUtil::TERUp; systdown = METUtil::TERDown; break;
  case METUtil::MERID: systup = METUtil::MERIDUp; systdown = METUtil::MERIDDown; break;
  case METUtil::MERMS: systup = METUtil::MERMSUp; systdown = METUtil::MERMSDown; break;
  case METUtil::MES: systup = METUtil::MESUp; systdown = METUtil::MESDown; break;
  case METUtil::CES: systup = METUtil::CESUp; systdown = METUtil::CESDown; break;
  case METUtil::CER: systup = METUtil::CERUp; systdown = METUtil::CERDown; break;
  case METUtil::TrkES: systup = METUtil::TrkESUp; systdown = METUtil::TrkESDown; break;
  case METUtil::TrkER: systup = METUtil::TrkERUp; systdown = METUtil::TrkERDown; break;
  case METUtil::ResoSoftTermsCo_ptHard: systup = METUtil::ResoSoftTermsUp_ptHard; systdown = METUtil::ResoSoftTermsDown_ptHard; break;
  case METUtil::ResoSoftTermsAnti_ptHard: systup = METUtil::ResoSoftTermsUpDown_ptHard; systdown = METUtil::ResoSoftTermsDownUp_ptHard; break;
  case METUtil::ScaleSoftTermsSyst_ptHard: systup = METUtil::ScaleSoftTermsUp_ptHard; systdown = METUtil::ScaleSoftTermsDown_ptHard; break;
  case METUtil::ResoSoftTermsSyst: systup = METUtil::ResoSoftTermsUp; systdown = METUtil::ResoSoftTermsDown; break;
  case METUtil::ScaleSoftTermsSyst: systup = METUtil::ScaleSoftTermsUp; systdown = METUtil::ScaleSoftTermsDown; break;
  case METUtil::All:
    cerr << "METUTILITY: Only individual uncertainties should be used with this method. 'METUtil::All' is only valid for deltaMissingET() and absDeltaMissingET()." << endl;
    return METObject(0.,0.,0.);
  default:
    cerr << "METUTILITY: This case " << systtype << " is currently unhandled in getMissingETDiff." << endl;
  } // end of switch

  float _etx = MissingETHelper(term).etx();
  if(_etx != 0) _etx = (fabs(_etx - MissingETHelper(term, systup).etx()) + fabs(_etx - MissingETHelper(term, systdown).etx()))/(2*_etx);
  float _ety = MissingETHelper(term).ety();
  if(_ety != 0) _ety = (fabs(_ety - MissingETHelper(term, systup).ety()) + fabs(_ety - MissingETHelper(term, systdown).ety()))/(2*_ety);
  float _sumet = MissingETHelper(term).sumet();
  if(_sumet != 0) _sumet = (fabs(_sumet - MissingETHelper(term, systup).sumet()) + fabs(_sumet - MissingETHelper(term, systdown).sumet()))/(2*_sumet);
  float _et = MissingETHelper(term).et();
  if(_et != 0) _et = (fabs(_et - MissingETHelper(term, systup).et()) + fabs(_et - MissingETHelper(term, systdown).et()))/(2*_et);

  METObject metDiff(_etx, _ety, _sumet);
  metDiff.setEt(_et); //substitutes our et diff for what would be summed in quadrature from etx/ety
  return metDiff;
  
}// end of getMissingETDiff


//most of the excitment happens here
METObject METUtility::MissingETHelper(METUtil::Terms term, METUtil::MultiSyst systematic)
{

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  if(m_verbose) {
    cout << "In MissingETHelper, with term = " << METUtil::getTermName(term)
	 << " and systematics = ";
    for(int iSyst=METUtil::None; iSyst<=METUtil::ScaleSoftTermsDown_ptHard; iSyst++)
      cout << METUtil::getSystName(METUtil::Systematics(iSyst)) << ", ";
    cout << endl;
  }

  // Default to loop over RefFinal
  METUtil::Terms firstTerm = METUtil::RefEle;
  METUtil::Terms lastTerm = METUtil::CellOutEflow;
  if(term == METUtil::SoftTerms) {
    firstTerm = METUtil::SoftJets;
    lastTerm = METUtil::CellOutEflow;
  } else if(term == METUtil::HardTerms) {
    firstTerm = METUtil::RefEle;
    lastTerm = METUtil::MuonTotal;
  } else if(term == METUtil::Truth) {
    return m_terms[term]; // this has no systematics!
  } else if(term != METUtil::RefFinal) {
    firstTerm = lastTerm = term;
  }

  // this keeps things consistent between SoftJets and CellOut
  // when doing the smearing of the resolution.
  m_rndGaus = gRandom->Gaus(0.,1.);
  // only loop if it's one of the individual, storeable terms
  for(int iTerm = firstTerm; iTerm <= lastTerm && iTerm <= METUtil::CellOutEflow; iTerm++) {
    if(m_doTerm[iTerm]) {

      if(m_verbose) {
	cout << "Current call:" << endl
	     << METUtil::getTermName(iTerm)
	     << "(";
	for(int iSyst=METUtil::None; iSyst<=METUtil::ScaleSoftTermsDown_ptHard; iSyst++)
	  cout << METUtil::getSystName(METUtil::Systematics(iSyst)) << ", ";
	cout << ")" << endl;
      }

      switch(iTerm) {
      case METUtil::RefEle:    m_termsScaled[iTerm] = RefEle(systematic); break;
      case METUtil::RefGamma:  m_termsScaled[iTerm] = RefGamma(systematic); break;
      case METUtil::RefJet:    m_termsScaled[iTerm] = RefJet(systematic); break;
      case METUtil::RefTau:    m_termsScaled[iTerm] = RefTau(systematic); break;
      case METUtil::RefMuon:   m_termsScaled[iTerm] = RefMuon(systematic); break;
      case METUtil::MuonTotal: m_termsScaled[iTerm] = MuonTotal(systematic); break;
      case METUtil::SoftJets:  m_termsScaled[iTerm] = SoftJets(systematic); break;
      case METUtil::CellOut:   m_termsScaled[iTerm] = CellOut(systematic); break;
      case METUtil::CellOutEflow: m_termsScaled[iTerm] = CellOutEflow(systematic); break;
      default:
	cerr << "METUTILITY: reached an illegal place in MissingETHelper! Returned zero MET." << endl;
	return METObject(0.,0.,0.);
      } // end of term switch statement

      _etx += m_termsScaled[iTerm].etx();
      _ety += m_termsScaled[iTerm].ety();

      if(iTerm != METUtil::MuonTotal)
	_sumet += m_termsScaled[iTerm].sumet();
      // Muons only contribute to sumET through RefMuon; it is a calorimetric quantity.
    } // end of doTerm check
  } // end of term loop

  // this term functions rather differently from the rest
  // If the "_ptHard" systematics are requested, use the SoftTerms_ptHard method.
  if( (term == METUtil::SoftTerms || term == METUtil::RefFinal) && 
      (
       systematic.hasSyst(METUtil::ResoSoftTermsUp_ptHard) ||
       systematic.hasSyst(METUtil::ResoSoftTermsUpDown_ptHard) ||
       systematic.hasSyst(METUtil::ResoSoftTermsDownUp_ptHard) ||
       systematic.hasSyst(METUtil::ResoSoftTermsDown_ptHard) ||
       systematic.hasSyst(METUtil::ScaleSoftTermsUp_ptHard) ||
       systematic.hasSyst(METUtil::ScaleSoftTermsDown_ptHard)
       )
      ) {
    METObject _softTerms = SoftTerms_ptHard(systematic);
    _etx += _softTerms.etx();
    _ety += _softTerms.ety();
    _sumet += _softTerms.sumet();
  }

  METObject metSum(_etx, _ety, _sumet);
  //refFinal.setEt(-1.0);
  float siggy = 0.0;
  if(m_doSignificance) {
    siggy = METSignificance(systematic);
  }// end of do sig
  
  metSum.setSignificance(siggy);
 
  return metSum;

}// end of MissingET main function


void METUtility::setObjectsHelper(METUtil::Objects type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, vector<vector<float> > wet, vector<vector<float> > wex, vector<vector<float> > wey, vector<vector<unsigned int> > statusWord)
{
  vector<PhysicsObject> _objectVector;

  // Check that pT and wet have the same size. Expect all kinematics and all weights should have the same length...
  if(pT.size() != wet.size()) {
    cerr << "METUtility: You have supplied pT with " << pT.size() << " elements, and wet with " << wet.size() << " elements." << endl;
    cerr << "METUtility: Objects not set. Tool will return rubbish." << endl;
    return;
  }

  //why do it like this? instead of giving an object multiple weights, I'll just do a different PhysicsObject and have multiple muons, 
  //such that I have combined muons, spectro muons, track muons, and refmuon muons and all have different weights, and the statusWord 
  //tells how to use them. It's sort of how it works when METPerformance uses the composition map at AOD level... 
  //Of course the filler for uncertainties and resolutions has to be adjusted to match...
  for(unsigned int i = 0; i < wet.size(); ++i) {
    for(unsigned int j = 0; j < wet.at(i).size(); ++j) {
      PhysicsObject _object;// = new PhysicsObject;
      _object.setMomenergy(pT.at(i), eta.at(i), phi.at(i), E.at(i));
      float _wex = wex.at(i).at(j);
      float _wey = wey.at(i).at(j);
      float _wet = wet.at(i).at(j);
      unsigned int _statusWord = statusWord.at(i).at(j);
      _object.setWeights(_wex, _wey, _wet);
      _object.setStatusCode(_statusWord);
      _object.setIndex(i);
      _objectVector.push_back(_object);
    }// end of inner weight loop
  }// end of outer weight loop
  
  if(type == METUtil::Electrons) {
    m_objects[METUtil::Electrons].clear();
    m_objects[METUtil::Electrons] = _objectVector;
    m_momentaSet[METUtil::Electrons] = true;
  }
  else if(type == METUtil::Photons) {
    m_objects[METUtil::Photons].clear();
    m_objects[METUtil::Photons] = _objectVector;
    m_momentaSet[METUtil::Photons] = true;
  }
  else if(type == METUtil::Taus) {
    m_objects[METUtil::Taus].clear();
    m_objects[METUtil::Taus] = _objectVector;
    m_momentaSet[METUtil::Taus] = true;
  }
  else if(type == METUtil::Jets) {
    m_objects[METUtil::Jets].clear();
    m_objects[METUtil::Jets] = _objectVector;
    m_momentaSet[METUtil::Jets] = true;
  }
  else if(type == METUtil::Muons) {
    m_objects[METUtil::Muons].clear();
    m_objects[METUtil::Muons] = _objectVector;
    m_momentaSet[METUtil::Muons] = true;
  }
  else if(type == METUtil::Clusters) {
    m_objects[METUtil::Clusters].clear();
    m_objects[METUtil::Clusters] = _objectVector;
    m_momentaSet[METUtil::Clusters] = true;
  }
  else if(type == METUtil::Tracks) {
    m_objects[METUtil::Tracks].clear();
    m_objects[METUtil::Tracks] = _objectVector;
    m_momentaSet[METUtil::Tracks] = true;
  } 

}// end of setobjecthelper


void METUtility::setObjectsHelper(METUtil::Objects type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, vector<float> wet, vector<float> wex, vector<float> wey, vector<unsigned short> statusWord)
{

  vector<vector<float> > newWex;// = new vector<vector<float> >;
  vector<vector<float> > newWey;// = new vector<vector<float> >;
  vector<vector<float> > newWet;// = new vector<vector<float> >;
  vector<vector<unsigned int> > newStatusWord;// = new vector<vector<unsigned int> >;

  for(unsigned int i = 0; i < statusWord.size(); ++i) {

    vector<float> wexHold;
    vector<float> weyHold;
    vector<float> wetHold;
    vector<unsigned int> statusWordHold;

    wexHold.push_back(wex.at(i));
    weyHold.push_back(wey.at(i));
    wetHold.push_back(wet.at(i));
    statusWordHold.push_back(static_cast<unsigned int>(statusWord.at(i)));

    newWex.push_back(wexHold);
    newWey.push_back(weyHold);
    newWet.push_back(wetHold);
    newStatusWord.push_back(statusWordHold);
   
  }// end of loop

  setObjectsHelper(type, (pT), (eta), (phi), (E), newWet, newWex, newWey, newStatusWord);

}// end of setObjectshelper 


void METUtility::setObjects(METUtil::Objects type, const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{

  setObjectsHelper(type, (*pT), (*eta), (*phi), (*E), (*wet), (*wex), (*wey), (*statusWord));
}//adapter to use vector<float> instead of vector<vector<float> >


void METUtility::setObjects(METUtil::Objects type, const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  
  setObjectsHelper(type, (*pT), (*eta), (*phi), (*E), (*wet), (*wex), (*wey), (*statusWord));

}// end of setObjects


void METUtility::setObjects(METUtil::Objects type, const vector<TLorentzVector> &momenta, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  vector<float> pT, eta, phi, E;
  for(vector<TLorentzVector>::const_iterator iobj=momenta.begin(); iobj != momenta.end(); iobj++) {
    pT.push_back(iobj->Pt());
    eta.push_back(iobj->Eta());
    phi.push_back(iobj->Phi());
    E.push_back(iobj->E());
  }
  setObjectsHelper(type, pT, eta, phi, E, (*wet), (*wex), (*wey), (*statusWord));

}// end of setObjects


void METUtility::setObjectMomenta(METUtil::Objects type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E)
{

  if(type == METUtil::SpectroMuons) { // For these we need to set the secondary momentum
    setObjectMomentaHelper(m_objects[METUtil::Muons], pT, eta, phi, E, true);
  } else if(type == METUtil::OriJets) { // For these we need to set the secondary momentum
    setObjectMomentaHelper(m_objects[METUtil::Jets], pT, eta, phi, E, true);
    m_oriJetsSet = true;
  } else if(type == METUtil::MuonsComboID || METUtil::MuonsComboMS) { // No need to set momenta for these.
    cerr << "METUTILITY: No need to set ComboID or ComboMS momenta separately! Nothing done." << endl;
  } else { // Everything else goes here
    setObjectMomentaHelper(m_objects[type], pT, eta, phi, E);
  }

}// end of setMomenta


void METUtility::setObjectMomentaHelper(vector<PhysicsObject> &object, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, bool isSecondary)
{

  for(unsigned int i = 0; i < object.size(); ++i) {
    object.at(i).setMomenergy(pT.at(object.at(i).index()), eta.at(object.at(i).index()), phi.at(object.at(i).index()), E.at(object.at(i).index()), isSecondary);
  }// end of object loops

}// end of setMomentaHelper


void METUtility::setObjectEnergyUncertainties(METUtil::Objects type, const vector<float> energyUp, const vector<float> energyDown, EnergyShift::ShiftType iShift)
{

  if(!m_momentaSet[type] && type <= METUtil::Tracks) {
    cerr << "METUTILITY: Set object four-momentum and weights before setting uncertainties!" << endl;
    return;
  }

  if(type == METUtil::SpectroMuons || type == METUtil::MuonsComboMS || type == METUtil::MuonsComboID) {
    cerr << "METUTILITY: Separate energy scale uncertainties are not defined for different muon types/components. Just use METUtil::Muons." << endl;
  } else {
    setObjectEnergyUncertaintiesHelper(m_objects[type], energyUp, energyDown, iShift);
    m_scaleShiftsSet[type] = true;
  }
  
}// end of setObjectUncertainties


void METUtility::setObjectEnergyUncertaintiesHelper(vector<PhysicsObject> &object, vector<float> energyUp, vector<float> energyDown, EnergyShift::ShiftType iShift)
{
  //don't assume number of objects matches number of uncertainties (muons and clusters won't), use index to match
  for(unsigned int i = 0; i < object.size(); ++i) {
    object.at(i).setEnergyUncertainty(energyUp.at(object.at(i).index()), energyDown.at(object.at(i).index()), iShift);
  }// end of object loops

}// end of setObjectUncertaintiesHelper


void METUtility::setObjectResolutionShift(METUtil::Objects type, const vector<float> resUp, const vector<float> resDown)
{

  if(!m_momentaSet[type] && type <= METUtil::Tracks) {
    cerr << "METUTILITY: Set object four-momentum and weights before setting resolution shifts!" << endl;
    return;
  }

  unsigned short restype;
  METUtil::Objects objtype;
  switch(type) {
  case METUtil::SpectroMuons:
    objtype = METUtil::Muons;
    restype = 2; break;
  case METUtil::MuonsComboMS:
    objtype = METUtil::Muons;
    restype = 1; break;
  case METUtil::MuonsComboID:
    objtype = METUtil::Muons;
    restype = 0; break;
  default:
    objtype = type;
    restype = 0;
  }
  setObjectResolutionShiftHelper(m_objects[objtype], resUp, resDown, restype);
  m_resShiftsSet[type] = true;
}// end of setObjectResolutions


void METUtility::setObjectResolutionShiftHelper(vector<PhysicsObject> &object, vector<float> resUp, vector<float> resDown, unsigned short type)
{
  
  //don't assume number of objects matches number of resolutions (muons and clusters won't), use index to match
  for(unsigned int i = 0; i < object.size(); ++i) {
    object.at(i).setResShift(resUp.at(object.at(i).index()), resDown.at(object.at(i).index()), type);
  }// end of object loops

}// end of setObjectUncertaintiesHelper
 

void METUtility::setObjectResolutions(METUtil::Objects type, const vector<float> resolutions)
{

  if(!m_momentaSet[type] && type <= METUtil::Tracks) {
    cerr << "METUTILITY: Set object four-momentum and weights before setting resolutions!" << endl;
    return;
  }
  
  if(type == METUtil::SpectroMuons || 
     type == METUtil::MuonsComboMS || 
     type == METUtil::MuonsComboID ||
     type == METUtil::OriJets) {
    cerr << "METUTILITY: Resolutions just pertain to the basic object types, e.g. Muons, Jets." << endl;
  }
  else {
    setObjectResolutionsHelper(m_objects[type], resolutions);
    m_resolutionsSet[type] = true;
  }
  
}// end of setObjectResolutions


void METUtility::setObjectResolutionsHelper(vector<PhysicsObject> &objects, vector<float> resolutions)
{

  for(unsigned int i = 0; i < objects.size(); ++i) {
    objects.at(i).setResolution(resolutions.at(objects.at(i).index()));
  }// end of object loops
  
}// end of setObjectResolutionsHelper


void METUtility::setMETTerm(METUtil::Terms term, float _etx, float _ety, float _sumet)
{

  m_terms[term].setEtx(_etx);
  m_terms[term].setEty(_ety);
  m_terms[term].setSumet(_sumet);
  m_terms[term].setIsValid(true);

}


void METUtility::defineMissingET(bool doRefEle, bool doRefGamma, bool doRefTau, bool doRefJet, bool doSoftJets, bool doRefMuon, bool doMuonTotal, bool doCellOut, bool doCellOutEflow)
{

  m_doTerm[METUtil::RefEle] = doRefEle;
  m_doTerm[METUtil::RefGamma] = doRefGamma;
  m_doTerm[METUtil::RefTau] = doRefTau;
  m_doTerm[METUtil::RefJet] = doRefJet;
  m_doTerm[METUtil::SoftJets] = doSoftJets;
  m_doTerm[METUtil::RefMuon] = doRefMuon;
  m_doTerm[METUtil::MuonTotal] = doMuonTotal;
  m_doTerm[METUtil::CellOut] = doCellOut;
  m_doTerm[METUtil::CellOutEflow] = doCellOutEflow;

  if(m_doTerm[METUtil::SoftJets] == false) m_softJetCut = 0.0; // so smeared jets don't get lost from the RefJet term

}


void METUtility::configMissingET(bool is2012, bool isSTVF)
{
  
  if(is2012 && isSTVF) {
    m_softTermsResUnc = 0.01;
    m_softTermsScaleUnc = 0.06;
  } else if(!is2012 && isSTVF) {
    m_softTermsResUnc = 0.077;
    m_softTermsScaleUnc = 0.066;
  } else if(is2012 && !isSTVF) {
    m_softTermsResUnc = 0.025;
    m_softTermsScaleUnc = 0.08;
  } else if(!is2012 && !isSTVF) {
    m_softTermsResUnc = 0.02;
    m_softTermsScaleUnc = 0.05;
  }

  m_is2012 = is2012;
  m_isSTVF = isSTVF;
}


METObject METUtility::RefEle(METUtil::MultiSyst systematic)
{

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_objects[METUtil::Electrons].size(); ++i) {
    scale = 1.0;
    PhysicsObject thisElectron = m_objects[METUtil::Electrons].at(i);
   
    if(thisElectron.Pt()>0 && (fabs(thisElectron.eta()) < m_etaCut[0] || fabs(thisElectron.eta()) > m_etaCut[1])) continue;
 
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::EERUp) && systematic.hasSyst(METUtil::EERDown)) {
      cerr << "METUtility: You have requested both EERUp and EERDown systematics -- please don't!" << endl;
    } else if(systematic.hasSyst(METUtil::EESUp) && systematic.hasSyst(METUtil::EESDown)) {
      cerr << "METUtility: You have requested both EESUp and EESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::EERUp) && m_resShiftsSet[METUtil::Electrons]) scale += thisElectron.resShift().first;
      else if(systematic.hasSyst(METUtil::EERDown) && m_resShiftsSet[METUtil::Electrons]) scale += thisElectron.resShift().second;
      if(systematic.hasSyst(METUtil::EESUp) && m_scaleShiftsSet[METUtil::Electrons]) scale += thisElectron.energyShift().first;
      else if(systematic.hasSyst(METUtil::EESDown) && m_scaleShiftsSet[METUtil::Electrons]) scale += thisElectron.energyShift().second;
    }
    
    _etx -= thisElectron.Px()*thisElectron.wex()*scale;
    _ety -= thisElectron.Py()*thisElectron.wey()*scale;
    _sumet += thisElectron.Pt()*thisElectron.wet()*scale;
   
  }// end of loop
 
  if(m_momentaSet[METUtil::Electrons] == false && m_terms[METUtil::RefEle].isValid()) {
    _etx += m_terms[METUtil::RefEle].etx();
    _ety += m_terms[METUtil::RefEle].ety();
    _sumet += m_terms[METUtil::RefEle].sumet();
  }

  METObject _refEle(_etx, _ety, _sumet);
  return _refEle;
}// end of RefEle term remaker


METObject METUtility::RefGamma(METUtil::MultiSyst systematic)
{

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_objects[METUtil::Photons].size(); ++i) {
    scale = 1.0;
    PhysicsObject thisPhoton = m_objects[METUtil::Photons].at(i);
  
    if(thisPhoton.Pt()>0 && (fabs(thisPhoton.eta()) < m_etaCut[0] || fabs(thisPhoton.eta()) > m_etaCut[1])) continue;
 
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::PERUp) && systematic.hasSyst(METUtil::PERDown)) {
      cerr << "METUtility: You have requested both PERUp and PERDown systematics -- please don't!" << endl;
    } else if(systematic.hasSyst(METUtil::PESUp) && systematic.hasSyst(METUtil::PESDown)) {
      cerr << "METUtility: You have requested both PESUp and PESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::PERUp) && m_resShiftsSet[METUtil::Photons]) scale += thisPhoton.resShift().first;
      else if(systematic.hasSyst(METUtil::PERDown) && m_resShiftsSet[METUtil::Photons]) scale += thisPhoton.resShift().second;
      if(systematic.hasSyst(METUtil::PESUp) && m_scaleShiftsSet[METUtil::Photons]) scale += thisPhoton.energyShift().first;
      else if(systematic.hasSyst(METUtil::PESDown) && m_scaleShiftsSet[METUtil::Photons]) scale += thisPhoton.energyShift().second;
    }
   
    _etx -= thisPhoton.Px()*thisPhoton.wex()*scale;
    _ety -= thisPhoton.Py()*thisPhoton.wey()*scale;
    _sumet += thisPhoton.Pt()*thisPhoton.wet()*scale;

  }// end of loop

  if(m_momentaSet[METUtil::Photons] == false  && m_terms[METUtil::RefGamma].isValid()) {
    _etx += m_terms[METUtil::RefGamma].etx();
    _ety += m_terms[METUtil::RefGamma].ety();
    _sumet += m_terms[METUtil::RefGamma].sumet();
  }

  METObject _refGamma(_etx, _ety, _sumet);
  return _refGamma;
}// end of RefGamma


METObject METUtility::RefTau(METUtil::MultiSyst systematic)
{

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_objects[METUtil::Taus].size(); ++i) {
    scale = 1.0;
    PhysicsObject thisTau = m_objects[METUtil::Taus].at(i);

    if(thisTau.Pt()>0 && (fabs(thisTau.eta()) < m_etaCut[0] || fabs(thisTau.eta()) > m_etaCut[1])) continue;
 
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::TERUp) && systematic.hasSyst(METUtil::TERDown)) {
      cerr << "METUtility: You have requested both TERUp and TERDown systematics -- please don't!" << endl;
    } else if(systematic.hasSyst(METUtil::TESUp) && systematic.hasSyst(METUtil::TESDown)) {
      cerr << "METUtility: You have requested both TESUp and TESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::TERUp) && m_resShiftsSet[METUtil::Taus]) scale += thisTau.resShift().first;
      else if(systematic.hasSyst(METUtil::TERDown) && m_resShiftsSet[METUtil::Taus]) scale += thisTau.resShift().second;
      if(systematic.hasSyst(METUtil::TESUp) && m_scaleShiftsSet[METUtil::Taus]) scale += thisTau.energyShift().first;
      else if(systematic.hasSyst(METUtil::TESDown) && m_scaleShiftsSet[METUtil::Taus]) scale += thisTau.energyShift().second;
    }

    _etx -= thisTau.Px()*thisTau.wex()*scale;
    _ety -= thisTau.Py()*thisTau.wey()*scale;
    _sumet += thisTau.Pt()*thisTau.wet()*scale;

  }// end of loop

  if(m_momentaSet[METUtil::Taus] == false  && m_terms[METUtil::RefTau].isValid()) {
    _etx += m_terms[METUtil::RefTau].etx();
    _ety += m_terms[METUtil::RefTau].ety();
    _sumet += m_terms[METUtil::RefTau].sumet();
  }

  METObject _refTau(_etx, _ety, _sumet);
  return _refTau;
}// end of RefTau


METObject METUtility::RefJet(METUtil::MultiSyst systematic)
{

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  if(m_momentaSet[METUtil::Jets] && !m_oriJetsSet)
    cerr << "METUTILITY: Warning -- you have not set the original jet pTs. RefJet will return 0." << endl;

  for(unsigned int i = 0; i < m_objects[METUtil::Jets].size(); ++i) {
    PhysicsObject thisJet = m_objects[METUtil::Jets].at(i);

    if(thisJet.Pt()>0 && (fabs(thisJet.eta()) < m_etaCut[0] || fabs(thisJet.eta()) > m_etaCut[1])) continue;
    if(!MissingETTags::usesCODE(MissingETTags::DEFAULT,thisJet.statusWord())) continue;
    if(!MissingETTags::isCODE(m_jet_pu_statusWord,thisJet.statusWord())) continue;

    float scale = 1.0;
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::JERUp) && systematic.hasSyst(METUtil::JERDown)) {
      cerr << "METUtility: You have requested both JERUp and JERDown systematics -- please don't!" << endl;
    } else if(systematic.hasSyst(METUtil::JESUp) && systematic.hasSyst(METUtil::JESDown)) {
      cerr << "METUtility: You have requested both JESUp and JESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::JERUp) && m_resShiftsSet[METUtil::Jets]) scale += thisJet.resShift().first;
      else if(systematic.hasSyst(METUtil::JERDown) && m_resShiftsSet[METUtil::Jets]) scale += thisJet.resShift().second;
      if(systematic.hasSyst(METUtil::JESUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift().first;
      else if(systematic.hasSyst(METUtil::JESDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift().second;
      // Expanded JES systematics
      if(systematic.hasSyst(METUtil::bJESUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::bJES).first;
      else if(systematic.hasSyst(METUtil::bJESDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::bJES).second;
      if(systematic.hasSyst(METUtil::JES1Up) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES1).first;
      else if(systematic.hasSyst(METUtil::JES1Down) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES1).second;
      if(systematic.hasSyst(METUtil::JES2Up) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES2).first;
      else if(systematic.hasSyst(METUtil::JES2Down) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES2).second;
      if(systematic.hasSyst(METUtil::JES3Up) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES3).first;
      else if(systematic.hasSyst(METUtil::JES3Down) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES3).second;
      if(systematic.hasSyst(METUtil::JES4Up) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES4).first;
      else if(systematic.hasSyst(METUtil::JES4Down) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES4).second;
      if(systematic.hasSyst(METUtil::JES5Up) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES5).first;
      else if(systematic.hasSyst(METUtil::JES5Down) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES5).second;
      if(systematic.hasSyst(METUtil::JES6Up) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES6).first;
      else if(systematic.hasSyst(METUtil::JES6Down) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JES6).second;
      if(systematic.hasSyst(METUtil::JESEtaIntercalibrationUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESEtaIntercalibration).first;
      else if(systematic.hasSyst(METUtil::JESEtaIntercalibrationDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESEtaIntercalibration).second;
      if(systematic.hasSyst(METUtil::JESHighPtUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESHighPt).first;
      else if(systematic.hasSyst(METUtil::JESHighPtDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESHighPt).second;
      if(systematic.hasSyst(METUtil::JESRelNonClosureUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESRelNonClosure).first;
      else if(systematic.hasSyst(METUtil::JESRelNonClosureDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESRelNonClosure).second;
      if(systematic.hasSyst(METUtil::JESNPVUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESNPV).first;
      else if(systematic.hasSyst(METUtil::JESNPVDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESNPV).second;
      if(systematic.hasSyst(METUtil::JESMuUp) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESMu).first;
      else if(systematic.hasSyst(METUtil::JESMuDown) && m_scaleShiftsSet[METUtil::Jets]) scale += thisJet.energyShift(0,EnergyShift::JESMu).second;
    }

    bool addjet = thisJet.Pt(true) > m_softJetCut && thisJet.wet()>0;

    if(addjet) {
//       cout << "Add jet with orig pt " << thisJet.Pt(true) << ", syst is " << systematic << endl;
      _etx   -= thisJet.Px()*thisJet.wex()*scale;
      _ety   -= thisJet.Py()*thisJet.wey()*scale;
      _sumet += thisJet.Pt()*thisJet.wet()*scale;
    }// end of softjet cut
  }// end of loop


  if(m_momentaSet[METUtil::Jets] == false  && m_terms[METUtil::RefJet].isValid()) {
    _etx += m_terms[METUtil::RefJet].etx();
    _ety += m_terms[METUtil::RefJet].ety();
    _sumet += m_terms[METUtil::RefJet].sumet();
  }

  METObject _refJet(_etx, _ety, _sumet);
  return _refJet;
}// end of jet term remaker


METObject METUtility::RefMuon(METUtil::MultiSyst systematic)
{
  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;
  
  if(m_terms[METUtil::RefMuon].isValid()) {
    _etx += m_terms[METUtil::RefMuon].etx();
    _ety += m_terms[METUtil::RefMuon].ety();
    _sumet += m_terms[METUtil::RefMuon].sumet();
  } else {
    float scale = 1.0;
    for(unsigned int i = 0; i < m_objects[METUtil::Muons].size(); ++i) {
      scale = 1.0;
      PhysicsObject thisMuon = m_objects[METUtil::Muons].at(i);

      if(systematic == METUtil::None) scale = 1.0;
      if(thisMuon.Pt()>0 && (fabs(thisMuon.eta()) < m_etaCut[0] || fabs(thisMuon.eta()) > m_etaCut[1])) continue;
 
      if(MissingETTags::usesCODE(MissingETTags::REFMUON,thisMuon.statusWord())) {
	_etx -= thisMuon.Px()*thisMuon.wex()*scale;
	_ety -= thisMuon.Py()*thisMuon.wey()*scale;
	_sumet += thisMuon.Pt()*thisMuon.wet()*scale;
      }// end of statusWord check

    }// end of loop
  }

  METObject _refMuon(_etx, _ety, _sumet);
  return _refMuon;
}// end of refmuon


METObject METUtility::MuonTotal(METUtil::MultiSyst systematic)
{

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_objects[METUtil::Muons].size(); ++i) {
    scale = 1.0;
    PhysicsObject thisMuon = m_objects[METUtil::Muons].at(i);

    if(MissingETTags::usesCODE(MissingETTags::REFMUON,thisMuon.statusWord())) continue;
    if(!m_isMuid && MissingETTags::usesCODE(MissingETTags::MUID,thisMuon.statusWord())) continue;
    if(m_isMuid && !(MissingETTags::usesCODE(MissingETTags::MUID,thisMuon.statusWord()))) continue;

    if(MissingETTags::usesCODE(MissingETTags::DEFAULT,thisMuon.statusWord()) ||
       MissingETTags::usesCODE(MissingETTags::TRACK,thisMuon.statusWord())) {

      if(systematic == METUtil::None) scale = 1.0;
      else if(systematic.hasSyst(METUtil::MERIDUp) && systematic.hasSyst(METUtil::MERIDDown)) {
	cerr << "METUtility: You have requested both MERIDUp and MERIDDown systematics -- please don't!" << endl;
      } else if(systematic.hasSyst(METUtil::MESUp) && systematic.hasSyst(METUtil::MESDown)) {
	cerr << "METUtility: You have requested both MESUp and MESDown systematics -- please don't!" << endl;
      } else {
	if(systematic.hasSyst(METUtil::MERIDUp) && m_resShiftsSet[METUtil::MuonsComboID])
	  scale += thisMuon.resShift(0).first;
	else if(systematic.hasSyst(METUtil::MERIDDown) && m_resShiftsSet[METUtil::MuonsComboID])
	  scale += thisMuon.resShift(0).second;
	if(systematic.hasSyst(METUtil::MERMSUp) && m_resShiftsSet[METUtil::MuonsComboMS])
	  scale += thisMuon.resShift(1).first;
	else if(systematic.hasSyst(METUtil::MERMSDown) && m_resShiftsSet[METUtil::MuonsComboMS])
	  scale += thisMuon.resShift(1).second;
	if(systematic.hasSyst(METUtil::MESUp) && m_scaleShiftsSet[METUtil::Muons])
	  scale += thisMuon.energyShift(0).first;
	else if(systematic.hasSyst(METUtil::MESDown) && m_scaleShiftsSet[METUtil::Muons])
	  scale += thisMuon.energyShift(0).second;
      }

      _etx -= thisMuon.Px()*thisMuon.wex()*scale;
      _ety -= thisMuon.Py()*thisMuon.wey()*scale;
      _sumet += thisMuon.Pt()*thisMuon.wet()*scale;
       
    }// end of default muons
    else if(MissingETTags::usesCODE(MissingETTags::SPECTRO,thisMuon.statusWord())) {
     
      if(systematic == METUtil::None) scale = 1.0;
      else if(systematic.hasSyst(METUtil::MERMSUp) && systematic.hasSyst(METUtil::MERMSDown)) {
	cerr << "METUtility: You have requested both MERMSUp and MERMSDown systematics -- please don't!" << endl;
      }  else if(systematic.hasSyst(METUtil::MESUp) && systematic.hasSyst(METUtil::MESDown)) {
	cerr << "METUtility: You have requested both MESUp and MESDown systematics -- please don't!" << endl;
      } else {
	if(systematic.hasSyst(METUtil::MERMSUp) && m_resShiftsSet[METUtil::SpectroMuons])
	  scale += thisMuon.resShift(2).first;
	else if(systematic.hasSyst(METUtil::MERMSDown) && m_resShiftsSet[METUtil::SpectroMuons])
	  scale += thisMuon.resShift(2).second;
	if(systematic.hasSyst(METUtil::MESUp) && m_scaleShiftsSet[METUtil::Muons])
	  scale += thisMuon.energyShift().first;
	else if(systematic.hasSyst(METUtil::MESDown) && m_scaleShiftsSet[METUtil::Muons])
	  scale += thisMuon.energyShift().second;
      }

      _etx -= thisMuon.Px(true)*thisMuon.wex()*scale;
      _ety -= thisMuon.Py(true)*thisMuon.wey()*scale;
      _sumet += thisMuon.Pt(true)*thisMuon.wet()*scale;
    }// end of spectro
    // else if(MissingETTags::usesCODE(MissingETTags::TRACK,thisMuon.statusWord())) {
    //        if(systematic == METUtil::None) scale = 1.0;
    //        else if(systematic.hasSyst(METUtil::MERIDUp) && systematic.hasSyst(METUtil::MERIDDown)) {
    //          cerr << "METUtility: You have requested both MERIDUp and MERIDDown systematics -- please don't!" << endl;
    //        } else if(systematic.hasSyst(METUtil::MESUp) && systematic.hasSyst(METUtil::MESDown)) {
    //          cerr << "METUtility: You have requested both MESUp and MESDown systematics -- please don't!" << endl;
    //        } else {
    //          if(systematic.hasSyst(METUtil::MERIDUp) && m_resShiftsSet[METUtil::MuonsComboID]) scale += m_objects[METUtil::Muon.resShift(0).first;
    //          else if(systematic.hasSyst(METUtil::MERIDDown) && m_resShiftsSet[METUtil::MuonsComboID]) scale += m_objects[METUtil::Muon.resShift(0).second;
    //          if(systematic.hasSyst(METUtil::MERMSUp) && m_resShiftsSet[METUtil::MuonsComboMS]) scale += m_objects[METUtil::Muon.resShift(1).first;
    //          else if(systematic.hasSyst(METUtil::MERMSDown) && m_resShiftsSet[METUtil::MuonsComboMS]) scale += m_objects[METUtil::Muon.resShift(1).second;
    //          if(systematic.hasSyst(METUtil::MESUp) && m_scaleShiftsSet[METUtil::Muons]) scale += m_objects[METUtil::Muon.energyShift().first;
    //          else if(systematic.hasSyst(METUtil::MESDown) && m_scaleShiftsSet[METUtil::Muons]]) scale += m_objects[METUtil::Muon.energyShift().second;
    //         }
       
    //        //this isn't any really different from default, but I keep them separate
    //       _etx -= thisMuon.Px()*thisMuon.wex()*scale;
    //       _ety -= thisMuon.Py()*thisMuon.wey()*scale;
    //       _sumet += thisMuon.Pt()*thisMuon.wet()*scale;
    //    }// end of track
    
  }// end of loop

  if(m_momentaSet[METUtil::Muons] == false && m_terms[METUtil::MuonTotal].isValid()) {
    _etx += m_terms[METUtil::MuonTotal].etx();
    _ety += m_terms[METUtil::MuonTotal].ety();
    _sumet += m_terms[METUtil::MuonTotal].sumet();
  }
  
  METObject _MuonTotal(_etx, _ety, _sumet);
  return _MuonTotal;
}// end of Muons


METObject METUtility::CellOut(METUtil::MultiSyst systematic)
{

  // If the "_ptHard" systematics are requested, don't compute.
  if( systematic.hasSyst(METUtil::ResoSoftTermsUp_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsUpDown_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsDownUp_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsDown_ptHard) ||
      systematic.hasSyst(METUtil::ScaleSoftTermsUp_ptHard) ||
      systematic.hasSyst(METUtil::ScaleSoftTermsDown_ptHard)
      ) return METObject(0.,0.,0.);

  float scale = 1.0;

  //   if(!m_momentaSet[METUtil::Clusters]) {
  //     if(m_verbose) cerr << "Clusters momenta and weights not loaded" << endl;
  //     return m_cellOut;
  //   }

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  short clusterEnergyCase = 0;
  short clusterResCase = 0;
  if(m_useStandardClusterRes) clusterResCase = 3;
  if(m_useStandardClusterEnergySigma) clusterEnergyCase = 3;

  bool objectsyst = false;

  for(unsigned int i = 0; i < m_objects[METUtil::Clusters].size(); ++i) {
    PhysicsObject thisCluster = m_objects[METUtil::Clusters].at(i);

    if(MissingETTags::usesCODE(MissingETTags::EFLOW_CLUSTER,thisCluster.statusWord())) {continue;}
    if(thisCluster.Pt()>0 && (fabs(thisCluster.eta()) < m_etaCut[0] || fabs(thisCluster.eta()) > m_etaCut[1])) continue;
 
    scale = 1.0;
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::CERUp) && systematic.hasSyst(METUtil::CERDown)) {
      cerr << "METUtility: You have requested both CERUp and CERDown systematics -- please don't!" << endl;
    } else if(systematic.hasSyst(METUtil::CESUp) && systematic.hasSyst(METUtil::CESDown)) {
      cerr << "METUtility: You have requested both CESUp and CESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::CERUp) && m_resShiftsSet[METUtil::Clusters]) scale += thisCluster.resShift(clusterResCase).first;
      else if(systematic.hasSyst(METUtil::CERDown) && m_resShiftsSet[METUtil::Clusters]) scale += thisCluster.resShift(clusterResCase).second;
      if(systematic.hasSyst(METUtil::CESUp) && m_scaleShiftsSet[METUtil::Clusters]) scale += thisCluster.energyShift(clusterEnergyCase).first;
      else if(systematic.hasSyst(METUtil::CESDown) && m_scaleShiftsSet[METUtil::Clusters]) scale += thisCluster.energyShift(clusterEnergyCase).second;
    }
    
    _etx -= thisCluster.Px()*thisCluster.wex()*scale;
    _ety -= thisCluster.Py()*thisCluster.wey()*scale;
    _sumet += thisCluster.Pt()*thisCluster.wet()*scale;

    objectsyst = true;
  }// end of loop

  // Done, no need to do overall systematics
  if(objectsyst) {
    METObject _cellOutEflow(_etx, _ety, _sumet);
    return _cellOutEflow;
  }
  
  if(m_momentaSet[METUtil::Clusters] == false && m_terms[METUtil::CellOut].isValid()) {
    _etx += m_terms[METUtil::CellOut].etx();
    _ety += m_terms[METUtil::CellOut].ety();
    _sumet += m_terms[METUtil::CellOut].sumet();
  }

  // Derived for CellOut_Eflow, but retained so CellOut has something
  if(systematic.hasSyst(METUtil::ScaleSoftTermsUp) && systematic.hasSyst(METUtil::ScaleSoftTermsDown)) {
    cerr << "METUtility: You have requested ScaleSoftTermsUp and ScaleSoftTermsDown systematics -- please don't!" << endl;
  } else if(systematic.hasSyst(METUtil::ScaleSoftTermsUp)) scale += m_softTermsScaleUnc;
  else if(systematic.hasSyst(METUtil::ScaleSoftTermsDown)) scale -= m_softTermsScaleUnc;
  // end of Milano SoftTerms scale systematic (April 2012)

  if(systematic.hasSyst(METUtil::ResoSoftTermsUp) && systematic.hasSyst(METUtil::ResoSoftTermsDown)) {
    cerr << "METUtility: You have requested ResoSoftTermsUp and ResoSoftTermsDown systematics -- please don't!" << endl;
  } else if(systematic.hasSyst(METUtil::ResoSoftTermsUp) || systematic.hasSyst(METUtil::ResoSoftTermsDown)) {
    METObject _SoftJets;
    if(m_doTerm[METUtil::SoftJets]) _SoftJets = SoftJets(METUtil::None);
    METObject _SoftTerms(_etx+_SoftJets.etx(), _ety+_SoftJets.ety(), _sumet+_SoftJets.sumet());
    float sumET_SoftTerms = fabs(_SoftTerms.sumet()/1000);   // in GeV
    double Met_sigma = 0.7 * sqrt(sumET_SoftTerms);    // now get the uncertainty
    double Met_smearing_sigma =  sqrt(pow(Met_sigma * (1 + m_softTermsResUnc),2) - pow(Met_sigma,2));

    // this keeps things consistent between SoftJets and CellOut
    float shift = _SoftTerms.et()==0. ? 0. : m_rndGaus*Met_smearing_sigma * 1000 / _SoftTerms.et();
    if(_SoftTerms.et()==0. )
      cerr << "METUTILITY: Tried to get SoftTerms uncertainty but SoftTerms et == 0!" << endl;
    if(systematic.hasSyst(METUtil::ResoSoftTermsUp)) scale += shift;
    else if(systematic.hasSyst(METUtil::ResoSoftTermsDown)) scale -= shift;
  }// end of Milano SoftTerms resolution systematic (April 2012)

  METObject _cellOut(_etx*scale, _ety*scale, _sumet);
  return _cellOut;

}// end of CellOut


METObject METUtility::CellOutEflow(METUtil::MultiSyst systematic)
{

  // If the "_ptHard" systematics are requested, don't compute.
  if( systematic.hasSyst(METUtil::ResoSoftTermsUp_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsUpDown_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsDownUp_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsDown_ptHard) ||
      systematic.hasSyst(METUtil::ScaleSoftTermsUp_ptHard) ||
      systematic.hasSyst(METUtil::ScaleSoftTermsDown_ptHard)
      ) return METObject(0.,0.,0.);

  float scale = 1.0;

  //if(!m_momentaSet[METUtil::Clusters] || !m_trackmomentaSet) {
  // if(m_verbose) cerr << "Clusters momenta and weights not loaded" << endl;
  //return m_cellOutEflow;
  //}

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  short clusterEnergyCase = 0;
  short clusterResCase = 0;
  if(m_useStandardClusterRes) clusterResCase = 3;
  if(m_useStandardClusterEnergySigma) clusterEnergyCase = 3;

  bool objectsyst = false;

  for(unsigned int i = 0; i < m_objects[METUtil::Clusters].size(); ++i) {
    PhysicsObject thisCluster = m_objects[METUtil::Clusters].at(i);

    if(!(MissingETTags::usesCODE(MissingETTags::EFLOW_CLUSTER,m_objects[METUtil::Clusters].at(i).statusWord()))) {continue;}

    if(fabs(thisCluster.eta()) < m_etaCut[0] || fabs(thisCluster.eta()) > m_etaCut[1]) continue;

    scale = 1.0;
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::CERUp) && systematic.hasSyst(METUtil::CERDown)) {
	cerr << "METUtility: You have requested both CERUp and CERDown systematics -- please don't!" << endl;
    }  else if(systematic.hasSyst(METUtil::CESUp) && systematic.hasSyst(METUtil::CESDown)) {
	cerr << "METUtility: You have requested both CESUp and CESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::CERUp) && m_resShiftsSet[METUtil::Clusters]) scale += thisCluster.resShift(clusterResCase).first;
      else if(systematic.hasSyst(METUtil::CERDown) && m_resShiftsSet[METUtil::Clusters]) scale += thisCluster.resShift(clusterResCase).second;
      if(systematic.hasSyst(METUtil::CESUp) && m_scaleShiftsSet[METUtil::Clusters]) scale += thisCluster.energyShift(clusterEnergyCase).first;
      else if(systematic.hasSyst(METUtil::CESDown) && m_scaleShiftsSet[METUtil::Clusters]) scale += thisCluster.energyShift(clusterEnergyCase).second;
    }
    
    _etx -= thisCluster.Px()*thisCluster.wex()*scale;
    _ety -= thisCluster.Py()*thisCluster.wey()*scale;
    _sumet += thisCluster.Pt()*thisCluster.wet()*scale;
    
    objectsyst = true;
  }// end of loop

  for(unsigned int i = 0; i < m_objects[METUtil::Tracks].size(); ++i) {
    PhysicsObject thisTrack = m_objects[METUtil::Tracks].at(i);

    scale = 1.0;
    if(systematic == METUtil::None) scale = 1.0;
    else if(systematic.hasSyst(METUtil::TrkERUp) && systematic.hasSyst(METUtil::TrkERDown)) {
	cerr << "METUtility: You have requested both TrkERUp and TrkERDown systematics -- please don't!" << endl;
    } else if(systematic.hasSyst(METUtil::TrkESUp) && systematic.hasSyst(METUtil::TrkESDown)) {
	cerr << "METUtility: You have requested both TrkESUp and TrkESDown systematics -- please don't!" << endl;
    } else {
      if(systematic.hasSyst(METUtil::TrkERUp) && m_resShiftsSet[METUtil::Tracks]) scale += thisTrack.resShift().first;
      else if(systematic.hasSyst(METUtil::TrkERDown) && m_resShiftsSet[METUtil::Tracks]) scale += thisTrack.resShift().second;
      if(systematic.hasSyst(METUtil::TrkESUp) && m_scaleShiftsSet[METUtil::Tracks]) scale += thisTrack.energyShift().first;
      else if(systematic.hasSyst(METUtil::TrkESDown) && m_scaleShiftsSet[METUtil::Tracks]) scale += thisTrack.energyShift().second;
    }
    _etx -= thisTrack.Px()*thisTrack.wex()*scale;
    _ety -= thisTrack.Py()*thisTrack.wey()*scale;
    _sumet += thisTrack.Pt()*thisTrack.wet()*scale;

    objectsyst = true;
  }// end of track loop

  // Done, no need to do overall systematics
  if(objectsyst) {
    METObject _cellOutEflow(_etx, _ety, _sumet);
    return _cellOutEflow;
  }

  if((m_momentaSet[METUtil::Tracks] == false || m_momentaSet[METUtil::Clusters] == false) && m_terms[METUtil::CellOutEflow].isValid()) {
    if(m_verbose) cerr << "Clusters momenta and weights not loaded" << endl;
    _etx = m_terms[METUtil::CellOutEflow].etx();
    _ety = m_terms[METUtil::CellOutEflow].ety();
    _sumet = m_terms[METUtil::CellOutEflow].sumet();
  }

  if(systematic.hasSyst(METUtil::ScaleSoftTermsUp) && systematic.hasSyst(METUtil::ScaleSoftTermsDown)) {
    cerr << "METUtility: You have requested ScaleSoftTermsUp and ScaleSoftTermsDown systematics -- please don't!" << endl;
  } else if(systematic.hasSyst(METUtil::ScaleSoftTermsUp)) scale += m_softTermsScaleUnc;
  else if(systematic.hasSyst(METUtil::ScaleSoftTermsDown)) scale -= m_softTermsScaleUnc;
  // end of Milano SoftTerms scale systematic (April 2012)

  if(systematic.hasSyst(METUtil::ResoSoftTermsUp) && systematic.hasSyst(METUtil::ResoSoftTermsDown)) {
    cerr << "METUtility: You have requested ResoSoftTermsUp and ResoSoftTermsDown systematics -- please don't!" << endl;
  } else if(systematic.hasSyst(METUtil::ResoSoftTermsUp) || systematic.hasSyst(METUtil::ResoSoftTermsDown)) {
    METObject _SoftJets;
    if(m_doTerm[METUtil::SoftJets]) _SoftJets = SoftJets(METUtil::None);
    METObject _SoftTerms(_etx+_SoftJets.etx(), _ety+_SoftJets.ety(), _sumet+_SoftJets.sumet());
    float sumET_SoftTerms = fabs(_SoftTerms.sumet()/1000);   // in GeV
    double Met_sigma = 0.7 * sqrt(sumET_SoftTerms);    // now get the uncertainty
    double Met_smearing_sigma =  sqrt(pow(Met_sigma * (1 + m_softTermsResUnc),2) - pow(Met_sigma,2));

    // this keeps things consistent between SoftJets and CellOut
    float shift =_SoftTerms.et()==0. ? 0. : m_rndGaus*Met_smearing_sigma * 1000 / _SoftTerms.et();
    if(_SoftTerms.et()==0. )
      cerr << "METUTILITY: Tried to get SoftTerms uncertainty but SoftTerms et == 0!" << endl;
    if(systematic.hasSyst(METUtil::ResoSoftTermsUp)) scale += shift;
    else if(systematic.hasSyst(METUtil::ResoSoftTermsDown)) scale -= shift;
    // _sumet *= scale;
  }// end of Milano SoftTerms resolution systematic (April 2012)

  METObject _cellOutEflow(_etx*scale, _ety*scale, _sumet);
  return _cellOutEflow;
  
}// end of CellOutEflow


METObject METUtility::SoftJets(METUtil::MultiSyst systematic)
{

  // If the "_ptHard" systematics are requested, don't compute.
  if( systematic.hasSyst(METUtil::ResoSoftTermsUp_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsUpDown_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsDownUp_ptHard) ||
      systematic.hasSyst(METUtil::ResoSoftTermsDown_ptHard) ||
      systematic.hasSyst(METUtil::ScaleSoftTermsUp_ptHard) ||
      systematic.hasSyst(METUtil::ScaleSoftTermsDown_ptHard)
      ) return METObject(0.,0.,0.);

  float scale = 1.0;

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  if( m_terms[METUtil::SoftJets].isValid()) {
    _etx = m_terms[METUtil::SoftJets].etx();
    _ety = m_terms[METUtil::SoftJets].ety();
    _sumet = m_terms[METUtil::SoftJets].sumet();
  }

  if(systematic.hasSyst(METUtil::ScaleSoftTermsUp) && systematic.hasSyst(METUtil::ScaleSoftTermsDown)) {
    cerr << "METUtility: You have requested ScaleSoftTermsUp and ScaleSoftTermsDown systematics -- please don't!" << endl;
  } else if(systematic.hasSyst(METUtil::ScaleSoftTermsUp)) scale += m_softTermsScaleUnc;
  else if(systematic.hasSyst(METUtil::ScaleSoftTermsDown)) scale -= m_softTermsScaleUnc;
  // end of Milano SoftTerms scale systematic (April 2012)

  if(systematic.hasSyst(METUtil::ResoSoftTermsUp) && systematic.hasSyst(METUtil::ResoSoftTermsDown)) {
    cerr << "METUtility: You have requested ResoSoftTermsUp and ResoSoftTermsDown systematics -- please don't!" << endl;
  } else if(systematic.hasSyst(METUtil::ResoSoftTermsUp) || systematic.hasSyst(METUtil::ResoSoftTermsDown)) {
    METObject _CellOut;
    if(m_doTerm[METUtil::CellOutEflow]) _CellOut = CellOutEflow(METUtil::None);
    else if(m_doTerm[METUtil::CellOut]) _CellOut = CellOut(METUtil::None);
    METObject _SoftTerms(_etx+_CellOut.etx(), _ety+_CellOut.ety(), _sumet+_CellOut.sumet());
    float sumET_SoftTerms = fabs(_SoftTerms.sumet()/1000);   // in GeV
    double Met_sigma = 0.7 * sqrt(sumET_SoftTerms);    // now get the uncertainty
    double Met_smearing_sigma =  sqrt(pow(Met_sigma * (1 + m_softTermsResUnc),2) - pow(Met_sigma,2));

    // this keeps things consistent between SoftJets and CellOut
    float shift = _SoftTerms.et()==0. ? 0. : m_rndGaus*Met_smearing_sigma * 1000 / _SoftTerms.et();
    if(_SoftTerms.et()==0. )
      cerr << "METUTILITY: Tried to get SoftTerms uncertainty but SoftTerms et == 0!" << endl;
    if(systematic.hasSyst(METUtil::ResoSoftTermsUp)) scale += shift;
    else if(systematic.hasSyst(METUtil::ResoSoftTermsDown)) scale -= shift;
  }// end of Milano SoftTerms resolution systematic (April 2012)
  
  METObject _softJets(_etx*scale, _ety*scale, _sumet);
  return _softJets;

}// end of SoftJets


METObject METUtility::SoftTerms_ptHard(METUtil::MultiSyst systematic)
{
  
  float _etx = 0;
  float _ety = 0;
  METObject _softTerms = MissingETHelper(METUtil::SoftTerms,METUtil::None);
  float _sumet = _softTerms.sumet();
  TVector2 ptHard= calculate_ptHard();
  if(ptHard.Mod()<1e-9) return _softTerms;

  double ux = ptHard.Px()/ptHard.Mod();
  double uy = ptHard.Py()/ptHard.Mod();
  double celloutl = _softTerms.etx()*ux + _softTerms.ety()*uy;
  double celloutt = _softTerms.etx()*uy - _softTerms.ety()*ux;

  double meanlmc =  0.0801085 + 0.507246*ptHard.Mod();
  if ( meanlmc > 52) meanlmc = 52;

  if(!m_is2012 && !m_isSTVF) {
    meanlmc = 0.518935 + 0.648222*ptHard.Mod();
    if (meanlmc>80) meanlmc = 80;
  }

  meanlmc = meanlmc*1000.;
  double reso_scalel = 1.;
  double reso_scalet = 1.;

  double newcelloutl = reso_scalel*(celloutl-meanlmc)+meanlmc;
  double newcelloutt = reso_scalet*(celloutt);

  // This is not used right now... comment out to avoid build warnings
  //  double ratiosigmalvspthard_nvx1_bin[3] = {30.,60.,200.};

  // Set up the master pointers, which we'll reassign to the appropriate parameters
  double zeroes[3] = {0.,0.,0.};
  double ones[3] = {1.,1.,1.};
  double* diffmeanlvspthard_mu1 = zeroes;
  double* diffmeanlvspthard_mu2 = zeroes;
  double* diffmeanlvspthard_mu3 = zeroes;
  double* ratiosigmalvspthard_mu1 = ones;
  double* ratiosigmalvspthard_mu2 = ones;
  double* ratiosigmalvspthard_mu3 = ones;
  double* ratiosigmatvspthard_mu1 = ones;
  double* ratiosigmatvspthard_mu2 = ones;
  double* ratiosigmatvspthard_mu3 = ones;

  // 2012, no STVF
  double diffmeanlvspthard_mu1_2012noSTVF[3] = {0.16197,0.596695,0.747181};
  double diffmeanlvspthard_mu2_2012noSTVF[3] = {0.133905,0.40287,0.429138};
  double diffmeanlvspthard_mu3_2012noSTVF[3] = {0.207049,0.18041,-0.175014};
  //
  double ratiosigmalvspthard_mu1_2012noSTVF[3] = {0.970996,0.96426,1.00098};
  double ratiosigmalvspthard_mu2_2012noSTVF[3] = {0.965591,0.976173,0.978536};
  double ratiosigmalvspthard_mu3_2012noSTVF[3] = {0.9635,0.967231,0.994153 };
  //
  double ratiosigmatvspthard_mu1_2012noSTVF[3] = {0.967715,0.957746,1.01353 };
  double ratiosigmatvspthard_mu2_2012noSTVF[3] = {0.966874,0.974869,0.999589};
  double ratiosigmatvspthard_mu3_2012noSTVF[3] = {0.963548,0.975741,0.982291};

  // 2012, STVF
  double diffmeanlvspthard_mu1_2012STVF[3] = {-0.00308496,0.629223,1.12903};
  double diffmeanlvspthard_mu2_2012STVF[3] = {0.104673,0.649614,1.12842};
  double diffmeanlvspthard_mu3_2012STVF[3] = {0.133939,0.574716,1.01711};
  //
  double ratiosigmalvspthard_mu1_2012STVF[3] ={0.93269,0.967472,1.08598 };
  double ratiosigmalvspthard_mu2_2012STVF[3] ={0.933112,1.00706,1.19131  };
  double ratiosigmalvspthard_mu3_2012STVF[3] ={0.975512,1.09078,1.14746 };
  //
  double ratiosigmatvspthard_mu1_2012STVF[3] ={0.93576,0.985547,0.945715 };
  double ratiosigmatvspthard_mu2_2012STVF[3] ={0.93754,0.994635,0.973336 };
  double ratiosigmatvspthard_mu3_2012STVF[3] ={0.974219,0.987653,1.11326 };

  if(m_is2012 && !m_isSTVF) {
    diffmeanlvspthard_mu1 = diffmeanlvspthard_mu1_2012noSTVF;
    diffmeanlvspthard_mu2 = diffmeanlvspthard_mu2_2012noSTVF;
    diffmeanlvspthard_mu3 = diffmeanlvspthard_mu3_2012noSTVF;
    ratiosigmalvspthard_mu1 = ratiosigmalvspthard_mu1_2012noSTVF;
    ratiosigmalvspthard_mu2 = ratiosigmalvspthard_mu2_2012noSTVF;
    ratiosigmalvspthard_mu3 = ratiosigmalvspthard_mu3_2012noSTVF;
    ratiosigmatvspthard_mu1 = ratiosigmatvspthard_mu1_2012noSTVF;
    ratiosigmatvspthard_mu2 = ratiosigmatvspthard_mu2_2012noSTVF;
    ratiosigmatvspthard_mu3 = ratiosigmatvspthard_mu3_2012noSTVF;
  } else if(m_is2012 && m_isSTVF) {
    diffmeanlvspthard_mu1 = diffmeanlvspthard_mu1_2012STVF;
    diffmeanlvspthard_mu2 = diffmeanlvspthard_mu2_2012STVF;
    diffmeanlvspthard_mu3 = diffmeanlvspthard_mu3_2012STVF;
    ratiosigmalvspthard_mu1 = ratiosigmalvspthard_mu1_2012STVF;
    ratiosigmalvspthard_mu2 = ratiosigmalvspthard_mu2_2012STVF;
    ratiosigmalvspthard_mu3 = ratiosigmalvspthard_mu3_2012STVF;
    ratiosigmatvspthard_mu1 = ratiosigmatvspthard_mu1_2012STVF;
    ratiosigmatvspthard_mu2 = ratiosigmatvspthard_mu2_2012STVF;
    ratiosigmatvspthard_mu3 = ratiosigmatvspthard_mu3_2012STVF;
  }

  //generate the numbers
  double rsig_l_pthard =-1;
  double rsig_t_pthard =-1;
  double dmean_l_pthard=-1;

  if(m_is2012 || m_isSTVF) {

    //at the moment all same binning ... calculate the right bin
    int the_index=-1;

    if (ptHard.Mod()<30.) the_index = 0;
    if (ptHard.Mod()>=30. && ptHard.Mod()<60. ) the_index = 1;
    if (ptHard.Mod()>=60. ) the_index = 2;

    if (m_averageIntPerXing<16){
      rsig_l_pthard  = ratiosigmalvspthard_mu1[the_index];
      rsig_t_pthard  = ratiosigmatvspthard_mu1[the_index];
      dmean_l_pthard = diffmeanlvspthard_mu1[the_index];
      if(m_is2012 && m_isSTVF) {
	meanlmc =   0.0456348 + 0.0676502 *ptHard.Mod();
	if (meanlmc>2) meanlmc = 2;
      }
    }else if (m_averageIntPerXing>=16&&m_averageIntPerXing<22){
      rsig_l_pthard  = ratiosigmalvspthard_mu2[the_index];
      rsig_t_pthard  = ratiosigmatvspthard_mu2[the_index];
      dmean_l_pthard = diffmeanlvspthard_mu2[the_index];
      if(m_is2012 && m_isSTVF) {
	meanlmc = 0.117416 +    0.0421375*ptHard.Mod();
	if (meanlmc>1.2) meanlmc = 1.2;
      }
    }else if (m_averageIntPerXing>=22){
      rsig_l_pthard  = ratiosigmalvspthard_mu3[the_index];
      rsig_t_pthard  = ratiosigmatvspthard_mu3[the_index];
      dmean_l_pthard = diffmeanlvspthard_mu3[the_index];
      if(m_is2012 && m_isSTVF) {
	meanlmc =      0.137743 +  0.0295313*ptHard.Mod();
	if (meanlmc>0.9) meanlmc = 0.9;
      }
    }
  }

  double systl_topo = fabs(rsig_l_pthard-1);
  double systt_topo = fabs(rsig_t_pthard-1);

  if(!m_is2012 && !m_isSTVF) {
    double systl_temp = fabs(0.984+0.001356*ptHard.Mod()-1);
    double systt_temp = fabs(0.9858+0.00117*ptHard.Mod()-1);
    double pileup_paral = abs(0.97 + 0.0053*m_nvtx-1); 
    double pileup_parat = abs(0.97 + 0.0048*m_nvtx-1);
    systl_topo = pow(0.02*0.02+pileup_paral*pileup_paral+systl_temp*systl_temp,0.5);
    systt_topo = pow(0.02*0.02+pileup_parat*pileup_parat+systt_temp*systt_temp,0.5);
    dmean_l_pthard = 1.;
  }

  if(systematic.hasSyst(METUtil::ResoSoftTermsUp_ptHard)) {

    reso_scalel = 1+systl_topo; 
    reso_scalet = 1+systt_topo;

    newcelloutl = reso_scalel*(celloutl-meanlmc)+meanlmc;
    newcelloutt = reso_scalet*(celloutt);
    _etx = (ux*newcelloutl + uy*newcelloutt);
    _ety = (uy*newcelloutl - ux*newcelloutt);

  } else if(systematic.hasSyst(METUtil::ResoSoftTermsDown_ptHard)) {

    reso_scalel = 1-systl_topo; 
    reso_scalet = 1-systt_topo;

    newcelloutl = reso_scalel*(celloutl-meanlmc)+meanlmc;
    newcelloutt = reso_scalet*(celloutt);
    _etx = (ux*newcelloutl + uy*newcelloutt);
    _ety = (uy*newcelloutl - ux*newcelloutt);    

  } else if(systematic.hasSyst(METUtil::ResoSoftTermsUpDown_ptHard)) {

    reso_scalel = 1+systl_topo; 
    reso_scalet = 1-systt_topo;

    newcelloutl = reso_scalel*(celloutl-meanlmc)+meanlmc;
    newcelloutt = reso_scalet*(celloutt);
    _etx = (ux*newcelloutl + uy*newcelloutt);
    _ety = (uy*newcelloutl - ux*newcelloutt);

  }  else if(systematic.hasSyst(METUtil::ResoSoftTermsDownUp_ptHard)) {

    reso_scalel = 1-systl_topo; 
    reso_scalet = 1+systt_topo;

    newcelloutl = reso_scalel*(celloutl-meanlmc)+meanlmc;
    newcelloutt = reso_scalet*(celloutt);
    _etx = (ux*newcelloutl + uy*newcelloutt);
    _ety = (uy*newcelloutl - ux*newcelloutt);    

  } else if(systematic.hasSyst(METUtil::ScaleSoftTermsUp_ptHard)) {

    reso_scalel = 1.;
    reso_scalet = 1.;
    double shift = dmean_l_pthard*1000.;
    newcelloutl = reso_scalel*(celloutl-meanlmc) + meanlmc + shift;
    newcelloutt = reso_scalet*(celloutt);
    _etx = (ux*newcelloutl + uy*newcelloutt);
    _ety = (uy*newcelloutl - ux*newcelloutt);    

  } else if(systematic.hasSyst(METUtil::ScaleSoftTermsDown_ptHard)) {

    reso_scalel = 1.;
    reso_scalet = 1.;
    double shift = -1*dmean_l_pthard*1000.;
    newcelloutl = reso_scalel*(celloutl-meanlmc) + meanlmc + shift;
    newcelloutt = reso_scalet*(celloutt);
    _etx = (ux*newcelloutl + uy*newcelloutt);
    _ety = (uy*newcelloutl - ux*newcelloutt);    

  } else {
    cerr << "METUTILITY: Requested SoftTerms_ptHard with an incompatible or no systematic! This returns nonsense." << endl;
  }

  METObject _softTerms_ptHard(_etx, _ety, _sumet);
  return _softTerms_ptHard;
}

TVector2 METUtility::calculate_ptHard()
{

  TVector2 ptHard;

  METObject _hardTerms = MissingETHelper(METUtil::HardTerms,METUtil::None);

  double dirx = -1. * _hardTerms.etx();
  double diry = -1. * _hardTerms.ety();
  
  double dirx_mettruth = m_terms[METUtil::Truth].etx();
  double diry_mettruth = m_terms[METUtil::Truth].ety();

  double dirx_tot = dirx + dirx_mettruth;
  double diry_tot = diry + diry_mettruth;
  
  ptHard.Set(dirx_tot/1000.,diry_tot/1000.);
  return ptHard;
}


float METUtility::METSignificance(METUtil::MultiSyst systematic)
{

  float denominator = 0;
  float shift = 0.0;
  float deltaPhi = 0;
  if(m_doSignificance) {
    
    if(m_doTerm[METUtil::RefEle]) {
      if(m_resolutionsSet[METUtil::Electrons]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Electrons].size(); ++i) {
	  shift = 0.0;
	  if(systematic == METUtil::None) shift = 0.0;
	  else if(systematic.hasSyst(METUtil::EERUp) && systematic.hasSyst(METUtil::EERDown)) {
	    cerr << "METUtility: You have requested both EERUp and EERDown systematics -- please don't!" << endl;
	  } else {
	    if(systematic.hasSyst(METUtil::EERUp) && m_resShiftsSet[METUtil::Electrons]) shift = m_objects[METUtil::Electrons].at(i).resShift().first;
	    else if(systematic.hasSyst(METUtil::EERDown) && m_resShiftsSet[METUtil::Electrons]) shift = m_objects[METUtil::Electrons].at(i).resShift().second;
	  }
	  deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Electrons].at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_objects[METUtil::Electrons].at(i).resolution()*(m_objects[METUtil::Electrons].at(i).resolution() + shift);
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::RefEle].sumet();
      }
    }// end of if doRefEle


    if(m_doTerm[METUtil::RefGamma]) {
      if(m_resolutionsSet[METUtil::Photons]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Photons].size(); ++i) {
	  shift = 0.0;
	  if(systematic == METUtil::None) shift = 0.0;
	 
	  deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Photons].at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_objects[METUtil::Photons].at(i).resolution()*(m_objects[METUtil::Photons].at(i).resolution() + shift);
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::RefGamma].sumet();
      }
    }// end of if doRefGamma

    if(m_doTerm[METUtil::RefTau]) {
      if(m_resolutionsSet[METUtil::Taus]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Taus].size(); ++i) {
	  shift = 0.0;
	  if(systematic == METUtil::None) shift = 0.0;
	  
	  deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Taus].at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_objects[METUtil::Taus].at(i).resolution()*(m_objects[METUtil::Taus].at(i).resolution() + shift);
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::RefTau].sumet();
      }
    }// end of if doRefTau

    if(m_doTerm[METUtil::RefJet]) {
      if(m_resolutionsSet[METUtil::Jets]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Jets].size(); ++i) {
	  if(m_objects[METUtil::Jets].at(i).Pt() > m_softJetCut) {
	    shift = 0.0;
	    if(systematic == METUtil::None) shift = 0.0;
	    
	    deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Jets].at(i).phi();
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	    deltaPhi = TMath::Cos(deltaPhi);
	    
	    denominator += deltaPhi*deltaPhi*m_objects[METUtil::Jets].at(i).resolution()*(m_objects[METUtil::Jets].at(i).resolution() + shift);
	  }//soft jet check	
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::RefJet].sumet();
      }
    }// end of if doRefJet
    
    if(m_doTerm[METUtil::MuonTotal]) {
      if(m_resolutionsSet[METUtil::Muons]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Muons].size(); ++i) {
	  shift = 0.0;

	  if(MissingETTags::usesCODE(MissingETTags::DEFAULT,m_objects[METUtil::Muons].at(i).statusWord())) {
	    if(systematic == METUtil::None) shift = 0.0;
	      
	  }// end of default muons
	  else if(MissingETTags::usesCODE(MissingETTags::SPECTRO,m_objects[METUtil::Muons].at(i).statusWord())) {
	    if(systematic == METUtil::None) shift = 0.0;
	     
	  }// end of spectro
	  else if(MissingETTags::usesCODE(MissingETTags::TRACK,m_objects[METUtil::Muons].at(i).statusWord())) {
	    if(systematic == METUtil::None) shift = 0.0;
	      
	  }// end of spectro
	    
	  deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Muons].at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	    
	  denominator += deltaPhi*deltaPhi*m_objects[METUtil::Muons].at(i).resolution()*(m_objects[METUtil::Muons].at(i).resolution() + shift);
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::MuonTotal].sumet();
      }
    }// end of if doMuonTotal
    
    if(m_doTerm[METUtil::RefMuon]) denominator += 0.5*0.5*m_termsScaled[METUtil::RefMuon].sumet();

    if(m_doTerm[METUtil::SoftJets]) {
      if(m_resolutionsSet[METUtil::Jets]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Jets].size(); ++i) {
	  if(m_objects[METUtil::Jets].at(i).Pt() <= m_softJetCut) {
	    shift = 0.0;
	    if(systematic == METUtil::None) shift = 0.0;
	    
	    deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Jets].at(i).phi();
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	    deltaPhi = TMath::Cos(deltaPhi);
	    
	    denominator += deltaPhi*deltaPhi*m_objects[METUtil::Jets].at(i).resolution()*(m_objects[METUtil::Jets].at(i).resolution() + shift);
	  }//soft jet check	
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::SoftJets].sumet();
      }
    }// end of if doSoftJets
    
    if(m_doTerm[METUtil::CellOut]) {
      if(m_resolutionsSet[METUtil::Clusters]) {
	for(unsigned int i = 0; i < m_objects[METUtil::Clusters].size(); ++i) {
	  if(MissingETTags::usesCODE(MissingETTags::EFLOW_CLUSTER,m_objects[METUtil::Clusters].at(i).statusWord())) {continue;}
	  
	  short clusterResCase = 0;
	  if(m_useStandardClusterRes) clusterResCase = 3;
	  
	  shift = 0.0;
	  if(systematic == METUtil::None) shift = 0.0;
	  
	  deltaPhi = m_terms[METUtil::RefFinal].phi() - m_objects[METUtil::Clusters].at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_objects[METUtil::Clusters].at(i).resolution()*(m_objects[METUtil::Clusters].at(i).resolution() + shift);
	}// end of loop
      }// end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_termsScaled[METUtil::CellOut].sumet();
      }
    }// end of if doCellOut

  }// end of if doSignificance

  return m_terms[METUtil::RefFinal].et()/TMath::Sqrt(denominator);
}// end of METSignificance


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Further helper function to deal with the fact that not every Performance group dumps the same variables///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void METUtility::setJetParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  setObjects(METUtil::Jets, pT, eta, phi, E, wet, wex, wey, statusWord);
}// end of setJetParameters


void METUtility::setElectronParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  vector<float> E;// = new const vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.51099891);
    E.push_back(_hlv.E());
  }
  setObjectsHelper(METUtil::Electrons, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setElectronParameters


void METUtility::setPhotonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.0);
    E.push_back(_hlv.E());
  }
  setObjectsHelper(METUtil::Photons, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setPhotonParameters


void METUtility::setTauParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{

  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 1776.8);
    E.push_back(_hlv.E());
  }
  
  setObjectsHelper(METUtil::Taus, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setTauParameters


void METUtility::setClusterParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  setObjects(METUtil::Clusters, pT, eta, phi, E, wet, wex, wey, statusWord);
}// end of setClusterParameters


void METUtility::setTrackParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{
  setObjects(METUtil::Tracks, pT, eta, phi, E, wet, wex, wey, statusWord);
}// end of setTrackParameters


void METUtility::setMuonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord)
{

  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 105.658367);
    E.push_back(_hlv.E());
  }
  
  setObjectsHelper(METUtil::Muons, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setMuonParameters


void METUtility::setExtraMuonParameters(const vector<float> *qOverPSpectro, const vector<float> *thetaSpectro, const vector<float> *phiSpectro, const vector<float> *charge)
{

  vector<float> mu_staco_ms_eta;// = new vector<float>;
  vector<float> mu_staco_ms_pt;// = new vector<float>;
  vector<float> mu_staco_ms_E;// = new vector<float>;
  
  for(unsigned int iMu = 0; iMu < phiSpectro->size(); ++iMu) {
    float p = 0.0;
    if(qOverPSpectro->at(iMu) != 0.0) p = fabs(charge->at(iMu)/qOverPSpectro->at(iMu));
    
    float E = TMath::Sqrt(105.658367*105.658367 + p*p);
    mu_staco_ms_E.push_back(E);
    
    float p_T = p*TMath::Sin(thetaSpectro->at(iMu));
    mu_staco_ms_pt.push_back(p_T);
    
    float eta = -1.0*TMath::Log(TMath::Tan(thetaSpectro->at(iMu))/2.0);
    mu_staco_ms_eta.push_back(eta);
  }// end of spectro muon loop

  setObjectMomenta(METUtil::SpectroMuons, mu_staco_ms_pt, mu_staco_ms_eta, (*phiSpectro), mu_staco_ms_E);
  
}// end of setExtraMuonParameters


void METUtility::setExtraMuonParameters(const vector<float> *mu_staco_ms_pt, const vector<float> *thetaSpectro, const vector<float> *phiSpectro)
{

  vector<float> mu_staco_ms_eta;// = new vector<float>;
  vector<float> mu_staco_ms_E;// = new vector<float>;
    
  float E = 0.0;
  float eta = 0.0;
   
  for(unsigned int iMu = 0; iMu < phiSpectro->size(); ++iMu) {
    float p = 0.0;
    if(TMath::Sin(thetaSpectro->at(iMu)) != 0) p = mu_staco_ms_pt->at(iMu)/TMath::Sin(thetaSpectro->at(iMu));  
    
    E = TMath::Sqrt(105.658367*105.658367 + p*p);
    mu_staco_ms_E.push_back(E);
    
    //float eta = -1.0*TMath::Log(TMath::Tan(thetaSpectro->at(iMu))/2.0);
    eta = -1.0*TMath::Log(fabs(TMath::Tan(thetaSpectro->at(iMu))/2.0))*TMath::Tan(thetaSpectro->at(iMu))/fabs(TMath::Tan(thetaSpectro->at(iMu)));
    mu_staco_ms_eta.push_back(eta);
    
  }// end of spectro muon loop

  setObjectMomenta(METUtil::SpectroMuons, (*mu_staco_ms_pt), mu_staco_ms_eta, (*phiSpectro), mu_staco_ms_E);
 
}// end of setExtraMuonParameters


void METUtility::setOriJetParameters(const vector<float> *pT)
{
  
  if(!m_momentaSet[METUtil::Jets])
    cerr << "Set jet momenta first!" << endl;
  vector<PhysicsObject>& jets = m_objects[METUtil::Jets];
  for(unsigned int i = 0; i < jets.size(); ++i) {
    PhysicsObject& jet = jets.at(i);
    float mass = jet.mom().M();
    TLorentzVector orivec;
    orivec.SetPtEtaPhiM(pT->at(jet.index()), jet.eta(), jet.phi(), mass);
    jet.setMomenergy(orivec, true);
  }// end of jet loop

  m_oriJetsSet = true;

}// end of setOriJetParameters


// void METUtility::setExtraJetParameters(const vector<float> *moment, const vector<float> */*mass*/, const vector<float> *eta, const vector<float> *phi)
// {

//   vector<float> lc_pt;
//   vector<float> lc_E;
   
//   float scale = 1.0;
//   float oldE = 0.0;
//   float oldPt = 0.0;
//   for(unsigned int i = 0; i < moment->size(); ++i) {
//     scale = 1.0;
//     if(moment->at(i) > 0) scale = 1.0/moment->at(i);
//     //scale = 1.0;
//     for(unsigned int j = 0; j < m_objects[METUtil::Jets].size(); ++j) {
//       if(i == m_objects[METUtil::Jets].at(j).index()) {
// 	oldE = m_objects[METUtil::Jets].at(j).E();
// 	oldPt = m_objects[METUtil::Jets].at(j).Pt();
// 	break;
//       }//index matching if
//     }//m_objects[METUtil::Jets] loop

//     float newE = oldE*scale;
//     //     float newMass = mass->at(i)*scale;
//     //     float theta = 2.0*tanh(exp(-1.0*eta->at(i)));
//     float newPt = oldPt*scale;//sqrt(newE*newE - newMass*newMass)*sin(theta);
          
//     lc_pt.push_back(newPt);
//     lc_E.push_back(newE);
     
//   }// end of loop

//   setObjectMomenta(METUtil::OriJets, lc_pt, (*eta), (*phi),lc_E);
   
// }// end of function
 

void METUtility::setJetParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{
  setObjects(METUtil::Jets, pT, eta, phi, E, wet, wex, wey, statusWord);
}// end of setJetParameters


void METUtility::setElectronParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.51099891);
    E.push_back(_hlv.E());
  }
  setObjectsHelper(METUtil::Electrons, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setElectronParameters


void METUtility::setPhotonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.0);
    E.push_back(_hlv.E());
  }
  setObjectsHelper(METUtil::Photons, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setPhotonParameters


void METUtility::setTauParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{

  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 1776.8);
    E.push_back(_hlv.E());
  }
  
  setObjectsHelper(METUtil::Taus, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setTauParameters


void METUtility::setClusterParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{
  setObjects(METUtil::Clusters, pT, eta, phi, E, wet, wex, wey, statusWord);
}// end of setClusterParameters


void METUtility::setTrackParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{
  setObjects(METUtil::Tracks, pT, eta, phi, E, wet, wex, wey, statusWord);
}// end of setTrackParameters


void METUtility::setMuonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord)
{
  
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 105.658367);
    E.push_back(_hlv.E());
  }
  setObjectsHelper(METUtil::Muons, (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}// end of setJetParameters


bool METUtility::checkConsistency(METObject met_test, const int metterm, METObject &diff)
{

  METObject met_util = getMissingET(METUtil::Terms(metterm));
  diff.setBase(met_util.etx() - met_test.etx(),
	       met_util.ety() - met_test.ety(),
	       met_util.sumet() - met_test.sumet());
  METObject sum(met_util.etx() + met_test.etx(),
		met_util.ety() + met_test.ety(),
		met_util.sumet() + met_test.sumet());

  if(diff.et()>1.0e3 && diff.et() / sum.et() > 1e-9) {
    if(m_verbose)
      cout << "CONSISTENCY CHECK: " << METUtil::getTermName(metterm) << " et differs by " << diff.et() << " MeV = "
	   << diff.et() / sum.et() * 100 << "%!" << endl;
    return false;
  }
  if(diff.sumet()>1.0e3 && diff.sumet() / sum.sumet() > 1e-9) {
    if(m_verbose)
      cout << "CONSISTENCY CHECK: " << METUtil::getTermName(metterm) << " sumet differs by " << diff.sumet() << " MeV = "
	   << diff.et() / sum.et() * 100 << "%!" << endl;
    return false;
  }

  return true;
}// end of checkConsistency


bool METUtility::checkConsistency(const vector<pair<int,METObject> > & testvector, METObject &diff)
{

  for(unsigned int i = 0; i<testvector.size(); i++) {
    METUtil::Terms metterm = METUtil::Terms(testvector.at(i).first);
    METObject met_test = testvector.at(i).second;
    bool result = checkConsistency(met_test, metterm, diff);
    if(!result) return result;
  }

  return true;
}// end of checkConsistency
