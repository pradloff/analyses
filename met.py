from common.functions import event_function
from common.particle import particle
from common.external import load
import os

class correct_missing_energy(event_function):

	def __init__(self):
		event_function.__init__(self)
		self.names = []
		self.initialize_tools()

	def __call__(self,event):

		jet_attributes = vectorAttributes(event.jets,{'pt','pt_corrected':'float','eta':'float','phi':'float','e_corrected':'float'})
		muon_attributes
		electron_attributes

		self.METUtility.reset()
	  	self.METUtility.setJetPUcode(0x3300)
		self.METUtility.setObjects(
			ROOT.METUtil.Jets,
			jet_attributes['pt_correctedd']
			jet_attributes['eta']
			jet_attributes['phi']
			jet_attributes['e_corrected']
			list_to_vector(event.jet_AntiKt4LCTopo_MET_wet,'float'),
			list_to_vector(event.jet_AntiKt4LCTopo_MET_wpx,'float'),
			list_to_vector(event.jet_AntiKt4LCTopo_MET_wpy,'float'),
			list_to_vector(event.jet_AntiKt4LCTopo_MET_statusWord,'float'),
			) 
		self.METUtility.setOriJetParameters(
			jet_attributes['pt'],
			)
		self.METUtility.setMuonParameters(
			new_pt_muon_forMET,
			mu_staco_eta,
			mu_staco_phi,
			mu_staco_MET_wet,
			mu_staco_MET_wpx,
			mu_staco_MET_wpy,
			mu_staco_MET_statusWord,
			)
		self.METUtility.setExtraMuonParameters(
			new_pt_muonMS_forMET,
			mu_staco_ms_theta,
			mu_staco_ms_phi,
			)
		self.METUtility.setElectronParameters(
			new_pt_elec_forMET,
			new_eta_elec_forMET,
			new_phi_elec_forMET,
			el_MET_wet,
			el_MET_wpx,
			el_MET_wpy,
			el_MET_statusWord,
			)
		self.METUtility.setMETTerm(
			ROOT.METUtil.RefTau,
			MET_RefTau_etx,
			MET_RefTau_ety,
			MET_RefTau_sumet,
			)
		self.METUtility.setMETTerm(
			ROOT.METUtil.RefGamma,
			MET_RefGamma_etx,
			MET_RefGamma_ety,
			MET_RefGamma_sumet,
			)
		self.METUtility.setMETTerm(
			ROOT.METUtil.CellOutEflow,
			MET_CellOut_Eflow_STVF_etx,
			MET_CellOut_Eflow_STVF_ety,
			MET_CellOut_Eflow_STVF_sumet,
			)

		pxMiss = self.METUtility.getMissingET(self.namespaces['METUtil'].RefFinal).etx()
		pyMiss = self.METUtility.getMissingET(self.namespaces['METUtil'].RefFinal).ety()
		ptMiss = self.METUtility.getMissingET(self.namespaces['METUtil'].RefFinal).et()
		sumET  = self.METUtility.getMissingET(self.namespaces['METUtil'].RefFinal).sumet()




	def initialize_tools(self):

		load('MissingETUtility')

		self.met_utility = ROOT.METUtility()
		self.met_utility.configMissingET(True,True)

['CellOut_Eflow_STVF_et!', 'CellOut_Eflow_STVF_etx', 'CellOut_Eflow_STVF_ety', 'CellOut_Eflow_STVF_sumet', 'RefFinal_STVF', 'RefGamma_etx', 'RefGamma_ety', 'RefGamma_sumet', 'RefTau_etx', 'RefTau_ety', 'RefTau_sumet', 'corr_etx', 'corr_ety', 'etx', 'ety', 'sigma', 'statusWord', 'sumet', 'wet', 'wpx', 'wpy']


"""
  // initialize:
  m_metUtil = new METUtility();
  m_metUtil->configMissingET(true,true); //is2012, isSTVF

  // in event loop:
  m_metUtil->reset();
  m_metUtil->setJetPUcode(MissingETTags::JPU_JET_JVFCUT);//this must be put before setting the jets parameters
  m_metUtil->setObjects(METUtil::Jets,new_pt_jet_forMET,jet_eta,jet_phi,new_e_jet_forMET,jet_AntiKt4LCTopo_MET_wet,jet_AntiKt4LCTopo_MET_wpx,jet_AntiKt4LCTopo_MET_wpy,jet_AntiKt4LCTopo_MET_statusWord); 
  m_metUtil->setOriJetParameters(jet_pt);
  m_metUtil->setMuonParameters(new_pt_muon_forMET,mu_staco_eta,mu_staco_phi,mu_staco_MET_wet,mu_staco_MET_wpx,mu_staco_MET_wpy,mu_staco_MET_statusWord);
  m_metUtil->setExtraMuonParameters(new_pt_muonMS_forMET,mu_staco_ms_theta,mu_staco_ms_phi);
  m_metUtil->setElectronParameters(new_pt_elec_forMET,new_eta_elec_forMET,new_phi_elec_forMET,el_MET_wet,el_MET_wpx,el_MET_wpy,el_MET_statusWord);
  m_metUtil->setMETTerm(METUtil::RefTau,MET_RefTau_etx,MET_RefTau_ety,MET_RefTau_sumet);
  m_metUtil->setMETTerm(METUtil::RefGamma,MET_RefGamma_etx,MET_RefGamma_ety,MET_RefGamma_sumet);
  m_metUtil->setMETTerm(METUtil::CellOutEflow,MET_CellOut_Eflow_STVF_etx,MET_CellOut_Eflow_STVF_ety,MET_CellOut_Eflow_STVF_sumet);

  // get it
  pxMiss = m_metUtil->getMissingET(METUtil::RefFinal).etx();
  pyMiss = m_metUtil->getMissingET(METUtil::RefFinal).ety();
  ptMiss = m_metUtil->getMissingET(METUtil::RefFinal).et();
  sumET  = m_metUtil->getMissingET(METUtil::RefFinal).sumet();
"""
