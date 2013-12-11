from common.functions import event_function
from common.external import load
import ROOT
from misc import vector_attributes

class correct_missing_energy(event_function):

	def __init__(self):
		event_function.__init__(self)
		self.required_branches = [
			'jet_pt',
			'jet_AntiKt4LCTopo_MET_wet',
			'jet_AntiKt4LCTopo_MET_wpx',
			'jet_AntiKt4LCTopo_MET_wpy',
			'jet_AntiKt4LCTopo_MET_statusWord',
			'mu_staco_eta',
			'mu_staco_phi',
			'mu_staco_MET_wet',
			'mu_staco_MET_wpx',
			'mu_staco_MET_wpy',
			'mu_staco_MET_statusWord',
			'mu_staco_ms_theta',
			'mu_staco_ms_phi',
			'el_MET_wet',
			'el_MET_wpx',
			'el_MET_wpy',
			'el_MET_statusWord',
			'MET_RefTau_etx',
			'MET_RefTau_ety',
			'MET_RefTau_sumet',
			'MET_RefGamma_etx',
			'MET_RefGamma_ety',
			'MET_RefGamma_sumet',
			'MET_CellOut_Eflow_STVF_etx',
			'MET_CellOut_Eflow_STVF_ety',
			'MET_CellOut_Eflow_STVF_sumet',
			'muons',
			'electrons',
			'jets',			
			]
		self.initialize_tools()

	def __call__(self,event):

		jet_attributes = vector_attributes(event.jets,{'pt_corrected':'float','eta':'float','phi':'float','e_corrected':'float'})
		muon_attributes = vector_attributes(event.muons,{'pt_corrected':'float','pt_corrected_MS':'float'})
		electron_attributes = vector_attributes(event.electrons,{'pt_corrected':'float','eta':'float','phi':'float'})

		self.met_utility.reset()
	  	self.met_utility.setJetPUcode(0x3300)
		self.met_utility.setObjects(
			ROOT.METUtil.Jets,
			jet_attributes['pt_corrected'],
			jet_attributes['eta'],
			jet_attributes['phi'],
			jet_attributes['e_corrected'],
			event.jet_AntiKt4LCTopo_MET_wet,
			event.jet_AntiKt4LCTopo_MET_wpx,
			event.jet_AntiKt4LCTopo_MET_wpy,
			event.jet_AntiKt4LCTopo_MET_statusWord,
			) 
		self.met_utility.setOriJetParameters(
			event.jet_pt,
			)
		self.met_utility.setMuonParameters(
			muon_attributes['pt_corrected'],
			event.mu_staco_eta,
			event.mu_staco_phi,
			event.mu_staco_MET_wet,
			event.mu_staco_MET_wpx,
			event.mu_staco_MET_wpy,
			event.mu_staco_MET_statusWord,
			)
		self.met_utility.setExtraMuonParameters(
			muon_attributes['pt_corrected_MS'],
			event.mu_staco_ms_theta,
			event.mu_staco_ms_phi,
			)
		self.met_utility.setElectronParameters(
			electron_attributes['pt_corrected'],
			electron_attributes['eta'],
			electron_attributes['phi'],
			event.el_MET_wet,
			event.el_MET_wpx,
			event.el_MET_wpy,
			event.el_MET_statusWord,
			)
		self.met_utility.setMETTerm(
			ROOT.METUtil.RefTau,
			event.MET_RefTau_etx,
			event.MET_RefTau_ety,
			event.MET_RefTau_sumet,
			)
		self.met_utility.setMETTerm(
			ROOT.METUtil.RefGamma,
			event.MET_RefGamma_etx,
			event.MET_RefGamma_ety,
			event.MET_RefGamma_sumet,
			)
		self.met_utility.setMETTerm(
			ROOT.METUtil.CellOutEflow,
			event.MET_CellOut_Eflow_STVF_etx,
			event.MET_CellOut_Eflow_STVF_ety,
			event.MET_CellOut_Eflow_STVF_sumet,
			)

		pxMiss = self.met_utility.getMissingET(ROOT.METUtil.RefFinal).etx()
		pyMiss = self.met_utility.getMissingET(ROOT.METUtil.RefFinal).ety()
		ptMiss = self.met_utility.getMissingET(ROOT.METUtil.RefFinal).et()
		sumET  = self.met_utility.getMissingET(ROOT.METUtil.RefFinal).sumet()

	def initialize_tools(self):

		load('MissingETUtility')

		self.met_utility = ROOT.METUtility()
		self.met_utility.configMissingET(True,True)

