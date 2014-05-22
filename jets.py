from common.functions import event_function
from common.particle import particle
from common.external import load
import os
import ROOT
from copy import copy
from math import sin


class collect_tracks(event_function):

	def __init__(self,collection_name='trk_'):
		event_function.__init__(self)
		
		self.collection_name = collection_name

		self.names = [
			'pt',
			'eta',
			#'phi',
			'phi_wrtPV',
			'd0',
			'z0_wrtPV',
			'nPixHits',
			'nBLHits',
			]

		self.required_branches += [self.collection_name+name for name in self.names]
		self.required_branches += [self.collection_name+'n']


	def __call__(self,event):
		
		event.trks = {}
		for trk in range(event.__dict__[self.collection_name+'n']):
			event.trks[trk] = particle(\
				**dict((name,event.__dict__[self.collection_name+name][trk]) for name in self.names)
				)

		for trk in event.trks.values():
			trk.jet_owner = None
			trk.jet_owner_dR = 10000.
			trk.set_pt_eta_phi_m(
				trk.pt,
				trk.eta,
				trk.phi_wrtPV,
				#trk.phi,
				0.,
				)
			trk.passed_b_selection = all([
				trk.pt>1000.,
				abs(trk.d0)<1.,
				abs(trk.z0_wrtPV)*sin(trk().Theta())<1.5,
				trk.nPixHits+trk.nBLHits>=2,
				trk.nBLHits>=1,
				])

class collect_jets(event_function):

	def __init__(self,collection_name='jet_'):
		event_function.__init__(self)
		self.collection_name = collection_name
		self.names = [
			'constscale_E',
			'constscale_eta',
			'constscale_phi',
			'constscale_m',
			'ActiveAreaPx',
			'ActiveAreaPy',
			'ActiveAreaPz',
			'ActiveAreaE',
			'flavor_weight_MV1',
			'isBadLooseMinus',
			'jvtxf',
			]
		
		self.required_branches += [self.collection_name+name for name in self.names]
		self.required_branches += [self.collection_name+'n']
		self.required_branches += [
			'Eventshape_rhoKt4LC',
			'averageIntPerXing',
			'nPV_2trks',
			'EventNumber',
			]

		#Exists for MC not data
		self.create_branches.update(dict((name,branch_type) for name,branch_type in [
			('jet_antikt4truth_n',None),
			('jet_antikt4truth_pt',None),
			('jet_antikt4truth_eta',None),
			('jet_antikt4truth_phi',None),
			('jet_antikt4truth_E',None),
			('jet_flavor_truth_label',None),
			]))

		self.create_branches.update(dict((name,branch_type) for name,branch_type in [
			('jet_n','int'),
			('jet_passed_b_preselection','std.vector.bool'),
			#('jet_b_preselection_pt','std.vector.float'),
			#('jet_b_preselection_eta','std.vector.float'),
			#('jet_b_preselection_phi','std.vector.float'),
			#('jet_b_preselection_E','std.vector.float'),
			('jet_pt','std.vector.float'),
			('jet_eta','std.vector.float'),
			('jet_phi','std.vector.float'),
			('jet_E','std.vector.float'),
			('jet_jvf','std.vector.float'),
			('jet_jvf_up_cut','std.vector.float'),
			('jet_jvf_down_cut','std.vector.float'),
			('jet_jes_Error_Baseline','std.vector.float'),
			('jet_jes_Error_Pileup','std.vector.float'),
			('jet_jes_Error_Closeby','std.vector.float'),
			('jet_jes_Error_FlvCmp','std.vector.float'),
			('jet_jes_Error_FlvRsp','std.vector.float'),
			('jet_jes_Error_Bjet','std.vector.float'),
			('jet_bJet_scale_factor','std.vector.float'),
			('jet_bJet_scale_factor_error','std.vector.float'),
			('jet_flavor_weight_MV1','std.vector.float'),
			('jets',None),	
			]))

		#Load jet correction tools
		self.initialize_tools()

	def __call__(self,event):

		#Collect jets
		event.jets = {}
		for jet in range(event.__dict__[self.collection_name+'n']):
			event.jets[jet] = particle(\
				**dict((name,event.__dict__[self.collection_name+name][jet]) for name in self.names)
				)

		#Apply jet corrections
		self.apply_corrections(event)
		self.apply_btag_corrections(event)

		#Categorize jet as b-jet preselection
		self.preselect_bjets(event)

		#Define selections
		for jet in event.jets.values():

			jet.passed_preselection = all([
				jet.pt>20000.,
				abs(jet.eta)<4.5,
				])
			jet.passed_selection = all([
				jet.passed_preselection,
				not ((abs(jet.eta)<2.4 and jet.pt<50000.) and not ((jet.jvtxf)>0.5))
				])
			jet.passed_bselection = all([
				jet.passed_selection,
				abs(jet.eta)<2.4,
				jet.flavor_weight_MV1> 0.7892,
				])

		event.original_jet_pt = copy(event.jet_pt)
		#saves
		event.jet_n = 0
		event.jet_pt = []
		event.jet_eta = []
		event.jet_phi = []
		event.jet_E = []

		event.jet_passed_b_preselection = []
		"""
		event.jet_b_preselection_pt = []
		event.jet_b_preselection_eta = []
		event.jet_b_preselection_phi = []
		event.jet_b_preselection_E = []
		"""

		event.jet_jvf = []
		event.jet_jvf_up_cut = []
		event.jet_jvf_down_cut = []

		event.jet_jes_Error_Baseline = []
		event.jet_jes_Error_Pileup = []
		event.jet_jes_Error_Closeby = []

		event.jet_jes_Error_FlvCmp = []
		event.jet_jes_Error_FlvRsp = []
		event.jet_jes_Error_Bjet = []

		event.jet_bJet_scale_factor = []
		event.jet_bJet_scale_factor_error = []	
		event.jet_flavor_weight_MV1 = []

		for jet in event.jets.values():
			if not jet.passed_preselection: continue
			event.jet_n+=1
			event.jet_pt.append(jet.pt)
			event.jet_eta.append(jet.eta)
			event.jet_phi.append(jet.phi)
			event.jet_E.append(jet.E)

			event.jet_passed_b_preselection.append(jet.passed_b_preselection)
			"""
			event.jet_b_preselection_pt.append(jet.b_preselection_pt)
			event.jet_b_preselection_eta.append(jet.b_preselection_eta)
			event.jet_b_preselection_phi.append(jet.b_preselection_phi)
			event.jet_b_preselection_E.append(jet.b_preselection_E)
			"""

			event.jet_jvf.append(jet.jvtxf)
			event.jet_jvf_up_cut.append(jet.jvf_up_cut)
			event.jet_jvf_down_cut.append(jet.jvf_down_cut)

			event.jet_jes_Error_Baseline.append(jet.jesErrorBaseline)
			event.jet_jes_Error_Pileup.append(jet.jesErrorPileup)
			event.jet_jes_Error_Closeby.append(jet.jesErrorCloseby)

			event.jet_jes_Error_FlvCmp.append(jet.jesErrorFlvCmp)
			event.jet_jes_Error_FlvRsp.append(jet.jesErrorFlvRsp)
			event.jet_jes_Error_Bjet.append(jet.jesErrorBjet)

			event.jet_bJet_scale_factor.append(jet.bJetScaleFactor)
			event.jet_bJet_scale_factor_error.append(jet.bJetScaleFactorError)
			event.jet_flavor_weight_MV1.append(jet.flavor_weight_MV1)
		return

	def preselect_bjets(self,event):
		cone_size = 0.4

		for jet in event.jets.values():
			jet.btrack_selection_ntrks = 0

		for trk in event.trks.values():
			if not trk.passed_b_selection: continue
			for jetN,jet in event.jets.items():
				if jet.pt<20000.: continue
				dR = trk().DeltaR(jet())
				if dR<cone_size and dR<trk.jet_owner_dR:
					trk.jet_owner = jetN
					trk.jet_owner_dR = dR
			if trk.jet_owner is None: continue
			event.jets[trk.jet_owner].btrack_selection_ntrks += 1

		for jet in event.jets.values():
			jet.passed_b_preselection = all([
				jet.btrack_selection_ntrks>=1,
				])

		"""
		for jet in event.jets.values():

			jet.b_preselection_pt = jet.btrack_selection_vectors.Pt()
			jet.b_preselection_eta = jet.btrack_selection_vectors.Eta()
			jet.b_preselection_phi = jet.btrack_selection_vectors.Phi()
			jet.b_preselection_E = jet.btrack_selection_vectors.E()

			jet.passed_b_preselection = all([
				jet.b_preselection_pt>15000.,
				jet.b_preselection_eta<2.5,
				])
		"""

	def apply_corrections(self,event):

		#get truth jet info if mc, set dummies otherwise
		if event.is_mc:
			truth_jets = ROOT.std.vector('TLorentzVector')()
			for jet_n in range(event.jet_antikt4truth_n):
				if event.jet_antikt4truth_pt[jet_n] < 10000.: continue
				jet = ROOT.TLorentzVector()
				jet.SetPtEtaPhiE(
					event.jet_antikt4truth_pt[jet_n],
					event.jet_antikt4truth_eta[jet_n],
					event.jet_antikt4truth_phi[jet_n],
					event.jet_antikt4truth_E[jet_n],
					)
				truth_jets.push_back(jet)
		else:
			event.jet_antikt4truth_n = 0
			for name in [
				'jet_antikt4truth_pt',
				'jet_antikt4truth_eta',
				'jet_antikt4truth_phi',
				'jet_antikt4truth_E',
				'jet_flavor_truth_label',
				]: event.__dict__[name] = []

		self.jet_smearing_tool.SetSeed(event.EventNumber)

		for jetN,jet in event.jets.items():
			#apply jet calibration
			if event.is_mc:
				jet.set_particle(self.jet_calibration_tool_mc.ApplyJetAreaOffsetEtaJES(
					jet.constscale_E,
					jet.constscale_eta,
					jet.constscale_phi,
					jet.constscale_m,
					jet.ActiveAreaPx,
					jet.ActiveAreaPy,
					jet.ActiveAreaPz,
					jet.ActiveAreaE,
					event.Eventshape_rhoKt4LC,
					event.averageIntPerXing,
					event.nPV_2trks
					))
			else:
				jet.set_particle(self.jet_calibration_tool_data.ApplyJetAreaOffsetEtaJES(
					jet.constscale_E,
					jet.constscale_eta,
					jet.constscale_phi,
					jet.constscale_m,
					jet.ActiveAreaPx,
					jet.ActiveAreaPy,
					jet.ActiveAreaPz,
					jet.ActiveAreaE,
					event.Eventshape_rhoKt4LC,
					event.averageIntPerXing,
					event.nPV_2trks
					))

		for jetN,jet in event.jets.items():
			if not event.is_mc:
				jet.is_pileup_jet = False
				jet.jvf_up_cut = jet.jvtxf
				jet.jvf_down_cut = jet.jvtxf

				jet.jesErrorBaseline = 0.
				jet.jesErrorPileup = 0.
				jet.jesErrorCloseby = 0.

				jet.jesErrorFlvCmp = 0.
				jet.jesErrorFlvRsp = 0.
				jet.jesErrorBjet = 0.

			if event.is_mc:
				#get jvf uncertainty
				jet.is_pileup_jet = self.jvf_uncertainty_tool.isPileUpJet(jet(),truth_jets)
				if jet().Pt()<50000. and abs(jet().Eta())<2.4:
					jet.jvf_up_cut = self.jvf_uncertainty_tool.getJVFcut(0.5,jet.is_pileup_jet,jet().Pt(),jet().Eta(),True)
					jet.jvf_down_cut = self.jvf_uncertainty_tool.getJVFcut(0.5,jet.is_pileup_jet,jet().Pt(),jet().Eta(),False)
				else:
					jet.jvf_up_cut=0.
					jet.jvf_down_cut=0.

				#get jes uncertainty
				close_by_fraction = 0.
				for jet_neighborN,jet_neighbor in event.jets.items():
					if any([
						jet_neighborN==jetN,
						jet_neighbor().Pt()<12000.,
						jet().DeltaR(jet_neighbor())>1.1
						]): continue
					close_by_fraction += jet().Vect().Dot(jet_neighbor().Vect())/(jet().P()**2.)

  				#jet.jesErrorBaseline = self.jes_uncertainty_provider.getRelUncert(jet().Pt(),jet().Eta(), close_by_fraction, jet.is_pileup_jet, event.jet_flavor_truth_label[jetN]==5,len(event.jets))
				#jet.jesErrorBaseline = self.jes_uncertainty_provider.getRelUncert(jet().Pt(),jet().Eta())
				jet.jesErrorBaseline = abs(self.jes_uncertainty_provider._JESUncertaintyProvider__getRelUncertUncorr(jet().Pt(),jet().Eta()))+\
					abs(self.jes_uncertainty_provider._JESUncertaintyProvider__getRelUncertBias(jet().Pt(),jet().Eta()))
				jet.jesErrorPileup = self.jes_uncertainty_provider.getRelOffsetUncert(jet().Pt(),jet().Eta(),event.nPV_2trks,event.averageIntPerXing)
				jet.jesErrorCloseby = self.jes_uncertainty_provider.getRelClosebyUncert(jet().Pt(),close_by_fraction)

				jet.jesErrorFlvCmp = 0.
				jet.jesErrorFlvRsp = 0.
				jet.jesErrorBjet = 0.

				if event.jet_flavor_truth_label[jetN]==5:
					jet.jesErrorBjet = self.jes_uncertainty_provider.getRelBJESUncert(jet().Pt(),jet().Eta())
				else:
					jet.jesErrorFlvCmp = self.jes_uncertainty_provider.getRelFlavorCompUncert(jet().Pt(),jet().Eta(),True);
					jet.jesErrorFlvRsp = self.jes_uncertainty_provider.getRelFlavorResponseUncert(jet().Pt(),jet().Eta());
				
				#smear jet
				jet_tlv = jet()
				jet_tlv *= self.jet_smearing_tool.GetRandomSmearingFactorSyst(jet().Pt(),jet().Eta())

			jet.pt = jet().Pt()
			jet.eta = jet().Eta()
			jet.phi = jet().Phi()
			jet.E = jet().E()

	def apply_btag_corrections(self,event):
		for jetN,jet in event.jets.items():

			if event.is_mc:
				calibration_jet = ROOT.Analysis.CalibrationDataVariables()
				if abs(jet().Eta())<2.4 and jet().Pt()<50000. and abs(jet.jvtxf)>0.5: calibration_jet.jetAuthor = 'AntiKt4TopoLCJVF0_5'
				else: calibration_jet.jetAuthor = 'AntiKt4TopoLCnoJVF'
				calibration_jet.jetPt = jet().Pt()
				calibration_jet.jetEta = jet().Pt()
				if event.jet_flavor_truth_label[jetN] == 4: jet_type = ROOT.std.string('C')
				elif event.jet_flavor_truth_label[jetN] == 5: jet_type = ROOT.std.string('B')
				elif event.jet_flavor_truth_label[jetN] == 15: jet_type = ROOT.std.string('T')
				else: jet_type = ROOT.std.string('Light')

				if jet.flavor_weight_MV1> 0.7892: result = self.btag_calibration_tool.getScaleFactor(calibration_jet,jet_type,"0_7892",1)
				else: result = self.btag_calibration_tool.getInefficiencyScaleFactor(calibration_jet,jet_type,"0_7892",1)

				bJetScaleFactor = result.first
				bJetScaleFactorError = result.second
			else:
        			bJetScaleFactor = 1.
        			bJetScaleFactorError = 0.

			jet.bJetScaleFactor = bJetScaleFactor
			jet.bJetScaleFactorError = bJetScaleFactorError

	def initialize_tools(self):
		analysis_home = os.getenv('ANALYSISHOME')

		load('ApplyJetCalibration')
		self.jet_calibration_tool_mc = ROOT.JetCalibrationTool(
			"AntiKt4LCTopo",
			"{0}/external/ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_Jan13.config".format(analysis_home),
			False,
			)
		self.jet_calibration_tool_data = ROOT.JetCalibrationTool(
			"AntiKt4LCTopo",
			"{0}/external/ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_Jan13.config".format(analysis_home),
			True,
			)

		load('JVFUncertaintyTool')
		self.jvf_uncertainty_tool = ROOT.JVFUncertaintyTool("AntiKt4LCTopo")

		load('JetUncertainties')
		self.jes_uncertainty_provider = ROOT.MultijetJESUncertaintyProvider(
			"JES_2012/Moriond2013/MultijetJES_2012.config",
			"JES_2012/Moriond2013/InsituJES2012_AllNuisanceParameters.config",
			"AntiKt4LCTopo",
			"MC12a",
			"{0}/external/JetUncertainties/share".format(analysis_home),
			)

		load('ApplyJetResolutionSmearing')
		self.jet_smearing_tool = ROOT.JetSmearingTool(
			"AntiKt4LCTopo",
			"{0}/external/JetResolution/share/JERProviderPlots_2012.root".format(analysis_home)
			)
		self.jet_smearing_tool.init()

		load('CalibrationDataInterface')
		self.btag_calibration_tool = ROOT.Analysis.CalibrationDataInterfaceROOT(
			'MV1',
			os.path.relpath('{0}/external/CalibrationDataInterface/share/BTagCalibration.env'.format(analysis_home)),
			'{0}/external/CalibrationDataInterface/share/'.format(analysis_home)
			)

