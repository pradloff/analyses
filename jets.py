from common.functions import event_function
from common.particle import particle
from common.external import load
import os
import ROOT

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
			'isBadLoose',
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

		self.create_branches.update(dict((name,branch_type) for name,branch_type in [
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


		#saves
		event.jet_pt = []
		event.jet_eta = []
		event.jet_phi = []
		event.jet_E = []
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
		
		for jet in event.jets.values():
			if not jet.passed_preselection: continue
			event.jet_pt.append(jet.pt)
			event.jet_eta.append(jet.eta)
			event.jet_phi.append(jet.phi)
			event.jet_E.append(jet.E)
			event.jet_jvf.append(jet.jvtxf)
			event.jet_jvf_up.append(jet.jvf_up_cut)
			event.jet_jvf_down.append(jet.jvf_down_cut)

			event.jet_jes_Error_Baseline.append(jet.jesErrorBaseline)
			event.jet_jes_Error_Pileup.append(jet.jesErrorPileup)
			event.jet_jes_Error_Closeby.append(jet.jesErrorCloseby)

			event.jet_jes_Error_FlvCmp.append(jet.jesErrorFlvCmp)
			event.jet_jes_Error_FlvRsp.append(jet.jesErrorFlvRsp)
			event.jet_jes_Error_Bjet.append(jet.jesErrorBjet)

			event.jet_bJet_scale_factor.append(jet.bJetScaleFactor)
			event.jet_bJet_scale_factor_error.append(jet.bJetScaleFactorError)

		return

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

