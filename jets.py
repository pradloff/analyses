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

		self.create_branches.update(dict((branch_name,branch_type) for branch_name,branch_type in [
			('jet_antikt4truth_n',None),
			('jet_antikt4truth_pt',None),
			('jet_antikt4truth_eta',None),
			('jet_antikt4truth_phi',None),
			('jet_antikt4truth_E',None),
			('jet_flavor_truth_label',None),
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

			jet.pt_corrected = jet().Pt()
			jet.eta = jet().Eta()
			jet.phi = jet().Phi()
			jet.e_corrected = jet().E()

			jet.passed_preselection_jets = all([
				jet().Pt()>20000.,
				abs(jet().Eta())<4.5,
				])
			jet.passed_preselection = jet.passed_preselection_jets
			jet.passed_selection = all([
				jet.passed_preselection,
				not ((abs(jet().Eta())<2.4 and jet().Pt()<50000.) and not ((jet.jvtxf)>0.5))
				])
			jet.passed_bselection = all([
				jet.passed_selection,
				abs(jet().Eta())<2.4,
				jet.flavor_weight_MV1> 0.7892,
				])

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
				jet.jvf_up_cut = 0.
				jet.jvf_down_cut = 0.

				jet.jesErrorBasline = 0.
				jet.jesErrorPileup = 0.
				jet.jesErrorCloseby = 0.

				jet.jesErrorFlvCmp = 0.
				jet.jesErrorFlvRsp = 0.
				jet.jesErrorBjet = 0.

				jet.smearFactor = 1.

			if event.is_mc:
				#get jvf uncertainty
				jet.is_pileup_jet = self.jvf_uncertainty_tool.isPileUpJet(jet(),truth_jets)
				if jet().Pt()<50000. and abs(jet().Eta())<2.4:
					jet.jvf_up_cut = self.jvf_uncertainty_tool.getJVFcut(0.5,jet.is_pileup_jet,jet().Pt(),jet().Eta(),True)
					jet.jvf_down_cut = self.jvf_uncertainty_tool.getJVFcut(0.5,jet.is_pileup_jet,jet().Pt(),jet().Eta(),False)
				else:
					jet.jvf_up_cut=0.
					jet.jvf_up_cut=0.
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
				jet.smearFactor = self.jet_smearing_tool.GetRandomSmearingFactorSyst(jet().Pt(),jet().Eta())

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


"""
jet_ActiveArea
jet_ActiveAreaE
jet_ActiveAreaPx
jet_ActiveAreaPy
jet_ActiveAreaPz
jet_AntiKt4LCTopo_MET_BDTMedium_n
jet_AntiKt4LCTopo_MET_BDTMedium_statusWord
jet_AntiKt4LCTopo_MET_BDTMedium_wet
jet_AntiKt4LCTopo_MET_BDTMedium_wpx
jet_AntiKt4LCTopo_MET_BDTMedium_wpy
jet_AntiKt4LCTopo_MET_n
jet_AntiKt4LCTopo_MET_statusWord
jet_AntiKt4LCTopo_MET_wet
jet_AntiKt4LCTopo_MET_wpx
jet_AntiKt4LCTopo_MET_wpy
jet_AntiKt4TopoEM_ActiveArea
jet_AntiKt4TopoEM_ActiveAreaE
jet_AntiKt4TopoEM_ActiveAreaPx
jet_AntiKt4TopoEM_ActiveAreaPy
jet_AntiKt4TopoEM_ActiveAreaPz
jet_AntiKt4TopoEM_AverageLArQF
jet_AntiKt4TopoEM_BCH_CORR_CELL
jet_AntiKt4TopoEM_BCH_CORR_DOTX
jet_AntiKt4TopoEM_BCH_CORR_JET
jet_AntiKt4TopoEM_E
jet_AntiKt4TopoEM_EMJES
jet_AntiKt4TopoEM_EMJES_EtaCorr
jet_AntiKt4TopoEM_EMJESnooffset
jet_AntiKt4TopoEM_EtaOrigin
jet_AntiKt4TopoEM_HECQuality
jet_AntiKt4TopoEM_LArQuality
jet_AntiKt4TopoEM_MOrigin
jet_AntiKt4TopoEM_NegativeE
jet_AntiKt4TopoEM_Offset
jet_AntiKt4TopoEM_PhiOrigin
jet_AntiKt4TopoEM_SamplingMax
jet_AntiKt4TopoEM_Timing
jet_AntiKt4TopoEM_emfrac
jet_AntiKt4TopoEM_emscale_E
jet_AntiKt4TopoEM_emscale_eta
jet_AntiKt4TopoEM_emscale_m
jet_AntiKt4TopoEM_emscale_phi
jet_AntiKt4TopoEM_emscale_pt
jet_AntiKt4TopoEM_eta
jet_AntiKt4TopoEM_flavor_weight_Comb
jet_AntiKt4TopoEM_flavor_weight_GbbNN
jet_AntiKt4TopoEM_flavor_weight_IP2D
jet_AntiKt4TopoEM_flavor_weight_IP3D
jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN
jet_AntiKt4TopoEM_flavor_weight_JetFitterCharm
jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN
jet_AntiKt4TopoEM_flavor_weight_MV1
jet_AntiKt4TopoEM_flavor_weight_MV2
jet_AntiKt4TopoEM_flavor_weight_SV0
jet_AntiKt4TopoEM_flavor_weight_SV1
jet_AntiKt4TopoEM_flavor_weight_SV2
jet_AntiKt4TopoEM_flavor_weight_SecondSoftMuonTagChi2
jet_AntiKt4TopoEM_flavor_weight_SoftMuonTagChi2
jet_AntiKt4TopoEM_fracSamplingMax
jet_AntiKt4TopoEM_hecf
jet_AntiKt4TopoEM_isBadLoose
jet_AntiKt4TopoEM_isBadLooseMinus
jet_AntiKt4TopoEM_isBadMedium
jet_AntiKt4TopoEM_isBadTight
jet_AntiKt4TopoEM_isUgly
jet_AntiKt4TopoEM_jvtx_x
jet_AntiKt4TopoEM_jvtx_y
jet_AntiKt4TopoEM_jvtx_z
jet_AntiKt4TopoEM_jvtxf
jet_AntiKt4TopoEM_jvtxfFull
jet_AntiKt4TopoEM_m
jet_AntiKt4TopoEM_n
jet_AntiKt4TopoEM_phi
jet_AntiKt4TopoEM_pt
jet_AntiKt4TopoEM_sumPtTrk
jet_AverageLArQF
jet_BCH_CORR_CELL
jet_BCH_CORR_DOTX
jet_BCH_CORR_JET
jet_E
jet_EMJES
jet_EMJES_EtaCorr
jet_EMJESnooffset
jet_EtaOrigin
jet_HECQuality
jet_LArQuality
jet_MOrigin
jet_NegativeE
jet_Offset
jet_PhiOrigin
jet_SamplingMax
jet_Timing
jet_constscale_E
jet_constscale_eta
jet_constscale_m
jet_constscale_phi
jet_constscale_pt
jet_emfrac
jet_emscale_E
jet_emscale_eta
jet_emscale_m
jet_emscale_phi
jet_emscale_pt
jet_eta
jet_flavor_weight_Comb
jet_flavor_weight_GbbNN
jet_flavor_weight_IP2D
jet_flavor_weight_IP3D
jet_flavor_weight_JetFitterCOMBNN
jet_flavor_weight_JetFitterCharm
jet_flavor_weight_JetFitterTagNN
jet_flavor_weight_MV1
jet_flavor_weight_MV2
jet_flavor_weight_SV0
jet_flavor_weight_SV1
jet_flavor_weight_SV2
jet_flavor_weight_SecondSoftMuonTagChi2
jet_flavor_weight_SoftMuonTagChi2
jet_fracSamplingMax
jet_hecf
jet_isBadLoose
jet_isBadLooseMinus
jet_isBadMedium
jet_isBadTight
jet_isUgly
jet_jvtx_x
jet_jvtx_y
jet_jvtx_z
jet_jvtxf
jet_jvtxfFull
jet_m
jet_n
jet_phi
jet_pt
jet_sumPtTrk
"""

