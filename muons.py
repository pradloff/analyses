from common.functions import event_function,arg
from common.particle import particle
from common.external import load
import ROOT
from math import sin,cos,log,tan,pi
from misc import list_attributes
import os

class collect_muons(event_function):

	def __init__(
		self,
		min_pT = arg(10000.,help='minimum pT [MeV]'),
		collection_name='mu_staco_'
		):
		event_function.__init__(self)

		self.collection_name = collection_name
		self.min_pT = min_pT

		self.names = [
			'E',
			'MET_statusWord', 
			'MET_wet',
			'MET_wpx',
			'MET_wpy',
			'charge',
			'eta',
			'etcone20',
			'id_qoverp',
			'id_theta',
			'id_z0_exPV',
			'isCombinedMuon',
			'isSegmentTaggedMuon',
			'isStandAloneMuon',
			'me_qoverp',
			'me_theta',
			'ms_phi',
			'ms_qoverp',
			'ms_theta',
			'phi',
			'pt',
			'ptcone40',
			'px',
			'py',
			'pz',
			'loose',
			'tight',
			'expectBLayerHit',
			'nBLHits',
			'nPixHits',
			'nPixelDeadSensors',
			'nSCTDeadSensors',
			'nSCTHits',
			'nPixHoles',
			'nSCTHoles',
			'nTRTHits',
			'nTRTOutliers',
			]

		"""
		self.new_collection_names = dict((name,branch_type) for name,branch_type in [
			('etcone20_corrected','std.vector.float'),
			('pt_corrected','std.vector.float'),
			('scaleFactorReco','std.vector.float'),
			('scaleFactorRecoError','std.vector.float'),
			('smearFactor','std.vector.float'),
			('passed_preselection_embedding','std.vector.bool'),
			('passed_preselection','std.vector.bool'),
			('passed_selection','std.vector.bool'),
			])
		"""
		self.required_branches += [self.collection_name+name for name in self.names]
		self.required_branches += [self.collection_name+'n']
		self.create_branches['muons'] = None
		#self.create_branches.update(dict((branch_name,branch_type) for branch_name,branch_type in [
		#	('muons',None),
		#	]+[(self.collection_name+name,branch_type) for name,branch_type in self.new_collection_names.items()]))

		#self.keep_branches += [self.collection_name+name for name in self.names]
		#self.keep_branches += [self.collection_name+'n']
		
		self.initialize_tools()

	def __call__(self,event):

		#Collect muons
		event.muons = {}
		for mu in range(event.__dict__[self.collection_name+'n']):
			event.muons[mu] = particle(\
				**dict((name,event.__dict__[self.collection_name+name][mu]) for name in self.names)
				)

		#Apply corrections
		self.apply_corrections(event)

		#Define selections
		for muon in event.muons.values():
			muon.passed_preselection_taus = all([
				muon.id_z0_exPV<10.,
				muon.pt_corrected>4000.,
				abs(muon.eta)<2.5,
				muon.loose,
				muon.expectBLayerHit==0 or muon.nBLHits>0,
				(muon.nPixHits+muon.nPixelDeadSensors)>1,
				(muon.nSCTHits+muon.nSCTDeadSensors)>5,
				(muon.nPixHoles+muon.nSCTHoles)<3,
				not((0.1<abs(muon.eta)<1.9) and not ((muon.nTRTHits+muon.nTRTOutliers)>5 and muon.nTRTOutliers<0.9*(muon.nTRTHits+muon.nTRTOutliers))),
				not(((0.1>abs(muon.eta) or abs(muon.eta)>1.9) and ((muon.nTRTHits+muon.nTRTOutliers)>5)) and not (muon.nTRTOutliers<0.9*(muon.nTRTHits+muon.nTRTOutliers))),			
				])

			muon.passed_preselection_embedding = all([
				muon.passed_preselection_taus,
				muon.pt_corrected>self.min_pT,
				])

			muon.passed_preselection = all([
				muon.passed_preselection_taus,
				muon.tight,
				muon.pt_corrected>self.min_pT,
				])

			muon.passed_selection = all([
				muon.passed_preselection,
				muon.isCombinedMuon,
				muon.etcone20_corrected/muon.pt_corrected<0.09,
				muon.ptcone40/muon.pt_corrected<0.18,	
				])

		#event.__dict__.update(list_attributes(event.muons,self.new_collection_names.keys(),self.collection_name))

	def apply_corrections(self,event):

		for muonN,muon in event.muons.items():

			muon.set_px_py_pz_e(muon.px,muon.py,muon.pz,muon.E)

			muon.me_pt = sin(muon.me_theta)/abs(muon.me_qoverp)
			muon.ms_pt = sin(muon.ms_theta)/abs(muon.ms_qoverp)
			muon.id_pt = sin(muon.id_theta)/abs(muon.id_qoverp)
			muon.me_eta = -log(tan(muon.me_theta/2.)+pi)

			if event.is_mc and muon().Pt()>4000.:
				#Compute scale factors
				scaleFactor = self.muon_scalefactor_tool.scaleFactor(int(muon.charge),muon())
 				scaleFactorError = self.muon_scalefactor_tool.scaleFactorUncertainty(int(muon.charge),muon())+\
					self.muon_scalefactor_tool.scaleFactorSystematicUncertainty(int(muon.charge),muon())
				#Smear energy
				self.mcp_smear.SetSeed(event.EventNumber,muonN)

				if muon.isCombinedMuon:
					self.mcp_smear.Event(muon.me_pt,muon.id_pt,muon.pt,muon.eta,int(muon.charge));
					smearFactor = self.mcp_smear.pTCB()/muon.pt;
					smearFactorMS = self.mcp_smear.pTMS()/muon.me_pt
				elif muon.isSegmentTaggedMuon:
					self.mcp_smear.Event(muon.pt,muon.eta,"ID",int(muon.charge));
					smearFactor = self.mcp_smear.pTID()/muon.pt;
					smearFactorMS = 1.
				elif muon.isStandAloneMuon:
					self.mcp_smear.Event(muon.me_pt,muon.me_eta,"MS",1 if muon.me_qoverp>0. else -1)
					smearFactor = self.mcp_smear.pTMS()/muon.me_pt
					smearFactorMS = smearFactor
				else:
					smearFactor = 1.
					smearFactorMS = 1.		
			else:
				scaleFactor = 1.
				scaleFactorError = 0.
				smearFactor = 1.
				smearFactorMS = 1.

			muon.scaleFactorReco = scaleFactor
			muon.scaleFactorRecoError = scaleFactorError
			muon.smearFactor = smearFactor
			muon.smearFactorMS = smearFactorMS
			muon.E_corrected = muon.E*smearFactor
			muon.pt_corrected = muon.pt*smearFactor
			muon.pt_corrected_MS = muon.ms_pt*smearFactorMS
			muon.set_particle(muon()*smearFactor)

			muon.scale_factor = muon.scaleFactorReco
			muon.scale_factor_error = muon.scaleFactorRecoError

			#Correct etcone values
			muon.etcone20_corrected = self.muon_isolation_tool.CorrectEtCone(
				muon.etcone20,
				event.nPV_3trks,
				muon.eta,
				"cone20Comb"
				)
		return

	def initialize_tools(self):
		load('MuonIsolationCorrection')
		self.muon_isolation_tool = ROOT.CorrectCaloIso()

		load('MuonEfficiencyCorrections')
		self.muon_scalefactor_tool = ROOT.Analysis.AnalysisMuonConfigurableScaleFactors(
			"{0}/external/MuonEfficiencyCorrections/share/".format(os.getenv('ANALYSISHOME')),
			"STACO_CB_plus_ST_2012_SF.txt",
			"MeV",
			ROOT.Analysis.AnalysisMuonConfigurableScaleFactors.AverageOverRuns
			)
		self.muon_scalefactor_tool.Initialise()

		load('MuonMomentumCorrections')
		self.mcp_smear = ROOT.MuonSmear.SmearingClass(
			"Data12",
			"staco",
			"pT",
			"Rel17.2Repro",
			"{0}/external/MuonMomentumCorrections/share/".format(os.getenv('ANALYSISHOME')),
			)
		self.mcp_smear.UseScale(1)
		self.mcp_smear.UseImprovedCombine()

