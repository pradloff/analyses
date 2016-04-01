from common.functions import event_function
from common.particle import particle
from math import cosh,sqrt
from common.external import load
import ROOT
import os
from misc import list_attributes
from common.commandline import commandline,arg
from common.branches import auto_branch,branch
import code
class collect_electrons(event_function):

	@commandline(
	    'collect_electrons',
        cut_crack = arg('--cc',action='store_true',help='Cut crack region'),
	    min_pT = arg('--min_pt',type=float,help='Minimum el pT'),
	    )
	def __init__(
		self,
		min_pT = 15000.,
		cut_crack = True,
		collection_name='el_'
		):
		super(collect_electrons,self).__init__()

		self.cut_crack = cut_crack
		self.min_pT = min_pT
		self.collection_name = collection_name

		self.names = [
			'MET_statusWord',
			'MET_wet',
			'MET_wpx',
			'MET_wpy',
			'OQ',
			'cl_E',
			'cl_pt',
			'cl_eta',
			'cl_phi',
			'tracketa',
			'trackphi',
			'etap',
			'etas2',
			'nPixHits',
			'nSCTHits',
			'Etcone20',
			'ptcone40',
			'author',
			'charge',
			'loosePP',
			'mediumPP',
			'tightPP',
			]

		self.new_collection_names = dict((name,branch_type) for name,branch_type in [
			('pt_corrected','std.vector.float'),
			('E_corrected','std.vector.float'),
			('etcone20_corrected','std.vector.float'),
			('mediumScaleFactor','std.vector.float'),
			('mediumScaleFactorError','std.vector.float'),
			('tightScaleFactor','std.vector.float'),
			('tightScaleFactorError','std.vector.float'),
			('is_tightPP','std.vector.bool'),
			('is_mediumPP','std.vector.bool'),
			('is_loosePP','std.vector.bool'),
			('passed_preselection','std.vector.bool'),
			('passed_selection','std.vector.bool'),
			])

		for name in self.names:
			self.branches.append(branch(self.collection_name+name,'r'))

		self.branches.append(branch(self.collection_name+'n','r'))
		#self.branches.append(branch('random_RunNumber','r'))

		self.initialize_tools()

	def __call__(self,event):
		super(collect_electrons,self).__call__(event)
		code.interact(local=locals())

		#Collect electrons
		event.electrons = {}
		for el in range(event.__dict__[self.collection_name+'n']):
			event.electrons[el] = particle(\
				**dict((name,event.__dict__[self.collection_name+name][el]) for name in self.names)
				)

		#set electron direction
		for electron in event.electrons.values():
			if (electron.nPixHits+electron.nSCTHits)<4:
				eta = electron.cl_eta
				phi = electron.cl_phi
				pt = electron.cl_pt
				E = electron.cl_E
			else:
				eta = electron.tracketa
				phi = electron.trackphi
				try:
					pt = electron.cl_E/cosh(electron.tracketa)
				except ArithmeticError:
					pt = -999.
				E = electron.cl_E
			electron.pt = pt
			electron.eta = eta
			electron.phi = phi
			electron.E = E

		#Apply electron collections
		self.apply_corrections(event)

		#Define selections
		for electron in event.electrons.values():
			electron.passed_preselection_taus = all([
				electron.is_loosePP,
				electron.pt_corrected>self.min_pT,
				abs(electron.cl_eta)<1.37 or (1.52<abs(electron.cl_eta)<2.47) if self.cut_crack else abs(electron.cl_eta)<2.47,
				electron.author==1 or electron.author==3,
				(electron.OQ&1446)==0
				])
			electron.passed_preselection_embedding = electron.passed_preselection_taus

			electron.passed_preselection = all([
				electron.passed_preselection_taus,
				electron.is_mediumPP,
				])
			electron.passed_selection = all([
				electron.passed_preselection,
				electron.etcone20_corrected/electron.pt_corrected<0.09,
				electron.ptcone40/electron.pt_corrected<0.17,
				True
				])

		event.__dict__.update(list_attributes(event.electrons,self.new_collection_names.keys(),self.collection_name))

	def apply_corrections(self,event):

		for electronN,electron in event.electrons.items():

			electron.is_tightPP = electron.tightPP
			electron.is_mediumPP = electron.mediumPP
			electron.is_loosePP = electron.loosePP

			"""

			#EMPlusPlus computation
			electron.is_tightPP = self.is_tightPP(
				electron.etas2,
        		electron.pt,
        		electron.f3,
        		electron.Ethad/electron.pt,
        		electron.Ethad1/electron.pt,
        		electron.reta,
        		electron.weta2,
        		electron.f1,
        		electron.wstot,
				(electron.emaxs1-electron.Emax2)/(electron.emaxs1+electron.Emax2) if (electron.emaxs1+electron.Emax2)>0. else 0.,
        		electron.deltaeta1,
        		electron.trackd0_physics,
        		electron.TRTHighTOutliersRatio,
        		electron.nTRTHits,
        		electron.nTRTOutliers,
        		electron.nSiHits,
        		electron.nSCTOutliers+electron.nPixelOutliers,
        		electron.nPixHits,
        		electron.nPixelOutliers,
        		electron.nBLHits,
        		electron.nBLayerOutliers,
        		electron.expectHitInBLayer>0.,
        		electron.cl_E*abs(electron.trackqoverp),
        		electron.deltaphi2,
        		(electron.isEM & 2),
				ROOT.egammaMenu.eg2012,
        		False,
        		False,
				)

			electron.is_mediumPP = self.is_mediumPP(
				electron.etas2,
				electron.pt,
				electron.f3,
				electron.Ethad/electron.pt,
				electron.Ethad1/electron.pt,
				electron.reta,
				electron.weta2,
				electron.f1,
				electron.wstot,
				(electron.emaxs1-electron.Emax2)/(electron.emaxs1+electron.Emax2) if (electron.emaxs1+electron.Emax2)>0. else 0.,
				electron.deltaeta1,
				electron.trackd0_physics,
				electron.TRTHighTOutliersRatio,
				electron.nTRTHits,
				electron.nTRTOutliers,
				electron.nSiHits,
				electron.nSCTOutliers+electron.nPixelOutliers,
				electron.nPixHits,
				electron.nPixelOutliers,
				electron.nBLHits,
				electron.nBLayerOutliers,
	        	electron.expectHitInBLayer>0.,
				ROOT.egammaMenu.eg2012,
				False,
				False,
				)

			electron.is_loosePP = self.is_loosePP(
				electron.etas2,
				electron.pt,
				electron.Ethad/electron.pt,
				electron.Ethad1/electron.pt,
				electron.reta,
				electron.weta2,
				electron.f1,
				electron.wstot,
				(electron.emaxs1-electron.Emax2)/(electron.emaxs1+electron.Emax2) if (electron.emaxs1+electron.Emax2)>0. else 0.,
				electron.deltaeta1,
				electron.nSiHits,
				electron.nSCTOutliers+electron.nPixelOutliers,
				electron.nPixHits,
				electron.nPixelOutliers,
				ROOT.egammaMenu.eg2012,
				False,
				False,
				)
			"""
			#Energy scaling
			if event.is_mc:
				self.energy_rescaler.SetRandomSeed(event.EventNumber*100+electronN)
				energyScale = self.energy_rescaler.getSmearingCorrection(
					electron.cl_eta,
					electron.cl_E,
					ROOT.egRescaler.EnergyRescalerUpgrade.Nominal,
					False,
					)
			else:
				energyScale = self.energy_rescaler.applyEnergyCorrection(
					electron.cl_eta,
					electron.cl_E,
					ROOT.egRescaler.EnergyRescalerUpgrade.Electron,
					ROOT.egRescaler.EnergyRescalerUpgrade.Nominal,
					)/electron.cl_E

			#Apply corrections

			electron.pt_corrected = electron.pt*energyScale
			electron.E_corrected = electron.cl_E*energyScale
			electron.set_pt_eta_phi_e(
				electron.pt_corrected,
				electron.eta,
				electron.phi,
				electron.E_corrected,
				)

			#Scale factors
			if event.is_mc:
				recoScaleFactorResult = self.reco_scalefactor_tool.calculate(
					1,
					event.random_RunNumber,
					electron.cl_eta,
					electron.pt_corrected
					)
				mediumIDScaleFactorResult = self.mediumID_scalefactor_tool.calculate(
					1,
					event.random_RunNumber,
					electron.cl_eta,
					electron.pt_corrected
					)
				tightIDScaleFactorResult = self.tightID_scalefactor_tool.calculate(
					1,
					event.random_RunNumber,
					electron.cl_eta,
					electron.pt_corrected
					)
				mediumScaleFactor = recoScaleFactorResult.getScaleFactor()*mediumIDScaleFactorResult.getScaleFactor()
				mediumScaleFactorError = sqrt(
					recoScaleFactorResult.getTotalUncertainty()**2.+\
					mediumIDScaleFactorResult.getTotalUncertainty()**2.
					)
				tightScaleFactor = recoScaleFactorResult.getScaleFactor()*tightIDScaleFactorResult.getScaleFactor()
				tightScaleFactorError = sqrt(
					recoScaleFactorResult.getTotalUncertainty()**2.+\
					tightIDScaleFactorResult.getTotalUncertainty()**2.
					)
			else:
				mediumScaleFactor = 1.
				mediumScaleFactorError = 0.
				tightScaleFactor = 1.
				tightScaleFactorError = 0.

			#Apply scale factor to electron objects
			electron.mediumScaleFactor = mediumScaleFactor
			electron.mediumScaleFactorError = mediumScaleFactorError
			electron.tightScaleFactor = tightScaleFactor
			electron.tightScaleFactorError = tightScaleFactorError

			electron.scale_factor = mediumScaleFactor
			electron.scale_factor_error = mediumScaleFactorError

			#Correct etcone isolation
			electron.etcone20_corrected = self.electron_isolation_tool.GetPtNPVCorrectedIsolation(
				event.nPV_2trks,
				electron.E_corrected,
				electron.etas2,
				electron.etap,
				electron.cl_eta,
				30,
				event.is_mc,
				electron.Etcone20,
				False,
				ROOT.CaloIsoCorrection.ELECTRON
				)

	def initialize_tools(self):

		analysis_home = os.getenv('ANALYSISHOME')

		load('egammaAnalysisUtils')
		#etcone correction tools
		self.electron_isolation_tool = ROOT.CaloIsoCorrection()
		#EMPlusPlus correction tools
		self.is_tightPP = ROOT.isTightPlusPlus
		self.is_mediumPP = ROOT.isMediumPlusPlus
		self.is_loosePP = ROOT.isLoosePlusPlus
		#energy rescaling tools
		self.energy_rescaler = ROOT.egRescaler.EnergyRescalerUpgrade()
		self.energy_rescaler.Init(
			"{0}/external/egammaAnalysisUtils/share/EnergyRescalerData.root".format(analysis_home),
			"2012",
			"es2010"
			)

		load('ElectronEfficiencyCorrection')
		#reco scalefactor
		self.reco_scalefactor_tool = ROOT.Root.TElectronEfficiencyCorrectionTool()
		self.reco_scalefactor_tool.addFileName("{0}/external/ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.v02.root".format(analysis_home))
		self.reco_scalefactor_tool.initialize()
		#medium ID scalefactor
		self.mediumID_scalefactor_tool = ROOT.Root.TElectronEfficiencyCorrectionTool()
		self.mediumID_scalefactor_tool.addFileName("{0}/external/ElectronEfficiencyCorrection/data/efficiencySF.offline.Medium.2012.8TeV.rel17p2.v02.root".format(analysis_home))
		self.mediumID_scalefactor_tool.initialize()
		#tight ID scalefactor
		self.tightID_scalefactor_tool = ROOT.Root.TElectronEfficiencyCorrectionTool()
		self.tightID_scalefactor_tool.addFileName("{0}/external/ElectronEfficiencyCorrection/data/efficiencySF.offline.Tight.2012.8TeV.rel17p2.v02.root".format(analysis_home))
		self.tightID_scalefactor_tool.initialize()
