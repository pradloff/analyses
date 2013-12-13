from common.analysis import analysis
from common.functions import event_function
from common.external import load

from misc import count_primary_vertices
from pileup import pileup_weighting
from muons import collect_muons
from electrons import collect_electrons
from selection import trigger, preselection
from taus import collect_taus
from jets import collect_jets
from overlap import remove_overlap
from met import correct_missing_energy
from metadata import lumi

import ROOT
import os
from math import sqrt

class make_preselection(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			count_primary_vertices(),
			pileup_weighting(),
			collect_muons(),
			collect_electrons(),
			collect_taus(),
			collect_jets(),
			remove_overlap(),
			correct_missing_energy(),
			preselection(),
			trigger(),
			)

		self.add_result_function(
			)

		self.add_meta_result_function(
			lumi()
			)

class trigger(event_function):
	
	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'EF_mu24i_tight',
			'EF_mu36_tight',
			'EF_e12Tvh_medium1_mu8',
			'EF_e24vhi_medium1',
			'EF_e60_medium1',
			'random_RunNumber',
			]

		self.periods_runnumbers = {
			"A_":(200804,201556),
			"B_":(202660,205113),
			"CtoE_":(206248,210308),
			"G_":(211522,212272),
			"HtoL_":(212619,215643),
			}

		self.create_branches['trigger_scale_factor'] = 'float'
		self.create_branches['trigger_scale_factor_error'] = 'float'

		self.initialize_tools()

	def __call__(self,event):

		if event.lepton_class == 0:
			if not any([
				event.EF_e24vhi_medium1 and event.l1_pt>25000.,
				event.EF_e60_medium1 and event.l1_pt>65000.,
				]):
				event.__break__ = True
				return

		if event.lepton_class == 1:
			if not any([
				event.EF_mu24i_tight and event.l1_pt>25000.,
				event.EF_mu36_tight and event.l1_pt>40000.,
				]):
				event.__break__ = True
				return
		
		if event.lepton_class == 2:
			if not event.EF_e12Tvh_medium1_mu8:
				event.__break__ = True
				return

		if event.is_mc: self.apply_corrections(event)
		else: 
			event.trigger_scale_factor = 1.
			event.trigger_scale_factor_error = 0.

	def apply_corrections(self,event):

		#Update configs
		run = event.random_RunNumber
		period = None
		for period_,(run_begin,run_end) in self.periods_runnumbers.items():
			if run_begin<=run<=run_end:
				period = period_
				break

		if period is None: 
			print 'Period could not be found for run {0}'.format(run)
			event.__break__=True
			return

		self.config_muon_trigger_mu24i_tight.runNumber = run
		self.config_muon_trigger_mu24i_tight.period = period
		self.config_muon_trigger_mu8.runNumber = run
		self.config_muon_trigger_mu8.period = period

		#ee
		if event.lepton_class == 0:
		  	result = self.electron_trigger_e24vhi_medium.calculate(1,run,event.l1.cl_eta,event.l1.pt)
		  	event.trigger_scale_factor = result.getScaleFactor()
			event.trigger_scale_factor_error = result.getTotalUncertainty()

		#mumu
		if event.lepton_class == 1:
			muon = event.l1()
			#get data efficiency
			self.config_muon_trigger_mu24i_tight.isData = True
			result_data = self.muon_trigger_mu24i_tight.getMuonEfficiency(
				self.config_muon_trigger_mu24i_tight,
				muon,
				1,
				self.config_muon_trigger_mu24i_tight.trigger
				)
			#get mc efficiency
			self.config_muon_trigger_mu24i_tight.isData = False
			result_mc = self.muon_trigger_mu24i_tight.getMuonEfficiency(
				self.config_muon_trigger_mu24i_tight,
				muon,
				1,
				self.config_muon_trigger_mu24i_tight.trigger
				)

			event.trigger_scale_factor = result_data.first/result_mc.first
			event.trigger_scale_factor_error = max([
				abs(event.trigger_scale_factor-((result_data.first+result_data.second)/(result_mc.first-result_mc.second))),
				abs(event.trigger_scale_factor-((result_data.first-result_data.second)/(result_mc.first+result_mc.second)))
				])		

		#emu
		if event.lepton_class == 2:
			#get electron scale factor
			electron = event.l1()
         		electron_scale_factor = self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",0)
			electron_scale_factor_error = max([
				abs(electron_scale_factor-self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",1)),
				abs(electron_scale_factor-self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",-1))
				])
			#get muon scale factor
			muon = event.l2()
			#get data efficiency
			self.config_muon_trigger_mu8.isData = True
			result_data = self.muon_trigger_mu8.getMuonEfficiency(
				self.config_muon_trigger_mu8,
				muon,
				1,
				self.config_muon_trigger_mu8.trigger
				)
			#get mc efficiency
			self.config_muon_trigger_mu8.isData = False
			result_mc = self.muon_trigger_mu8.getMuonEfficiency(
				self.config_muon_trigger_mu8,
				muon,
				1,
				self.config_muon_trigger_mu8.trigger
				)

			muon_scale_factor = result_data.first/result_mc.first
			muon_scale_factor_error = max([
				abs(muon_scale_factor-((result_data.first+result_data.second)/(result_mc.first-result_mc.second))),
				abs(muon_scale_factor-((result_data.first-result_data.second)/(result_mc.first+result_mc.second)))
				])			

			event.trigger_scale_factor = electron_scale_factor*muon_scale_factor
			event.trigger_scale_factor_error = sqrt((electron_scale_factor*muon_scale_factor_error)**2.+(muon_scale_factor*electron_scale_factor_error)**2.)

	def initialize_tools(self):
		analysis_home = os.getenv('ANALYSISHOME')

		load('TrigMuonEfficiency')
		#scale factor tool and config for EF_mu24i_tight
		self.muon_trigger_mu24i_tight = ROOT.LeptonTriggerSF(
			2012,
			'{0}/external/TrigMuonEfficiency/share'.format(analysis_home),
			'muon_trigger_sf_2012_AtoL.p1328.root',
			'{0}/external/ElectronEfficiencyCorrection/data'.format(analysis_home),
			"rel17p2.v02")
		self.config_muon_trigger_mu24i_tight = ROOT.TrigMuonEff.Configuration(True,False,False,False,-1,-1,0,'mu24i_tight','period_','fine')

		#scale factor tool and config for EF_mu8
		self.muon_trigger_mu8 = ROOT.LeptonTriggerSF(
			2012,
			'{0}/external/TrigMuonEfficiency/share'.format(analysis_home),
			'muon_trigger_sf_2012_AtoL.p1328.root',
			'{0}/external/ElectronEfficiencyCorrection/data'.format(analysis_home),
			"rel17p2.v02")
		self.config_muon_trigger_mu8 = ROOT.TrigMuonEff.Configuration(True,False,False,False,-1,-1,0,'mu8','period_','fine')

		load('ElectronEfficiencyCorrection')
		#scale factor tool for EF_e24vhi_medium
		self.electron_trigger_e24vhi_medium = ROOT.Root.TElectronEfficiencyCorrectionTool()
		self.electron_trigger_e24vhi_medium.addFileName("{0}/external/ElectronEfficiencyCorrection/data/efficiencySF.e24vhi_medium1_e60_medium1.Medium.2012.8TeV.rel17p2.v02.root".format(analysis_home))
		self.electron_trigger_e24vhi_medium.initialize()

		load('HSG4LepLepTriggerSF')
		#scale factor tool for EF_e12Tvh_medium1
		self.electron_trigger_e12Tvh_medium1 = ROOT.HSG4LepLepTriggerSF("{0}/external/HSG4LepLepTriggerSF/data/".format(analysis_home),False)

 
class preselection(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'electrons',
			'muons',
			'jets',
			'larError',
			'tileError',
			'coreFlags',
			'lbn',
			'random_RunNumber',
			'EventNumber',
			]
		self.create_branches.update(dict((key,value) for key,value in [
			('l1_eta','float'),
			('l1_phi','float'),
			('l1_pt','float'),
			('l1_E','float'),
			('l1_charge','float'),
			('l1_ptcone40','float'),
			('l1_etcone20','float'),
			('l1_passed_selection','bool'),
			('l1_scale_factor','float'),
			('l1_scale_factor_error','float'),
			]))

		self.create_branches.update(dict((key,value) for key,value in [
			('l2_eta','float'),
			('l2_phi','float'),
			('l2_pt','float'),
			('l2_E','float'),
			('l2_charge','float'),
			('l2_ptcone40','float'),
			('l2_etcone20','float'),
			('l2_passed_selection','bool'),
			('l2_scale_factor','float'),
			('l2_scale_factor_error','float'),
			]))

		self.create_branches.update(dict((key,value) for key,value in [
			('lepton_class','int'),
			]))

		self.initialize_tools()

	def __call__(self,event):

		#2 preselection leptons, no hadronic taus, at least one preselection jet
		if not all([
			sum(1 for lepton in event.electrons.values()+event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed)==2,
			sum(1 for tau in event.taus.values() if tau.passed_preselection and not tau.overlap_removed)==0,
			sum(1 for jet in event.jets.values() if jet.passed_preselection and not jet.overlap_removed)>0,
			sum(1 for jet in event.jets.values() if jet.passed_preselection and not jet.overlap_removed and jet.isBadLooseMinus)==0,
			event.larError!=2,
			event.tileError!=2,
			(event.coreFlags&0x40000)==0,
			self.tile_trip_reader.checkEvent(event.random_RunNumber,event.lbn,event.EventNumber),
			]):
			event.__break__=True
			return

		#ee
		if sum(1 for lepton in event.electrons.values() if lepton.passed_preselection and not lepton.overlap_removed)==2:
			event.l1,event.l2 = [lepton for lepton in event.electrons.values() if lepton.passed_preselection and not lepton.overlap_removed]
			if event.l2.pt>event.l1.pt: event.l1,event.l2 = event.l2,event.l1
			event.lepton_class = 0
		#mumu
		elif sum(1 for lepton in event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed)==2:
			event.l1,event.l2 = [lepton for lepton in event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed]
			if event.l2.pt>event.l1.pt: event.l1,event.l2 = event.l2,event.l1
			event.lepton_class = 1
		#emu
		else:
			event.l1,event.l2 = [lepton for lepton in event.electrons.values()+event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed]
			event.lepton_class = 2

		event.l1_eta = event.l1.eta
		event.l1_phi = event.l1.eta
		event.l1_pt = event.l1.pt_corrected
		event.l1_E = event.l1.E_corrected
		event.l1_charge = event.l1.charge
		event.l1_ptcone40 = event.l1.ptcone40
		event.l1_etcone20 = event.l1.etcone20_corrected
		event.l1_passed_selection = event.l1.passed_selection
		event.l1_scale_factor = event.l1.scale_factor
		event.l1_scale_factor_error = event.l1.scale_factor_error

		event.l2_eta = event.l2.eta
		event.l2_phi = event.l2.eta
		event.l2_pt = event.l2.pt_corrected
		event.l2_E = event.l2.E_corrected
		event.l2_charge = event.l2.charge
		event.l2_ptcone40 = event.l2.ptcone40
		event.l2_etcone20 = event.l2.etcone20_corrected
		event.l2_passed_selection = event.l2.passed_selection
		event.l2_scale_factor = event.l2.scale_factor
		event.l2_scale_factor_error = event.l2.scale_factor_error

	def initialize_tools(self):
		load('TileTripReader')
		self.tile_trip_reader = ROOT.Root.TTileTripReader()
		#comment
