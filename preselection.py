from common.analysis import analysis
from common.functions import event_function,EventBreak
from common.external import load

from misc import count_primary_vertices
from pileup import pileup_weighting
from muons import collect_muons
from electrons import collect_electrons
from taus import collect_taus
from jets import collect_jets, collect_tracks
from overlap import remove_overlap
#from met import correct_missing_energy
from metadata import lumi
from common.branches import auto_branch,branch

import ROOT
import os
from math import sqrt
import json

class make_preselection(analysis):
	def __init__(self):
		super(make_preselection,self).__init__()

		self.add_event_function(
			count_primary_vertices(),
			pileup_weighting(),
			#get_weight_pileup(),
			collect_muons(),
			collect_electrons(),
			collect_taus(),
			collect_tracks(),
			collect_jets(),
			remove_overlap(),
			#correct_missing_energy(),
			preselection(),
			trigger(),
			)

		self.add_result_function(
			)
		"""
		self.add_meta_result_function(
			lumi()
			)
		"""
class make_preselection_mumu_embedding(analysis):
	def __init__(self):
		super(make_preselection_mumu_embedding,self).__init__()

		self.add_event_function(
			count_primary_vertices(),
			pileup_weighting(),
			#get_weight_pileup(),
			collect_muons(),
			collect_electrons(),
			collect_taus(),
			collect_tracks(),
			collect_jets(),
			remove_overlap(),
			#correct_missing_energy(),
			preselection(),
			trigger_mumu_embed(),
			)

		self.add_result_function(
			)
		"""
		self.add_meta_result_function(
			lumi()
			)
		"""

class get_weight_pileup(event_function):
	def __init__(self):
		super(get_weight_pileup,self).__init__()
		self.branches.append(branch('weight_pileup', 'r'))

	def __call__(self,event):
		super(get_weight_pileup,self).__call__(event)

		for weight in [
			event.weight_pileup,
			]: event.__weight__*=weight

class get_weight(event_function):
	def __init__(self):
		super(get_weight,self).__init__()
		self.branches.append(branch('mc_channel_number', 'r'))

		self.initialize()

	def __call__(self,event):
		super(get_weight,self).__call__(event)

		if event.mc_channel_number == 0: lumi_event_weight = 1.
		else: lumi_event_weight = self.mc_lumi_info['lumi_event_weight'][str(event.mc_channel_number)] #= Lumi_data*(xsec*k_factor)/N_gen / 1 for data
		for weight in [
			lumi_event_weight,
			]: event.__weight__*=weight

	def initialize(self):
		analysis_home = os.getenv('ANALYSISHOME')
		mc_lumi_file = '{0}/data/mc_lumi.json'.format(analysis_home)
		with open(mc_lumi_file) as f: self.mc_lumi_info = json.loads(f.read())


class trigger_mumu_embed(event_function):

	class muon_event(EventBreak): pass
	class trigger(EventBreak): pass

	def __init__(self):
		super(trigger_mumu_embed,self).__init__()

		self.break_exceptions += [
			trigger_mumu_embed.muon_event,
			trigger_mumu_embed.trigger,
			]

		self.branches.append(branch('EF_mu18_tight_mu8_EFFS', 'r'))
		#self.branches.append(branch('random_RunNumber', 'r'))

		self.periods_runnumbers = {
			"A_":(200804,201556),
			"B_":(202660,205113),
			"CtoE_":(206248,210308),
			"G_":(211522,212272),
			"HtoL_":(212619,215643),
			}


		self.branches.append(auto_branch('trigger_scale_factor','w','Float_t'))
		self.branches.append(auto_branch('trigger_scale_factor','w','Float_t'))


	def __call__(self,event):
		super(trigger_mumu_embed,self).__call__(event)

		if event.lepton_class != 1: raise trigger_mumu_embed.muon_event()

		if not any([
			event.EF_mu18_tight_mu8_EFFS and event.l1_pt>20000. and event.l2_pt>10000.,
			]): raise trigger_mumu_embed.trigger()

		if event.is_mc:
			event.trigger_scale_factor = 1.
			event.trigger_scale_factor_error = 0.
			#raise RuntimeError('mumu embed trigger should not be run on MC')
		else:
			event.trigger_scale_factor = 1.
			event.trigger_scale_factor_error = 0.

class trigger(event_function):

	class trigger(EventBreak): pass

	def __init__(self):
		super(trigger,self).__init__()

		self.break_exceptions += [
			trigger.trigger,
			]

		for name in [
			'EF_mu24i_tight',
			'EF_mu36_tight',
			'EF_e12Tvh_medium1_mu8',
			'EF_e24vhi_medium1',
			'EF_e60_medium1',
			#'random_RunNumber',
			]:
			self.branches.append(branch(name, 'r'))

		self.branches.append(auto_branch('trigger_scale_factor', 'w','Float_t'))
		self.branches.append(auto_branch('trigger_scale_factor_error', 'w','Float_t'))

		self.periods_runnumbers = {
			"A_":(200804,201556),
			"B_":(202660,205113),
			"CtoE_":(206248,210308),
			"G_":(211522,212272),
			"HtoL_":(212619,215643),
			}

		self.initialize_tools()

	def __call__(self,event):
		super(trigger,self).__call__(event)

		if event.lepton_class == 0:
			if not any([
				event.EF_e24vhi_medium1 and event.l1_pt>25000.,
				event.EF_e60_medium1 and event.l1_pt>65000.,
				]):
				raise trigger.trigger()

		if event.lepton_class == 1:
			if not any([
				event.EF_mu24i_tight and event.l1_pt>25000.,
				event.EF_mu36_tight and event.l1_pt>40000.,
				]):
				raise trigger.trigger()

		if event.lepton_class == 2:
			if not event.EF_e12Tvh_medium1_mu8:
				raise trigger.trigger()

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

		if period is None: raise RuntimeError('Period could not be found for run {0}'.format(run))

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

	class two_leptons(EventBreak): pass
	class no_taus(EventBreak): pass
	class one_jet(EventBreak): pass
	class no_bad_jets(EventBreak): pass
	class no_bad_flags(EventBreak): pass
	class no_tile_trip(EventBreak): pass

	def __init__(self):
		super(preselection,self).__init__()

		self.break_exceptions += [
			preselection.two_leptons,
			preselection.no_taus,
			preselection.one_jet,
			preselection.no_bad_jets,
			preselection.no_bad_flags,
			preselection.no_tile_trip,
			]


		for name in [
			'larError',
			'tileError',
			'coreFlags',
			'lbn',
			#'random_RunNumber',
			'EventNumber'
			]:
			self.branches.append(branch(name, 'r'))

		for name,branch_type in [
			('l1_eta','Float_t'),
			('l1_phi','Float_t'),
			('l1_pt','Float_t'),
			('l1_E','Float_t'),
			('l1_charge','Float_t'),
			('l1_ptcone40','Float_t'),
			('l1_etcone20','Float_t'),
			('l1_passed_selection','Bool_t'),
			('l1_scale_factor','Float_t'),
			('l1_scale_factor_error','Float_t'),
			('l2_eta','Float_t'),
			('l2_phi','Float_t'),
			('l2_pt','Float_t'),
			('l2_E','Float_t'),
			('l2_charge','Float_t'),
			('l2_ptcone40','Float_t'),
			('l2_etcone20','Float_t'),
			('l2_passed_selection','Bool_t'),
			('l2_scale_factor','Float_t'),
			('l2_scale_factor_error','Float_t'),
			('lepton_class','Int_t'),
			('top_hfor_type','Int_t'),
			]:
			self.branches.append(auto_branch(name,'w',branch_type))

		self.initialize_tools()

	def __call__(self,event):
		super(preselection,self).__call__(event)

		event.top_hfor_type = getattr(event,'top_hfor_type',-1)

		#2 preselection leptons, no hadronic taus, at least one preselection jet
		for requirement,exception in [
			(sum(1 for lepton in event.electrons.values()+event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed)==2,preselection.two_leptons),
			(sum(1 for tau in event.taus.values() if tau.passed_preselection and not tau.overlap_removed)==0,preselection.no_taus),
			(sum(1 for jet in event.jets.values() if jet.passed_preselection and not jet.overlap_removed)>0,preselection.one_jet),
			(sum(1 for jet in event.jets.values() if jet.passed_preselection and not jet.overlap_removed and jet.isBadLooseMinus)==0,preselection.no_bad_jets),
			(event.larError!=2 and event.tileError!=2 and (event.coreFlags&0x40000)==0,preselection.no_bad_flags),
			(self.tile_trip_reader.checkEvent(event.random_RunNumber,event.lbn,event.EventNumber),preselection.no_tile_trip),
			]:
			if not requirement: raise exception()

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
		event.l1_phi = event.l1.phi
		event.l1_pt = event.l1.pt_corrected
		event.l1_E = event.l1.E_corrected
		event.l1_charge = event.l1.charge
		event.l1_ptcone40 = event.l1.ptcone40
		event.l1_etcone20 = event.l1.etcone20_corrected
		event.l1_passed_selection = event.l1.passed_selection
		event.l1_scale_factor = event.l1.scale_factor
		event.l1_scale_factor_error = event.l1.scale_factor_error

		event.l2_eta = event.l2.eta
		event.l2_phi = event.l2.phi
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
