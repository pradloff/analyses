from common.functions import event_function
from common.particle import particle
from math import cosh,sqrt
from common.external import load
import ROOT
import os
from misc import list_attributes
from common.node import node
from common.analysis import analysis
from math import cos,sin,atan,exp
import time

from mc import truth_tree

from common.analysis import analysis
from common.functions import event_function
from common.external import load

from misc import count_primary_vertices
from pileup import pileup_weighting
from muons import collect_muons
from electrons import collect_electrons
from taus import collect_taus
from jets import collect_jets, collect_tracks
from overlap import remove_overlap
from met import correct_missing_energy
from metadata import lumi

import ROOT
import os
from math import sqrt

class select_mumu(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(),
			identify_z_leptons(mode=0),
			count_primary_vertices(),
			pileup_weighting(),
			collect_muons(),
			collect_electrons(),
			collect_taus(),
			collect_tracks(),
			collect_jets(),
			remove_overlap(),
			correct_missing_energy(),
			match(0.1),
			trigger(),
			preselection(),
			)

		self.add_result_function(
			)

		self.add_meta_result_function(
			lumi()
			)

class select_tautau(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(),
			identify_z_leptons(mode=1),
			count_primary_vertices(),
			pileup_weighting(),
			collect_muons(),
			collect_electrons(),
			collect_taus(),
			collect_tracks(),
			collect_jets(),
			remove_overlap(),
			correct_missing_energy(),
			match(0.1),
			trigger(),
			preselection(),
			)

		self.add_result_function(
			)

		self.add_meta_result_function(
			lumi()
			)

from itertools import product

class identify_z_leptons(event_function):

	def __init__(self,mode):

		self.mode = mode

		event_function.__init__(self)

		self.names = [
			'pt',
			'eta',
			'phi',
			'm',
			'charge',
			]

		self.required_branches += ['truth']
		self.create_branches['l1'] = None
		self.create_branches['l2'] = None
		self.create_branches['mode'] = None
		for l,name in product(['l1','l2'],self.names):
			self.create_branches[l+'_'+name] = 'float'

	def __call__(self,event):

		event.mode = self.mode

		try: z = [p for p in event.truth.values() if p().pdgId==23 and p().status in [2,10902]][0]
		except IndexError:
			print 'Z could not be found!'
			event.__break__=True
			return

		if event.mode==0: #mumu
			try: event.l1,event.l2 = [c() for c in z.children if abs(c().pdgId)==13 and c().status==1]
			except IndexError:
				print 'Muons could not be found'
				event.__break__=True
				return

			if event.l1.pt<event.l2.pt:
				event.l1,event.l2 = event.l2,event.l1 #swap pt ordered

		elif event.mode==1: #tautau->emu
			try: tau1,tau2 = [p for p in z.children if abs(p().pdgId)==15 and p().status==2]
			except IndexError:
				print 'Taus could not be found'
				event.__break__=True
				return
			try:
				event.l1 = [c() for c in tau1.children+tau2.children if abs(c().pdgId)==11][0] #electron
				event.l2 = [c() for c in tau1.children+tau2.children if abs(c().pdgId)==13][0] #muon
			except IndexError:
				event.__break__=True
				return

		else:
			print 'Unrecognized mode'
			event.__break__=True
			return

		for l,name in product(['l1','l2'],self.names):
			event.__dict__[l+'_'+name] = event.__dict__[l].__dict__[name]

class match(event_function):
	
	def __init__(self,dR_max):

		self.dR_max = dR_max

		event_function.__init__(self)

		self.required_branches += [
			'electrons',
			'muons',
			'l1',
			'l2',
			'mode',
			]

		self.names = [
			('eta','float',0.),
			('phi','float',0.),
			('pt','float',0.),
			('E','float',0.),
			('charge','float',0.),
			('ptcone40','float',0.),
			('etcone20','float',0.),
			('passed_preselection','bool',False),
			('passed_selection','bool',False),
			('scale_factor','float',0.),
			('scale_factor_error','float',0.),
			]

		for l,(name,type_,_) in product(['l1','l2'],self.names):
			self.create_branches[l+'_offline_'+name] = type_

		self.dummy = particle(**dict((name,default) for name,_,default in self.names))

	def __call__(self,event):
		if event.mode==0:
			#find matching to leading muon
			try: event.l1.offline_match = sorted([muon for muon in event.muons.values() if muon().DeltaR(event.l1())<self.dR_max and muon.passed_preselection and muon.pt_corrected>20000.], key=lambda mu: mu().DeltaR(event.l1()))[0]
			except IndexError: event.l1.offline_match = self.dummy
			#find matching to subleading muon
			try: event.l2.offline_match = sorted([muon for muon in event.muons.values() if muon().DeltaR(event.l2())<self.dR_max and muon != event.l1.offline_match and muon.passed_preselection and muon], key=lambda mu: mu().DeltaR(event.l2()))[0]
			except IndexError: event.l2.offline_match = self.dummy
		elif event.mode==1:
			#find matching to electron
			try: event.l1.offline_match = sorted([electron for electron in event.electrons.values() if electron().DeltaR(event.l1())<self.dR_max and electron.passed_preselection], key=lambda el: el().DeltaR(event.l1()))[0]
			except IndexError: event.l1.offline_match = self.dummy
			#find matching muon
			try: event.l2.offline_match = sorted([muon for muon in event.muons.values() if muon().DeltaR(event.l2())<self.dR_max and muon.passed_preselection], key=lambda mu: mu().DeltaR(event.l2()))[0]
			except IndexError: event.l2.offline_match = self.dummy
		else:
			print 'Unrecognized mode'
			event.__break__=True
			return

		if event.l1.offline_match!=self.dummy: 
			event.l1.offline_match.pt = event.l1.offline_match.pt_corrected
			event.l1.offline_match.E = event.l1.offline_match.E_corrected


		if event.l2.offline_match!=self.dummy: 
			event.l2.offline_match.pt = event.l2.offline_match.pt_corrected
			event.l2.offline_match.E = event.l2.offline_match.E_corrected

		for l,(name,type_,_) in product(['l1','l2'],self.names):
			event.__dict__[l+'_offline_'+name] = event.__dict__[l].offline_match.__dict__[name]

class trigger(event_function):
	
	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'EF_mu18_tight_mu8_EFFS',
			'EF_e12Tvh_medium1_mu8',
			'random_RunNumber',
			]

		self.periods_runnumbers = {
			"A_":(200804,201556),
			"B_":(202660,205113),
			"CtoE_":(206248,210308),
			"G_":(211522,212272),
			"HtoL_":(212619,215643),
			}

		self.create_branches['triggered'] = 'bool'
		self.create_branches['trigger_scale_factor'] = 'float'
		self.create_branches['trigger_scale_factor_error'] = 'float'

		self.initialize_tools()

	def __call__(self,event):

		if event.mode == 0:
			event.triggered = event.EF_mu18_tight_mu8_EFFS
		
		if event.mode == 1:
			event.triggered = event.EF_e12Tvh_medium1_mu8

		if event.l1.offline_match.passed_preselection and event.l2.offline_match.passed_preselection: self.apply_corrections(event)
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

		self.config_muon_trigger_mu8.runNumber = run
		self.config_muon_trigger_mu8.period = period
		self.config_muon_trigger_mu18_tight_mu8_EFFS.runNumber = run
		self.config_muon_trigger_mu18_tight_mu8_EFFS.period = period

		#mumu
		if event.mode == 0:
			muon1 = event.l1.offline_match()
			muon2 = event.l1.offline_match()

			muons = ROOT.std.vector('TLorentzVector')()
			muons.push_back(muon1)
			muons.push_back(muon2)

			muons_quality = ROOT.std.vector('int')()
			muons_quality.push_back(1)
			muons_quality.push_back(1)

			#get data efficiency
			self.config_muon_trigger_mu18_tight_mu8_EFFS.isData = True
			result_data = self.muon_trigger_mu18_tight_mu8_EFFS.getDimuonEfficiency(
				self.config_muon_trigger_mu18_tight_mu8_EFFS,
				muons,
				muons_quality,
				self.config_muon_trigger_mu18_tight_mu8_EFFS.trigger
				)
			#get mc efficiency
			self.config_muon_trigger_mu18_tight_mu8_EFFS.isData = False
			result_mc = self.muon_trigger_mu18_tight_mu8_EFFS.getDimuonEfficiency(
				self.config_muon_trigger_mu18_tight_mu8_EFFS,
				muons,
				muons_quality,
				self.config_muon_trigger_mu18_tight_mu8_EFFS.trigger
				)

			event.trigger_scale_factor = result_data.first/result_mc.first
			event.trigger_scale_factor_error = max([
				abs(event.trigger_scale_factor-((result_data.first+result_data.second)/(result_mc.first-result_mc.second))),
				abs(event.trigger_scale_factor-((result_data.first-result_data.second)/(result_mc.first+result_mc.second)))
				])		

		#emu
		if event.mode == 1:
			#get electron scale factor
			electron = event.l1.offline_matched()
         		electron_scale_factor = self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",0)
			electron_scale_factor_error = max([
				abs(electron_scale_factor-self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",1)),
				abs(electron_scale_factor-self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",-1))
				])
			#get muon scale factor
			muon = event.l2.offline_matched()
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
		#scale factor tool and config for EF_mu18_tight_mu8_EFFS
		self.muon_trigger_mu18_tight_mu8_EFFS = ROOT.LeptonTriggerSF(
			2012,
			'{0}/external/TrigMuonEfficiency/share'.format(analysis_home),
			'muon_trigger_sf_2012_AtoL.p1328.root',
			'{0}/external/ElectronEfficiencyCorrection/data'.format(analysis_home),
			"rel17p2.v02")
		self.config_muon_trigger_mu18_tight_mu8_EFFS = ROOT.TrigMuonEff.Configuration(True,False,False,False,-1,-1,0,'mu18_tight_mu8_EFFS','period_','fine')

		#scale factor tool and config for EF_mu8
		self.muon_trigger_mu8 = ROOT.LeptonTriggerSF(
			2012,
			'{0}/external/TrigMuonEfficiency/share'.format(analysis_home),
			'muon_trigger_sf_2012_AtoL.p1328.root',
			'{0}/external/ElectronEfficiencyCorrection/data'.format(analysis_home),
			"rel17p2.v02")
		self.config_muon_trigger_mu8 = ROOT.TrigMuonEff.Configuration(True,False,False,False,-1,-1,0,'mu8','period_','fine')

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

		self.create_branches['top_hfor_type'] = 'int'

		self.initialize_tools()

	def __call__(self,event):

		event.top_hfor_type = getattr(event,'top_hfor_type',-1)

		#2 preselection leptons, no hadronic taus, at least one preselection jet
		if not all([
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

		load('TileTripReader')
		self.tile_trip_reader = ROOT.Root.TTileTripReader()
		#comment

	def initialize_tools(self):
		load('TileTripReader')
		self.tile_trip_reader = ROOT.Root.TTileTripReader()
		#comment
