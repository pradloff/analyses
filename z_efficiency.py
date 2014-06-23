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
from common.functions import event_function,result_function
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
import json
from itertools import product

truth_tree(pdgIds = [11,-11,23])

class select_ee(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(pdgIds = [11,-11,23]),
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

class select_mumu(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(pdgIds = [13,-13,23]),
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

class select_tautau(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(pdgIds = [15,-15,23]),
			identify_z_leptons(mode=2),
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

class leplep_efficiency(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			#build_events(),
			get_weight(),
			)

		self.add_result_function(
			efficiency(),
			)

		self.add_meta_result_function(
			)

class build_events(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'triggered'
			]

		self.lepton_names = [
			'eta',
			'pt',
			'offline_eta',
			'offline_passed_preselection',
			'offline_pt',
			'offline_scale_factor',
			'offline_scale_factor_error',
			]
		self.required_branches += ['l1_'+name for name in self.lepton_names]
		self.required_branches += ['l2_'+name for name in self.lepton_names]

		self.jet_names = [
			'E',
			'eta',
			'phi',
			'pt',
			'flavor_weight_MV1',
			'jvf',
			'bJet_scale_factor',
			]

		self.create_branches['top_hfor_type'] = None

		self.required_branches += ['jet_'+name for name in self.jet_names]
		self.required_branches += ['jet_n']

		self.create_branches['jets'] = None
		self.create_branches['l1'] = None
		self.create_branches['l2'] = None

	def __call__(self,event):

		#collect jets
		event.jets = {}		
		for jet in range(event.jet_n):
			if not all([
				not ((abs(event.jet_eta[jet])<2.4 and event.jet_pt[jet]<50000.) and not ((event.jet_jvf[jet])>0.5)),
				]): continue
			event.jets[jet] = particle(\
				**dict((name,event.__dict__['jet_'+name][jet]) for name in self.jet_names)
				)
			event.jets[jet].set_pt_eta_phi_e(
				event.jets[jet].pt,
				event.jets[jet].eta,
				event.jets[jet].phi,
				event.jets[jet].E,
				)

		event.bjets_preselected = {}
		for jet in range(event.jet_n):
			if not all([
				jet in event.jets,
				abs(event.jet_eta[jet])<2.4,
				]): continue
			event.bjets_preselected[jet] = particle(\
				**dict((name,event.__dict__['jet_'+name][jet]) for name in self.jet_names)
				)
			event.bjets_preselected[jet].set_pt_eta_phi_e(
				event.jets[jet].pt,
				event.jets[jet].eta,
				event.jets[jet].phi,
				event.jets[jet].E,
				)

		event.bjets = {}
		for jet in range(event.jet_n):
			if not all([
				jet in event.bjets_preselected,
				event.jet_flavor_weight_MV1[jet] > 0.7892,
				]): continue
			event.bjets[jet] = particle(\
				**dict((name,event.__dict__['jet_'+name][jet]) for name in self.jet_names)
				)
			event.bjets[jet].set_pt_eta_phi_e(
				event.jets[jet].pt,
				event.jets[jet].eta,
				event.jets[jet].phi,
				event.jets[jet].E,
				)

		#if not event.jets:
		#	event.__break__ = True
		#	return

		if getattr(event,'top_hfor_type',0)==4:
			event.__break__ = True
			return

class get_weight(event_function):
	def __init__(self):
		event_function.__init__(self)
		

		self.required_branches += [
			'triggered'
			]

		self.lepton_names = [
			'eta',
			'pt',
			'offline_eta',
			'offline_passed_preselection',
			'offline_pt',
			'offline_scale_factor',
			'offline_scale_factor_error',
			]
		self.required_branches += ['l1_'+name for name in self.lepton_names]
		self.required_branches += ['l2_'+name for name in self.lepton_names]
			
		self.required_branches += [
			'mc_channel_number',
			'trigger_scale_factor',
			'trigger_scale_factor_error',
			'weight_pileup',
			]

		self.initialize()

	def __call__(self,event):
		if event.mc_channel_number == 0: lumi_event_weight = 1.
		else: lumi_event_weight = self.mc_lumi_info['lumi_event_weight'][str(event.mc_channel_number)] #= Lumi_data*(xsec*k_factor)/N_gen / 1 for data
		for weight in [
			lumi_event_weight,
			#event.l1_offline_scale_factor if event.l1_offline_scale_factor>1./100000. else 1.,
			#event.l2_offline_scale_factor if event.l2_offline_scale_factor>1./100000. else 1.,
			#event.trigger_scale_factor if event.trigger_scale_factor>1./100000. else 1.,
			event.weight_pileup,
			]: event.__weight__*=weight
		#event.__weight__*=reduce(mul,[jet.bJet_scale_factor for jet in event.jets.values()],1)

	def initialize(self):
		analysis_home = os.getenv('ANALYSISHOME')
		mc_lumi_file = '{0}/data/mc_lumi_embedding.json'.format(analysis_home)
		with open(mc_lumi_file) as f: self.mc_lumi_info = json.loads(f.read())

import array

class efficiency(result_function):

	def __init__(self):
		result_function.__init__(self)

		etas = [
			0.,
			0.1,
			0.5,
			1.0,
			1.37,
			1.52,
			1.8,
			2.0,
			2.2,
			2.4,
			2.5,
			]

		pts = [
			10.,
			12.,
			14.,
			16.,
			18.,
			20.,
			22.,
			24.,
			26.,
			30.,
			34.,
			40.,
			50.,
			70.,
			90.,
			140.,
			200.,
			1000.,
			]

		self.eta_bins = array.array('d',etas)		
		self.pt_bins = array.array('d',[1000.*num for num in pts])

		self.results['eta_binning'] = ROOT.TH1F('eta_binning','eta_binning',25,0.,2.5)
		self.results['eta_binning'].GetXaxis().Set(len(self.eta_bins)-1,self.eta_bins)

		self.results['pt_binning'] = ROOT.TH1F('pt_binning','pt_binning',100,0.,200000.)
		self.results['pt_binning'].GetXaxis().Set(len(self.pt_bins)-1,self.pt_bins)
		
		for cut_level in [
			'total',
			'trigger',
			'reco_id',
			]:
			for name in [
				'l1',
				'l2',
				]:
				self.results['{0}_{1}_eta'.format(cut_level,name)] = ROOT.TH1F('{0}_{1}_eta'.format(cut_level,name),'{0}_{1}_eta'.format(cut_level,name),25,0.,2.5)
				self.results['{0}_{1}_pt'.format(cut_level,name)] = ROOT.TH1F('{0}_{1}_pt'.format(cut_level,name),'{0}_{1}_pt'.format(cut_level,name),100,0.,200000.)

		for name_ in [
			'total_counts_eta_{0}_{1}',
			'trigger_counts_eta_{0}_{1}',
			'reco_id_counts_eta_{0}_{1}',
			]:
			for eta1,eta2 in product(range(1,len(self.eta_bins)),range(1,len(self.eta_bins))):
				name = name_.format(eta1,eta2)
				self.results[name] = ROOT.TH2F(name,name,100,0.,200000,100,0.,200000.)
				self.results[name].GetXaxis().Set(len(self.pt_bins)-1,self.pt_bins)
				self.results[name].GetYaxis().Set(len(self.pt_bins)-1,self.pt_bins)

		for name_ in [
			'total_counts_pt_{0}_{1}',
			'trigger_counts_pt_{0}_{1}',
			'reco_id_counts_pt_{0}_{1}',
			]:
			for pt1,pt2 in product(range(1,len(self.pt_bins)),range(1,len(self.pt_bins))):
				name = name_.format(pt1,pt2)
				self.results[name] = ROOT.TH2F(name,name,25,0.,2.5,25,0.,2.5)
				self.results[name].GetXaxis().Set(len(self.eta_bins)-1,self.eta_bins)
				self.results[name].GetYaxis().Set(len(self.eta_bins)-1,self.eta_bins)

		"""
		self.lepton_resolution_name = '{lepton}_{bin}_{dist}_distribution{reversed_}'

		for lepton in [
			'l1',
			'l2',
			]:
			for reversed_ in ['_reversed','']:
				for dist in ['p
				for pt_bin,eta_bin in product(range(1,len(self.pt_bins)),range(1,len(self.eta_bins))):
					eta_low = self.eta_binning_resolution.GetBinLowEdge(eta_bin)
					eta_high = eta_low + self.eta_binning_resolution.GetBinWidth(eta_bin)
					pt_low = self.pt_binning_resolution.GetBinLowEdge(pt_bin)
					pt_high = pt_low + self.pt_binning_resolution.GetBinWidth(pt_bin)
					name = self.lepton_resolution_name.format(
						lepton=lepton,
						pt=pt_bin,
						eta=eta_bin,
						dist = 'pt',
						reversed_=reversed_,
						)
					self.results[name] = ROOT.TH1F(name,name,100,pt_low*.5,pt_high*1.5)
					name = self.lepton_resolution_name.format(
						lepton=lepton,
						pt=pt_bin,
						eta=eta_bin,
						dist = 'eta',
						reversed_=reversed_,
						)
					self.results[name] = ROOT.TH1F(name,name,100,eta_low*.5,eta_high*1.5)

		"""
		for lepton in ['l1','l2']:
			for dist in ['pt','eta']:
				for reversed_ in [True,False]:
					name = '_'.join([lepton,dist,'resolution']) + ('_reversed' if reversed_ else '')
					self.results[name] = ROOT.TProfile2D(name,name,50,-2.5,2.5,100,0,200000.)
					self.results[name].GetYaxis().Set(len(self.pt_bins)-1,self.pt_bins)
				
		"""
		for name in [
			'pt1_resolution',
			'pt2_resolution',
			'pt1_resolution_reversed',
			'pt2_resolution_reversed',
			]:
			self.results[name] = ROOT.TProfile2D(name,name,25,-2.5,2.5,100,0,200000.) #pt_truth-pt_off/pt_off:pt_off,eta_off
			self.results[name].GetYaxis().Set(len(self.pt_bins)-1,self.pt_bins)

		for name in [
			'pt1_resolution_',
			'pt2_resolution_',
			'pt1_resolution_reversed_',
			'pt2_resolution_reversed_',
			]:
			self.results[name] = ROOT.TProfile2D(name,name,25,-2.5,2.5,100,0,200000.) #pt_truth-pt_off/pt_off:pt_off,eta_off
			self.results[name].GetYaxis().Set(len(self.pt_bins)-1,self.pt_bins)

		for name in [
			'eta1_resolution',
			'eta2_resolution',
			'eta1_resolution_reversed',
			'pt2_resolution_reversed',
			]:
			self.results[name] = ROOT.TProfile2D(name,name,25,-2.5,2.5,100,0,200000.) #pt_truth-pt_off/pt_off:pt_off,eta_off
			self.results[name].GetYaxis().Set(len(self.pt_bins)-1,self.pt_bins)

		for name in [
			'pt1_resolution_',
			'pt2_resolution_',
			'pt1_resolution_reversed_',
			'pt2_resolution_reversed_',
			]:
			self.results[name] = ROOT.TProfile2D(name,name,25,-2.5,2.5,100,0,200000.) #pt_truth-pt_off/pt_off:pt_off,eta_off
			self.results[name].GetYaxis().Set(len(self.pt_bins)-1,self.pt_bins)
		"""

	def __call__(self,event):

		if event.__break__: return

		i = self.results['eta_binning'].FindBin(abs(event.l1_eta))
		j = self.results['eta_binning'].FindBin(abs(event.l2_eta))

		if not all([
			0<i<len(self.eta_bins),
			0<j<len(self.eta_bins),
			]): return

		total_counts_eta = self.results['total_counts_eta_{0}_{1}'.format(i,j)]
		trigger_counts_eta = self.results['trigger_counts_eta_{0}_{1}'.format(i,j)]
		reco_id_counts_eta = self.results['reco_id_counts_eta_{0}_{1}'.format(i,j)]

		i = self.results['pt_binning'].FindBin(event.l1_pt)
		j = self.results['pt_binning'].FindBin(event.l2_pt)

		if not all([
			0<i<len(self.pt_bins),
			0<j<len(self.pt_bins),
			]): return

		total_counts_pt = self.results['total_counts_pt_{0}_{1}'.format(i,j)]
		trigger_counts_pt = self.results['trigger_counts_pt_{0}_{1}'.format(i,j)]
		reco_id_counts_pt = self.results['reco_id_counts_pt_{0}_{1}'.format(i,j)]

		total_counts_eta.Fill(event.l1_pt,event.l2_pt,event.__weight__)
		total_counts_pt.Fill(event.l1_eta,event.l2_eta,event.__weight__)

		self.results['total_l1_pt'].Fill(event.l1_pt,event.__weight__)
		self.results['total_l1_eta'].Fill(event.l1_eta,event.__weight__)
		self.results['total_l2_pt'].Fill(event.l2_pt,event.__weight__)
		self.results['total_l2_eta'].Fill(event.l2_eta,event.__weight__)

		if not event.triggered: return
		#event.__weight__*= event.trigger_scale_factor
		
		trigger_counts_eta.Fill(event.l1_pt,event.l2_pt,event.__weight__)
		trigger_counts_pt.Fill(event.l1_eta,event.l2_eta,event.__weight__)

		self.results['trigger_l1_pt'].Fill(event.l1_pt,event.__weight__)
		self.results['trigger_l1_eta'].Fill(event.l1_eta,event.__weight__)
		self.results['trigger_l2_pt'].Fill(event.l2_pt,event.__weight__)
		self.results['trigger_l2_eta'].Fill(event.l2_eta,event.__weight__)

		if not all([
			event.l1_offline_passed_preselection,
			event.l2_offline_passed_preselection,
			]): return
			
		#event.__weight__*= event.l1_offline_scale_factor
		#event.__weight__*= event.l2_offline_scale_factor
		
		reco_id_counts_eta.Fill(event.l1_pt,event.l2_pt,event.__weight__)
		reco_id_counts_pt.Fill(event.l1_eta,event.l2_eta,event.__weight__)

		self.results['reco_id_l1_pt'].Fill(event.l1_pt,event.__weight__)
		self.results['reco_id_l1_eta'].Fill(event.l1_eta,event.__weight__)
		self.results['reco_id_l2_pt'].Fill(event.l2_pt,event.__weight__)
		self.results['reco_id_l2_eta'].Fill(event.l2_eta,event.__weight__)

		for lepton in ['l1','l2']:
			for dist in ['pt','eta']:
				for reversed_ in [True,False]:
					name = '_'.join([lepton,dist,'resolution']) + ('_reversed' if reversed_ else '')

					official_lepton = lepton+'_offline_' if reversed_ else lepton+'_'
					official_pt = getattr(event,official_lepton+'pt') 
					official_eta = getattr(event,official_lepton+'eta')
				
					match_lepton = lepton+'_offline_' if not reversed_ else lepton+'_' 
					match_pt = getattr(event,match_lepton+'pt') 
					match_eta = getattr(event,match_lepton+'eta')
						
					residual = (match_pt-official_pt)/official_pt if dist == 'pt' else (match_eta-official_eta)
					if dist == 'pt' and residual > 0.3: continue
					#print official_eta,official_pt,match_eta,match_pt,residual
					self.results[name].Fill(
						official_eta,
						official_pt,
						residual,
						event.__weight__,
						)
					
		"""
		for lepton in [
			'l1',
			'l2',
			]:
			for reversed_ in [True,False]:
				official_lepton = lepton+'_offline_' if reversed_ else lepton+'_'
				official_pt = getattr(event,official_lepton+'pt') 
				official_eta = getattr(event,official_lepton+'eta')
				
				match_lepton = lepton+'_offline_' if not reversed_ else lepton+'_' 
				match_pt = getattr(event,match_lepton+'pt') 
				match_eta = getattr(event,match_lepton+'eta')

				pt_bin = self.pt_binning_resolution.FindBin(official_pt)
				eta_bin = self.eta_binning_resolution.FindBin(official_eta)
				
				pt_dist = self.results.get(self.lepton_resolution_name.format(
					lepton=lepton,
					pt=pt_bin,
					eta=eta_bin,
					dist = 'pt',
					reversed_='reversed_' if reversed_ else '',				
					))
				if pt_dist is not None: pt_dist.Fill(match_pt,event.__weight__)
				eta_dist = self.results.get(self.lepton_resolution_name.format(
					lepton=lepton,
					pt=pt_bin,
					eta=eta_bin,
					dist = 'eta',
					reversed_='reversed_' if reversed_ else '',				
					))
				if eta_dist is not None: eta_dist.Fill(match_eta,event.__weight__)
		"""
		"""

		if abs(event.l1_pt-event.l1_offline_pt)/event.l1_pt<.3:
			self.results['pt1_resolution'].Fill(
				event.l1_eta,
				event.l1_pt,
				abs(event.l1_pt-event.l1_offline_pt)/event.l1_pt,
				event.__weight__
				)
		if abs(event.l2_pt-event.l2_offline_pt)/event.l2_pt<.3:
			self.results['pt2_resolution'].Fill(
				event.l2_eta,
				event.l2_pt,
				abs(event.l2_pt-event.l2_offline_pt)/event.l2_pt,
				event.__weight__
				)

		if abs(event.l1_pt-event.l1_offline_pt)/event.l1_offline_pt<.3:
			self.results['pt1_resolution_reversed'].Fill(
				event.l1_offline_eta,
				event.l1_offline_pt,
				abs(event.l1_pt-event.l1_offline_pt)/event.l1_offline_pt,
				event.__weight__
				)
		if abs(event.l2_pt-event.l2_offline_pt)/event.l2_offline_pt<.3:
			self.results['pt2_resolution_reversed'].Fill(
				event.l2_offline_eta,
				event.l2_offline_pt,
				abs(event.l2_pt-event.l2_offline_pt)/event.l2_offline_pt,
				event.__weight__
				)

		#dlfk
		if abs(event.l1_pt-event.l1_offline_pt)/event.l1_pt<.3:
			self.results['pt1_resolution_'].Fill(
				event.l1_eta,
				event.l1_pt,
				(event.l1_pt-event.l1_offline_pt)/event.l1_pt,
				event.__weight__
				)
		if abs(event.l2_pt-event.l2_offline_pt)/event.l2_pt<.3:
			self.results['pt2_resolution_'].Fill(
				event.l2_eta,
				event.l2_pt,
				(event.l2_pt-event.l2_offline_pt)/event.l2_pt,
				event.__weight__
				)

		if abs(event.l1_pt-event.l1_offline_pt)/event.l1_offline_pt<.3:
			self.results['pt1_resolution_reversed_'].Fill(
				event.l1_offline_eta,
				event.l1_offline_pt,
				(event.l1_pt-event.l1_offline_pt)/event.l1_offline_pt,
				event.__weight__
				)
		if abs(event.l2_pt-event.l2_offline_pt)/event.l2_offline_pt<.3:
			self.results['pt2_resolution_reversed_'].Fill(
				event.l2_offline_eta,
				event.l2_offline_pt,
				(event.l2_pt-event.l2_offline_pt)/event.l2_offline_pt,
				event.__weight__
				)
		"""
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

		try: z = [p for p in event.truth.values() if p().pdgId==23 and p().status in [2,155,10902]][0]
		except IndexError:
			print 'Z could not be found!'
			event.__break__=True
			return

		if event.mode==0: #ee
			try: event.l1,event.l2 = [c() for c in z.children if abs(c().pdgId)==11]
			except ValueError:
				print [p().pdgId for p in z.children]
				print 'Electrons could not be found'
				event.__break__=True
				return

			if event.l1.pt<event.l2.pt:
				event.l1,event.l2 = event.l2,event.l1 #swap pt ordered

		elif event.mode==1: #mumu
			try: event.l1,event.l2 = [c() for c in z.children if abs(c().pdgId)==13]
			except ValueError:
				print [p().pdgId for p in z.children]
				print 'Muons could not be found'
				event.__break__=True
				return

			if event.l1.pt<event.l2.pt:
				event.l1,event.l2 = event.l2,event.l1 #swap pt ordered

		elif event.mode==2: #tautau->emu
			try: tau1,tau2 = [p for p in z.children if abs(p().pdgId)==15]
			except ValueError:
				print [p().pdgId for p in z.children]
				print 'Taus could not be found'
				event.__break__=True
				return
			tau1 = ([p for p in tau1.children if abs(p().pdgId)==15]+[tau1])[0]
			tau2 = ([p for p in tau2.children if abs(p().pdgId)==15]+[tau2])[0]
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
			('scale_factor','float',1.),
			('scale_factor_error','float',0.),
			]

		for l,(name,type_,_) in product(['l1','l2'],self.names):
			self.create_branches[l+'_offline_'+name] = type_

		self.dummy = particle(**dict((name,default) for name,_,default in self.names))

	def __call__(self,event):
		if event.mode==0:
			#find matching to leading electron
			try: event.l1.offline_match = sorted([electron for electron in event.electrons.values() if electron().DeltaR(event.l1())<self.dR_max and electron.passed_preselection], key=lambda el: el().DeltaR(event.l1()))[0]
			except IndexError: event.l1.offline_match = self.dummy
			#find matching to subleading muon
			try: event.l2.offline_match = sorted([electron for electron in event.electrons.values() if electron().DeltaR(event.l2())<self.dR_max and electron != event.l1.offline_match and electron.passed_preselection], key=lambda el: el().DeltaR(event.l2()))[0]
			except IndexError: event.l2.offline_match = self.dummy
		elif event.mode==1:
			#find matching to leading muon
			try: event.l1.offline_match = sorted([muon for muon in event.muons.values() if muon().DeltaR(event.l1())<self.dR_max and muon.passed_preselection and muon.pt_corrected>20000.], key=lambda mu: mu().DeltaR(event.l1()))[0]
			except IndexError: event.l1.offline_match = self.dummy
			#find matching to subleading muon
			try: event.l2.offline_match = sorted([muon for muon in event.muons.values() if muon().DeltaR(event.l2())<self.dR_max and muon != event.l1.offline_match and muon.passed_preselection and muon], key=lambda mu: mu().DeltaR(event.l2()))[0]
			except IndexError: event.l2.offline_match = self.dummy
		elif event.mode==2:
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
			event.l1.offline_match.etcone20 = event.l1.offline_match.etcone20_corrected

		if event.l2.offline_match!=self.dummy: 
			event.l2.offline_match.pt = event.l2.offline_match.pt_corrected
			event.l2.offline_match.E = event.l2.offline_match.E_corrected
			event.l2.offline_match.etcone20 = event.l2.offline_match.etcone20_corrected

		for l,(name,type_,_) in product(['l1','l2'],self.names):
			event.__dict__[l+'_offline_'+name] = event.__dict__[l].offline_match.__dict__[name]

class trigger(event_function):
	
	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'EF_mu18_tight_mu8_EFFS',
			'EF_e12Tvh_medium1_mu8',
			'EF_e24vhi_medium1',
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
			event.triggered = event.EF_e24vhi_medium1
		if event.mode == 1:
			event.triggered = event.EF_mu18_tight_mu8_EFFS
		if event.mode == 2:
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
		#ee
		if event.mode == 0:
			result = self.electron_trigger_e24vhi_medium.calculate(1,run,event.l1.offline_match.cl_eta,event.l1.offline_match.pt)
			event.trigger_scale_factor = result.getScaleFactor()
			event.trigger_scale_factor_error = result.getTotalUncertainty()
		#mumu
		if event.mode == 1:
			muon1 = event.l1.offline_match()
			muon2 = event.l2.offline_match()

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
		if event.mode == 2:
			#get electron scale factor
			electron = event.l1.offline_match()
         		electron_scale_factor = self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",0)
			electron_scale_factor_error = max([
				abs(electron_scale_factor-self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",1)),
				abs(electron_scale_factor-self.electron_trigger_e12Tvh_medium1.getSFElec(electron,run,"e12Tvhm1",-1))
				])
			#get muon scale factor
			muon = event.l2.offline_match()
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

		self.create_branches['top_hfor_type'] = 'int'

		self.initialize_tools()

	def __call__(self,event):

		event.top_hfor_type = getattr(event,'top_hfor_type',-1)

		#2 preselection leptons, no hadronic taus, at least one preselection jet
		if not all([
			sum(1 for tau in event.taus.values() if tau.passed_preselection and not tau.overlap_removed)==0,
			#sum(1 for jet in event.jets.values() if jet.passed_preselection and not jet.overlap_removed)>0,
			sum(1 for jet in event.jets.values() if jet.passed_preselection and not jet.overlap_removed and jet.isBadLooseMinus)==0,
			event.larError!=2,
			event.tileError!=2,
			(event.coreFlags&0x40000)==0,
			self.tile_trip_reader.checkEvent(event.random_RunNumber,event.lbn,event.EventNumber),
			]):
			event.__break__=True
			return

	def initialize_tools(self):
		load('TileTripReader')
		self.tile_trip_reader = ROOT.Root.TTileTripReader()
		#comment
