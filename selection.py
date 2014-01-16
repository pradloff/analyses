from common.analysis import analysis
from common.functions import event_function,result_function
from common.external import load
from common.particle import particle

import ROOT
import os
from math import sqrt
import json

class make_selection_Z_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			get_weight(),
			select_Z_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class make_selection_Z_scaled_Z_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			get_weight(),
			Z_scale(),
			select_Z_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class make_selection_Z_scaled_tt_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			get_weight(),
			Z_scale(),
			select_tt_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class make_selection_tt_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			get_weight(),
			select_tt_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)


class make_selection_signal(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			get_weight(),
			select_signal_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class make_selection_Z_scaled_signal(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			get_weight(),
			Z_scale(),
			select_signal_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class get_weight(event_function):
	def __init__(self):
		event_function.__init__(self)
		self.required_branches += [
			'l1_scale_factor',
			'l1_scale_factor_error',
			'l2_scale_factor',
			'l2_scale_factor_error',
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
			event.l1_scale_factor,
			event.l2_scale_factor,
			event.trigger_scale_factor,
			event.weight_pileup,
			]: event.__weight__*=weight
		for jet in event.jets.values(): event.__weight__*=jet.bJet_scale_factor

	def initialize(self):
		analysis_home = os.getenv('ANALYSISHOME')
		mc_lumi_file = '{0}/data/mc_lumi.json'.format(analysis_home)
		with open(mc_lumi_file) as f: self.mc_lumi_info = json.loads(f.read())

class Z_scale(event_function):
	def __init__(self):
		event_function.__init__(self)
		self.initialize()

	def __call__(self,event):
		"""
		profile = getattr(self.Z_scale,'l2_pt_{0}_scale'.format(event.lepton_class))
		if event.l2_pt > profile.GetBinLowEdge(profile.GetNbinsX()+1): weight_bin = profile.GetNbinsX()
		else: weight_bin = profile.FindBin(event.l2_pt)
		weight = profile.GetBinContent(weight_bin)
		event.__weight__*=weight
		"""
		"""
		profile = getattr(self.Z_scale,'missing_energy_{0}_scale'.format(event.lepton_class))
		if event.missing_energy > profile.GetBinLowEdge(profile.GetNbinsX()+1): weight_bin = profile.GetNbinsX()
		elif event.missing_energy < profile.GetBinLowEdge(1): weight_bin = 1
		else: weight_bin = profile.FindBin(event.missing_energy)
		weight = profile.GetBinContent(weight_bin)
		event.__weight__*=weight
		"""
		profile = getattr(self.Z_scale,'jet_energy_{0}_scale'.format(event.lepton_class))
		if event.jet_energy > profile.GetBinLowEdge(profile.GetNbinsX()+1): weight_bin = profile.GetNbinsX()
		elif event.jet_energy < profile.GetBinLowEdge(1): weight_bin = 1
		else: weight_bin = profile.FindBin(event.missing_energy)
		weight = profile.GetBinContent(weight_bin)
		event.__weight__*=weight


	def initialize(self):
		analysis_home = os.getenv('ANALYSISHOME')
		Z_scale_file = '{0}/data/Z_scale.root'.format(analysis_home)
		self.Z_scale = ROOT.TFile(Z_scale_file)

class select_Z_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		if not all([
			#event.jet_energy < 100000.,
			#event.missing_energy < 50000.,
			event.l1.etcone20/event.l1.pt<0.09,
			event.l1.ptcone40/event.l1.pt<0.17,
			event.l2.etcone20/event.l2.pt<0.09,
			event.l2.ptcone40/event.l2.pt<0.17,
			any([
				event.lepton_class in [0,1] and 40000.<event.lepton_pair_mass<100000.,
				event.lepton_class == 2 and 35000.<event.lepton_pair_mass<65000.,
				]),
			len(event.bjets_preselected)>=1,
			len(event.bjets)==0,
			]):
			event.__break__=True
			return


class select_signal_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		if not all([
			event.jet_energy < 100000.,
			event.l1.etcone20/event.l1.pt<0.09,
			event.l1.ptcone40/event.l1.pt<0.17,
			event.l2.etcone20/event.l2.pt<0.09,
			event.l2.ptcone40/event.l2.pt<0.17,
			any([
				event.lepton_class in [0,1],
				event.lepton_class == 2 and 30000.<event.missing_energy<80000.,
				]),
			len(event.jets)>=1,
			len(event.bjets)==1,
			]):
			event.__break__=True
			return

class select_tt_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		if not all([
			event.jet_energy > 100000.,
			event.missing_energy > 50000.,
			event.l1.etcone20/event.l1.pt<0.09,
			event.l1.ptcone40/event.l1.pt<0.17,
			event.l2.etcone20/event.l2.pt<0.09,
			event.l2.ptcone40/event.l2.pt<0.17,
			len(event.jets)>=1,
			len(event.bjets)>=1,
			]):
			event.__break__=True
			return

class build_events(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.lepton_names = [
			'E',
			'charge',
			'eta',
			'etcone20',
			'phi',
			'pt',
			'ptcone40',
			]
		self.required_branches += ['l1_'+name for name in self.lepton_names]
		self.required_branches += ['l2_'+name for name in self.lepton_names]

		self.required_branches += ['lepton_class']

		self.jet_names = [
			'E',
			'eta',
			'phi',
			'pt',
			'flavor_weight_MV1',
			'jvf',
			'bJet_scale_factor',
			'passed_b_preselection',
			]
		self.required_branches += ['jet_'+name for name in self.jet_names]
		self.required_branches += ['jet_n']

		self.required_branches += [
			'phi_miss',
			'pt_miss',
			'px_miss',
			'py_miss',
			]

		self.create_branches['miss'] = None
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
				not ((abs(event.jet_eta[jet])<2.4 and event.jet_pt[jet]<50000.) and not ((event.jet_jvf[jet])>0.5)),
				abs(event.jet_eta[jet])<2.4,
				event.passed_b_preselection[jet],
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
				not ((abs(event.jet_eta[jet])<2.4 and event.jet_pt[jet]<50000.) and not ((event.jet_jvf[jet])>0.5)),
				abs(event.jet_eta[jet])<2.4,
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

		#create missing energy particle
		event.miss = particle()
		event.miss.set_px_py_pz_e(
			event.px_miss,
			event.py_miss,
			0.,
			event.pt_miss
			)

		#create lepton 1
		event.l1 = particle(\
			**dict((name,event.__dict__['l1_'+name]) for name in self.lepton_names)
			)
		event.l1.set_pt_eta_phi_e(
			event.l1.pt,
			event.l1.eta,
			event.l1.phi,
			event.l1.E,
			)

		#create lepton 2
		event.l2 = particle(\
			**dict((name,event.__dict__['l2_'+name]) for name in self.lepton_names)
			)
		event.l2.set_pt_eta_phi_e(
			event.l2.pt,
			event.l2.eta,
			event.l2.phi,
			event.l2.E,
			)

		#remove jets from electrons, muons
		for jetN,jet in event.jets.items():
			if jet().DeltaR(event.l2())<0.2:
				del event.jets[jetN]
				if jetN in event.bjets: del event.bjets[jetN]
				continue
			if jet().DeltaR(event.l1())<0.2:
				del event.jets[jetN]
				if jetN in event.bjets: del event.bjets[jetN]

		event.missing_energy = event.miss().Pt()
		event.lepton_pair_mass = (event.l1()+event.l2()).M()
		event.lepton_dR = event.l1().DeltaR(event.l2())
		event.same_sign = (event.l1.charge*event.l2.charge)>0.
		try: event.jet_energy = sum(jet.pt for jet in event.jets.values())
		except ValueError: event.jet_energy = 0.
		try: event.bjet_energy = sum(jet.pt for jet in event.bjets.values())
		except ValueError: event.bjet_energy = 0.

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low)) for name,binning,high,low in [
			('missing_energy',100,0.,100000.),
			('lepton_pair_mass',100,0.,150000.),
			('lepton_dR',100,0.,10.),
			('jet_energy',100,0.,200000.),
			('bjet_energy',100,0.,200000.),
			('l1_pt',100,0.,100000.),
			('l1_eta',24,-3.,3.),
			('l1_phi',32,-3.2,3.2),
			('l2_pt',100,0.,100000.),
			('l2_eta',24,-3.,3.),
			('l2_phi',32,-3.2,3.2),
			])

		for name_,(binning,high,low) in self.names.items():
			for lepton_class in [0,1,2]:
				name = name_+'_'+str(lepton_class)
				self.results[name] = ROOT.TH1F(name,name,binning,high,low)
				self.results[name].Sumw2()

	def __call__(self,event):
		if event.__break__: return

		weight = event.__weight__
		weight*= -1 if event.same_sign else 1.

		for name_ in self.names:
			name = name_+'_'+str(event.lepton_class)
			self.results[name].Fill(event.__dict__[name_],weight)

