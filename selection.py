from common.analysis import analysis
from common.functions import event_function,result_function
from common.external import load
from common.particle import particle

import ROOT
import os
from math import sqrt
import json
import array
import random

from operator import itemgetter, attrgetter, mul

class mutate_make_selection_Z_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_tautau(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			select_Z_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class mutate_make_selection_signal(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_tautau(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			select_signal_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class make_selection_Z_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			remove_overlapped_jets(),
			compute_kinematics(),
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
			remove_overlapped_jets(),
			compute_kinematics(),
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
			remove_overlapped_jets(),
			compute_kinematics(),
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
			remove_overlapped_jets(),
			compute_kinematics(),
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
			remove_overlapped_jets(),
			compute_kinematics(),
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
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			Z_scale(),
			select_signal_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

from tauola import tauola_

class mutate_mumu_to_tautau(event_function):
	def __init__(self):
		event_function.__init__(self)
		self.required_branches += [
			]
		self.electron_mass = 0.5/1000.
		self.muon_mass = 100.
		self.tau_mass = 1776.82

		self.electron_decay = array.array('d',[self.electron_mass,0.,0.])
		self.muon_decay = array.array('d',[self.muon_mass,0.,0.])
		self.tauola = tauola_()

	def __call__(self,event):
		if not event.lepton_class==1:
			event.__break__=True
			return

		if random.getrandbits(1): event.l1,event.l2 = event.l2,event.l1 #flip e<->mu decay
	
		tauola_call = []

		mother = event.l1()+event.l2()
		boost = mother.BoostVector()
		for muon in [event.l1(),event.l2()]:
			muon.Boost(-boost)
			try: scale = sqrt(muon.E()**2.-self.tau_mass**2.)/muon.P()
			except ValueError:
				event.__break__ == True
				return
			muon.SetPxPyPzE(muon.Px()*scale,muon.Py()*scale,muon.Pz()*scale,muon.E())
			muon.Boost(boost)
			
			tauola_call+=[muon.Px()/1000.,muon.Py()/1000.,muon.Pz()/1000.] #GEV for tauola

		result = self.tauola.leptonic_decay(*tauola_call)
		event.l1.set_px_py_pz_e(*[energy*1000. for energy in result[:4]])
		event.l2.set_px_py_pz_e(*[energy*1000. for energy in result[4:]])
				
		event.l1.pt = event.l1().Pt()
		event.l1.eta = event.l1().Eta()
		event.l1.phi = event.l1().Phi()
		event.l1.E = event.l1().E()
		event.l2.pt = event.l2().Pt()
		event.l2.eta = event.l2().Eta()
		event.l2.phi = event.l2().Phi()
		event.l2.E = event.l2().E()

		additional_missing_energy = mother-event.l1()-event.l2()
		event.miss.set_particle(event.miss()+additional_missing_energy)

		"""
		for muon,decay in zip([event.l1,event.l2],[decay1,decay2]):
			phase_space = ROOT.TGenPhaseSpace()
			tau = ROOT.TLorentzVector()
			pt = sqrt(muon.pt**2.-self.tau_mass**2+self.muon_mass**2.)
			tau.SetPtEtaPhiM(pt,muon.eta,muon.phi,self.tau_mass)
			if not phase_space.SetDecay(tau,3,decay):
				event.__break__ = True
				return
			while phase_space.Generate() < random.random(): pass
			#weight = phase_space.Generate()
			#event.__weight__*=weight

			additional_missing_energy = phase_space.GetDecay(1)+phase_space.GetDecay(2)
			event.miss.set_particle(event.miss()+additional_missing_energy)

			muon.set_particle(phase_space.GetDecay(0))
			muon.etcone20 = 0.
			muon.ptcone40 = 0.
			muon.pt = muon().Pt()
			muon.eta = muon().Eta()
			muon.phi = muon().Phi()
			muon.E = muon().E()
		"""

		

		#if flip: event.l1,event.l2 = event.l2,event.l1

		if not all([
			event.l1.pt>15000. and any([abs(event.l1.eta)<1.37 or 1.52<abs(event.l1.eta)<2.5]), #electron selection
			event.l2.pt>10000. and abs(event.l2.eta)<2.5, #muon selection
			]):
			event.__break__ = True
			return

		event.__weight__ *= 0.06197796 #tau branching ratio to emu
		event.lepton_class = 2 #now this is emu event

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
		if event.mass_range == 0: return

		profile = getattr(self.Z_scale,'lepton_pair_pT_1_{0}_scale'.format(event.lepton_class))
		if event.lepton_pair_pT > profile.GetBinLowEdge(profile.GetNbinsX()+1): weight_bin = profile.GetNbinsX()
		elif event.lepton_pair_pT < profile.GetBinLowEdge(1): weight_bin = 1
		else: weight_bin = profile.FindBin(event.lepton_pair_pT)
		weight = profile.GetBinContent(weight_bin)
		event.__weight__*=weight

		profile = getattr(self.Z_scale,'jet_energy_1_{0}_scale'.format(event.lepton_class))
		if event.jet_energy > profile.GetBinLowEdge(profile.GetNbinsX()+1): weight_bin = profile.GetNbinsX()
		elif event.jet_energy < profile.GetBinLowEdge(1): weight_bin = 1
		else: weight_bin = profile.FindBin(event.jet_energy)
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
	
		event.__weight__*=event.bveto_scale_factor #apply bveto scale factor

		if not all([
			#event.jet_energy < 100000.,
			#event.missing_energy < 50000.,
                        #event.lepton_dR<2.0,
			#not (event.lepton_pair_pT<10000. and event.lepton_class in [0,1]),
			#not (event.lepton_pair_pT<5000. and event.lepton_class ==2 ),
			event.l1.etcone20/event.l1.pt<0.09,
			event.l1.ptcone40/event.l1.pt<0.17,
			event.l2.etcone20/event.l2.pt<0.09,
			event.l2.ptcone40/event.l2.pt<0.17,
			#any([
			#	event.lepton_class in [0,1] and 40000.<event.lepton_pair_mass<100000.,
			#	event.lepton_class == 2 and 35000.<event.lepton_pair_mass<65000.,
			#	]),
			len(event.bjets_preselected)>=1,
			len(event.bjets)==0,
			]):
			event.__break__=True
			return

class select_signal_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		#apply both because simultaneous requirement of btag and all others bvetoed
		event.__weight__*=event.btag_scale_factor #apply btag scale factor
		event.__weight__*=event.bveto_scale_factor #apply bveto scale factor

		if not all([
			not (event.lepton_pair_pT<10000. and event.lepton_class in [0,1]),
			not (event.lepton_pair_pT<5000. and event.lepton_class ==2 ),
                        event.lepton_dR<2.5,
			#event.jet_energy < 120000.,
			event.l1.etcone20/event.l1.pt<0.09,
			event.l1.ptcone40/event.l1.pt<0.17,
			event.l2.etcone20/event.l2.pt<0.09,
			event.l2.ptcone40/event.l2.pt<0.17,
			#10000.<event.missing_energy<60000.,
			event.transverse_W_mass<40000.,
			len(event.bjets)==1,
			]):
			event.__break__=True
			return

class select_tt_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		event.__weight__*=event.btag_scale_factor #apply btag scale factor

		if not all([
			not (event.lepton_pair_pT<10000. and event.lepton_class in [0,1]),
			not (event.lepton_pair_pT<5000. and event.lepton_class ==2 ),
			event.jet_energy > 100000.,
			event.missing_energy > 50000.,
			event.l1.etcone20/event.l1.pt<0.09,
			event.l1.ptcone40/event.l1.pt<0.17,
			event.l2.etcone20/event.l2.pt<0.09,
			event.l2.ptcone40/event.l2.pt<0.17,
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
			#'passed_b_preselection',
			]
		self.create_branches['top_hfor_type'] = None

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

		if event.l1().DeltaR(event.l2())<0.4:
			event.l2.ptcone40-=event.l1.pt
			event.l1.ptcone40-=event.l2.pt

		if getattr(event,'top_hfor_type',0)==4:
			event.__break__ = True
			return

class remove_overlapped_jets(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		#remove jets from electrons, muons
		for jetN,jet in event.jets.items():
			if jet().DeltaR(event.l2())<0.2:
				del event.jets[jetN]
				if jetN in event.bjets_preselected: del event.bjets_preselected[jetN]
				if jetN in event.bjets: del event.bjets[jetN]
				continue
			if jet().DeltaR(event.l1())<0.2:
				del event.jets[jetN]
				if jetN in event.bjets_preselected: del event.bjets_preselected[jetN]
				if jetN in event.bjets: del event.bjets[jetN]

def collinear_mass(l1,l2,miss):
	m_frac_1 = ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())) / ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())+(l2.Py()*miss.Px())-(l2.Px()*miss.Py()))
	m_frac_2 = ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())) / ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())+(l1.Px()*miss.Py())-(l1.Py()*miss.Px()))

	if m_frac_1*m_frac_2 > 0.: return (l1+l2).M()/sqrt(m_frac_1*m_frac_2)
	return -1.

class compute_kinematics(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		sorted_jets = sorted(event.jets.values(),key=attrgetter('pt'), reverse=True) #jets sorted highest pt first
		lepton_pair = event.l1()+event.l2()

		event.missing_energy = event.miss().Pt()
		event.lepton_pair_miss_dPhi = abs(event.miss().DeltaPhi(lepton_pair))
		event.lepton_pair_pT = lepton_pair.Pt()
		event.lepton_pair_pT_diff = abs(event.l1.pt-event.l2.pt)
		event.lepton_pair_mass = lepton_pair.M()
		event.lepton_pair_mass_low = event.lepton_pair_mass

		
		try:
			event.transverse_W_mass = max([
				sqrt(2*(event.miss().Et()*event.l1().Et()-event.l1().Px()*event.miss().Px()-event.l1().Py()*event.miss().Py())),
				sqrt(2*(event.miss().Et()*event.l2().Et()-event.l2().Px()*event.miss().Px()-event.l2().Py()*event.miss().Py())),
				#(event.l1()+event.miss()).Mt(),
				#(event.l2()+event.miss()).Mt(),
				])
		except: event.transverse_W_mass = -1.

		if event.lepton_class==0:
			event.off_threshold = min([event.l1.pt-25000.,event.l2.pt-15000.])
			event.collinear_mass = -1.
		elif event.lepton_class==1:
			event.off_threshold = min([event.l1.pt-25000.,event.l2.pt-10000.])
			event.collinear_mass = -1.
		else:
			event.off_threshold = min([event.l1.pt-15000.,event.l2.pt-10000.])
			event.collinear_mass = collinear_mass(event.l1(),event.l2(),event.miss())
		if not event.lepton_pair_mass<150000.:
			event.__break__ = True
			return

		event.mass_range = 0 if event.lepton_pair_mass<5000. else 1

		for lepton,name in zip([event.l1,event.l2],['l1','l2']):
			for attr in ['pt','eta','phi']:
				setattr(event,name+'_'+attr,getattr(lepton,attr))
				

		event.lepton_dR = event.l1().DeltaR(event.l2())
		event.same_sign = (event.l1.charge*event.l2.charge)>0.
		try: event.jet_energy = sum(jet.pt for jet in event.jets.values())
		except ValueError: event.jet_energy = 0.
		try: event.bjet_energy = sum(jet.pt for jet in event.bjets.values())
		except ValueError: event.bjet_energy = 0.
		if len(sorted_jets)>=1: event.leading_jet_miss_dPhi = abs(event.miss().DeltaPhi(sorted_jets[0]()))
		else: event.leading_jet_miss_dPhi = -1.
		if len(sorted_jets)>=2: event.subleading_jet_miss_dPhi = abs(event.miss().DeltaPhi(sorted_jets[1]()))
		else: event.subleading_jet_miss_dPhi = -1.
		event.jet_n = len(event.jets)
		event.bjet_n = len(event.bjets)

		event.btag_scale_factor = reduce(mul,[jet.bJet_scale_factor for jet in event.bjets.values()],1)
		event.bveto_scale_factor = reduce(mul,[jet.bJet_scale_factor for jet in event.bjets_preselected.values() if not jet.flavor_weight_MV1>0.7892],1)

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low)) for name,binning,high,low in [
			('off_threshold',40,0.,40000.),
			('transverse_W_mass',50,0.,200000.),
			('missing_energy',100,0.,100000.),
			('collinear_mass',50,0.,150000.),
			('lepton_pair_mass',50,5000.,150000.),
			('lepton_pair_mass_low',180,0.,45000.),
			('lepton_dR',100,0.,10.),
			('jet_energy',100,0.,200000.),
			('bjet_energy',100,0.,200000.),
			('leading_jet_miss_dPhi',32,-3.2,3.2),
			('subleading_jet_miss_dPhi',32,-3.2,3.2),
			('lepton_pair_miss_dPhi',32,-3.2,3.2),
			('lepton_pair_pT',100,0.,100000.),
			('lepton_pair_pT_diff',100,0.,100000.),
			('l1_pt',100,0.,100000.),
			('l1_eta',24,-3.,3.),
			('l1_phi',32,-3.2,3.2),
			('l2_pt',100,0.,100000.),
			('l2_eta',24,-3.,3.),
			('l2_phi',32,-3.2,3.2),
			('jet_n',10,0,10),
			('bjet_n',10,0,10),
			])

		self.names_2d = [
			('lepton_pair_mass','collinear_mass'),
			('lepton_pair_mass','off_threshold'),
			('lepton_pair_mass','missing_energy'),
			('collinear_mass','off_threshold'),
			('missing_energy','collinear_mass'),
			('collinear_mass','transverse_W_mass'),
			]

		for name_,(binning,high,low) in self.names.items():
			for mass_range in [0,1]:
				for lepton_class in [0,1,2]:
					name = '{0}_{1}_{2}'.format(name_,mass_range,lepton_class)
					self.results[name] = ROOT.TH1F(name,name,binning,high,low)
					self.results[name].Sumw2()

		for name1,name2 in self.names_2d:
			binning1,high1,low1 = self.names[name1]
			binning2,high2,low2 = self.names[name2]
			for mass_range in [0,1]:
				for lepton_class in [0,1,2]:
					name = '{0}_{1}_{2}_{3}'.format(name1,name2,mass_range,lepton_class)
					self.results[name] = ROOT.TH2F(name,name,binning1,high1,low1,binning2,high2,low2)
					self.results[name].Sumw2()

	def __call__(self,event):
		if event.__break__: return

		weight = event.__weight__
		weight*= -1 if event.same_sign else 1.

		for name_ in self.names:
			name = '{0}_{1}_{2}'.format(name_,event.mass_range,event.lepton_class)
			self.results[name].Fill(event.__dict__[name_],weight)

		for name1,name2 in self.names_2d:
			name = '{0}_{1}_{2}_{3}'.format(name1,name2,event.mass_range,event.lepton_class)
			self.results[name].Fill(event.__dict__[name1],event.__dict__[name2],weight)
