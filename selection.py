from common.analysis import analysis
from common.functions import event_function,result_function
from common.external import load
from common.particle import particle

import ROOT
import os
from math import sqrt,pi,cos
import json
import array
import random

from operator import itemgetter, attrgetter, mul

class make_selection_preselection(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			preselection_events(), #just one jet
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class mutate_make_selection_preselection(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_tautau(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			preselection_events(), #just one jet
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

class make_selection_tt_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(b=True),
			select_tt_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)


class mutate_make_selection_tt_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_tautau(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(b=True),
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
			get_weight(b=True),
			select_signal_events()
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
			get_weight(b=True),
			select_signal_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

#----

from tauola import tauola_

class mutate_mumu_to_tautau(event_function):
	def __init__(self,min_mass=0.,max_mass=1000000000.):
		self.min_mass = min_mass
		self.max_mass = max_mass
		event_function.__init__(self)
		self.required_branches += [
			'random_RunNumber',
			]

		self.periods_runnumbers = {
			"A_":(200804,201556),
			"B_":(202660,205113),
			"CtoE_":(206248,210308),
			"G_":(211522,212272),
			"HtoL_":(212619,215643),
			}

		self.electron_mass = 0.5/1000.
		self.muon_mass = 100.
		self.tau_mass = 1776.82

		self.electron_decay = array.array('d',[self.electron_mass,0.,0.])
		self.muon_decay = array.array('d',[self.muon_mass,0.,0.])
		self.tauola = tauola_()

		self.initialize_tools()

	def __call__(self,event):
		if not event.lepton_class==1:
			event.__break__=True
			return

		#Remove offline muon effect on event
		event.miss.set_particle(event.miss()+(event.l1()+event.l2()))
		#event.miss().RotateZ(random.random()*2.*pi)
		event.sum_Et_miss-=event.l1.pt
		event.sum_Et_miss-=event.l2.pt
		if event.sum_Et_miss<0.: event.sum_Et_miss = 0.

		#get reverse smeared muons (now we have smeared truth muons)
		binx = self.mumu.pt1_resolution_reversed.GetXaxis().FindBin(event.l1.eta)
		biny = self.mumu.pt1_resolution_reversed.GetYaxis().FindBin(event.l1.pt)
		sigma = self.mumu.pt1_resolution_reversed.GetBinContent(binx,biny)
		smear1 = random.gauss(0.,sigma)
		event.l1.set_particle(event.l1()*(1+smear1))

		binx = self.mumu.pt2_resolution_reversed.GetXaxis().FindBin(event.l2.eta)
		biny = self.mumu.pt2_resolution_reversed.GetYaxis().FindBin(event.l2.pt)
		sigma = self.mumu.pt2_resolution_reversed.GetBinContent(binx,biny)
		smear2 = random.gauss(0.,sigma)
		event.l2.set_particle(event.l2()*(1+smear2))

		#get eta bins
		eta1 = self.mumu.eta_binning.FindBin(event.l1.eta)
		eta2 = self.mumu.eta_binning.FindBin(event.l2.eta)
		
		#get muon inefficiency
		binx = self.mumu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetXaxis().FindBin(event.l1().Pt())
		biny = self.mumu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetYaxis().FindBin(event.l2().Pt())
		total = self.mumu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetBinContent(binx,biny)
		selected = self.mumu.Get('reco_id_counts_{0}_{1}'.format(eta1,eta2)).GetBinContent(binx,biny)

		uncovered1 = False
		if total != 0.: inefficiency = float(selected)/total
		else: 
			uncovered1 = True
			inefficiency = 1.
		if inefficiency<0.1: inefficiency=0.1

		"""
		if uncovered1:
			print 'Uncovered: {0}'.format([round(num,2) for num in [
				eta1,
				eta2,
				binx,
				biny,
				event.l1_original.eta,
				event.l1_original.pt,
				event.l2_original.eta,
				event.l2_original.pt,
				event.l1().Eta(),
				event.l1().Pt(),
				event.l2().Eta(),
				event.l2().Pt(),
				]])
			event.__break__ = True
			return
		"""
		"""
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

		#divide by muon efficiency to recover full phase-space

		self.config_muon_trigger_mu18_tight_mu8_EFFS.runNumber = run
		self.config_muon_trigger_mu18_tight_mu8_EFFS.period = period

		muons = ROOT.std.vector('TLorentzVector')()
		muons.push_back(event.l1())
		muons.push_back(event.l2())

		muons_quality = ROOT.std.vector('int')()
		muons_quality.push_back(1)
		muons_quality.push_back(1)

		event.__weight__ /= self.muon_trigger_mu18_tight_mu8_EFFS.getDimuonEfficiency(
			self.config_muon_trigger_mu18_tight_mu8_EFFS,
			muons,
			muons_quality,
			self.config_muon_trigger_mu18_tight_mu8_EFFS.trigger
			).first
		"""
		if random.getrandbits(1): event.l1,event.l2 = event.l2,event.l1 #flip e<->mu decay
	
		tauola_call = []

		mother = event.l1()+event.l2()

		if not (self.min_mass<mother.M()<self.max_mass):
			event.__break__ == True
			return

		boost = mother.BoostVector()
		for muon in [event.l1(),event.l2()]:
			muon.Boost(-boost)
			try: scale = sqrt(muon.E()**2.-self.tau_mass**2.)/muon.P()
			except ValueError:
				event.__break__ = True
				return
			muon.SetPxPyPzE(muon.Px()*scale,muon.Py()*scale,muon.Pz()*scale,muon.E())
			muon.Boost(boost)
			
			tauola_call+=[muon.Px()/1000.,muon.Py()/1000.,muon.Pz()/1000.] #GEV for tauola
		tauola_call.append(23) #Z emulation

		#We now have truth electron and muon
		result = self.tauola.leptonic_decay(*tauola_call)
		event.l1.set_px_py_pz_e(*[energy*1000. for energy in result[:4]])
		event.l2.set_px_py_pz_e(*[energy*1000. for energy in result[4:]])
				
		#event.l1.etcone20 = 0.
		#event.l1.ptcone40 = 0.
		#event.l2.etcone20 = 0.
		#event.l2.ptcone40 = 0.

		#additional_missing_energy = mother-event.l1()-event.l2()
		#event.miss.set_particle(event.miss()+additional_missing_energy)


		#event.miss.set_particle(event.miss()-event.l1())

		#get eta bins
		eta1 = self.emu.eta_binning.FindBin(event.l1.eta)
		eta2 = self.emu.eta_binning.FindBin(event.l2.eta)
		
		#get electron/muon inefficiency
		binx = self.emu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetXaxis().FindBin(event.l1().Pt())
		biny = self.emu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetYaxis().FindBin(event.l2().Pt())
		total = self.emu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetBinContent(binx,biny)
		selected = self.emu.Get('reco_id_counts_{0}_{1}'.format(eta1,eta2)).GetBinContent(binx,biny)
		#total = self.emu.Get('total_counts_{0}_{1}'.format(eta1,eta2)).GetBinContent(event.l1().Pt(),event.l2().Pt())
		#selected = self.emu.Get('reco_id_counts_{0}_{1}'.format(eta1,eta2)).GetBinContent(event.l1().Pt(),event.l2().Pt())

		uncovered2 = False
		if total != 0.: efficiency = float(selected)/total
		else: 
			uncovered2 = True
			efficiency = 1.
			"""
			print 'Uncovered: {0}'.format([round(num,2) for num in [
				eta1,
				eta2,
				binx,
				biny,
				event.l1().Eta(),
				event.l1().Pt(),
				event.l2().Eta(),
				event.l2().Pt(),
				]])
			"""
		#get smeared electrons/muons
		binx = self.emu.pt1_resolution.GetXaxis().FindBin(event.l1().Eta())
		biny = self.emu.pt1_resolution.GetYaxis().FindBin(event.l1().Pt())
		sigma = self.emu.pt1_resolution.GetBinContent(binx,biny)
		smear1 = random.gauss(0.,sigma)
		event.l1.set_particle(event.l1()*(1+smear1))

		binx = self.emu.pt2_resolution.GetXaxis().FindBin(event.l2().Eta())
		biny = self.emu.pt2_resolution.GetYaxis().FindBin(event.l2().Pt())
		sigma = self.emu.pt2_resolution.GetBinContent(binx,biny)
		smear2 = random.gauss(0.,sigma)
		event.l2.set_particle(event.l2()*(1+smear2))

		event.l1.pt = event.l1().Pt()
		event.l1.eta = event.l1().Eta()
		event.l1.phi = event.l1().Phi()
		event.l1.E = event.l1().E()

		event.l2.pt = event.l2().Pt()
		event.l2.eta = event.l2().Eta()
		event.l2.phi = event.l2().Phi()
		event.l2.E = event.l2().E()

		if not all([
			event.l1.pt>15000. and any([abs(event.l1.eta)<1.37 or 1.52<abs(event.l1.eta)<2.5]), #electron selection
			event.l2.pt>10000. and abs(event.l2.eta)<2.5, #muon selection
			]):
			event.__break__ = True
			return

		if uncovered1 or uncovered2:

			event.__break__ = True
			return


		#Update sum energy information
		event.miss.set_particle(event.miss()-(event.l1()+event.l2()))
		event.sum_Et_miss += event.l1.pt
		event.sum_Et_miss += event.l2.pt

		event.__weight__/= inefficiency
		event.__weight__*= efficiency
		event.__weight__*= 0.06197796 #tautau branching ratio to emu
		event.lepton_class = 2 #now this is emu event

	def initialize_tools(self):

		analysis_home = os.getenv('ANALYSISHOME')
		mumu_file = '{0}/data/mumu_efficiency.root'.format(analysis_home)
		self.mumu = ROOT.TFile(mumu_file)
		emu_file = '{0}/data/emu_efficiency.root'.format(analysis_home)
		self.emu = ROOT.TFile(emu_file)
		"""
		load('TrigMuonEfficiency')
		#scale factor tool and config for EF_mu18_tight_mu8_EFFS
		self.muon_trigger_mu18_tight_mu8_EFFS = ROOT.LeptonTriggerSF(
			2012,
			'{0}/external/TrigMuonEfficiency/share'.format(analysis_home),
			'muon_trigger_sf_2012_AtoL.p1328.root',
			'{0}/external/ElectronEfficiencyCorrection/data'.format(analysis_home),
			"rel17p2.v02")
		self.config_muon_trigger_mu18_tight_mu8_EFFS = ROOT.TrigMuonEff.Configuration(True,False,False,False,-1,-1,0,'mu18_tight_mu8_EFFS','period_','fine')
		self.config_muon_trigger_mu18_tight_mu8_EFFS.isData = True
		"""

class mutation_scale(event_function):
	def __init__(self):
		event_function.__init__(self)
		self.initialize()

	def __call__(self,event):

		if event.mass_range == 0: return

		profile = self.Z_scale.lepton_dR_1_2_scale
		if event.lepton_dR > profile.GetBinLowEdge(profile.GetNbinsX()+1): weight_bin = profile.GetNbinsX()
		elif event.lepton_dR < profile.GetBinLowEdge(1): weight_bin = 1
		else: weight_bin = profile.FindBin(event.lepton_dR)
		weight = profile.GetBinContent(weight_bin)
		event.__weight__*=weight

	def initialize(self):
		analysis_home = os.getenv('ANALYSISHOME')
		Z_scale_file = '{0}/data/Z_scale.root'.format(analysis_home)
		self.Z_scale = ROOT.TFile(Z_scale_file)

class get_weight(event_function):
	def __init__(self,b=False):
		event_function.__init__(self)

		self.b = b
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
		if self.b: event.__weight__*=reduce(mul,[jet.bJet_scale_factor for jet in event.jets.values()],1)

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


class preselection_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
	
		if not all([
			event.jet_n>0,
			]):
			event.__break__=True
			return

class select_Z_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
	
		if not all([
			event.Mt1<75000.,
			event.Mt2<75000.,
			event.sum_Et_miss<175000.,
			abs(event.lepton_pair_miss_dPhi)<2.0,
			#abs(event.lepton_dPhi)<2.8,
			event.jet_n>0,
			]):
			event.__break__=True
			return

class select_signal_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		if not all([
			#not (event.lepton_pair_pT<10000. and event.lepton_class in [0,1]),
			#not (event.lepton_pair_pT<5000. and event.lepton_class ==2 ),
                        #event.lepton_dR<2.5,
			#event.jet_energy < 120000.,
			#10000.<event.missing_energy<60000.,
			#event.Mt1<40000.,
			#event.Mt2<40000.,
			event.Mt1<75000.,
			event.Mt2<75000.,
			event.sum_Et_miss<175000.,
			len(event.bjets)==1,
			]):
			event.__break__=True
			return

class select_tt_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		if not all([
			175000.<event.sum_Et_miss<250000.,
			len(event.bjets)>=1,
			event.Mt1>75000. or event.Mt2>75000.,
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
			'sum_Et_miss',
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

		#create lepton 1 (2 copies)
		event.l1 = particle(\
			**dict((name,event.__dict__['l1_'+name]) for name in self.lepton_names)
			)
		event.l1_original = particle(\
			**dict((name,event.__dict__['l1_'+name]) for name in self.lepton_names)
			)
		event.l1.set_pt_eta_phi_e(
			event.l1.pt,
			event.l1.eta,
			event.l1.phi,
			event.l1.E,
			)
		event.l1_original.set_pt_eta_phi_e(
			event.l1.pt,
			event.l1.eta,
			event.l1.phi,
			event.l1.E,
			)
		#create lepton 2
		event.l2 = particle(\
			**dict((name,event.__dict__['l2_'+name]) for name in self.lepton_names)
			)
		event.l2_original = particle(\
			**dict((name,event.__dict__['l2_'+name]) for name in self.lepton_names)
			)
		event.l2.set_pt_eta_phi_e(
			event.l2.pt,
			event.l2.eta,
			event.l2.phi,
			event.l2.E,
			)
		event.l2_original.set_pt_eta_phi_e(
			event.l2.pt,
			event.l2.eta,
			event.l2.phi,
			event.l2.E,
			)

		event.lepton_dR_original = event.l1().DeltaR(event.l2())
		event.lepton_pair_mass_low_original = (event.l1()+event.l2()).M()
		if event.l1().DeltaR(event.l2())<0.4:
			event.l2.ptcone40-=event.l1.pt
			event.l1.ptcone40-=event.l2.pt

		if event.l1.ptcone40<0.: event.l1.ptcone40=0.
		if event.l2.ptcone40<0.: event.l2.ptcone40=0.
		if event.l1.etcone20<0.: event.l1.etcone20=0.
		if event.l2.etcone20<0.: event.l2.etcone20=0.

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


		event.miss.set_particle(ROOT.TLorentzVector())
		event.sum_Et_miss = 0.

		for p in event.jets.values()+[event.l1,event.l2]:
			event.miss.set_particle(event.miss()-p())
			event.sum_Et_miss+= p().Et()

		sorted_jets = sorted(event.jets.values(),key=attrgetter('pt'), reverse=True) #jets sorted highest pt first
		lepton_pair = event.l1()+event.l2()

		event.missing_energy = -event.miss().Et()

		event.lepton_pair_pT = lepton_pair.Pt()
		event.lepton_pair_pT_diff = abs(event.l1.pt-event.l2.pt)
		event.lepton_pair_mass = lepton_pair.M()
		event.lepton_pair_mass_low = event.lepton_pair_mass

		event.lepton_pair_miss_dPhi = abs(event.miss().DeltaPhi(lepton_pair))
		event.l1_miss_dPhi = event.miss().DeltaPhi(event.l1())
		event.l2_miss_dPhi = event.miss().DeltaPhi(event.l2())

		try:
			#event.Mt1 = sqrt(2*(event.miss().Et()*event.l1().Et()-event.l1().Px()*event.miss().Px()-event.l1().Py()*event.miss().Py()))
			#event.Mt2 = sqrt(2*(event.miss().Et()*event.l2().Et()-event.l2().Px()*event.miss().Px()-event.l2().Py()*event.miss().Py()))
			#print event.miss().Et(),event.l1().Et(),(1-cos(event.l1_miss_dPhi))
			event.Mt1 = sqrt(-2*event.miss().Et()*event.l1().Et()*(1-cos(event.l1_miss_dPhi)))
			event.Mt2 = sqrt(-2*event.miss().Et()*event.l2().Et()*(1-cos(event.l2_miss_dPhi)))
		except:
			event.Mt1 = -1.
			event.Mt2 = -1.

		if event.lepton_class==0:
			event.off_threshold = min([event.l1.pt-25000.,event.l2.pt-15000.])
			event.collinear_mass = -1.
		elif event.lepton_class==1:
			event.off_threshold = min([event.l1.pt-25000.,event.l2.pt-10000.])
			event.collinear_mass = -1.
		else:
			event.off_threshold = min([event.l1.pt-15000.,event.l2.pt-10000.])
			try:
				event.collinear_mass = collinear_mass(event.l1(),event.l2(),event.miss())
			except ZeroDivisionError:
				event.collinear_mass = -1.
		if not event.lepton_pair_mass<150000.:
			event.__break__ = True
			return

		#event.mass_range = 0 if event.lepton_pair_mass<5000. else 1

		if event.lepton_pair_mass<5000.:
			event.__break__=True
			return

		for lepton,name in zip([event.l1,event.l2],['l1','l2']):
			for attr in ['pt','eta','phi']:
				setattr(event,name+'_'+attr,getattr(lepton,attr))
			setattr(event,name+'_etcone20_rat',lepton.etcone20/lepton.pt)
			setattr(event,name+'_ptcone40_rat',lepton.ptcone40/lepton.pt)

		event.l1.isolated = all([
			event.l1.etcone20/event.l1.pt<0.05,
			event.l1.ptcone40/event.l1.pt<0.08,
			])

		event.l2.isolated = all([
			event.l2.etcone20/event.l2.pt<0.05,
			event.l2.ptcone40/event.l2.pt<0.08,
			])
		
		if all([event.l1.isolated,event.l2.isolated]): event.isolated = 1
		elif not any([event.l1.isolated,event.l2.isolated]): event.isolated = 0
		else:
			event.__break__ = True
			return
		
		event.lepton_dR = abs(event.l1().DeltaR(event.l2()))
		event.lepton_dPhi = abs(event.l1().DeltaPhi(event.l2()))
		event.same_sign = 0 if (event.l1.charge*event.l2.charge)>0. else 1
		try: event.jet_energy = sum(jet.pt for jet in event.jets.values())
		except ValueError: event.jet_energy = 0.
		try: event.bjet_energy = sum(jet.pt for jet in event.bjets.values())
		except ValueError: event.bjet_energy = 0.
		if len(sorted_jets)>=1: 
			event.leading_jet_miss_dPhi = abs(event.miss().DeltaPhi(sorted_jets[0]()))
			event.l1_leading_jet_dR = abs(event.l1().DeltaR(sorted_jets[0]()))
			event.l2_leading_jet_dR = abs(event.l2().DeltaR(sorted_jets[0]()))
			event.lepton_pair_j1_dR = abs(sorted_jets[0]().DeltaR(lepton_pair))
		else: 
			event.leading_jet_miss_dPhi = -1.
			event.l1_leading_jet_dR = 0.
			event.l2_leading_jet_dR = 0.
		if len(sorted_jets)>=2: event.subleading_jet_miss_dPhi = abs(event.miss().DeltaPhi(sorted_jets[1]()))
		else: event.subleading_jet_miss_dPhi = -1.
		event.jet_n = len(event.jets)
		event.bjet_n = len(event.bjets)

from itertools import product

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low,xlabel)) for name,binning,high,low,xlabel in [
			('off_threshold',25,0.,25000.,"max(p_{T}^{l_{1}} - p_{T}^{off_{1}}, p_{T}^{l_{2}} - p_{T}^{off_{2}} [MeV]"),
			('sum_Et_miss',25,0.,250000.,"\Sigma E_{T} [MeV]"),
			('Mt1',26,-8000.,200000.,"M_{T}(l_{1}, MET) [MeV]"),
			('Mt2',26,-8000.,200000.,"M_{T}(l_{2}, MET) [MeV]"),
			('missing_energy',25,0.,100000.,"MET [MeV]"),
			('collinear_mass',20,0.,140000.,"M_{C}(l_{1}, l_{2}, MET) [MeV]"),
			('lepton_pair_mass',25,0.,150000.,"M(l_{1}, l_{2}) [MeV]"),
			('lepton_pair_mass_low',22,0.,45000.,"M(l_{1}, l_{2}) [MeV]"),
			#('lepton_pair_mass_low_original',22,0.,45000.,"M(\mu_{1}, \mu_{2}) [MeV]"),
			#('lepton_dR_original',60,0.,6.,"\Delta R(l_{1}, l_{2})"),
			('lepton_dR',15,0.,6.,"\DeltaR(l_{1}, l_{2})"),
			('lepton_pair_miss_dPhi',16,0.,3.2,"\Delta\phi(l_{1}+l_{2},MET)"),
			('lepton_pair_j1_dR',15,0.,6.,"\DeltaR(l_{1}+l_{2},j1)"),
			('l1_leading_jet_dR',15,0.,6.,"\DeltaR(l_{1}, j_{1})"),
			('l2_leading_jet_dR',15,0.,6.,"\DeltaR(l_{2}, j_{1})"),
			('lepton_dPhi',16,0.,3.2,"\Delta\phi(l_{1}, l_{2})"),
			('jet_energy',25,0.,200000.,"H_{T} [MeV]"),
			('bjet_energy',25,0.,200000.,"H_{T}^{b-tagged} [MeV]"),
			('leading_jet_miss_dPhi',21,-1,3.2,"\Delta\phi (j_{1},MET)"),
			('subleading_jet_miss_dPhi',21,-1,3.2,"\Delta\phi (j_{2},MET)"),
			('l1_miss_dPhi',16,0.,3.2,"\Delta\phi (l_{1},MET)"),
			('l2_miss_dPhi',16,0.,3.2,"\Delta\phi (l_{2},MET)"),
			('lepton_pair_pT',25,0.,100000.,"p_{T}^{l_{1} + l_{2}} [MeV]"),
			('lepton_pair_pT_diff',30,0.,60000.,"|p_{T}^{l_{1}} - p_{T}^{l_{2}}| [MeV]"),
			('l1_pt',20,0.,80000.,"p_{T}^{l_{1}} [MeV]"),
			('l1_ptcone40_rat',25,0.,0.2,"\Sigma^{\Delta R=0.4} p_{T}^{O}/p_{T}^{l_{1}}"),
			('l1_etcone20_rat',25,0.,0.2,"\Sigma^{\Delta R=0.2} E_{T}^{O}/p_{T}^{l_{1}}"),
			('l1_eta',24,-3.,3.,"\eta^{l_{1}}"),
			('l1_phi',32,-3.2,3.2,"\phi^{l_{1}}"),
			('l2_pt',15,0.,60000.,"p_{T}^{l_{2}} [MeV]"),
			('l2_ptcone40_rat',25,0.,0.2,"\Sigma^{\Delta R=0.4} p_{T}^{O}/p_{T}^{l_{2}}"),
			('l2_etcone20_rat',25,0.,0.2,"\Sigma^{\Delta R=0.2} E_{T}^{O}/p_{T}^{l_{2}}"),
			('l2_eta',24,-3.,3.,"\eta^{l_{2}}"),
			('l2_phi',32,-3.2,3.2,"\phi^{l_{2}}"),
			('jet_n',10,0,10,"# of jets"),
			('bjet_n',10,0,10,"# of b-tagged jets"),
			])

		self.names_2d = [
			('lepton_pair_miss_dPhi','lepton_pair_j1_dR'),
			('lepton_pair_miss_dPhi','lepton_pair_mass'),
			('lepton_pair_j1_dR','lepton_pair_mass'),
			#('lepton_pair_mass_low','lepton_pair_mass_low_original'),
			('Mt1','missing_energy'),
			('Mt2','missing_energy'),
			('Mt1','l1_miss_dPhi'),
			('Mt1','Mt2'),
			('l1_leading_jet_dR','l2_leading_jet_dR'),
			('Mt2','l2_miss_dPhi'),
			('sum_Et_miss','Mt1'),
			('sum_Et_miss','Mt2'),
			('sum_Et_miss','missing_energy'),
			('lepton_pair_mass','collinear_mass'),
			('lepton_pair_mass','off_threshold'),
			('lepton_pair_mass','missing_energy'),
			('collinear_mass','off_threshold'),
			('missing_energy','l1_pt'),
			('missing_energy','l2_pt'),
			('Mt1','l1_pt'),
			('Mt2','l2_pt'),
			('missing_energy','collinear_mass'),
			('collinear_mass','Mt1'),
			]

		for name_,(binning,high,low,xlabel) in self.names.items():
			for sign,isolated,lepton_class in product([0,1],[0,1],[0,1,2]):
				name = '{0}_{1}_{2}_{3}'.format(name_,sign,isolated,lepton_class)
				self.results[name] = ROOT.TH1F(name,name,binning,high,low)
				self.results[name].Sumw2()
				self.results[name].GetXaxis().SetTitle(xlabel)
				self.results[name].GetYaxis().SetTitle('Events')
				self.results[name].GetYaxis().CenterTitle()

		for name1,name2 in self.names_2d:
			binning1,high1,low1,xlabel = self.names[name1]
			binning2,high2,low2,ylabel = self.names[name2]
			for sign,isolated,lepton_class in product([0,1],[0,1],[0,1,2]):
				name = '{0}_{1}_{2}_{3}_{4}'.format(name1,name2,sign,isolated,lepton_class)
				self.results[name] = ROOT.TH2F(name,name,binning1,high1,low1,binning2,high2,low2)
				self.results[name].Sumw2()
				self.results[name].GetXaxis().SetTitle(xlabel)
				self.results[name].GetXaxis().CenterTitle()
				self.results[name].GetYaxis().SetTitle(ylabel)
				self.results[name].GetYaxis().CenterTitle()

	def __call__(self,event):
		if event.__break__: return

		weight = event.__weight__
		#weight*= -1 if event.same_sign else 1.

		for name_ in self.names:
			name = '{0}_{1}_{2}_{3}'.format(name_,event.same_sign,event.isolated,event.lepton_class)
			self.results[name].Fill(event.__dict__[name_],weight)

		for name1,name2 in self.names_2d:
			name = '{0}_{1}_{2}_{3}_{4}'.format(name1,name2,event.same_sign,event.isolated,event.lepton_class)
			self.results[name].Fill(event.__dict__[name1],event.__dict__[name2],weight)
