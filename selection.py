from common.analysis import analysis
from common.functions import event_function,result_function,arg,EventBreak
from common.external import load
from common.particle import particle

import ROOT
import os
from math import sqrt,pi,cos,sin
import json
import array
import random

from copy import deepcopy,copy
from operator import itemgetter, attrgetter, mul

class embedding(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_tautau(),
			)


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

class mutate_ee_make_selection_preselection(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_ee(),
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

class mutate_ee_make_selection_Z_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_ee(),
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

class make_selection_W_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			select_W_events()
			)

		self.add_result_function(
			plot_kinematics()
			)

		self.add_meta_result_function(
			)

class mutate_make_selection_W_control(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			mutate_mumu_to_tautau(),
			remove_overlapped_jets(),
			compute_kinematics(),
			get_weight(),
			select_W_events()
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

def get_mean_error_hist(hist,x,y):
		binx = hist.GetXaxis().FindBin(x)
		biny = hist.GetYaxis().FindBin(y)
		error = hist.GetBinError(binx,biny)
		mean = hist.GetBinContent(binx,biny)	
		return mean,error

def get_smear_hist(hist_file,particle,lepton,dist='pt'):

		if dist == 'pt':
			i = hist_file.pt_binning_resolution.FindBin(particle.pt)
			j = hist_file.eta_binning_resolution.FindBin(abs(particle.eta))
		
			if any([
				not 0<i<=hist_file.pt_binning_resolution.GetNbinsX(),
				not 0<j<=hist_file.eta_binning_resolution.GetNbinsX(),
				]):
				#print 'uncovered',lepton,i,j,particle.pt,particle.eta
				return None

			name = '{0}_pt_resolution_{1}_{2}'.format(lepton,i,j)
			resolution_histogram = getattr(hist_file,name)

		if dist == 'E':
			i = hist_file.pt_binning_resolution.FindBin(particle().E())
			j = hist_file.eta_binning_resolution.FindBin(abs(particle.eta))
		
			if any([
				not 0<i<=hist_file.pt_binning_resolution.GetNbinsX(),
				not 0<j<=hist_file.eta_binning_resolution.GetNbinsX(),
				]):
				#print 'uncovered',lepton,i,j,particle().E(),particle.eta
				return None

			name = '{0}_E_resolution_{1}_{2}'.format(lepton,i,j)

			resolution_histogram = getattr(hist_file,name)

		return resolution_histogram
	
def smear_particle_pt(hist_file,particle,lepton,dist='pt'):

		if dist == 'pt':
			i = hist_file.pt_binning_resolution.FindBin(particle.pt)
			j = hist_file.eta_binning_resolution.FindBin(abs(particle.eta))
		
			if any([
				not 0<i<=hist_file.pt_binning_resolution.GetNbinsX(),
				not 0<j<=hist_file.eta_binning_resolution.GetNbinsX(),
				]):
				print 'resolution uncovered',lepton,i,j,particle.pt,particle.eta
				return None

			name = '{0}_pt_resolution_{1}_{2}'.format(lepton,i,j)
			resolution_histogram = getattr(hist_file,name)

			smear = resolution_histogram.GetRandom()
			
			if smear==0.: print 'empty',lepton,i,j,particle.pt,particle.eta
		if dist == 'E':
			i = hist_file.pt_binning_resolution.FindBin(particle().E())
			j = hist_file.eta_binning_resolution.FindBin(abs(particle.eta))
		
			if any([
				not 0<i<=hist_file.pt_binning_resolution.GetNbinsX(),
				not 0<j<=hist_file.eta_binning_resolution.GetNbinsX(),
				]):
				print 'resolution uncovered',lepton,i,j,particle().E(),particle.eta
				return None

			name = '{0}_E_resolution_{1}_{2}'.format(lepton,i,j)

			resolution_histogram = getattr(hist_file,name)

			smear = resolution_histogram.GetRandom()

			if smear==0.: print 'empty',lepton,i,j,particle().E(),particle.eta
		
		#if abs(smear)<0.0000000000000001: print lepton,i,j,particle.pt,particle.eta
		#print lepton,i,j,particle.pt,particle.eta,smear

		particle.set_particle(particle()*(1+smear))
		particle.pt = particle().Pt()
		particle.eta = particle().Eta()
		particle.phi = particle().Phi()
		particle.E = particle().E()
		return smear
		
def smear_particle_eta(particle,smear):
		particle().SetEta(particle.eta+smear)
		particle.eta = particle().Eta()

def get_reco_efficiency(hist_file,l1_eta,l2_eta,l1_pt,l2_pt,debug=False):
		eta1 = hist_file.eta_binning.FindBin(abs(l1_eta))
		eta2 = hist_file.eta_binning.FindBin(abs(l2_eta))
		total_hist = hist_file.Get('reco_counts_eta_{0}_{1}'.format(eta1,eta2))
		selected_hist = hist_file.Get('trigger_counts_eta_{0}_{1}'.format(eta1,eta2))
		if debug: print eta1,eta2,total_hist,selected_hist
		if total_hist and selected_hist:
			binx = total_hist.GetXaxis().FindBin(l1_pt)
			biny = total_hist.GetYaxis().FindBin(l2_pt)
			total = total_hist.GetBinContent(binx,biny)
			selected = selected_hist.GetBinContent(binx,biny)
			if debug: print total,selected
			if total>0.: efficiency = selected/total 
			else: efficiency = -1.
		else: efficiency = -1.
		return efficiency
		
def get_selection_efficiency(hist_file,l1_eta,l2_eta,l1_pt,l2_pt,debug=False):
		eta1 = hist_file.eta_binning.FindBin(abs(l1_eta))
		eta2 = hist_file.eta_binning.FindBin(abs(l2_eta))
		total_hist = hist_file.Get('total_counts_eta_{0}_{1}'.format(eta1,eta2))
		selected_hist = hist_file.Get('id_counts_eta_{0}_{1}'.format(eta1,eta2))
		if debug: print eta1,eta2,total_hist,selected_hist
		if total_hist and selected_hist:
			binx = total_hist.GetXaxis().FindBin(l1_pt)
			biny = total_hist.GetYaxis().FindBin(l2_pt)
			total = total_hist.GetBinContent(binx,biny)
			selected = selected_hist.GetBinContent(binx,biny)
			if debug: print total,selected
			if total>0.: efficiency = selected/total 
			else: efficiency = -1.
		else: efficiency = -1.
		return efficiency

class mutate_mumu_to_tautau(event_function):

	class mumu_event(EventBreak): pass
	class min_inefficiency(EventBreak): pass
	class min_efficiency(EventBreak): pass
	class kinematic_cuts(EventBreak): pass
	
	def __init__(self):
		event_function.__init__(self)

		self.break_exceptions += [
			mutate_mumu_to_tautau.mumu_event,
			mutate_mumu_to_tautau.min_inefficiency,
			mutate_mumu_to_tautau.min_efficiency,
			mutate_mumu_to_tautau.kinematic_cuts,
			]
	
		from tauola import tauola_
		self.tauola = tauola_()

		self.electron_mass = 0.5/1000.
		self.muon_mass = 100.
		self.tau_mass = 1776.82

		self.lepton_names = [
			'eta',
			'phi',
			'pt',
			'E',
			]

		for name in self.lepton_names:
			for lepton in ['l1','l2']:
				self.create_branches[lepton+'_'+name] = 'float'

		self.create_branches['lepton_class'] = 'int'
		self.create_branches['mutation_weight'] = 'float'
	
		self.initialize_tools()

	def __call__(self,event):

		if not event.lepton_class == 1: raise mutate_mumu_to_tautau.mumu_event()

		if event.l1.pt < event.l2.pt: 
			event.l1,event.l2 = event.l2,event.l1

		inefficiency = get_selection_efficiency(self.inefficiency_file,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)
		if inefficiency <= 0.01: raise mutate_mumu_to_tautau.min_inefficiency()

		if random.getrandbits(1): event.l1,event.l2 = event.l2,event.l1 #flip e<->mu decay
		
		#do tauola decay
		tauola_call = []

		mother = event.l1()+event.l2()

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
		
		for particle in [
			event.l1,
			event.l2,
			]:
			particle.pt = particle().Pt()
			particle.eta = particle().Eta()
			particle.phi = particle().Phi()
			particle.E = particle().E()
			
		efficiency = get_selection_efficiency(self.efficiency_file,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)
		if efficiency < 0.: raise chain_weight.min_efficiency()

		if not all([
			event.l1().DeltaR(event.l2()) > 0.2,
			event.l1.pt>15000. and abs(event.l1.eta)<2.47, #electron selection
			event.l2.pt>10000. and abs(event.l2.eta)<2.5, #muon selection
			]):
			raise mutate_mumu_to_tautau.kinematic_cuts()

		event.mutation_weight = efficiency*0.0619779
		event.lepton_class = 2


		"""
		if not event.lepton_class==1: raise mutate_mumu_to_tautau.muons_event()

		etx = event.px_miss
		ety = event.py_miss
	
		for p in [event.l1,event.l2]:
			etx -= p().Et()*cos(p().Phi())
			ety -= p().Et()*sin(p().Phi())
			event.sum_Et_miss-= p().Et()

		#get reverse smeared muons (now we have smeared truth muons)
		for particle,hist in [
			(event.l1,self.mumu.l1_pt_resolution_reversed),
			(event.l2,self.mumu.l2_pt_resolution_reversed),
			]:
			smear = random.gauss(*get_mean_error_hist(hist,particle.eta,particle.pt))
			smear_particle_pt(particle,smear)
			#smear_particle_eta(event.l1,get_mean_error_hist(self.mumu.l1_eta_resolution_reversed,event.l1.eta,event.l1.pt)[0])
		

		inefficiency = get_efficiency(self.mumu,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)

		if random.getrandbits(1): event.l1,event.l2 = event.l2,event.l1 #flip e<->mu decay
	
		tauola_call = []

		mother = event.l1()+event.l2()

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
		smear_particle_pt(event.l1,1.)
		smear_particle_pt(event.l2,1.)
		
		efficiency = get_efficiency(self.emu,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)

		for particle,hist in [
			(event.l1,self.emu.l1_pt_resolution),
			(event.l2,self.emu.l2_pt_resolution),
			]:
			smear = random.gauss(*get_mean_error_hist(hist,particle.eta,particle.pt))
			smear_particle_pt(particle,smear)

		for name in self.lepton_names:
			for lepton in ['l1','l2']:
				overwrite_name = lepton+'_'+name
				new_value = getattr(getattr(event,lepton),name)
				setattr(event,overwrite_name,new_value)

		if event.l1().DeltaR(event.l2()) < 0.2:
			event.__break__=True
			return

		if not all([
			event.l1.pt>15000. and abs(event.l1.eta)<2.47, #electron selection
			event.l2.pt>10000. and abs(event.l2.eta)<2.5, #muon selection
			]):
			event.__break__ = True
			return

		if inefficiency < 0. or efficiency < 0.:
			event.__break__ = True
			return
			
		if inefficiency < 0.01: inefficiency = 0.01
		
		for p in [event.l1,event.l2]:
			etx += p().Et()*cos(p().Phi())
			ety += p().Et()*sin(p().Phi())
			event.sum_Et_miss += p().Et()

		event.miss.set_px_py_pz_e(-etx,-ety,0.,sqrt(etx**2.+ety**2.))
		#Update sum energy information
		#event.miss.set_particle(event.miss()-(event.l1()+event.l2()))
		#event.sum_Et_miss += event.l1.pt
		#event.sum_Et_miss += event.l2.pt

		event.__weight__/= inefficiency
		event.__weight__*= efficiency
		event.__weight__*= 0.06197796 #tautau branching ratio to emu
		event.lepton_class = 2 #now this is emu event


		"""
	def initialize_tools(self):

		analysis_home = os.getenv('ANALYSISHOME')
		file_name = '{0}/data/mumu_efficiency_alt2.root'.format(analysis_home)
		self.inefficiency_file = ROOT.TFile(file_name)
		if not self.inefficiency_file: raise RuntimeError('Unknown file {0}'.format(file_name))
		file_name = '{0}/data/ee_efficiency_alt2.root'.format(analysis_home)
		self.efficiency_file = ROOT.TFile(file_name)
		if not self.efficiency_file: raise RuntimeError('Unknown file {0}'.format(file_name))



"""
	
class mutate_mumu_to_tautau(event_function):

	class muons_event(EventBreak): pass

	def __init__(self):

		self.break_exceptions += [
			mutate_mumu_to_tautau.muons_event,
			]
	
		from tauola import tauola_
		self.min_mass = min_mass
		self.max_mass = max_mass
		event_function.__init__(self)

		self.electron_mass = 0.5/1000.
		self.muon_mass = 100.
		self.tau_mass = 1776.82

		self.electron_decay = array.array('d',[self.electron_mass,0.,0.])
		self.muon_decay = array.array('d',[self.muon_mass,0.,0.])
		self.tauola = tauola_()

		self.lepton_names = [
			'eta',
			'phi',
			'pt',
			'E',
			]

		for name in self.lepton_names:
			for lepton in ['l1','l2']:
				self.create_branches[lepton+'_'+name] = 'float'

		self.create_branches['lepton_class'] = 'int'
		
		self.initialize_tools()

	def __call__(self,event):

		if not event.lepton_class==1: raise mutate_mumu_to_tautau.muons_event()

		etx = event.px_miss
		ety = event.py_miss
	
		for p in [event.l1,event.l2]:
			etx -= p().Et()*cos(p().Phi())
			ety -= p().Et()*sin(p().Phi())
			event.sum_Et_miss-= p().Et()

		#get reverse smeared muons (now we have smeared truth muons)
		for particle,hist in [
			(event.l1,self.mumu.l1_pt_resolution_reversed),
			(event.l2,self.mumu.l2_pt_resolution_reversed),
			]:
			smear = random.gauss(*get_mean_error_hist(hist,particle.eta,particle.pt))
			smear_particle_pt(particle,smear)
			#smear_particle_eta(event.l1,get_mean_error_hist(self.mumu.l1_eta_resolution_reversed,event.l1.eta,event.l1.pt)[0])
		

		inefficiency = get_efficiency(self.mumu,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)

		if random.getrandbits(1): event.l1,event.l2 = event.l2,event.l1 #flip e<->mu decay
	
		tauola_call = []

		mother = event.l1()+event.l2()

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
		smear_particle_pt(event.l1,1.)
		smear_particle_pt(event.l2,1.)
		
		efficiency = get_efficiency(self.emu,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)

		for particle,hist in [
			(event.l1,self.emu.l1_pt_resolution),
			(event.l2,self.emu.l2_pt_resolution),
			]:
			smear = random.gauss(*get_mean_error_hist(hist,particle.eta,particle.pt))
			smear_particle_pt(particle,smear)

		for name in self.lepton_names:
			for lepton in ['l1','l2']:
				overwrite_name = lepton+'_'+name
				new_value = getattr(getattr(event,lepton),name)
				setattr(event,overwrite_name,new_value)

		if event.l1().DeltaR(event.l2()) < 0.2:
			event.__break__=True
			return

		if not all([
			event.l1.pt>15000. and abs(event.l1.eta)<2.47, #electron selection
			event.l2.pt>10000. and abs(event.l2.eta)<2.5, #muon selection
			]):
			event.__break__ = True
			return

		if inefficiency < 0. or efficiency < 0.:
			event.__break__ = True
			return
			
		if inefficiency < 0.01: inefficiency = 0.01
		
		for p in [event.l1,event.l2]:
			etx += p().Et()*cos(p().Phi())
			ety += p().Et()*sin(p().Phi())
			event.sum_Et_miss += p().Et()

		event.miss.set_px_py_pz_e(-etx,-ety,0.,sqrt(etx**2.+ety**2.))
		#Update sum energy information
		#event.miss.set_particle(event.miss()-(event.l1()+event.l2()))
		#event.sum_Et_miss += event.l1.pt
		#event.sum_Et_miss += event.l2.pt

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

class mutate_mumu_to_ee(event_function):
	def __init__(self,min_mass=0.,max_mass=1000000000.):
		self.min_mass = min_mass
		self.max_mass = max_mass
		event_function.__init__(self)


		self.electron_mass = 0.5/1000.
		self.muon_mass = 100.
		self.tau_mass = 1776.82

		self.electron_decay = array.array('d',[self.electron_mass,0.,0.])
		self.muon_decay = array.array('d',[self.muon_mass,0.,0.])

		self.lepton_names = [
			'eta',
			'phi',
			'pt',
			'E',
			'etcone20',
			'ptcone40',
			]

		for name in self.lepton_names:
			for lepton in ['l1','l2']:
				self.create_branches[lepton+'_'+name] = 'float'

		self.create_branches['lepton_class'] = 'int'

		self.initialize_tools()

	def __call__(self,event):

		if not all([
			event.l1.pt>20000., #first muon selection
			event.l2.pt>15000., #second muon selection
			]):
			event.__break__ = True
			return

		if not event.lepton_class==1:
			event.__break__=True
			return

		etx = event.px_miss
		ety = event.py_miss
	
		for p in [event.l1,event.l2]:
			etx -= p().Et()*cos(p().Phi())
			ety -= p().Et()*sin(p().Phi())
			event.sum_Et_miss-= p().Et()

		#get reverse smeared muons (now we have smeared truth muons)
		for particle,hist in [
			(event.l1,self.mumu.l1_pt_resolution_reversed),
			(event.l2,self.mumu.l2_pt_resolution_reversed),
			]:
			smear = random.gauss(*get_mean_error_hist(hist,particle.eta,particle.pt))
			smear_particle_pt(particle,smear)
			#smear_particle_eta(event.l1,get_mean_error_hist(self.mumu.l1_eta_resolution_reversed,event.l1.eta,event.l1.pt)[0])
		
		if event.l1.pt<event.l2.pt: event.l1,event.l2 = event.l2,event.l1
		
		mass = (event.l1()+event.l2()).M()
		
		inefficiency = get_efficiency(self.mumu,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)
		if inefficiency<0.: print 'inefficiency uncovered:', map(round,[event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt],[2]*4),mass
		efficiency = get_efficiency(self.ee,event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt)
		if efficiency<0.: print 'efficiency uncovered:', map(round,[event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt],[2]*4),mass

		if 0. <= inefficiency < 0.01: 
			print 'low inefficiency:', map(round,[event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt,mass],[2]*5),round(inefficiency,5)
			inefficiency = 0.01
			
		#leading electron passes isolation so should be excempt from any other cuts
		event.l1.etcone20 = 0.
		event.l2.ptcone40 = 0.

		for particle,hist in [
			(event.l1,self.ee.l1_pt_resolution),
			(event.l2,self.ee.l2_pt_resolution),
			]:
			smear = random.gauss(*get_mean_error_hist(hist,particle.eta,particle.pt))
			smear_particle_pt(particle,smear)

		if event.l1.pt<event.l2.pt: event.l1,event.l2 = event.l2,event.l1

		new_mass = (event.l1()+event.l2()).M()

		for name in self.lepton_names:
			for lepton in ['l1','l2']:
				overwrite_name = lepton+'_'+name
				new_value = getattr(getattr(event,lepton),name)
				setattr(event,overwrite_name,new_value)

		if not all([
			event.l1.pt>30000. and abs(event.l1.eta)<2.47, #electron selection
			event.l2.pt>20000. and abs(event.l2.eta)<2.47, #muon selection
			]):
			event.__break__ = True
			return

		if inefficiency < 0. or efficiency < 0.:
			event.__break__ = True
			return

		#if inefficiency < 0.01: 
		#	print 'low inefficiency:', map(round,[event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt],[2]*4),mass,inefficiency
		#	#inefficiency = 0.01

		if 80000. < new_mass < 100000. or 80000. < mass < 100000.:
			print 'mass window:', map(round,[event.l1.eta,event.l2.eta,event.l1.pt,event.l2.pt,mass,new_mass],[2]*6),round(inefficiency,4),round(efficiency,4),round(efficiency/inefficiency,4)
			
		for p in [event.l1,event.l2]:
			etx += p().Et()*cos(p().Phi())
			ety += p().Et()*sin(p().Phi())
			event.sum_Et_miss += p().Et()

		event.miss.set_px_py_pz_e(-etx,-ety,0.,sqrt(etx**2.+ety**2.))
		#Update sum energy information
		#event.miss.set_particle(event.miss()-(event.l1()+event.l2()))
		#event.sum_Et_miss += event.l1.pt
		#event.sum_Et_miss += event.l2.pt

		event.__weight__/= inefficiency
		event.__weight__*= efficiency
		event.lepton_class = 0 #now this is ee event
		
	def initialize_tools(self):

		analysis_home = os.getenv('ANALYSISHOME')
		mumu_file = '{0}/data/mumu_efficiency.root'.format(analysis_home)
		self.mumu = ROOT.TFile(mumu_file)
		ee_file = '{0}/data/ee_efficiency.root'.format(analysis_home)
		self.ee = ROOT.TFile(ee_file)

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
"""
class get_weight(event_function):
	def __init__(
		self,
		b=False,
		l1_fluctuation=arg(0.,help='Fluctation on l1 scale factor'),
		l2_fluctuation=arg(0.,help='Fluctation on l1 scale factor'),
		trigger_fluctuation=arg(0.,help='Fluctation on trigger scale factor'),
		bjet_fluctuation=arg(0.,help='Fluctation on b-jet ID scale factor'),
		):
		event_function.__init__(self)

		self.b = b
		self.l1_fluctuation = l1_fluctuation
		self.l2_fluctuation = l2_fluctuation
		self.trigger_fluctuation = trigger_fluctuation
		self.bjet_fluctuation = bjet_fluctuation

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
			event.l1_scale_factor+self.l1_fluctuation*event.l1_scale_factor_error,
			event.l2_scale_factor+self.l2_fluctuation*event.l2_scale_factor_error,
			event.trigger_scale_factor+self.trigger_fluctuation*event.trigger_scale_factor_error,
			event.weight_pileup,
			]: event.__weight__*=weight
		if self.b: event.__weight__*=reduce(mul,[jet.bJet_scale_factor+jet.bJet_scale_factor_error*self.bjet_fluctuation if jet_n in event.bjets else jet.bJet_scale_factor-jet.bJet_scale_factor_error*self.bjet_fluctuation for jet_n,jet in event.jets.items()],1)

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

	class lepton_pt(EventBreak): pass
	class one_jet(EventBreak): pass

	def __init__(self):
		event_function.__init__(self)

		self.break_exceptions += [
			preselection_events.lepton_pt,
			preselection_events.one_jet,
			]

	def __call__(self,event):
	
		for requirement,exception in [
			(any([
				event.lepton_class == 0 and all([
					event.l1.pt>30000.,
					event.l2.pt>20000.,
					]),
				event.lepton_class == 1 and all([
					event.l1.pt>30000.,
					event.l2.pt>20000.,
					]),
				event.lepton_class == 2 and all([
					event.l1.pt>15000.,
					event.l2.pt>10000.,
					]),				
				]),preselection_events.lepton_pt),
			(event.jet_n>0,preselection_events.one_jet),
			]:
			if not requirement: raise exception()

"""
class select_Z_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
	
		if not all([
			#event.Mt1<75000.,
			#event.Mt2<75000.,
			event.sum_Mt<70000.,
			#event.sum_Et_miss<175000.,
			event.sum_Et_miss<175000.-35./24.*event.sum_Mt,
			event.miss_direction_lepton_pair>event.lepton_pair_pT-40000.,
			event.subleading_jet_pT<30000.,
			#not (event.lepton_pair_miss_dPhi>pi/2 and event.lepton_pair_pT>30000.),
			#abs(event.l2_fraction-event.l1_fraction)<0.15,
			#event.l1_fraction*event.l2_fraction>0.,
			#not (event.lepton_pair_mass<20000. and event.missing_energy>
			#abs(event.lepton_dPhi)<2.8,
			event.jet_n>0,
			]):
			event.__break__=True
			return

class select_W_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
	
		if not all([
			#event.Mt1<75000.,
			#event.Mt2<75000.,
			event.sum_Et_miss<150000.,
			event.sum_Mt>100000.,
			#event.miss_direction_lepton_pair<(4./5.*event.lepton_pair_pT-40000.),
			#not (event.lepton_pair_miss_dPhi>pi/2 and event.lepton_pair_pT>30000.),
			#abs(event.l2_fraction-event.l1_fraction)<0.15,
			#event.l1_fraction*event.l2_fraction>0.,
			#not (event.lepton_pair_mass<20000. and event.missing_energy>
			#abs(event.lepton_dPhi)<2.8,
			event.jet_n>0,
			]):
			event.__break__=True
			return
"""

class select_signal_events(event_function):

	class sum_Mt(EventBreak): pass
	class miss_direction(EventBreak): pass
	#class subleading_jet(EventBreak): pass
	class one_bjet(EventBreak): pass
	
	def __init__(self):
		event_function.__init__(self)


	def __init__(self):
		event_function.__init__(self)

		self.break_exceptions += [
			select_signal_events.sum_Mt,
			select_signal_events.miss_direction,
			#select_signal_events.subleading_jet,
			select_signal_events.one_bjet,
			]


	def __call__(self,event):

		for requirement,exception in [
			#(event.sum_Mt<70000.,select_signal_events.sum_Mt),
			#(event.miss_direction_lepton_pair>event.lepton_pair_pT-40000.,select_signal_events.miss_direction),
			#(event.subleading_jet_pT<30000.,select_signal_events.subleading_jet),
			(len(event.bjets)==1,select_signal_events.one_bjet),
			]:
			if not requirement: raise exception()
			
"""
class select_tt_events(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		if not all([
			175000.<event.sum_Et_miss<250000.,
			#event.sum_Mt>70000.,
			#event.sum_Et_miss>175000.-35./24.*event.sum_Mt,
			len(event.bjets)>=1,
			event.Mt1>75000. or event.Mt2>75000.,
			]):
			event.__break__=True
			return

"""

class build_events(event_function):

	class heavy_flavor_removal(EventBreak): pass

	def __init__(
		self,
		jvf_fluctuation=arg(0,help='Fluctuation jvf cut choose between [-1,0,1] to fluctuate cut down, nominal and up respectively'),
		):

		event_function.__init__(self)
		
		self.break_exceptions += [
			build_events.heavy_flavor_removal,
			]

		self.jvf_fluctuation = jvf_fluctuation
		#print jvf_fluctuation,type(self.jvf_fluctuation)
		#print '__deferred_init__'
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
			'bJet_scale_factor_error',
			'jvf_down_cut',
			'jvf_up_cut',
			#'passed_b_preselection',
			]
		self.create_branches['top_hfor_type'] = None

		self.required_branches += ['jet_'+name for name in self.jet_names]
		self.required_branches += ['jet_n']

		self.create_branches['jets'] = None
		self.create_branches['l1'] = None
		self.create_branches['l2'] = None

	def __call__(self,event):

		event.inefficiency_weight = 1.
		event.efficiency_weight = 1.
		event.total_efficiency_weight = 1.

		#collect jets
		event.jets = {}		
		for jet in range(event.jet_n):
			if not all([
				#requirements
				not ((abs(event.jet_eta[jet])<2.4 and event.jet_pt[jet]<50000.) and not any([
					event.jet_jvf[jet]>event.jet_jvf_down_cut[jet] and self.jvf_fluctuation == -1,
					event.jet_jvf[jet]>0.5 and self.jvf_fluctuation == 0,
					event.jet_jvf[jet]>event.jet_jvf_up_cut[jet] and self.jvf_fluctuation == 1,
					])),
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

		event.l1 = particle(\
			**dict((name,event.__dict__['l1_'+name]) for name in self.lepton_names)
			)

		event.l1.set_pt_eta_phi_e(
			event.l1.pt,
			event.l1.eta,
			event.l1.phi,
			event.l1.E,
			)

		event.l2 = particle(\
			**dict((name,event.__dict__['l2_'+name]) for name in self.lepton_names)
			)

		event.l2.set_pt_eta_phi_e(
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

		if getattr(event,'top_hfor_type',0)==4: raise build_events.heavy_flavor_removal()

class remove_overlapped_jets(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		#remove jets from electrons, muons
		for jetN,jet in event.jets.items():
			if jet().DeltaR(event.l2())<0.4:
				del event.jets[jetN]
				if jetN in event.bjets_preselected: del event.bjets_preselected[jetN]
				if jetN in event.bjets: del event.bjets[jetN]
				continue
			if jet().DeltaR(event.l1())<0.4:
				del event.jets[jetN]
				if jetN in event.bjets_preselected: del event.bjets_preselected[jetN]
				if jetN in event.bjets: del event.bjets[jetN]

def collinear_mass(l1,l2,miss):
	m_frac_1 = ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())) / ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())+(l2.Py()*miss.Px())-(l2.Px()*miss.Py()))
	m_frac_2 = ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())) / ((l1.Px()*l2.Py())-(l1.Py()*l2.Px())+(l1.Px()*miss.Py())-(l1.Py()*miss.Px()))

	if m_frac_1*m_frac_2 > 0.: return (l1+l2).M()/sqrt(m_frac_1*m_frac_2)
	return -1.

class compute_kinematics(event_function):

	class sign_requirement(EventBreak): pass
	class lepton_class(EventBreak): pass
	class partial_isolation_requirement(EventBreak): pass
	class isolation_requirement(EventBreak): pass
	class lepton_pair_mass_window(EventBreak): pass

	def __init__(
		self,
		opposite_sign=arg(1,help='Sign required of leptons {0:same-sign,1:opposite-sign}'),
		lepton_class=arg(2,help='Sign of leptons {0:ee,1:mumu,2:emu}'),
		lower_mass_window=arg(0.,help='lower cut on mass window [GeV]'),
		upper_mass_window=arg(100.,help='upper cut on mass window [GeV]'),
		l1_isolated=arg(1,help='l1 is isolated {0:False,1:True}'),
		l2_isolated=arg(1,help='l2 is isolated {0:False,1:True}'),
		):
		event_function.__init__(self)

		self.break_exceptions += [
			compute_kinematics.sign_requirement,
			compute_kinematics.lepton_class,
			compute_kinematics.partial_isolation_requirement,
			compute_kinematics.isolation_requirement,
			compute_kinematics.lepton_pair_mass_window,
			]

		self.opposite_sign = bool(opposite_sign)
		self.lepton_class = lepton_class
		self.l1_isolated = bool(l1_isolated)
		self.l2_isolated = bool(l2_isolated)
		self.lower_mass_window = lower_mass_window*1000.
		self.upper_mass_window = upper_mass_window*1000.
		
	def __call__(self,event):
	
		event.opposite_sign = event.l1.charge*event.l2.charge<0.
		if event.opposite_sign is not self.opposite_sign: raise compute_kinematics.sign_requirement()
		if event.lepton_class != self.lepton_class: raise compute_kinematics.lepton_class()

		event.l1.partially_isolated = all([
			event.l1.etcone20/event.l1.pt<0.15,
			event.l1.ptcone40/event.l1.pt<0.3,
			])

		event.l2.partially_isolated = all([
			event.l2.etcone20/event.l2.pt<0.15,
			event.l2.ptcone40/event.l2.pt<0.3,
			])
			
		if not all([event.l1.partially_isolated,event.l2.partially_isolated]): raise compute_kinematics.partial_isolation_requirement()

		event.l1.isolated = all([
			event.l1.etcone20/event.l1.pt<0.06,
			event.l1.ptcone40/event.l1.pt<0.16,
			])

		event.l2.isolated = all([
			event.l2.etcone20/event.l2.pt<0.06,
			event.l2.ptcone40/event.l2.pt<0.16,
			])


		if not all([
			event.l1.isolated is self.l1_isolated,
			event.l2.isolated is self.l2_isolated,
			]): raise compute_kinematics.isolation_requirement()

		#compute missing energy/sum Et
		event.sum_Et_miss = 0.
		etx = 0.
		ety = 0.
		for p in event.jets.values()+[event.l1,event.l2]:
			etx += p().Et()*cos(p().Phi())
			ety += p().Et()*sin(p().Phi())
			event.sum_Et_miss+= p().Et()
		event.miss = particle()
		event.miss.set_px_py_pz_e(-etx,-ety,0.,sqrt(etx**2.+ety**2.))

		lepton_pair = event.l1()+event.l2()

		event.miss_phi = event.miss().Phi()
		event.missing_energy = event.miss().Et()

		event.lepton_pair_pT = lepton_pair.Pt()
		event.lepton_pair_pT_diff = abs(event.l1.pt-event.l2.pt)
		event.lepton_pair_mass = lepton_pair.M()
		event.lepton_pair_mass_low = event.lepton_pair_mass
		event.lepton_pair_mass_high = event.lepton_pair_mass
		event.lepton_dR = abs(event.l1().DeltaR(event.l2()))
		event.lepton_dPhi = abs(event.l1().DeltaPhi(event.l2()))
				
		if not self.lower_mass_window<event.lepton_pair_mass<self.upper_mass_window: raise compute_kinematics.lepton_pair_mass_window()
		
		event.lepton_pair_miss_dPhi = abs(event.miss().DeltaPhi(lepton_pair))
		event.miss_direction_lepton_pair = event.missing_energy*cos(event.lepton_pair_miss_dPhi)
		event.lepton_pair_pT_direction_miss = event.lepton_pair_pT*cos(event.lepton_pair_miss_dPhi)
		event.l1_miss_dPhi = abs(event.miss().DeltaPhi(event.l1()))
		event.l2_miss_dPhi = abs(event.miss().DeltaPhi(event.l2()))
		event.sum_l1_miss_dPhi_l2_miss_dPhi = event.l1_miss_dPhi+event.l2_miss_dPhi

		try:
			event.Mt1 = sqrt(2*event.miss().Et()*event.l1().Et()*(1-cos(event.l1_miss_dPhi)))
			event.Mt2 = sqrt(2*event.miss().Et()*event.l2().Et()*(1-cos(event.l2_miss_dPhi)))
		except:
			event.Mt1 = -1.
			event.Mt2 = -1.

		event.sum_Mt = event.Mt1+event.Mt2

		metx = event.miss().Px()
		mety = event.miss().Py()

		try: event.l1_fraction = ((event.l1().Px()*event.l2().Py())-(event.l1().Py()*event.l2().Px())) / ((event.l1().Px()*event.l2().Py())-(event.l1().Py()*event.l2().Px())+(event.l2().Py()*metx)-(event.l2().Px()*mety))
		except ZeroDivisionError: event.l1_fraction = -4.

		try: event.l2_fraction = ((event.l1().Px()*event.l2().Py())-(event.l1().Py()*event.l2().Px())) / ((event.l1().Px()*event.l2().Py())-(event.l1().Py()*event.l2().Px())+(event.l1().Px()*mety)-(event.l1().Py()*metx))
		except ZeroDivisionError: event.l2_fraction = -4.

		if event.l1_fraction*event.l2_fraction > 0. and event.l2_fraction>-4. and event.l2_fraction>-4.: event.collinear_mass = event.lepton_pair_mass/sqrt(event.l1_fraction*event.l2_fraction)
		else: event.collinear_mass = -1.

		if event.lepton_class==0:
			event.off_threshold = min([event.l1.pt-30000.,event.l2.pt-20000.])
			event.collinear_mass = -1.
		elif event.lepton_class==1:
			event.off_threshold = min([event.l1.pt-25000.,event.l2.pt-10000.])
			event.collinear_mass = -1.
		else:
			event.off_threshold = min([event.l1.pt-15000.,event.l2.pt-10000.])

		for lepton,name in zip([event.l1,event.l2],['l1','l2']):
			for attr in ['pt','eta','phi']:
				setattr(event,name+'_'+attr,getattr(lepton,attr))
			setattr(event,name+'_etcone20_rat',lepton.etcone20/lepton.pt)
			setattr(event,name+'_ptcone40_rat',lepton.ptcone40/lepton.pt)

		sorted_jet_keys = sorted(event.jets.keys(), key = lambda index: event.jets[index].pt, reverse=True)
		sorted_jets = sorted(event.jets.values(),key=attrgetter('pt'), reverse=True) #jets sorted highest pt first

		try: event.jet_energy = sum(jet.pt for jet in event.jets.values())
		except ValueError: event.jet_energy = 0.
		try: event.bjet_energy = sum(jet.pt for jet in event.bjets.values())
		except ValueError: event.bjet_energy = 0.

		if len(sorted_jet_keys)>1:
			if sorted_jet_keys[1] in event.bjets and sorted_jet_keys[0] not in event.bjets:
				sorted_jets[0:2]=reversed(sorted_jets[0:2])

		if len(sorted_jets)>=1: 
			event.leading_jet_miss_dPhi = abs(event.miss().DeltaPhi(sorted_jets[0]()))
			event.l1_leading_jet_dR = abs(event.l1().DeltaR(sorted_jets[0]()))
			event.l2_leading_jet_dR = abs(event.l2().DeltaR(sorted_jets[0]()))
			event.lepton_pair_j1_dR = abs(sorted_jets[0]().DeltaR(lepton_pair))
			event.lepton_pair_jet_mass = (sorted_jets[0]()+lepton_pair).M()
			event.leading_jet_pT = sorted_jets[0]().Pt()
		else: 
			event.leading_jet_miss_dPhi = -1.
			event.l1_leading_jet_dR = 0.
			event.l2_leading_jet_dR = 0.
			event.lepton_pair_j1_dR = 0.
			event.lepton_pair_jet_mass = -1.
			event.leading_jet_pT = -1.
		if len(sorted_jets)>=2: 
			event.subleading_jet_miss_dPhi = abs(event.miss().DeltaPhi(sorted_jets[1]()))
			event.lepton_pair_2jet_mass = (sorted_jets[0]()+sorted_jets[1]()+lepton_pair).M()
			event.subleading_jet_pT = sorted_jets[1]().Pt()
		else: 
			event.subleading_jet_miss_dPhi = -1.
			event.lepton_pair_2jet_mass = -1.
			event.subleading_jet_pT = 0.

		event.jet_n = len(event.jets)
		event.bjet_n = len(event.bjets)

		if self.sign_requirement==1:
			if event.l1_charge<0.: l = copy(event.l1())
			else: l = copy(event.l2())
		else:
			if random.getrandbits(1): l = copy(event.l1())
			else: l = copy(event.l2())

		b = lepton_pair.BoostVector()
		l.Boost(-b)
		a = l.Angle(b)
		#if a>pi/2.: a=pi-a
		event.cos_helicity_angle = abs(cos(a))

		l1 = ROOT.TLorentzVector()
		l1.SetPxPyPzE(event.l1().Px(),event.l1().Py(),0.,sqrt(event.l1().Px()**2.+event.l1().Py()**2.))
		l2 = ROOT.TLorentzVector()
		l2.SetPxPyPzE(event.l2().Px(),event.l2().Py(),0.,sqrt(event.l2().Px()**2.+event.l2().Py()**2.))
		m = ROOT.TLorentzVector()
		m.SetPxPyPzE(event.miss().Px(),event.miss().Py(),0.,sqrt(event.miss().Px()**2.+event.miss().Py()**2.))
		b = (l1+l2+m).BoostVector()
		l1.Boost(-b)
		l2.Boost(-b)
		m.Boost(-b)
		event.transverse_com_l1_l2_dPhi = abs(l1.DeltaPhi(l2))
		event.transverse_com_l1_miss_dPhi = abs(l1.DeltaPhi(m))
		event.transverse_com_l2_miss_dPhi = abs(l2.DeltaPhi(m))
		event.transverse_com_mass = (l1+l2+m).M()

from itertools import product

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low,xlabel)) for name,binning,high,low,xlabel in [
			('efficiency_weight',22,0.,1.1,"Efficiency Weight"),
			('inefficiency_weight',22,0.,1.1,"Inefficiency Weight"),
			('total_efficiency_weight',20,0.,5.0,"Total Efficiency Weight"),
			('transverse_com_l1_l2_dPhi',16,0.,3.2,"\Delta\phi(l_{1},l_{2})"),
			('transverse_com_l1_miss_dPhi',16,0.,3.2,"\Delta\phi(l_{1},MET)"),
			('transverse_com_l2_miss_dPhi',16,0.,3.2,"\Delta\phi(l_{2},MET)"),
			('transverse_com_mass',25,0.,150000.,"M_{T}(l_{1},l_{2},MET) [MeV]"),
			('off_threshold',25,0.,25000.,"max(p_{T}^{l_{1}} - p_{T}^{off_{1}}, p_{T}^{l_{2}} - p_{T}^{off_{2}} [MeV]"),
			('lepton_pair_pT_direction_miss',20,-50000.,50000.,r"p_{T}^{l_{1}+l_{2}} \times cos(\phi^{MET}-\phi^{l_{1}+l_{2}}) [MeV]"),
			('miss_direction_lepton_pair',24,-50000.,70000.,r"MET \times cos(\phi^{MET}-\phi^{l_{1}+l_{2}}) [MeV]"),
			('sum_Et_miss',25,0.,250000.,"\Sigma E_{T} [MeV]"),
			('sum_Mt',25,0.,200000.,"M_{T}(l_{1},MET) + M_{T}(l_{2},MET) [MeV]"),
			('Mt1',25,0.,200000.,"M_{T}(l_{1},MET) [MeV]"),
			('Mt2',25,0.,200000.,"M_{T}(l_{2},MET) [MeV]"),
			#('miss_miss_original_dPhi',16,0.,3.2,"\Delta\phi(MET,MET_{0})"
			('miss_phi',32,-3.2,3.2,"\phi^{MET}"),
			('missing_energy',25,0.,100000.,"MET [MeV]"),
			('collinear_mass',21,-7000.,140000.,"M_{C}(l_{1},l_{2},MET) [MeV]"),
			('lepton_pair_mass',20,0.,100000.,"M(l_{1},l_{2}) [MeV]"),
			('lepton_pair_mass_low',20,0.,40000.,"M(l_{1},l_{2}) [MeV]"),
			('lepton_pair_mass_high',20,60000.,100000.,"M(l_{1},l_{2}) [MeV]"),
			#('lepton_pair_mass_low_original',22,0.,45000.,"M(\mu_{1}, \mu_{2}) [MeV]"),
			#('lepton_dR_original',60,0.,6.,"\Delta R(l_{1}, l_{2})"),
			('lepton_pair_jet_mass',20,0.,200000.,"M(l_{1},l_{2},j_{1}) [MeV]"),
			('lepton_pair_2jet_mass',20,0.,200000.,"M(l_{1},l_{2},j_{1},j_{2}) [MeV]"),
			('lepton_dR',15,0.,6.,"\DeltaR(l_{1}, l_{2})"),
			('lepton_pair_miss_dPhi',16,0.,3.2,"\Delta\phi(l_{1}+l_{2},MET)"),
			('lepton_pair_j1_dR',15,0.,6.,"\DeltaR(l_{1}+l_{2},j1)"),
			('l1_leading_jet_dR',15,0.,6.,"\DeltaR(l_{1},j_{1})"),
			('l2_leading_jet_dR',15,0.,6.,"\DeltaR(l_{2},j_{1})"),
			('lepton_dPhi',16,0.,3.2,"\Delta\phi(l_{1},l_{2})"),
			('jet_energy',15,0.,150000.,"H_{T} [MeV]"),
			('leading_jet_pT',20,0.,80000.,"p_{T}^{j_{1}} [MeV]"),
			('subleading_jet_pT',20,0.,80000.,"p_{T}^{j_{2}} [MeV]"),
			('bjet_energy',20,0.,80000.,"H_{T}^{b-tagged} [MeV]"),
			('leading_jet_miss_dPhi',21,-1,3.2,"\Delta\phi(j_{1},MET)"),
			('subleading_jet_miss_dPhi',21,-1,3.2,"\Delta\phi(j_{2},MET)"),
			('cos_helicity_angle',10,0,1.,r"Cos(\theta^{*})"),
			('l1_miss_dPhi',16,0.,3.2,"\Delta\phi (l_{1},MET)"),
			('l2_miss_dPhi',16,0.,3.2,"\Delta\phi (l_{2},MET)"),
			('sum_l1_miss_dPhi_l2_miss_dPhi',16,0.,6.4,"\Delta\phi (l_{1},MET)+\Delta\phi (l_{2},MET)"),
			('lepton_pair_pT',25,0.,100000.,"p_{T}^{l_{1}+l_{2}} [MeV]"),
			('lepton_pair_pT_diff',30,0.,60000.,"|p_{T}^{l_{1}}-p_{T}^{l_{2}}| [MeV]"),
			('l1_pt',14,0.,70000.,"p_{T}^{l_{1}} [MeV]"),
			('l1_ptcone40_rat',17,0.,0.34,"\Sigma^{\Delta R=0.4} p_{T}^{O}/p_{T}^{l_{1}}"),
			('l1_etcone20_rat',10,0.,0.2,"\Sigma^{\Delta R=0.2} E_{T}^{O}/p_{T}^{l_{1}}"),
			('l1_eta',24,-3.,3.,"\eta^{l_{1}}"),
			('l1_phi',32,-3.2,3.2,"\phi^{l_{1}}"),
			('l1_fraction',120,-4.,4.,"l_{1} energy fraction"),
			('l2_pt',14,0.,70000.,"p_{T}^{l_{2}} [MeV]"),
			('l2_ptcone40_rat',17,0.,0.34,"\Sigma^{\Delta R=0.4} p_{T}^{O}/p_{T}^{l_{2}}"),
			('l2_etcone20_rat',10,0.,0.2,"\Sigma^{\Delta R=0.2} E_{T}^{O}/p_{T}^{l_{2}}"),
			('l2_eta',24,-3.,3.,"\eta^{l_{2}}"),
			('l2_phi',32,-3.2,3.2,"\phi^{l_{2}}"),
			('l2_fraction',120,-4.,4.,"l_{2} energy fraction"),
			('jet_n',4,0,4,"# of jets"),
			('bjet_n',4,0,4,"# of b-tagged jets"),
			])

		self.names_2d = [
			('l1_ptcone40_rat','l1_etcone20_rat'),
			('l2_ptcone40_rat','l2_etcone20_rat'),
			('efficiency_weight','lepton_pair_mass'),
			('inefficiency_weight','lepton_pair_mass'),
			('total_efficiency_weight','lepton_pair_mass'),
			('sum_Mt','collinear_mass'),
			('sum_Mt','lepton_pair_mass'),
			('sum_Mt','sum_Et_miss'),
			('sum_Mt','cos_helicity_angle'),
			('sum_Mt','jet_energy'),
			('sum_Mt','leading_jet_pT'),
			('sum_Mt','lepton_pair_pT_direction_miss'),
			('leading_jet_pT','miss_direction_lepton_pair'),
			('jet_energy','miss_direction_lepton_pair'),
			('lepton_pair_mass','l1_leading_jet_dR'),
			('lepton_pair_pT_diff','lepton_pair_miss_dPhi'),
			('lepton_pair_pT_diff','miss_direction_lepton_pair'),
			('jet_energy','missing_energy'),
			('lepton_pair_mass','transverse_com_l1_l2_dPhi'),
			('lepton_pair_mass','transverse_com_l1_miss_dPhi'),
			('lepton_pair_mass','transverse_com_l2_miss_dPhi'),
			('lepton_pair_mass','transverse_com_mass'),
			('lepton_pair_mass','cos_helicity_angle'),
			('lepton_pair_mass','lepton_dPhi'),
			('lepton_pair_mass','lepton_dR'),
			('l1_fraction','l2_fraction'),
			('lepton_pair_mass','leading_jet_pT'),
			('leading_jet_pT','missing_energy'),
			('leading_jet_pT','lepton_pair_pT'),
			('leading_jet_pT','lepton_pair_miss_dPhi'),
			('lepton_pair_mass','l1_etcone20_rat'),
			('lepton_pair_mass','l2_etcone20_rat'),
			('lepton_pair_pT','lepton_pair_miss_dPhi'),
			('lepton_pair_pT','miss_direction_lepton_pair'),
			('lepton_pair_mass','miss_direction_lepton_pair'),
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
			('lepton_pair_mass','lepton_pair_pT'),
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

		for name,(binning,high,low,xlabel) in self.names.items():
			self.results[name] = ROOT.TH1F(name,name,binning,high,low)
			self.results[name].Sumw2()
			self.results[name].GetXaxis().SetTitle(xlabel)
			self.results[name].GetYaxis().SetTitle('Events')
			self.results[name].GetYaxis().CenterTitle()
		
		for name1,name2 in self.names_2d:
			binning1,high1,low1,xlabel = self.names[name1]
			binning2,high2,low2,ylabel = self.names[name2]
			name = '{0}_{1}'.format(name1,name2)
			self.results[name] = ROOT.TH2F(name,name,binning1,high1,low1,binning2,high2,low2)
			self.results[name].Sumw2()
			self.results[name].GetXaxis().SetTitle(xlabel)
			self.results[name].GetXaxis().CenterTitle()
			self.results[name].GetYaxis().SetTitle(ylabel)
			self.results[name].GetYaxis().CenterTitle()

	def __call__(self,event):
		if event.__break__: return

		for name in self.names:
			self.results[name].Fill(event.__dict__[name],event.__weight__)
		for name1,name2 in self.names_2d:
			name = '{0}_{1}'.format(name1,name2)
			self.results[name].Fill(event.__dict__[name1],event.__dict__[name2],event.__weight__)



