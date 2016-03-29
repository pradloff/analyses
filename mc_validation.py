from common.functions import event_function, EventBreak
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
from unweight_mcfm import decay_fermions_as_taus
import ROOT
import os
from math import sqrt
import json
from itertools import products

class truth_analysis_sherpa(analysis):
	def __init__(self):
		analysis.__init__(self)

		self.add_event_function(
			truth_tree(pdgIds = [5,-5,15,-15,11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,25]),
			identify_sherpa_truth(),
			collect_truth_jets(),
			match_truth_jets(),
			select_emu_events(),
			build_events(),
			)

		self.add_result_function(
			plot_kinematics(),
			)

		self.add_meta_result_function(
			)

class truth_analysis_pythia(analysis):
	def __init__(self):
		analysis.__init__(self)

		self.add_event_function(
			truth_tree(pdgIds = [5,-5,15,-15,11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,25]),
			identify_pythia_truth(),
			decay_fermions_as_taus(),
			select_emu_events(),
			build_events(),
			)

		self.add_result_function(
			plot_kinematics(),
			)

		self.add_meta_result_function(
			)

class identify_pythia_truth(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += ['truth']

		self.create_branches.update(dict((particle_type+'_'+particle_value,'float') for particle_type,particle_value in product(['b1','b2','b3','b4','A','f1','f2'],['pt','eta','phi','m'])))

	def __call__(self,event):

		event.A = [p for p in event.truth.values() if abs(p().pdgId)==25 and p().status==3][0]

		try: event.f1,event.f2 = [p() for p in event.truth.values() if abs(p().pdgId) in [5,15] and p in event.A.children and p().status==3]
		except ValueError:
			print [p().pdgId for p in event.truth.values() if abs(p().pdgId) in [5,15] and p in event.A.children and p().status==3]
			raise
		event.A = event.A()

		bs = [p for p in event.truth.values() if abs(p().pdgId)==5 and p().status==2 and not p.parents]

		if len(bs)>4:
			event.b1,event.b2,event.b3,event.b4 = tuple([b() for b in sorted(bs,key=lambda b: b()().Pt(),reverse=True)[:4]])

		else:
			event.b1,event.b2,event.b3,event.b4 = tuple([b() for b in sorted(bs,key=lambda b: b()().Pt(),reverse=True)]+[
				particle(\
				pt = -10.,
				eta = -10.,
				phi = -10.,
				m = -10.,
				charge = -10.
        	                ) for i in range(4-len(bs))])

		for item,name in [\
			(event.b1,'b1_'),
			(event.b2,'b2_'),
			(event.b3,'b3_'),
			(event.b4,'b4_'),
			(event.A,'A_'),
			(event.f1,'f1_'),
			(event.f2,'f2_'),
			]:

			event.__dict__[name+'pt'] = item().Pt()
			event.__dict__[name+'eta'] = item().Eta()
			event.__dict__[name+'phi'] = item().Phi()
			event.__dict__[name+'m'] = item().M()

	def __init__(self):
		event_function.__init__(self)



class identify_sherpa_truth(event_function):

	class invalid_configuration(EventBreak): pass
	class select_emu(EventBreak): pass

	def __init__(self):
		event_function.__init__(self)

		self.break_exceptions += [
			identify_sherpa_truth.invalid_configuration,
			identify_sherpa_truth.select_emu,
			]

		self.required_branches += ['truth']

		self.create_branches.update(dict((particle_type+'_'+particle_value,'float') for particle_type,particle_value in product(
			['b1','b2','b3','b4','A','tau1','tau2','l1','l2'],
			['pt','eta','phi','m'],
			)))

	def __call__(self,event):

		bs = [p for p in event.truth.values() if abs(p().pdgId)==5 and p().status==11]

		taus = [p for p in event.truth.values() if abs(p().pdgId)==15 and p().status==2 and 15 in [parent().pdgId for parent in p.parents]]

		if len(taus)!=2: raise identify_sherpa_truth.invalid_configuration()

		if len(bs)>4:
			b1,b2,b3,b4 = tuple([b() for b in sorted(bs,key=lambda b: b()().Pt(),reverse=True)[:4]])

		else:
			b1,b2,b3,b4 = tuple([b() for b in sorted(bs,key=lambda b: b()().Pt(),reverse=True)]+[
				particle(\
					pt = -10.,
					eta = -10.,
					phi = -10.,
					m = -10.,
					charge = -10.
					) for i in range(4-len(bs))])

		tau1,tau2 = taus

		try:
			e = [item for item in tau1.children+tau2.children if abs(item().pdgId) == 11][0]
			nu_e = [item for item in tau1.children+tau2.children if abs(item().pdgId) == 12][0]
			mu = [item for item in tau1.children+tau2.children if abs(item().pdgId) == 13][0]
			nu_mu = [item for item in tau1.children+tau2.children if abs(item().pdgId) == 14][0]
			nu_tau1,nu_tau2 = [item for item in tau1.children+tau2.children if abs(item().pdgId) == 16][0:2]

		except IndexError: raise identify_sherpa_truth.select_emu()

		event.b1 = b1
		event.b2 = b2
		event.b3 = b3
		event.b4 = b4

		if mu in tau1.children: tau1,tau2 = tau2,tau1

		event.tau1 = tau1()
		event.tau2 = tau2()

		event.l2 = mu()
		event.l1 = e()
		A = tau1()()+tau2()()
		event.A = particle(); event.A.set_particle(A)

		for item,name in [\
			(event.b1,'b1_'),
			(event.b2,'b2_'),
			(event.b3,'b3_'),
			(event.b4,'b4_'),
			(event.A,'A_'),
			(event.tau1,'tau1_'),
			(event.tau2,'tau2_'),
			(event.l1,'l1_'),
			(event.l2,'l2_'),
			]:

			event.__dict__[name+'pt'] = item().Pt()
			event.__dict__[name+'eta'] = item().Eta()
			event.__dict__[name+'phi'] = item().Phi()
			event.__dict__[name+'m'] = item().M()

		return

class collect_truth_jets(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.collection_name = 'jet_antikt4truth_'

		self.names = [
			'pt',
			'eta',
			'phi',
			'E',
			]

		self.required_branches += [self.collection_name + name for name in self.names]
		self.required_branches += [self.collection_name+'n']
	def __call__(self,event):

		event.truth_jets = {}

		for jet_n in range(event.__dict__[self.collection_name+'n']):

			jet = particle(\
				**dict((name,event.__dict__[self.collection_name+name][jet_n]) for name in self.names)
				)
			jet.set_pt_eta_phi_e(jet.pt,jet.eta,jet.phi,jet.E)
			event.truth_jets[jet_n] = jet

class match_truth_jets(event_function):
	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		event.matched_truth_jets = []

		for jet in event.truth_jets.values():
			jet.matched = False
			for b in [
				event.b1,
				event.b2,
				event.b3,
				event.b4,
				]:
				if b().DeltaR(jet())<0.2: jet.matched = True

		for jet in event.truth_jets.values():
			if not jet.matched: continue
			event.matched_truth_jets.append(jet)

		event.matched_truth_jets.sort(key=lambda b: b().Pt(),reverse=True)


class select_emu_events(event_function):

	class leptons(EventBreak): pass
	class jets(EventBreak): pass

    """@commandline(
        'select_emu_events',
        pt_low = arg('--pt_low',type=float,help='Minimum Phi pT'),
        pt_high = arg('--pt_high',type=float,help='Maximum Phi pT'),
        )"""
	def __init__(self,pt_low=0., pt_high=20000000.):
		event_function.__init__(self)

		self.pt_low = pt_low
		self.pt_high = pt_high

		self.break_exceptions += [
			select_emu_events.leptons,
			select_emu_events.jets,
			]

	def __call__(self,event):

		for requirement,exception in [
			(all([
				event.l1_pt>12000.,
				abs(event.l1_eta)<3.0,
				event.l2_pt>8000.,
				abs(event.l2_eta)<3.0,
				event.A_pt>self.pt_low,
				event.A_pt<self.pt_high,
				]),select_emu_events.leptons),
			(len([1 for jet in event.truth_jets.values() if jet.pt>15000. and abs(jet.eta)<3.0 and jet.matched])>0,select_emu_events.jets),
			]:
			if not requirement: raise exception()


class build_events(event_function):
	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		event.lepton_pair_mass = (event.l1()+event.l2()).M()
		event.lepton_pair_mass_low = event.lepton_pair_mass
		event.lepton_dR = event.l1().DeltaR(event.l2())
		event.lepton_dPhi = event.l1().DeltaPhi(event.l2())

		event.matched_truth_jets = [jet for jet in event.matched_truth_jets if jet.pt>15000. and abs(jet.eta)<3.0]

		if event.matched_truth_jets:
			event.j1_pt = event.matched_truth_jets[0]().Pt()
			event.j1_eta = event.matched_truth_jets[0]().Eta()
		else:
			event.j1_pt = 0.
			event.j1_eta = -5.

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low,xlabel)) for name,binning,high,low,xlabel in [
			('lepton_pair_mass',25,0.,150000.,"M(l_{1},l_{2}) [MeV]"),
			('lepton_pair_mass_low',22,0.,45000.,"M(l_{1},l_{2}) [MeV]"),
			('lepton_dR',15,0.,6.,"\DeltaR(l_{1}, l_{2})"),
			('lepton_dPhi',16,0.,3.2,"\Delta\phi(l_{1},l_{2})"),
			('b1_pt',30,0.,120000.,"p_{T}^{b_{1}} [MeV]"),
			('b1_eta',40,-5.,5.,"\eta^{b_{1}}"),
			('j1_pt',30,0.,120000.,"p_{T}^{j_{1}} [MeV]"),
			('j1_eta',40,-5.,5.,"\eta^{j_{1}}"),
			('l1_pt',20,0.,80000.,"p_{T}^{l_{1}} [MeV]"),
			('l1_eta',24,-3.,3.,"\eta^{l_{1}}"),
			('l1_phi',32,-3.2,3.2,"\phi^{l_{1}}"),
			('l2_pt',15,0.,60000.,"p_{T}^{l_{2}} [MeV]"),
			('l2_eta',24,-3.,3.,"\eta^{l_{2}}"),
			('l2_phi',32,-3.2,3.2,"\phi^{l_{2}}"),
			('tau1_pt',20,0.,80000.,"p_{T}^{l_{1}} [MeV]"),
			('tau1_eta',24,-3.,3.,"\eta^{l_{1}}"),
			('tau1_phi',32,-3.2,3.2,"\phi^{l_{1}}"),
			('tau2_pt',15,0.,60000.,"p_{T}^{l_{2}} [MeV]"),
			('tau2_eta',24,-3.,3.,"\eta^{l_{2}}"),
			('tau2_phi',32,-3.2,3.2,"\phi^{l_{2}}"),
			])

		self.names_2d = [
			]

		for name_,(binning,high,low,xlabel) in self.names.items():
			name = name_
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

		weight = event.__weight__

		for name_ in self.names:
			name = name_
			self.results[name].Fill(event.__dict__[name_],weight)

		for name1,name2 in self.names_2d:
			name = '{0}_{1}'.format(name1,name2)
			self.results[name].Fill(event.__dict__[name1],event.__dict__[name2],weight)
