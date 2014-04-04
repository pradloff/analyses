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
from unweight_mcfm import decay_fermions_as_taus
from mc import identify as identify_sherpa_truth
import ROOT
import os
from math import sqrt
import json
from itertools import product

class truth_analysis_sherpa(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(),
			identify_sherpa_truth(),
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
			truth_tree(),
			identify_pythia_truth(),
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

		event.A = [p for p in event.truth.values() if abs(p().pdgId)==25 and p().status==2][0]

		event.f1,event.f2 = [p() for p in event.truth.values() if abs(p().pdgId) in [5,15] and p in event.A.children and p().status==2]
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


class identify_sherpa_truth(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += ['truth']

		self.create_branches.update(dict((particle_type+'_'+particle_value,'float') for particle_type,particle_value in product(['b1','b2','b3','b4','A','l1','l2'],['pt','eta','phi','m'])))

	def __call__(self,event):

		bs = [p for p in event.truth.values() if abs(p().pdgId)==5 and p().status==11]

		taus = [p for p in event.truth.values() if abs(p().pdgId)==15 and p().status==2 and 15 in [parent().pdgId for parent in p.parents]]

		if len(taus)!=2:
			event.__break__ = True
			return
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

		except IndexError:
			event.__break__ = True
			return
		
		event.b1 = b1
		event.b2 = b2
		event.b3 = b3
		event.b4 = b4

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
			(event.l1,'l1_'),
			(event.l2,'l2_'),
			]:
			
			event.__dict__[name+'pt'] = item().Pt()
			event.__dict__[name+'eta'] = item().Eta()
			event.__dict__[name+'phi'] = item().Phi()
			event.__dict__[name+'m'] = item().M()

		return

class select_emu_events(event_function):
	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		if not all([
			event.l1_pt>15000.,
			abs(event.l1_eta)<2.5,
			event.l2_pt>15000.,
			abs(event.l2_eta)<2.5,
			]):
			event.__break__= True
			return
	
class build_events(event_function):
	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):
		event.lepton_pair_mass = (event.l1()+event.l2()).M()
		event.lepton_pair_mass_low = event.lepton_pair_mass
		event.lepton_dR = event.l1().DeltaR(event.l2())
		event.lepton_dPhi = event.l1().DeltaPhi(event.l2())

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low,xlabel)) for name,binning,high,low,xlabel in [
			('lepton_pair_mass',25,0.,150000.,"M(l_{1},l_{2}) [MeV]"),
			('lepton_pair_mass_low',22,0.,45000.,"M(l_{1},l_{2}) [MeV]"),
			('lepton_dR',15,0.,6.,"\DeltaR(l_{1}, l_{2})"),
			('lepton_dPhi',16,0.,3.2,"\Delta\phi(l_{1},l_{2})"),
			('b1_pT',30,0.,120000.,"p_{T}^{j_{1}} [MeV]"),
			('b1_eta',40,-5.,5.,"\eta^{j_{1}}"),
			('l1_pt',20,0.,80000.,"p_{T}^{l_{1}} [MeV]"),
			('l1_eta',24,-3.,3.,"\eta^{l_{1}}"),
			('l1_phi',32,-3.2,3.2,"\phi^{l_{1}}"),
			('l2_pt',15,0.,60000.,"p_{T}^{l_{2}} [MeV]"),
			('l2_eta',24,-3.,3.,"\eta^{l_{2}}"),
			('l2_phi',32,-3.2,3.2,"\phi^{l_{2}}"),
			])

		self.names_2d = [
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
