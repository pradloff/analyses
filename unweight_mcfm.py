import random
from tauola import tauola_
from common.analysis import analysis
from common.functions import event_function,result_function
from common.particle import particle
import ROOT
from math import sqrt

class unweight_mcfm(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			get_weight(),
			)

		self.add_result_function(
			)

		self.add_meta_result_function(
			)

class mcfm(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			build_events(),
			decay_fermions_as_taus(),
			compute_kinematics(),
			)

		self.add_result_function(
			plot_kinematics(),
			)

		self.add_meta_result_function(
			)


class get_weight(event_function):
	def __init__(self):
		event_function.__init__(self)
		self.required_branches += [
			'wt_ALL',
			]
		self.maximum = 0.02432399056851864

	def __call__(self,event):
		if event.wt_ALL/self.maximum < random.random():
			event.__break__ = True		

class build_events(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.lepton_names = [
			'px',
			'py',
			'pz',
			'E_',
			]

		for name in self.lepton_names:
			for f in [3,4,5]:
				self.required_branches.append(name+str(f))

	def __call__(self,event):

		f1 = particle(
			px=event.px3*1000.,
			py=event.py3*1000.,
			pz=event.pz3*1000.,
			E=event.E_3*1000.
			)
		f1.set_px_py_pz_e(f1.px,f1.py,f1.pz,f1.E)

		f2 = particle(
			px=event.px4*1000.,
			py=event.py4*1000.,
			pz=event.pz4*1000.,
			E=event.E_4*1000.,
			)
		f2.set_px_py_pz_e(f2.px,f2.py,f2.pz,f2.E)

		b = particle(
			px=event.px5*1000.,
			py=event.py5*1000.,
			pz=event.pz5*1000.,
			E=event.E_5*1000.,
			)
		b.set_px_py_pz_e(b.px,b.py,b.pz,b.E)

		event.f1 = f1
		event.f2 = f2
		event.b = b


class decay_fermions_as_taus(event_function):
	def __init__(self):
		event_function.__init__(self)

		self.tau_mass = 1776.82
		self.tauola = tauola_()

	def __call__(self,event):

		if random.getrandbits(1): event.f1,event.f2 = event.f2,event.f1 #flip e<->mu decay
	
		tauola_call = []

		mother = event.f1()+event.f2()
		boost = mother.BoostVector()
		for fermion in [event.f1(),event.f2()]:
			fermion.Boost(-boost)
			try: scale = sqrt(fermion.E()**2.-self.tau_mass**2.)/fermion.P()
			except ValueError:
				print 'Failed rescaling fermion of energy {0} mass {1} momentum {2}'.format(fermion.E(),fermion.M(),fermion.P())
				event.__break__ = True
				return
			fermion.SetPxPyPzE(fermion.Px()*scale,fermion.Py()*scale,fermion.Pz()*scale,fermion.E())
			fermion.Boost(boost)
			
			tauola_call+=[fermion.Px()/1000.,fermion.Py()/1000.,fermion.Pz()/1000.] #GEV for tauola

		result = self.tauola.leptonic_decay(*tauola_call)

		event.l1 = particle()
		event.l2 = particle()

		event.l1.set_px_py_pz_e(*[energy*1000. for energy in result[:4]])
		event.l2.set_px_py_pz_e(*[energy*1000. for energy in result[4:]])

		event.l1_pt = event.l1().Pt()
		event.l1_eta = event.l1().Eta()
		event.l1_phi = event.l1().Phi()
		event.l1_E = event.l1().E()
		event.l2_pt = event.l2().Pt()
		event.l2_eta = event.l2().Eta()
		event.l2_phi = event.l2().Phi()
		event.l2_E = event.l2().E()

class compute_kinematics(event_function):

	def __init__(self):
		event_function.__init__(self)

	def __call__(self,event):

		higgs = event.f1()+event.f2()
		lepton_pair = event.l1()+event.l2()

		event.higgs_mass = higgs.M()
		event.higgs_pt = higgs.Pt()
		event.lepton_pair_mass = lepton_pair.M()
		event.lepton_pair_pT = lepton_pair.Pt()
		event.lepton_pair_dR = event.l1().DeltaR(event.l2())

		event.f1_eta = event.f1().Eta()
		event.f1_pt = event.f1().Pt()
		event.f2_eta = event.f2().Eta()
		event.f2_pt = event.f2().Pt()
		event.b_eta = event.b().Eta()
		event.b_pt = event.b().Pt()

class plot_kinematics(result_function):
	def __init__(self):
		result_function.__init__(self)
		self.names = dict((name,(binning,high,low)) for name,binning,high,low in [
			('higgs_mass',50,0.,150000.),
			('higgs_pt',50,0.,150000.),
			('lepton_pair_mass',50,0.,150000.),
			('lepton_pair_pT',100,0.,100000.),
			('lepton_pair_dR',100,0.,10.),
			('f1_pt',50,0.,50000.),
			('f1_eta',24,-3.,3.),
			('f2_pt',50,0.,50000.),
			('f2_eta',24,-3.,3.),
			('l1_pt',50,0.,50000.),
			('l1_eta',24,-3.,3.),
			('l2_pt',50,0.,50000.),
			('l2_eta',24,-3.,3.),
			('b_pt',50,0.,50000.),
			('b_eta',24,-3.,3.),
			])

		self.names_2d = [
			]

		for name,(binning,high,low) in self.names.items():
			self.results[name] = ROOT.TH1F(name,name,binning,high,low)
			self.results[name].Sumw2()

		for name1,name2 in self.names_2d:
			binning1,high1,low1 = self.names[name1]
			binning2,high2,low2 = self.names[name2]
			name = '{0}_{1}'.format(name1,name2)
			self.results[name] = ROOT.TH2F(name,name,binning1,high1,low1,binning2,high2,low2)
			self.results[name].Sumw2()

	def __call__(self,event):
		if event.__break__: return

		weight = event.__weight__

		for name in self.names:
			self.results[name].Fill(event.__dict__[name],weight)

		for name1,name2 in self.names_2d:
			name = '{0}_{1}'.format(name1,name2)
			self.results[name].Fill(event.__dict__[name1],event.__dict__[name2],weight)



