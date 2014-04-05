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

class skim_truth(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			truth_tree(),
			identify(),
			collinear(),
			)

		self.add_result_function(
			)

		self.add_meta_result_function(
			)


#--------------------------------------------------------------------------------------------------------------

class truth_tree(event_function):

	def __init__(self,pdgIds=None):
		event_function.__init__(self)
		self.names = [
			'n',
			'pt',
			'eta',
			'phi',
			'm',
			'charge',
			'pdgId',
			'parent_index',
			'status',
			]
		self.pdgIds = pdgIds

		self.required_branches += ['mc_'+name for name in self.names]
		self.create_branches['truth'] = None

	def __call__(self,event):
		event.truth = build_truth_tree(*[event.__dict__['mc_'+name] for name in self.names],self.pdgIds)
		return

#--------------------------------------------------------------------------------------------------------------

def build_truth_tree(\
	mc_n,
	mc_pt,
	mc_eta,
	mc_phi,
	mc_m,
	mc_charge,
	mc_pdgId,
	mc_parent_index,
	mc_status,
	pdgIds = None,
	):
	truth = {}
	if pdgId is not None: indices = [n for n in range(mc_n) if mc_pdgId[n] in pdgIds] 
	else: indices = range(mc_n)
	for n in indices:
		truth[n] = node(\
			[],
			[],
			particle(\
				pt = mc_pt[n],
				eta = mc_eta[n],
				phi = mc_phi[n],
				m = mc_m[n],
				charge = mc_charge[n],
				pdgId = mc_pdgId[n],
				status = mc_status[n],
				))
	for index,p in truth.items(): p.add_parents(*[truth[parent] for parent in mc_parent_index[index] and parent in truth])
	return truth
	
#--------------------------------------------------------------------------------------------------------------

from itertools import product

class identify(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += ['truth']

		self.create_branches.update(dict((particle_type+'_'+particle_value,'float') for particle_type,particle_value in product(['b1','b2','b3','b4','A','mu','e'],['pt','eta','phi','m'])))

		self.create_branches['MET'] = 'float'
		self.create_branches['METx'] = 'float'
		self.create_branches['METy'] = 'float'
		self.create_branches['A_m_visible'] = 'float'

		self.create_branches['mu'] = None
		self.create_branches['e'] = None

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

		event.mu = mu()
		event.e = e()
		A = tau1()()+tau2()()
		event.A = particle(); event.A.set_particle(A)

		event.A_m_visible = (event.mu()+event.e()).M()

		METLV = (
			nu_tau1()()+
			nu_tau2()()+
			nu_mu()()+
			nu_e()()
			)

		event.MET = METLV.Pt()
		event.METx = METLV.Px()
		event.METy = METLV.Py()

		for item,name in [\
			(event.b1,'b1_'),
			(event.b2,'b2_'),
			(event.b3,'b3_'),
			(event.b4,'b4_'),
			(event.A,'A_'),
			(event.mu,'mu_'),
			(event.e,'e_'),
			]:
			
			event.__dict__[name+'pt'] = item().Pt()
			event.__dict__[name+'eta'] = item().Eta()
			event.__dict__[name+'phi'] = item().Phi()
			event.__dict__[name+'m'] = item().M()

		return
			
#--------------------------------------------------------------------------------------------------------------

def eta_theta(eta):
	return 2*atan(exp(-eta))

class collinear(event_function):

	def __init__(self):
		event_function.__init__(self)
		
		self.required_branches += [
			'mu',
			'e',
			'METx',
			'METy',
			'MET',
			]

		self.create_branches['A_m_collinear'] = 'float'
		self.create_branches['A_m_super_collinear'] = 'float'

	def __call__(self,event):
		muon = event.mu
		electron = event.e
		metx = event.METx
		mety = event.METy

		m_frac_1 = ((muon().Px()*electron().Py())-(muon().Py()*electron().Px())) / ((muon().Px()*electron().Py())-(muon().Py()*electron().Px())+(electron().Py()*metx)-(electron().Px()*mety))
		m_frac_2 = ((muon().Px()*electron().Py())-(muon().Py()*electron().Px())) / ((muon().Px()*electron().Py())-(muon().Py()*electron().Px())+(muon().Px()*mety)-(muon().Py()*metx))

		x1,y1,z1 = (cos(electron.phi)*sin(eta_theta(electron.eta)),sin(electron.phi)*sin(eta_theta(electron.eta)),cos(eta_theta(electron.eta)))
		x2,y2,z2 = (cos(muon.phi)*sin(eta_theta(muon.eta)),sin(muon.phi)*sin(eta_theta(muon.eta)),cos(eta_theta(muon.eta)))

		if m_frac_1*m_frac_2 > 0.: event.A_m_collinear = (muon()+electron()).M()/sqrt(m_frac_1*m_frac_2)
		else: event.A_m_collinear = 0.

		p1 = ROOT.TLorentzVector()
		p2 = ROOT.TLorentzVector()

		E1 = event.MET/sin(eta_theta(electron.eta))
		E2 = event.MET/sin(eta_theta(muon.eta))

		p1.SetPxPyPzE(E1*x1,E1*y1,E1*z1,E1)
		p2.SetPxPyPzE(E2*x2,E2*y2,E2*z2,E2)

		if (x1*y2-x2*y1)<0.1: event.A_m_super_collinear = max([(electron()+muon()+p1).M(),(electron()+muon()+p2).M()])
		else: event.A_m_super_collinear = 0.
