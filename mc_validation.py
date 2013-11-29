from common.analysis import analysis
from common.EventFunction import EventFunction, ResultFunction
from functions.eventFunctions.EventFunctions import *
from itertools import product
from common.Particle import Particle
from math import sqrt,cos,sin
import ROOT
import sys

class skim_truth(analysis):
	def __init__(self,*args,**kwargs):
		analysis.__init__(self,*args,**kwargs)
		
		self.AddEventFunction(
			truthTree(),
			identify(),
			collinear(),
			)

		self.AddResultFunction(
			)

		self.AddMetaResultFunction(
			)

		for item,itemType in [
			]: self.AddItem(item,itemType)

		for items,itemType in [		
			]: self.AddItems(itemType,*items)
			
#--------------------------------------------------------------------------------------------------------------

from common.node import node
from common.Particle import Particle

class truthTree(EventFunction):

	def __init__(self):
		EventFunction.__init__(self)
		self.addItems()

	def addItems(self):
		self.addItem('mc_n',0)
		self.addItem('mc_parent_index',0)
		self.addItem('mc_pt',0)
		self.addItem('mc_eta',0)
		self.addItem('mc_phi',0)
		self.addItem('mc_m',0)
		self.addItem('mc_charge',0)
		self.addItem('mc_pdgId',0)
		self.addItem('mc_status',0)	

	def __call__(self,event):
		event.truth = buildTruthTree(\
			event.mc_pt,
			event.mc_eta,
			event.mc_phi,
			event.mc_m,
			event.mc_charge,
			event.mc_pdgId,
			event.mc_parent_index,
			event.mc_status,
			)
		return

#--------------------------------------------------------------------------------------------------------------

def buildTruthTree(\
	mc_pt,
	mc_eta,
	mc_phi,
	mc_m,
	mc_charge,
	mc_pdgId,
	mc_parent_index,
	mc_status,
	):
	truth = {}
	for particle in range(len(mc_pt)): 
		truth[particle] = node(\
			[],
			[],
			Particle(\
				pt = mc_pt[particle],
				eta = mc_eta[particle],
				phi = mc_phi[particle],
				m = mc_m[particle],
				charge = mc_charge[particle],
				pdgId = mc_pdgId[particle],
				status = mc_status[particle],
				))
	for index,particle in truth.items(): particle.addParents(*[truth[parent] for parent in mc_parent_index[index]])
	return truth
	
#--------------------------------------------------------------------------------------------------------------

class identify(EventFunction):

	def __init__(self):
		EventFunction.__init__(self)
		self.addItems()

	def addItems(self):
		for item in [\
			'b1_',
			'b2_',
			'b3_',
			'b4_',
			'A_',
			'mu_',
			'e_',
			]: 
			for sub in [\
				'pt',
				'eta',
				'phi',
				'm',
				]: self.addItem(item+sub,2,type_='float')
		self.addItem('A_m_visible',2,type_='float')
		self.addItem('MET',2,type_='float')

	def __call__(self,event):

		bs = [particle for particle in event.truth.values() if abs(particle.__item__.pdgId)==5 and particle.__item__.status==11]
		#bs = [particle for particle in bs if [abs(c.__item__.pdgId) for c in particle.__children__].count(15)==2]

		#if len(bs)!=2: 
		#	event.__break__==True
		#	return

		#taus = list(set(c for c in bs[0].__children__ if abs(c.__item__.pdgId)==15)&set(c for c in bs[1].__children__ if abs(c.__item__.pdgId)==15))
		#taus = list(set(c for c in taus[0].__children__ if abs(c.__item__.pdgId)==15)&set(c for c in taus[1].__children__ if abs(c.__item__.pdgId)==15))
		#taus = list(set(c for c in taus[0].__children__ if abs(c.__item__.pdgId)==15)&set(c for c in taus[1].__children__ if abs(c.__item__.pdgId)==15))

		taus = [particle for particle in event.truth.values() if abs(particle.__item__.pdgId)==15 and particle.__item__.status==2 and 15 in [p.__item__.pdgId for p in particle.__parents__]]

		if len(taus)!=2:
			event.__break__ = True
			return
		if len(bs)>4:
			b1,b2,b3,b4 = tuple([b.__item__ for b in sorted(bs,key=lambda b: b.__item__.__particle__.Pt(),reverse=True)[:4]])

		else:
			b1,b2,b3,b4 = tuple([b.__item__ for b in sorted(bs,key=lambda b: b.__item__.__particle__.Pt(),reverse=True)]+[
				Particle(\
				pt = -10.,
				eta = -10.,
				phi = -10.,
				m = -10.,
				charge = -10.
        	                ) for i in range(4-len(bs))])

		tau1,tau2 = taus

		try: 
			e = [item for item in tau1.__children__+tau2.__children__ if abs(item.__item__.pdgId) == 11][0]
			nu_e = [item for item in tau1.__children__+tau2.__children__ if abs(item.__item__.pdgId) == 12][0]
			mu = [item for item in tau1.__children__+tau2.__children__ if abs(item.__item__.pdgId) == 13][0]
			nu_mu = [item for item in tau1.__children__+tau2.__children__ if abs(item.__item__.pdgId) == 14][0]
			nu_tau1,nu_tau2 = [item for item in tau1.__children__+tau2.__children__ if abs(item.__item__.pdgId) == 16][0:2]

		except IndexError:
			event.__break__ = True
			return
			#raise
		
		event.b1 = b1
		event.b2 = b2
		event.b3 = b3
		event.b4 = b4

		event.mu = mu.__item__
		event.e = e.__item__
		event.A = tau1.__item__.__particle__+tau2.__item__.__particle__

		event.A_m_visible = (event.mu.__particle__+event.e.__particle__).M()

		METLV = (
			nu_tau1.__item__.__particle__+
			nu_tau2.__item__.__particle__+
			nu_mu.__item__.__particle__+
			nu_e.__item__.__particle__
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
			
			event.__dict__[name+'pt'] = item.Pt()
			event.__dict__[name+'eta'] = item.Eta()
			event.__dict__[name+'phi'] = item.Phi()
			event.__dict__[name+'m'] = item.M()

		return
			
#--------------------------------------------------------------------------------------------------------------

class collinear(EventFunction):

	def __init__(self):
		EventFunction.__init__(self)
		self.addItems()
		
	def addItems(self):
		for item,type_ in [
			('A_m_collinear','float')
			]: self.addItem(item,2,type_=type_)
		
	def __call__(self,event):
		muon = event.mu
		electron = event.e
		metx = event.METx
		mety = event.METy
		
		x1 = ((muon.Px()*electron.Py())-(muon.Py()*electron.Px())) / ((muon.Px()*electron.Py())-(muon.Py()*electron.Px())+(electron.Py()*metx)-(electron.Px()*mety))
		x2 = ((muon.Px()*electron.Py())-(muon.Py()*electron.Px())) / ((muon.Px()*electron.Py())-(muon.Py()*electron.Px())+(muon.Px()*mety)-(muon.Py()*metx))
		if x1*x2 > 0.: event.A_m_collinear = (muon.__particle__+electron.__particle__).M()/sqrt(x1*x2)
		else: event.A_m_collinear = 0.
