from common.functions import event_function
from common.particle import particle
from math import cosh,sqrt
from common.external import load
import ROOT
import os
from misc import list_attributes
from common.commandline import commandline,arg
from common.branches import auto_branch,branch
import code
class collect_taus(event_function):

	def __init__(self,collection_name='tau_'):
		super(collect_taus,self).__init__()

		self.collection_name = collection_name

		self.names = [
			'author',
			'JetBDTSigMedium',
			'eta',
			'phi',
			'track_n',
			'numTrack',
			#'track_eta',
			'Et',
			'charge',
			'EleBDTMedium',
			'muonVeto',
			]

		for name in self.names:
			self.branches.append(branch(self.collection_name+name,'r'))

		self.branches.append(branch(self.collection_name+'n','r'))

	def __call__(self,event):
		super(collect_taus,self).__call__(event)

		#Collect electrons
		event.taus = {}
		for tau in range(event.__dict__[self.collection_name+'n']):
			event.taus[tau] = particle(\
				**dict((name,event.__dict__[self.collection_name+name][tau]) for name in self.names)
				)

		for tau in event.taus.values():
			tau.set_pt_eta_phi_m(tau.Et,tau.eta,tau.phi,1777.)

		#Define selections
		for tau in event.taus.values():

			tau.passed_preselection = all([
				tau.Et>15000.,
				tau.author==1 or tau.author==3,
				tau.JetBDTSigMedium==1,
				abs(tau.eta)<2.47 and tau.track_n>0,# and abs(tau.track_eta[0])<2.47,
				tau.numTrack==1 or tau.numTrack==3,
				abs(tau.charge-1.)<0.1,
				not (tau.numTrack==1 and tau.EleBDTMedium==1),
				tau.muonVeto==0
				])

			tau.passed_selection = tau.passed_preselection

		return
