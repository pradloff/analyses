from common.analysis import analysis

from misc import count_primary_vertices
from muons import collect_muons
from electrons import collect_electrons
from selection import trigger, preselection
from pileup import pileup_weighting

from metadata import lumi

class test(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			count_primary_vertices(),
			pileup_weighting(),
			collect_muons(),
			collect_electrons(),
			collect_jets(),
			trigger(),
			selection(),
			)

		self.add_result_function(
			)

		self.add_meta_result_function(
			lumi()
			)

class trigger(event_function):
	
	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'EF_mu18_tight_mu8_EFFS',
			'EF_e12Tvh_medium1_mu8',
			'EF_e24vhi_medium1',
			'EF_2e12Tvh_loose1',
			]

		self.initialize_tools()

	def __call__(self,event):
	
		if not any([
			event.EF_mu18_tight_mu8_EFFS,
			event.EF_e12Tvh_medium1_mu8,
			event.EF_e24vhi_medium1,
			event.EF_2e12Tvh_loose1,			
			]):
			:
			event.__break__ = True
			return

		self.apply_corrections(event)

	def apply_corrections(self,event):
		pass

	def initialize_tools(self):
		pass

class selection(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'electrons',
			'muons',
			'jets',
			]

	def __call__(self,event):
		if any([
			sum(1 for lepton in event.electrons.values()+event.muons.values() if lepton.passed_preselection)<2,
			sum(1 for jet in event.jets.values() if jet.passed_preselection),
			]):
			event.__break__=True
			return

		
