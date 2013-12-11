from common.analysis import analysis
from common.functions import event_function

from misc import count_primary_vertices
from pileup import pileup_weighting
from muons import collect_muons
from electrons import collect_electrons
from selection import trigger, preselection
from taus import collect_taus
from jets import collect_jets
from overlap import remove_overlap
from met import correct_missing_energy
from metadata import lumi

class make_preselection(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			count_primary_vertices(),
			pileup_weighting(),
			collect_muons(),
			collect_electrons(),
			collect_taus(),
			collect_jets(),
			remove_overlap(),
			correct_missing_energy(),
			#trigger(),
			#preselection(),
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
			'EF_mu24i_tight',
			'EF_mu36_tight',
			'EF_e12Tvh_medium1_mu8',
			'EF_e24vhi_medium1',
			'EF_e60_medium1',
			]

		self.create_branches['trigger_scale_factor'] = 'float'

		self.initialize_tools()

	def __call__(self,event):
	
		if event.lepton_class == 0:
			if not any([
				event.EF_e24vhi_medium1 and event.l1_pt>25000.,
				event.EF_e60_medium1 and event.l1_pt>65000.,
				]):
				event.__break__ = True
				return

		if event.lepton_class == 1:
			if not any([
				event.EF_mu24i_tight and event.l1_pt>25000.,
				event.EF_mu36_tight and event.l1_pt>40000.,
				]):
				event.__break__ = True
				return
		
		if event.lepton_class == 2:
			if not event.EF_e12Tvh_medium1_mu8:
				event.__break__ = True
				return

	def initialize_tools(self):
		pass

class preselection(event_function):

	def __init__(self):
		event_function.__init__(self)

		self.required_branches += [
			'electrons',
			'muons',
			'jets',
			]
		self.create_branches.update(dict((key,value) for key,value in [
			('l1_eta','float'),
			('l1_phi','float'),
			('l1_pt','float'),
			('l1_ptcone40','float'),
			('l1_etcone20','float'),
			('l1_passed_selection','bool'),
			('l1_scale_factor','float'),
			]))

		self.create_branches.update(dict((key,value) for key,value in [
			('l2_eta','float'),
			('l2_phi','float'),
			('l2_pt','float'),
			('l2_ptcone40','float'),
			('l2_etcone20','float'),
			('l2_passed_selection','bool'),
			('l2_scale_factor','float'),
			]))

		self.create_branches.update(dict((key,value) for key,value in [
			('lepton_class','int'),
			]))

	def __call__(self,event):

		#2 preselection leptons, no hadronic taus, at least one preselection jet
		if any([
			sum(1 for lepton in event.electrons.values()+event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed)!=2,
			sum(1 for tau in event.taus.values() if tau.passed_preselection and not tau.overlap_removed)>0,
			sum(1 for jet in event.jets.values() if jet.passed_preselection)<1,
			]):
			event.__break__=True
			return
		#ee
		if sum(1 for lepton in event.electrons.values() if lepton.passed_preselection and not lepton.overlap_removed)==2:
			l1,l2 = [lepton for lepton in event.electrons.values() if lepton.passed_preselection and not lepton.overlap_removed]
			if l2.pt>l1.pt: l1,l2 = l2,l1
			event.lepton_class = 0
		#mumu
		elif sum(1 for lepton in event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed)==2:
			l1,l2 = [lepton for lepton in event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed]
			if l2.pt>l1.pt: l1,l2 = l2,l1
			event.lepton_class = 1
		#emu
		else:
			l1,l2 = [lepton for lepton in event.electrons.values()+event.muons.values() if lepton.passed_preselection and not lepton.overlap_removed]
			event.lepton_class = 2

		event.l1_eta = l1.eta
		event.l1_phi = l1.eta
		event.l1_pt = l1.pt_corrected
		event.l1_E = l1.E_corrected
		event.l1_ptcone40 = l1.ptcone40
		event.l1_etcone20 = l1.etcone20_corrected
		event.l1_scale_factor = l1.scale_factor
		event.l1_scale_factor_error = l1.scale_factor_error

		event.l2_eta = l2.eta
		event.l2_phi = l2.eta
		event.l2_pt = l2.pt_corrected
		event.l2_E = l2.E_corrected
		event.l2_ptcone40 = l2.ptcone40
		event.l2_etcone20 = l2.etcone20_corrected
		event.l2_scale_factor = l2.scale_factor
		event.l2_scale_factor_error = l2.scale_factor_error



