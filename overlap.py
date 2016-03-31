from common.functions import event_function

def overlapped(imposter,collection,dr):
	for particle in collection:
		if particle().DeltaR(imposter())<dr:
			return True
	return False


class remove_overlap(event_function):

	def __init__(self,dR=0.2):
		super(remove_overlap,self).__init__()
		self.dR = dR
		event_function.__init__(self)

	def __call__(self,event):
		super(remove_overlap,self).__call__(event)
		#Remove jets from taus,electrons and muons
		for jet in event.jets.values():
			jet.overlap_removed = any([overlapped(jet,collection,self.dR) for collection in [[particle for particle in collection_ if particle.passed_preselection] for collection_ in [
				event.taus.values(),
				event.electrons.values(),
				event.muons.values(),
				]]])

		#Remove taus from electrons and muons
		for tau in event.taus.values():
			tau.overlap_removed = any([overlapped(tau,collection,self.dR) for collection in [[particle for particle in collection_ if particle.passed_preselection_taus] for collection_ in [
				event.electrons.values(),
				event.muons.values(),
				]]])

		#Remove electrons from muons
		for electron in event.electrons.values():
			electron.overlap_removed = overlapped(electron,[particle for particle in event.muons.values() if particle.passed_preselection],self.dR)

		for muon in event.muons.values(): muon.overlap_removed = False
