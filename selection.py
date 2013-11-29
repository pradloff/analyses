from common.functions import event_function

class trigger(event_function):
	def __init__(self):
		event_function.__init__(self)
		self.create_branches['EF_e24vhi_medium1'] = 'bool'
	
	def __call__(self,event):
		if 'EF_e24vhi_medium1' not in event: event.EF_e24vhi_medium1 = False

class preselection(event_function):

	def __init__(self):
		event_function.__init__(self)
		self.required_branches += [
			'muons',
			'electrons',
			]

	def __call__(self,event):
		if sum(1 for muon in event.muons.values() if muon.passed_selection):
			event.__break__ = True
			return

		if sum(1 for electron in event.electrons.values() if electron.passed_selection) != 2:
			event.__break__ = True
			return

		electron1,electron2 = sorted([electron_ for electron_ in event.electrons.values() if electron_.passed_selection],key=lambda electron: electron.pt_corrected,reverse=True)

		if not electron1.pt_corrected>26000.:
			event.__break__ = True
			return

		if not event.EF_e24vhi_medium1:
			event.__break__ = True
			return
