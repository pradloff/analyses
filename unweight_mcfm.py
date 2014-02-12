import random

from common.analysis import analysis
from common.functions import event_function

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

