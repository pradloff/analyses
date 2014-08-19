from common.analysis import analysis
from common.functions import event_function,result_function,arg,EventBreak

import ROOT

class example1(analysis):
	def __init__(self):
		analysis.__init__(self)
		
		self.add_event_function(
			cut_on_l1(),
			)

		self.add_result_function(
			plot_l1_pt(),
			)

		self.add_meta_result_function(
			)
			

class cut_on_l1(event_function):

	class l1_pt_cut(EventBreak): pass

	def __init__(
		self,
		):

		event_function.__init__(self)
		
		self.break_exceptions += [
			cut_on_l1.l1_pt_cut
			]

		self.required_branches += [
			'l1_pt',
			]

	def __call__(self,event):
		if event.l1_pt < 15000.: raise cut_on_l1.l1_pt_cut()
	
class plot_l1_pt(result_function):
	def __init__(self):
		result_function.__init__(self)

		self.names = dict((name,(binning,high,low,xlabel)) for name,binning,high,low,xlabel in [
			('l1_pt',14,0.,70000.,"p_{T}^{l_{1}} [MeV]"),
			])

		for name,(binning,high,low,xlabel) in self.names.items():
			self.results[name] = ROOT.TH1F(name,name,binning,high,low)
			self.results[name].Sumw2()
			self.results[name].GetXaxis().SetTitle(xlabel)
			self.results[name].GetYaxis().SetTitle('Events')
			self.results[name].GetYaxis().CenterTitle()
		
	def __call__(self,event):
		if event.__break__: return

		for name in self.names:
			self.results[name].Fill(event.__dict__[name],event.__weight__)

