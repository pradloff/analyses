from common.functions import event_function
from common.external import load
import os
import ROOT

class pileup_weighting(event_function):

	def __init__(self):
		event_function.__init__(self)
		
		self.required_branches += [
			'averageIntPerXing',
			'RunNumber',
			'EventNumber',
			]
		self.create_branches['random_RunNumber'] = 'int'
		self.create_branches['weight_pileup'] = 'float'
		self.create_branches['mc_channel_number'] = None
		self.initialize_tools()
		
	def __call__(self,event):
		if event.is_mc:
			self.pileup_reweighting_tool.SetRandomSeed(event.mc_channel_number+event.EventNumber)
			event.random_RunNumber = self.pileup_reweighting_tool.GetRandomRunNumber(event.RunNumber)

			event.weight_pileup = self.pileup_reweighting_tool.GetCombinedWeight(
				int(event.RunNumber),
				int(event.mc_channel_number),
				float(event.averageIntPerXing),
				)
		else:
			event.random_RunNumber = event.RunNumber
			event.weight_pileup = 1.

	def initialize_tools(self):
		analysis_home = os.getenv('ANALYSISHOME')
		load('PileupReweighting')
		self.pileup_reweighting_tool = ROOT.Root.TPileupReweighting('pileup_reweighting_tool')
		self.pileup_reweighting_tool.SetDataScaleFactors(1/1.11)
		#self.pileup_reweighting_tool.MergeMCRunNumbers(195847,195848)
		self.pileup_reweighting_tool.AddConfigFile("{0}/external/PileupReweighting/share/mc12ab_defaults.prw.root".format(analysis_home))
		#self.pileup_reweighting_tool.AddConfigFile("{0}/external/PileupReweighting/share/my_PileUpReweighting_mc12a.root".format(analysis_home))
		self.pileup_reweighting_tool.AddLumiCalcFile("{0}/external/PileupReweighting/share/ilumicalc_histograms_None_200842-215643.root".format(analysis_home))
		self.pileup_reweighting_tool.SetUnrepresentedDataAction(2)
		self.pileup_reweighting_tool.Initialize()

