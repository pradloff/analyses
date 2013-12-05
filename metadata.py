import ROOT
from common.functions import meta_result_function
import os
import time

class lumi(meta_result_function):

	def __init__(self):
		meta_result_function.__init__(self)
		#ROOT.gDirectory.pwd()
		self.tree = ROOT.TTree('lumi','lumi')
		#ROOT.SetOwnership(self.tree,False)
		self.results['lumi'] = self.tree
		self.string = ROOT.TString()
		self.tree.Branch('grl','TString',ROOT.AddressOf(self.string),1600,0)
		
	def __call__(self,files):
		for f in files:
			try:
				f = ROOT.TFile.Open(f)
				if not f: raise OSError('File {0} not found or could not be opened'.format(f))
				#Look for Lumi TDirectory(ies)
				lumi_directories = [f.Get(key.GetName()+';'+str(key.GetCycle())) for key in f.GetListOfKeys() if isinstance(f.Get(key.GetName()+';'+str(key.GetCycle())),ROOT.TDirectory) and key.GetName()=='Lumi']
				for lumi_directory in lumi_directories:
					for key in lumi_directory.GetListOfKeys():
						self.string.Resize(0)
						self.string.Append(lumi_directory.Get(key.GetName()+';'+str(key.GetCycle())).String())
						self.tree.Fill()
				#Look for lumi trees
				lumi_trees = [f.Get(key.GetName()+';'+str(key.GetCycle())) for key in f.GetListOfKeys() if isinstance(f.Get(key.GetName()+';'+str(key.GetCycle())),ROOT.TTree) and key.GetName()=='lumi']
				for lumi_tree in lumi_trees:
					if 'grl' not in [b.GetName() for b in lumi_tree.GetListOfBranches()]: continue
					for entry in range(lumi_tree.GetEntries()):
						lumi_tree.GetEntry(entry)
						self.string.Resize(0)
						self.string.Append(lumi_tree.grl)
						self.tree.Fill()
			except Exception as error:
				print 'Problem getting lumi data in file'.format(f)
				raise error
