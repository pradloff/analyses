from common.external import load

if __name__=='__main__':
	ignores = []
	for package in [
		'ApplyJetCalibration',
		'ApplyJetResolutionSmearing',
		'CalibrationDataInterface',
		'ElectronEfficiencyCorrection',
		'egammaAnalysisUtils',
		'HSG4LepLepTriggerSF',
		'JetResolution',
		'JetUncertainties',
		'JVFUncertaintyTool',
		'MissingETUtility',
		'MuonIsolationCorrection',
		'MuonEfficiencyCorrections',
		'MuonMomentumCorrections',
		'PileupReweighting',
		'TileTripReader',
		'TrigMuonEfficiency',
		]:
		ignores = load(package,verbose=True,clean=True,ignores=ignores)
