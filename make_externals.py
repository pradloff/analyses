from common.external import load

if __name__=='__main__':

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
		load(package,verbose=True,clean=True)
