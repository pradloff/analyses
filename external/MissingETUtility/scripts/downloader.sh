#!/bin/bash
cd $TestArea/Reconstruction/MissingETUtility/run

svn co $SVNROOT/Reconstruction/Jet/JetUncertainties/tags/JetUncertainties-00-03-05-01 JetUncertainties
svn co $SVNROOT/Reconstruction/Jet/JetResolution/tags/JetResolution-01-00-00 JetResolution
svn co $SVNROOT/Reconstruction/egamma/egammaAnalysis/egammaAnalysisUtils/tags/egammaAnalysisUtils-00-02-76 egammaAnalysisUtils
svn co $SVNROOT/PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections/tags/MuonMomentumCorrections-00-05-03 MuonMomentumCorrections
