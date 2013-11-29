theApp.EvtMax = 1000

# Necessary and sufficient to read Athena pool files!
import AthenaPoolCnvSvc.ReadAthenaPool
import os
sourcedir = 'sample_AOD/mc11_7TeV.105200.T1_McAtNlo_Jimmy.merge.AOD.e835_s1272_s1274_r2730_r2700_tid517944_00/'
filelist = os.listdir(sourcedir)
for i in filelist:
    svcMgr.EventSelector.InputCollections.append(sourcedir+i)
#svcMgr.EventSelector.InputCollections = ['AOD.pool.root']
svcMgr.EventSelector.SkipEvents = 0

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from GaudiSvc.GaudiSvcConf import THistSvc
svcMgr += THistSvc()

from MissingETUtility.MissingETUtilityConf import MissingETUtility__METUtilAlg as METUtilAlg
from MissingETUtility.MissingETUtilityConf import MissingETUtility__METUtilityAthTool as METUtilityAthTool

svcMgr.THistSvc.Output += ["outfile DATAFILE='metUtilTest.root' OPT='RECREATE'"]
svcMgr.THistSvc.PrintAll = True

# Set up the algorithm
alg = METUtilAlg( 'METUtilTest' )
alg.OutputLevel = DEBUG

# Setup for the MET Utility tool
ToolSvc = Service( 'ToolSvc' )
ToolSvc += METUtilityAthTool( 'MyMETUtil' )
ToolSvc.MyMETUtil.OutputLevel = DEBUG

# Assign the tool to the algorithm
alg.MetUtil = ToolSvc.MyMETUtil

# Add the algorithm to the sequence
job += alg

THistSvc.Dump = True
