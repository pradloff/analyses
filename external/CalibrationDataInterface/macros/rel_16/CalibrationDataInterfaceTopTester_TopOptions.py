#-----------------------------------------------------------------------------
# Athena imports
#-----------------------------------------------------------------------------

#from AthenaCommon.Constants import *
from AthenaCommon.AppMgr import theApp
from AthenaCommon.AppMgr import ToolSvc
from AthenaCommon.AppMgr import ServiceMgr
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AlgSequence import AlgSequence
theJob = AlgSequence()


#-----------------------------------------------------------------------------
# Message Service
#-----------------------------------------------------------------------------
# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
ServiceMgr.MessageSvc.OutputLevel = WARNING
import AthenaServices
AthenaServices.AthenaServicesConf.AthenaEventLoopMgr.OutputLevel = WARNING

# the POOL converters
#include( "ParticleBuilderOptions/AOD_PoolCnv_jobOptions.py" )
#include( "ParticleBuilderOptions/McAOD_PoolCnv_jobOptions.py" )

#-----------------------------------------------------------------------------
# Input Datasets
#-----------------------------------------------------------------------------
ServiceMgr.EventSelector.InputCollections = ['/data/atlas3/users/filthaut/mc09_7TeV.106201.TTbar_McAtNlo_Jimmy_170GeV.merge.AOD.e522_s765_s767_r1302_r1306_tid139662_00/AOD.139662._000012.pool.root.1']
theApp.EvtMax = 10 # -1 means all events

#-----------------------------------------------------------------------------
# Algorithms
#-----------------------------------------------------------------------------

if not 'BTaggingFlags' in dir():
    from BTagging.BTaggingFlags import BTaggingFlags
BTaggingFlags.OutputLevel = VERBOSE
BTaggingFlags.Jets = ['AntiKt4H1TowerAODJets']
BTaggingFlags.CalibrationFromLocalReplica = True
BTaggingFlags.CalibrationTag = 'MyCalibV1'

include('CalibrationDataInterface/BTagPerformance_jobOptions.py')

from CalibrationDataInterface.CalibrationDataInterfaceConf import Analysis__CalibrationDataInterfaceTester as Tester

theJob += Tester(name = 'PerfCalibTester',
                 JetCollection = 'AntiKt4H1TowerAODJets',
                 Tagger = 'SV0',
                 OperatingPoint = '5_72',
                 CalibrationInterface = JetTagPerformanceSV0,
                 CalibrationUncertainty = 'Syst',
                 OutputLevel = VERBOSE)

print theJob
