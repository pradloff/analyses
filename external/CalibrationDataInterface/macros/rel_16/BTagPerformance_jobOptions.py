# This is to create a separate instance of the CalibrationBroker (we cannot use the
# same instance as the taggers themselves, as different root folders are used etc.).

# main settings:

PerfCalibrationSingleFolder = True
SingleFolderName = 'Cuts'
# analogue of BTaggingFlags.CalibrationTag
#CalibrationTag = 'MyCalibV1'
CalibrationTag = 'TopCalibV1'
# analogue of BTaggingFlags.CalibrationFromLocalReplica
CalibrationFromLocalReplica = True

# Setup aliases
# aliases = ['AntiKt4H1TowerAODJets->Antikt4H1Tower', 'AntiKt4H1Tower->Antikt4H1Tower', 'Kt4H1Tower->Antikt4H1Tower' ]
aliases = ['AntiKt4H1TowerAODJets->AntiKt4Topo', 'AntiKt4H1Tower->AntiKt4Topo', 'Kt4H1Tower->AntiKt4Topo' ]

PerfCalibrationRoot = '/GLOBAL/BTagPerfCalib/'
myfolders = [ 'SV0' ]
PerfCalibrationFolders = []
for folder in myfolders:
    PerfCalibrationFolders.append(PerfCalibrationRoot+folder)

from JetTagCalibration.JetTagCalibrationConf import Analysis__CalibrationBroker
JetTagPerformanceCalibrationBroker = Analysis__CalibrationBroker(
    name = "PerformanceBroker",
    folderRoot = PerfCalibrationRoot,
    folders = PerfCalibrationFolders,
    channels = ['AntiKt4H1TowerAODJets'],
    channelAliases = aliases,
    shadowFoldersAndChannels = PerfCalibrationSingleFolder,
    singleFolderName = SingleFolderName,
    OutputLevel = VERBOSE
    )
ToolSvc += JetTagPerformanceCalibrationBroker

# This creates a particular instance of the performance evaluation tool.

# OPs = [ 'Loose', 'Medium', 'Tight' ]

from CalibrationDataInterface.CalibrationDataInterfaceConf import Analysis__CalibrationDataInterfaceTool
JetTagPerformanceSV0 = Analysis__CalibrationDataInterfaceTool(
    name = "BTagPerformanceSV0",
    taggerName = "SV0",
#    operatingPoints = [ 'Loose', 'Medium', 'Tight' ],
    operatingPoints = [ '5_72' ],
    PerformanceBroker = JetTagPerformanceCalibrationBroker,
    OutputLevel = BTaggingFlags.OutputLevel
    )
ToolSvc += JetTagPerformanceSV0

from IOVDbSvc.CondDB import conddb 

if PerfCalibrationSingleFolder:
    PerfCalibrationFolders = []
    PerfCalibrationFolders.append(PerfCalibrationRoot+SingleFolderName)

# Note that the code below uses local variables rather than BTaggingFlags

for folder in PerfCalibrationFolders:
    if CalibrationFromLocalReplica:
        connection = "<dbConnection>sqlite://X;schema=mycool.db;dbname=OFLP200</dbConnection> "+folder+" <tag>"+CalibrationTag+"</tag>"
        conddb.addFolder("", connection)
        print 'connection string: ' + connection
    else:
        if CalibrationFromCERN:
            conddb.addFolder("","<dbConnection>oracle://ATLAS_COOLPROD;schema=ATLAS_COOLOFL_GLOBAL;dbname=OFLP200</dbConnection> "+folder+" <tag>"+CalibrationTag+"</tag>") 
        else: 
            if CalibrationTag == "":
                conddb.addFolder('GLOBAL_OFL',folder) 
            else:
                conddb.addFolder('GLOBAL_OFL',folder+" <tag>"+CalibrationTag+"</tag>") 

if BTaggingFlags.OutputLevel < 3: 
  print JetTagPerformanceCalibrationBroker
  print JetTagPerformanceSV0

