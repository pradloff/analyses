#!/usr/bin/python

import os, sys

run_script_name = '_run_fits.sh'
run_script = open(run_script_name, 'w')

taggername = [ 'SV0',
               'JetProb',
               'SV1IP3D', 'SV1IP3D', 'SV1IP3D',
               'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN',
               'MV1','MV1','MV1','MV1',
               'SV0',
               'JetProb', 
               'SV1IP3D', 'SV1IP3D', 'SV1IP3D',
               'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN',
               'MV1','MV1','MV1','MV1'
               ]

opname = [ '5_65',
           '2_65',
           '4_55', '1_70', '-0_80',
           '2_20', '1_80', '0_35', '-1_25',
           '0_910','0_614','0_416','0_080',
           '5_65',
           '2_65',
           '4_55', '1_70', '-0_80',
           '2_20', '1_80', '0_35', '-1_25',
           '0_910','0_614','0_416','0_080'
           ]

# to merge forward and central SF mistag files"
EffName = {'sv0_5.65': 'SV0_50',
           'jetp_2.65': 'JetProb_50',
           'cmb_4.55' : 'SV1IP3D_60', 'cmb_1.70' : 'SV1IP3D_70', 'cmb_-0.80' : 'SV1IP3D_80',
           'jfitCOMBNN_2.20' : 'JetFitterCOMBNN_57', 'jfitCOMBNN_1.80' : 'JetFitterCOMBNN_60',
           'jfitCOMBNN_0.35' : 'JetFitterCOMBNN_70', 'jfitCOMBNN_-1.25' : 'JetFitterCOMBNN_80',
           'MV1_0.910' : 'MV1_60','MV1_0.614' : 'MV1_70','MV1_0.416' : 'MV1_75','MV1_0.080' : 'MV1_85'
           }

negtagdirname = 'mistag_negtag/' # 'negtag_smear/'
for eff in EffName:
    cname = negtagdirname + 'out_' + eff + '_central'
    fname = negtagdirname + 'out_' + eff + '_forward'
    mergename = negtagdirname + 'scaleFactors_' + EffName[eff] + '_mistags.txt'
    print 'Merging %s %s to %s' % (cname, fname, mergename)
    # first concat the Central and Forward SF files:
    os.system('rm -f ' + mergename)
    os.system('echo abseta 0. >> ' + mergename) 
    os.system('echo abseta 1.2 >> ' + mergename) 
    os.system('cat ' + cname +  ' >> ' + mergename)
    os.system('echo abseta 2.5 >> ' + mergename) 
    os.system('cat ' + fname +  ' >> ' + mergename)
# end of merge

#sys.exit(0)

ptreldirname = 'eff_ptrel_new/' # 'new_ptrel/'

files = [ ptreldirname + 'scaleFactors_SV0_50_ptrel.txt', 
          ptreldirname + 'scaleFactors_JetProb_50_ptrel.txt', 
          ptreldirname + 'scaleFactors_SV1IP3D_60_ptrel.txt',
          ptreldirname + 'scaleFactors_SV1IP3D_70_ptrel.txt',
          ptreldirname + 'scaleFactors_SV1IP3D_80_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_57_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_60_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_70_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_80_ptrel.txt',
          ptreldirname + 'scaleFactors_MV1_60_ptrel.txt',
          ptreldirname + 'scaleFactors_MV1_70_ptrel.txt',
          ptreldirname + 'scaleFactors_MV1_75_ptrel.txt',
          ptreldirname + 'scaleFactors_MV1_85_ptrel.txt',
          
          negtagdirname + 'scaleFactors_SV0_50_mistags.txt',
          negtagdirname + 'scaleFactors_JetProb_50_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D_60_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D_70_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D_80_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN_57_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN_60_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN_70_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN_80_mistags.txt',
          negtagdirname + 'scaleFactors_MV1_60_mistags.txt',
          negtagdirname + 'scaleFactors_MV1_70_mistags.txt',
          negtagdirname + 'scaleFactors_MV1_75_mistags.txt',
          negtagdirname + 'scaleFactors_MV1_85_mistags.txt'
          ]

flavour = [ 'Heavy',
            'Heavy', 
            'Heavy', 'Heavy', 'Heavy',
            'Heavy', 'Heavy', 'Heavy', 'Heavy',
            'Heavy', 'Heavy', 'Heavy', 'Heavy',
            'Light',
            'Light',
            'Light', 'Light', 'Light',
            'Light', 'Light', 'Light', 'Light',
            'Light', 'Light', 'Light', 'Light'
            ]
writeopt = [ 'recreate',
             'update',
             'update', 'update', 'update',
             'update', 'update', 'update', 'update',
             'update', 'update', 'update', 'update',
             '',
             '',
             '', '', '',
             '', '', '', '',
             '', '', '', ''
             ]

script=[ 'ParseAscii_pt.py',
         'ParseAscii_pt.py',
         'ParseAscii_pt.py', 'ParseAscii_pt.py', 'ParseAscii_pt.py',
         'ParseAscii_pt.py', 'ParseAscii_pt.py', 'ParseAscii_pt.py', 'ParseAscii_pt.py',
         'ParseAscii_pt.py', 'ParseAscii_pt.py', 'ParseAscii_pt.py', 'ParseAscii_pt.py',
         'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py'
         ]

calibrootfilename = 'TopCalibrations_rel17_MC11a_new.root'
calibEffrootfilename = 'TopCalibrations_rel17_MC11a_Eff_new.root'

# first prepare the hpp file for light SF:
# so let's go from the back;)
ils=range(-1, -13, -1)
for il in ils:
    hppfile='SF_' + opname[il] + '_' + flavour[il] + '.hpp'
    os.system('python ' + script[il] + ' ' + files[il] + ' ' + hppfile + ' ' + flavour[il])


hpplink = 'SF.hpp'
hpplinkLight = 'SF_Light.hpp'

# now make hpp files for each tagger, OP and create the calib root file
# by appending the results:
for i in range(0, len(opname)-len(ils)):
    if writeopt[i] == 'SKIP':
        continue
    hppfile = 'SF_' + opname[i] + '_' + flavour[i] + '.hpp'
    hppfileLight = 'SF_' + opname[i] + '_Light.hpp'
    # make the hpp file:
    os.system('python ' + script[i] + ' ' + files[i] + ' ' + hppfile + ' ' + flavour[i])

    # use the hpp file to append the tagger operation point:
    FileBaseName='fitTopHists'
    FileBaseNameEff='fitEffTopHists'
    run_script.write('#############################\n')
    run_script.write('rm -f ' + hpplink + ' ' + hpplinkLight + '\n')
    run_script.write('ln -s ' + hppfile + ' ' + hpplink + '\n')
    run_script.write('ln -s ' + hppfileLight + ' ' + hpplinkLight + '\n')
    run_script.write('root -l -b -q \'' + FileBaseName + '.C("' + taggername[i] + '", "' + opname[i] + '", "' + calibrootfilename + '", "' + writeopt[i] + '")\' > log_' + taggername[i] + '_' + opname[i] + '.txt\n')
    run_script.write('root -l -b -q \'' + FileBaseNameEff + '.C("' + taggername[i] + '", "' + opname[i] + '", "' + calibEffrootfilename + '", "' + writeopt[i] + '")\' > log_Eff_' + taggername[i] + '_' + opname[i] + '.txt\n')

os.system('chmod +x ' + run_script_name)
print 'Run: '
print '     ./%s' % (run_script_name,)
print '     ./ParseLog.sh > tables.tex'
print '     python GenEffSFTexCode.py'
run_script.close()
