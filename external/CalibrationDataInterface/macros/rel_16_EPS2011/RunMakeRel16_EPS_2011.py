#!/usr/bin/python

import os, sys

run_script_name = '_run_fits.sh'
run_script = open(run_script_name, 'w')

taggername = [ 'SV0',
               'JetProb', 'JetProb', 'JetProb', 'JetProb',
               'SV1IP3D', 'SV1IP3D', 'SV1IP3D', 'SV1IP3D',
               'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN',
               'SV0', 'JetProb', 'JetProb', 'JetProb', 'JetProb', 
               'SV1IP3D', 'SV1IP3D', 'SV1IP3D', 'SV1IP3D',
               'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN'
               ]

opname = [ '5_85',
           '3_25', '2_05', '1_40', '0_60',
           '7_60', '4_50', '1_55', '-0_85',
           '3_00', '2_00', '0_35', '-1_25', '2_40',
           '5_85',
           '3_25', '2_05', '1_40', '0_60',
           '7_60', '4_50', '1_55', '-0_85',
           '3_00', '2_00', '0_35', '-1_25', '2_40'
           ]

# to merge forward and central SF mistag files"
EffName = {'sv0_5.85': 'SV0_50',
           'jetp_3.25': 'JetProb50', 'jetp_2.05': 'JetProb70', 'jetp_1.4': 'JetProb80', 'jetp_0.6': 'JetProb90',
           'cmb_7.6' : 'SV1IP3D50', 'cmb_4.5' : 'SV1IP3D60', 'cmb_1.55' : 'SV1IP3D70', 'cmb_-0.85' : 'SV1IP3D80',
           'jfitCOMBNN_3' : 'JetFitterCOMBNN50', 'jfitCOMBNN_2' : 'JetFitterCOMBNN60',
           'jfitCOMBNN_0.35' : 'JetFitterCOMBNN70', 'jfitCOMBNN_-1.25' : 'JetFitterCOMBNN80', 'jfitCOMBNN_2.40' : 'JetFitterCOMBNN57'
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
          ptreldirname + 'scaleFactors_JetProb_70_ptrel.txt', 
          ptreldirname + 'scaleFactors_JetProb_80_ptrel.txt', 
          ptreldirname + 'scaleFactors_JetProb_90_ptrel.txt',
          ptreldirname + 'scaleFactors_SV1IP3D_50_ptrel.txt',
          ptreldirname + 'scaleFactors_SV1IP3D_60_ptrel.txt',
          ptreldirname + 'scaleFactors_SV1IP3D_70_ptrel.txt',
          ptreldirname + 'scaleFactors_SV1IP3D_80_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_50_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_60_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_70_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_80_ptrel.txt',
          ptreldirname + 'scaleFactors_JetFitterCOMBNN_57_ptrel.txt',
          
          negtagdirname + 'scaleFactors_SV0_50_mistags.txt',
          negtagdirname + 'scaleFactors_JetProb50_mistags.txt',
          negtagdirname + 'scaleFactors_JetProb70_mistags.txt',
          negtagdirname + 'scaleFactors_JetProb80_mistags.txt',
          negtagdirname + 'scaleFactors_JetProb90_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D50_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D60_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D70_mistags.txt',
          negtagdirname + 'scaleFactors_SV1IP3D80_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN50_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN60_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN70_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN80_mistags.txt',
          negtagdirname + 'scaleFactors_JetFitterCOMBNN57_mistags.txt'
          ]

flavour = [ 'Heavy',
            'Heavy', 'Heavy', 'Heavy', 'Heavy',
            'Heavy', 'Heavy', 'Heavy', 'Heavy',
            'Heavy', 'Heavy', 'Heavy', 'Heavy', 'Heavy',
            'Light',
            'Light', 'Light', 'Light', 'Light',
            'Light', 'Light', 'Light', 'Light',
            'Light', 'Light', 'Light', 'Light', 'Light'
            ]
writeopt = [ 'recreate',
             'update', 'update', 'update', 'update',
             'SKIP', 'update', 'update', 'update',
             'SKIP', 'update', 'update', 'update', 'update',
             '',
             '', '', '', '',
             '', '', '', '',
             '', '', '', '', ''
             ]

script=[ 'ParseAscii_pt.py',
         'ParseAscii_pt.py',  'ParseAscii_pt.py', 'ParseAscii_pt.py',  'ParseAscii_pt.py',
         'ParseAscii_pt.py',  'ParseAscii_pt.py', 'ParseAscii_pt.py',  'ParseAscii_pt.py',
         'ParseAscii_pt.py',  'ParseAscii_pt.py', 'ParseAscii_pt.py',  'ParseAscii_pt.py',  'ParseAscii_pt.py',
         'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py'
         ]

calibrootfilename = 'TopCalibrations_EPS_2011_v1.root'

# first prepare the hpp file for light SF:
# so let's go from the back;)
ils=range(-1, -14, -1)
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
    run_script.write('#############################\n')
    run_script.write('rm ' + hpplink + ' ' + hpplinkLight + '\n')
    run_script.write('ln -s ' + hppfile + ' ' + hpplink + '\n')
    run_script.write('ln -s ' + hppfileLight + ' ' + hpplinkLight + '\n')
    run_script.write('rm ' + FileBaseName + '_C.so\n')
    run_script.write('root -l -b -q \'' + FileBaseName + '.C+("' + taggername[i] + '", "' + opname[i] + '", "' + calibrootfilename + '", "' + writeopt[i] + '")\' >& log_' + taggername[i] + '_' + opname[i] + '.txt\n')

os.system('chmod +x ' + run_script_name)
print 'Run: '
print '     ./%s' % (run_script_name,)
print '     ./ParseLog.sh > tables.tex'
print '     python GenEffSFTexCode.py'
run_script.close()
