#!/usr/bin/python

import os, sys

run_script_name = '_run_fits.sh'
run_script = open(run_script_name, 'w')

taggername = [ 'SV0', 'JetProb', 'JetProb', 'JetProb', 'JetProb', 'SV0', 'JetProb', 'JetProb', 'JetProb', 'JetProb']
opname = [ '5_85', '3_25', '2_05', '1_40', '0_60', '5_85', '3_25', '2_05', '1_40', '0_60']

EffName = ['SV050', 'JetProb50', 'JetProb70', 'JetProb80', 'JetProb90']

for eff in EffName:
    cname = 'scaleFactors_' + eff + '_C_mistags.txt'
    fname = 'scaleFactors_' + eff + '_F_mistags.txt'
    mergename = 'scaleFactors_' + eff + '_mistags.txt'

    # first concat the Central and Forward SF files:
    os.system('rm -f ' + mergename)
    os.system('echo abseta 0. >> ' + mergename) 
    os.system('echo abseta 1.2 >> ' + mergename) 
    os.system('cat ' + cname +  ' >> ' + mergename)
    os.system('echo abseta 2.5 >> ' + mergename) 
    os.system('cat ' + fname +  ' >> ' + mergename)
# end of merge

# sys.exit(0)
    
files = [ 'scaleFactors_SV050_ptrel.txt', 
          'scaleFactors_JetProb_50_ptrel.txt', 
          'scaleFactors_JetProb_70_ptrel.txt', 
          'scaleFactors_JetProb_80_ptrel.txt', 
          'scaleFactors_JetProb_90_ptrel.txt',

          'scaleFactors_SV050_mistags.txt',
          'scaleFactors_JetProb50_mistags.txt',
          'scaleFactors_JetProb70_mistags.txt',
          'scaleFactors_JetProb80_mistags.txt',
          'scaleFactors_JetProb90_mistags.txt',

          ]

flavour = [ 'Heavy', 'Heavy', 'Heavy', 'Heavy', 'Heavy',
            'Light', 'Light', 'Light', 'Light', 'Light']
writeopt = [ 'recreate', 'update', 'update', 'update', 'update',
             '', '', '', '', '']

script=[ 'ParseAscii_pt.py', 'ParseAscii_pt.py',  'ParseAscii_pt.py', 'ParseAscii_pt.py',  'ParseAscii_pt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py']

calibrootfilename = 'TopCalibrations_PLHC_2011.root'

# first prepare the hpp file for light SF:
ils=[-1, -2, -3, -4, -5]
for il in ils:
    hppfile='SF_' + opname[il] + '_' + flavour[il] + '.hpp'
    os.system('python ' + script[il] + ' ' + files[il] + ' ' + hppfile + ' ' + flavour[il])


hpplink = 'SF.hpp'
hpplinkLight = 'SF_Light.hpp'

# now make hpp files for each tagger, OP and create the calib root file
# by appending the results:
for i in range(0, len(opname)-len(ils)):
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
print 'Run: ./%s' % (run_script_name,)
run_script.close()
