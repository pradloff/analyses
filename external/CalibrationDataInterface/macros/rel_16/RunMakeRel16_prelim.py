#!/usr/bin/python

import os, sys

run_script_name = '_run_fits.sh'
run_script = open(run_script_name, 'w')

taggername = [ 'SV0', 'JetProb', 'JetProb', 'SV0', 'JetProb', 'JetProb']
opname = [ '5_85', '3_25', '2_05', '5_85', '3_25', '2_05']
files = [ 'SV050_PreliminaryNumbers_Efficiencies_2011_02_15.txt',
          'JetProb50_PreliminaryNumbers_Efficiencies_2011_02_15.txt',
          'JetProb70_PreliminaryNumbers_Efficiencies_2011_02_15.txt',
          'SV050_PreliminaryNumbers_Mistags_2011_02_15.txt',
          'JetProb50_PreliminaryNumbers_Mistags_2011_02_15.txt',
          'JetProb70_PreliminaryNumbers_Mistags_2011_02_15.txt'
          ]
flavour = [ 'Heavy', 'Heavy', 'Heavy', 'Light', 'Light', 'Light']
writeopt = [ 'recreate', 'update', 'update', '', '', '']

script=[ 'ParseAscii_pt.py', 'ParseAscii_pt.py',  'ParseAscii_pt.py',
         'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py', 'ParseAscii_absetapt.py']

calibrootfilename = 'TopCalibrations_rel16_prelim.root'

# first prepare the hpp file for light SF:
ils=[-1, -2, -3]
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
    FileBaseName='fitTopHists_rel16'
    run_script.write('#############################\n')
    run_script.write('rm ' + hpplink + ' ' + hpplinkLight + '\n')
    run_script.write('ln -s ' + hppfile + ' ' + hpplink + '\n')
    run_script.write('ln -s ' + hppfileLight + ' ' + hpplinkLight + '\n')
    run_script.write('rm ' + FileBaseName + '_C.so\n')
    run_script.write('root -l -b -q \'' + FileBaseName + '.C+("' + taggername[i] + '", "' + opname[i] + '", "' + calibrootfilename + '", "' + writeopt[i] + '")\' >& log_' + taggername[i] + '_' + opname[i] + '.txt\n')

os.system('chmod +x ' + run_script_name)
print 'Run: ./%s' % (run_script_name,)
run_script.close()
