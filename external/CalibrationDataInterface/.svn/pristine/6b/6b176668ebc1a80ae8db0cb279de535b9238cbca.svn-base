#!/usr/bin/python

import sys, os

flavours = [ 'B', 'C', 'Light']
SFs = ['true', 'false']
vartype = { SFs[0]: 'TF1', SFs[1]: 'TF2'}


for flav in flavours:
    for SF in SFs:
        hpp=open('MyTypeDef.h', 'w');
        mytype = vartype[SF]
        if flav == 'Light':
            mytype = 'TF2'
        hpp.write('typedef ' + mytype + ' MyType;\n')
        hpp.close()
        os.system('rm SimpleDraw_C.so')
        os.system('root -l \'SimpleDraw.C+("' + flav+ '", ' + SF+ ')\'')
