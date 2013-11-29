#!/usr/env python
# jiri kvita
# Wed Jun 22 13:37:21 CEST 2011


import sys
import os

TaggersOP = [ 'SV050',
            'JetProb50',
            'SV1IP3D60', 'SV1IP3D70', 'SV1IP3D80',
            'JetFitterCOMBNN57', 'JetFitterCOMBNN60', 'JetFitterCOMBNN70', 'JetFitterCOMBNN80'
            ]


Taggers = [ 'SV0',
            'JetProb',
            'SV1IP3D', 'SV1IP3D', 'SV1IP3D',
            'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN', 'JetFitterCOMBNN'
            ]


OPs = [ '5_65',
        '2_65',
        '4_55', '1_70', '-0_80',
        '2_20', '1_80', '0_35', '-1_25'
        ]

Draws = [ 1,
          1,
          1, 1, 1,
          1, 1, 1, 1
          ]

Flavours = { 'b' : '$b$', 'c':'$c$', 'light':'Light'}
outfile = open('EffFits_new.tex', 'w')
for flavour in Flavours:
    for tagop,tag,op,draw in zip(TaggersOP,Taggers,OPs,Draws):
        if draw == 0:
            continue
        cutval = op.replace('_', '.')
        outfile.write('% _____________________________________________________________________ %\n')
        outfile.write('\\frame{\n')
        outfile.write('\\frametitle{' + tagop + ' (' + cutval + ') ' + Flavours[flavour] + '-jets Efficiency}\n')
        outfile.write('\\vskip-0.30cm\n')
        outfile.write('%\\begin{itemize}\n')
        outfile.write('%  \\item \n')
        outfile.write('%\\end{itemize}\n')
        outfile.write('\\begin{tabular}{cc}\n')
        outfile.write('\\includegraphics[width=0.46\\textwidth]{eps/c' + flavour + '' + tag + '_' + op + '.eps} & \n')
        outfile.write('\\includegraphics[width=0.46\\textwidth]{eps/c' + flavour + 'Parm' + tag + '_' + op + '.eps} \\\\ \n')
        outfile.write('\\includegraphics[width=0.46\\textwidth]{eps/c' + flavour + 'Check' + tag + '_' + op + '.eps} & \n')
        # use either of the two lines here:
        outfile.write('\\includegraphics[width=0.46\\textwidth]{eps/c' + flavour + 'Check' + tag + '_' + op + '_zoom.eps}\\\\ \n')
        #outfile.write('\\includegraphics[width=0.46\\textwidth]{eps/c' + flavour + 'Check' + tag + '_' + op + '_box.eps}\\\\ \n')
        outfile.write('Histo/Fit &  zoomed to $\\pm$ 10\\% \\\\\n')
        outfile.write('\\end{tabular}\n')
        outfile.write('}\n')
        outfile.write('\n')
    #
#

outfile.close()

outfile = open('SFs_new.tex', 'w')
for flavour in Flavours:
    outfile.write('% _____________________________________________________________________ % \n')
    for tagop,tag,op in zip(TaggersOP,Taggers,OPs):
        cutval = op.replace('_', '.')
        outfile.write('% _____________________________________________________________________ % \n')
        outfile.write('\\frame{ \n')
        outfile.write('\\frametitle{' + tagop + ' (' + cutval + ') ' + Flavours[flavour] + '-jets SF and Error} \n')
        outfile.write('\\vskip-0.30cm \n')
        outfile.write('%\\begin{itemize} \n')
        outfile.write('%  \\item  \n')
        outfile.write('%\\end{itemize} \n')
        outfile.write('\\begin{tabular}{cc} \n')
        outfile.write('\\includegraphics[width=0.48\\textwidth]{eps/c' + flavour + '_SF_' + tag + '_' + op + '_1.eps} & \n')
        outfile.write('\\includegraphics[width=0.48\\textwidth]{eps/c' + flavour + '_SF_' + tag + '_' + op + '_2.eps} \\\\ \n')
        outfile.write('\\end{tabular}\n')
        outfile.write('}\n')
    #
#

