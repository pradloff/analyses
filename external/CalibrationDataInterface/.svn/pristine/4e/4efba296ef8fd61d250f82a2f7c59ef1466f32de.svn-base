#!/usr/bin/python

import sys, os

#print "Called with arguments:"
print sys.argv

if len(sys.argv) > 1:
    infilename = sys.argv[1]
else:
    infilename = 'test_ascii.txt'
print 'Will open infile %s' % (infilename,)
asciifile = open(infilename, 'r')

if len(sys.argv) > 2:
    outfilename = sys.argv[2]
else:
    outfilename = 'test_out.hpp'
print 'Will open outfile %s' % (outfilename,)
outfile = open(outfilename, 'w')


# Heavy or Light?
if len(sys.argv) > 3:
    flavour = sys.argv[3]
else:
    flavour = 'Heavy'
    
iline = 0
etabins = []
nptbins = []
ptbins = []
SF = []
SFerr = []
npt = 0
etaptbins = []
for line in asciifile:
    print "Processing line %i" % (iline,)
    print line
    elements = line.split()
    if len(elements) > 0:
        if str(elements[0])[0:2] == '//':
            continue
        if str(elements[0]) == 'abseta':
            if len(elements) > 1:
                etabins.append(elements[1])
            if len(etabins) > 2:
                print 'Appending npt=' + str(npt)
                ptbins.append(etaptbins)
                etaptbins = []
        else:
            etaptbins.append(elements[0])
            if iline > 0:
                if len(elements) > 1:
                    SF.append(elements[1])
                if len(elements) > 2:
                    SFerr.append(elements[2])
        #
    else:
        print "Error: empty line!"
    
            
    # elements
    iline = iline + 1
    #
#
ptbins.append(etaptbins)
etaptbins = []
                
print "Read:"
print len(etabins)
print etabins
print len(ptbins)
print ptbins
print len(SF)
print SF
print len(SFerr)
print SFerr

####
###
##
#

# print SF code:
outfile.write("  // _________________________________________________________\n")
outfile.write("  double " + flavour + "FlavourSF(Double_t* xx, Double_t* par) {\n")
outfile.write("    double abseta = xx[1];\n")
outfile.write("    double pt = xx[0];\n")
counteta = 0
countpt = 0
count = 0
Neta=len(etabins)-1
for j in range(0, Neta):
    etastr = str(etabins[j+1])
    if counteta == 0:
        outfile.write("    if (abseta < " + etastr + ") {\n")
    else:
        outfile.write("    else if (abseta < " + etastr + ") {\n")

    Npt=len(ptbins[counteta])
    for i in range(0, Npt-1):
        ptstr = str(ptbins[counteta][i+1])
        sfstr = str(SF[count])
        if i == 0:
            outfile.write("      if (pt < " + ptstr + ")\n")
            outfile.write("        return " + sfstr + ";\n")
        elif i < Npt-2:
            outfile.write("      else if (pt < " + ptstr + ")\n")
            outfile.write("        return " + sfstr + ";\n")
        else:
            outfile.write("      else\n")
            outfile.write("        return " + sfstr + ";\n")
        #
        countpt = countpt + 1
        if i < Npt-1:
            count = count + 1

    outfile.write("    }\n")
    counteta = counteta + 1
#
outfile.write("    else return 1.;\n")
outfile.write("  }\n")

# print SFerr code:
# scale parameter serves to scale the error by an additional factor
# this has been used to scale the ctag SF errors by 2 w.r.t. btag results

outfile.write("  // _________________________________________________________\n")
outfile.write("  double " + flavour + "FlavourSFSyst(Double_t* xx, Double_t* par) {\n")
outfile.write("    double abseta = xx[1];\n")
outfile.write("    double pt = xx[0];\n")
outfile.write("    double scale = par[0];\n")
outfile.write("    double stat = 0.;\n")
outfile.write("    double syst = 0.;\n")
counteta = 0
countpt = 0
count = 0
Neta=len(etabins)-1
for j in range(0, Neta):
    etastr = str(etabins[j+1])
    if counteta == 0:
        outfile.write("    if (abseta < " + etastr + ") {\n")
    else:
        outfile.write("    else if (abseta < " + etastr + ") {\n")

    Npt=len(ptbins[counteta])
    for i in range(0, Npt-1):
        ptstr = str(ptbins[counteta][i+1])
        sfstr = str(SFerr[count])
        if i == 0:
            outfile.write("      if (pt < " + ptstr + ")\n")
            outfile.write("        syst = " + sfstr + ";\n")
        elif i < Npt-2:
            outfile.write("      else if (pt < " + ptstr + ")\n")
            outfile.write("        syst = " + sfstr + ";\n")
        else:
            outfile.write("      else\n")
            outfile.write("        syst = " + sfstr + ";\n")
        #
        countpt = countpt + 1
        if i < Npt-1:
            count = count + 1

    outfile.write("    }\n")
    counteta = counteta + 1
#

outfile.write("    // return TMath::Sqrt(stat*stat+syst*syst)*scale;\n")
outfile.write("    return syst;\n")
outfile.write("  }\n")


asciifile.close()
outfile.close()


