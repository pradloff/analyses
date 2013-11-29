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
ptbins = []
SF = []
SFerr = []
for line in asciifile:
    print "Processing line %i" % (iline,)
    elements = line.split()
    if len(elements) > 0:
        if str(elements[0])[0:2] == '//':
            continue
        ptbins.append(elements[0])
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

print "Read:"
print ptbins
print SF
print SFerr


# print SF code:
outfile.write("  // _________________________________________________________\n")
outfile.write("  double " + flavour + "FlavourSF(Double_t* xx, Double_t* par) {\n")
outfile.write("    double pt = xx[0];\n")
count = 0
N=len(ptbins)-1
for i in range(0, N):
    ptstr = str(ptbins[i+1])
    sfstr = str(SF[i])
    if count == 0:
        outfile.write("    if (pt < " + ptstr + ")\n")
        outfile.write("      return " + sfstr + ";\n")
    elif count < N-1:
        outfile.write("    else if (pt < " + ptstr + ")\n")
        outfile.write("      return " + sfstr + ";\n")
    else:
        outfile.write("    else\n")
        outfile.write("      return " + sfstr + ";\n")
    #
    count = count + 1
#
outfile.write("  }")


# scale parameter serves to scale the error by an additional factor
# this has been used to scale the ctag SF errors by 2 w.r.t. btag results

outfile.write("  // _________________________________________________________\n")
outfile.write("  double " + flavour + "FlavourSFSyst(Double_t* xx, Double_t* par) {\n")
outfile.write("    double pt = xx[0];\n")
outfile.write("    double scale = par[0];\n")
outfile.write("    double stat = 0.;\n")
outfile.write("    double syst = 0.;\n")
count = 0
N=len(ptbins)-1
for i in range(0, N):
    ptstr = str(ptbins[i+1])
    sfstr = str(SFerr[i])
    if count == 0:
        outfile.write("    if (pt < " + ptstr + ")\n")
        outfile.write("      syst = " + sfstr + ";\n")
    elif count < N-1:
        outfile.write("    else if (pt < " + ptstr + ")\n")
        outfile.write("      syst = " + sfstr + ";\n")
    else:
        outfile.write("    else\n")
        outfile.write("      syst = " + sfstr + ";\n")
    #
    count = count + 1
#
outfile.write("    return TMath::Sqrt(stat*stat+syst*syst)*scale;\n")
outfile.write("  }\n")


asciifile.close()
outfile.close()


