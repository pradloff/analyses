


Instructions on how to prepare the btag calibration ROOT file
Jiri Kvita, 18.2.2011

=== Rel.16 ===

Needed inputs: 

1) SF and SF errors
(.txt files obtained from the flavour tagging group 
and manually slit into separate files)

SV050_PreliminaryNumbers_Efficiencies_2011_02_15.txt
JetProb50_PreliminaryNumbers_Efficiencies_2011_02_15.txt
JetProb70_PreliminaryNumbers_Efficiencies_2011_02_15.txt
SV050_PreliminaryNumbers_Mistags_2011_02_15.txt
JetProb50_PreliminaryNumbers_Mistags_2011_02_15.txt
JetProb70_PreliminaryNumbers_Mistags_2011_02_15.txt

2) ttbar7.root
 - provided by Top Group (Martin zur Nedden)
 - includes the b,c,l(=uds) efficiencies

3) Scripting code:
RunMakeRel16_prelim.py
 - prepares header files which hold the SF and SFerr TF1,TF2-generating code.
 - this parsing is done by ParseAscii_pt.py for b,c SF in pT
                           ParseAscii_absetapt.py for light SF in abseta, pt
 - prepares script _run_fits.sh, which runs the fitting and makes the calibration file

4) The real code:
fitTopHists_rel16.C
+additional methods in Tools.C

*** So the running if all inputs are present looks as follows ***
./RunMakeRel16_prelim.py
./_run_fits.sh

5) Auxiliary:
 i) Displaying results:
  
./RunSimpleDraw.py
which calls root SimpleDraw.C
This needed some C++ typedef tweaks as some SFs are TF1 while other TF2.
