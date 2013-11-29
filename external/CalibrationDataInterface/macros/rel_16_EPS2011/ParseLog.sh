#!/bin/sh
for i in `ls log*.txt` ; do
  echo "%%% $i %%%"
  name=`echo $i | sed "s/log_//g" | sed "s/\.txt//g" | sed "s/JetProb_/JetProb\ /g" | sed "s/SV0_/SV0\ /g" | sed "s/SV1IP3D_/SV1IP3D\ /g" | sed "s/JetFitterCOMBNN_/JetFitterCOMBNN\ /g" | sed "s/_/\./g"`
  echo "  % _____________________________________________________________________ %"
  echo "\\frame{"
  echo "\\frametitle{Tables ${name}}"
#  echo "\\vskip-0.35cm"
  echo "\\begin{itemize}"
  echo ""
  cat $i | grep chi2 | sed "s/chi2/$\\\chi^2$/g" | sed "s/\ /,\ /g"| awk '{printf("  \\item "); print $ALL;}'
  cat $i | grep Nneg | cut -d " " -f 4,5,6,7 | awk '{printf("  \\item "); print $ALL;}'
  echo "  \\bigskip"
  echo "  \\item \$N_{\\rm bad}\$: (histo-fit) / histo \$>\$ 0.10"
  echo "  \\item \$N_{\\rm sig}\$: (histo-fit) / error \$>\$ 2.5"
  echo "\\end{itemize}"
  echo "}"
  echo ""

done
