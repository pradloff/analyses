#include "egammaAnalysisUtils/EnergyRescaler.h"

#include <iostream>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TROOT.h"

void plotAFCorrection(string year) {

   gROOT->ProcessLine(".L EnergyRescaler.cxx");
  
   eg2011::EnergyRescaler  ers;
   ers.useDefaultCalibConstants(year);

   double etaMin = -2.5;
   double etaMax = +2.5;
   double etaStep = .01;
   double eta = etaMin;

   TGraph* gAF = new TGraph();
   int i=0;

   while( eta<=etaMax ) {

     gAF->SetPoint( i, eta, ers.applyAFtoG4(eta) );

     i++;
     eta += etaStep;

   }
 
   TCanvas * canvas = new TCanvas();

   gAF->SetMaximum(1.05);
   gAF->SetMinimum(0.95);

   char grafname[99];
   sprintf(grafname,"AF -> G4 Correction factor, %s", year.c_str());
   gAF->SetTitle(grafname);
   
   gAF->SetLineWidth(2);
   gAF->SetLineStyle(2);
   gAF->SetLineColor(4);
   gAF->GetXaxis()->SetTitleOffset(1.5);
   gAF->GetXaxis()->SetTitle("#eta");
   gAF->Draw("AL");
   
   char plotname[99];
   sprintf(plotname,"aftog4_%s.png", year.c_str());
   canvas->SaveAs(plotname);
   
}
