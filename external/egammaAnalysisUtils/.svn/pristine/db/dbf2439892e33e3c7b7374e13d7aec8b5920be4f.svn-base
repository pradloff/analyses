// eta : from cluster
// emin, emax : energies (NOT transverse), in GeV
// ptype : 0 (electron), 1 (unconverted photon), 2 (converted photon)
// logx : 0/1

TCanvas* plotLinearity_ERUpgrade(double eta, double etmin /* GeV */, double etmax /* GeV */, int ptype, string year, 
				 bool debug=false, bool logx=false, bool subtractAverage=false,
				 TGraph* overlay=0 ) {

  double GeV = 1000;

   gROOT->ProcessLine(".L EnergyRescalerUpgrade.cxx");
  
   egRescaler::EnergyRescalerUpgrade ers;
   if( year=="2011" )
     ers.Init("../share/EnergyRescalerData.root","2011","es2011a");
   else if( year=="2012" )
     ers.Init("../share/EnergyRescalerData.root","2012","es2012");

   ers.SetDebugFlag( debug );

   double eMin = etmin * cosh(eta) * GeV;
   double eMax = etmax * cosh(eta) * GeV;

   int nStep = 100;
   double eStep;
   if( logx )
     eStep = (log(eMax)-log(eMin))/(double)nStep;
   else
     eStep = (eMax-eMin)/(double)nStep;
   double e = eMin;

   TGraphErrors* gZee = new TGraphErrors();  // Z scale
   TGraph* gNom = new TGraph();        // nominal scale correction
   TGraph* gZeeAllUp = new TGraph();   // Zee scale uncertainty, Up
   TGraph* gZeeAllDown = new TGraph(); //                      , Down
   TGraph* gZeeUp = new TGraph();      // Zee scale uncertainty, Up
   TGraph* gZeeDown = new TGraph();    //                      , Down
   TGraph* gZeeGenUp = new TGraph();   // Zee scale uncertainty, Up
   TGraph* gZeeGenDown = new TGraph(); //                      , Down
   TGraph* gMATUp = new TGraph();      //       MAT uncertainty, Up
   TGraph* gMATDown = new TGraph();    //                      , Down
   TGraph* gPSUp = new TGraph();       //        ps uncertainty, Up
   TGraph* gPSDown = new TGraph();     //                      , Down
   TGraph* gLowUp = new TGraph();      //    low-pt uncertainty, Up
   TGraph* gLowDown = new TGraph();    //                      , Down

   int i=0;

   while( e<=eMax ) {

     double alpha          = ers.getAlphaZee(eta, egRescaler::EnergyRescalerUpgrade::Nominal);
     double alphaUp        = ers.getAlphaZee(eta, egRescaler::EnergyRescalerUpgrade::ZeeStatUp);

     double eNominal       = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::Nominal);

     double eVarZeeAllUp   = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::ZeeAllUp);
     double eVarZeeAllDown = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::ZeeAllDown);

     double eVarZeeUp      = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::ZeeMethodUp);
     double eVarZeeDown    = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::ZeeMethodDown);

     double eVarGenUp      = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::ZeeGenUp);
     double eVarGenDown    = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::ZeeGenDown);

     double eVarMATUp      = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::R12StatUp);
     double eVarMATDown    = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::R12StatDown);

     double eVarPSUp       = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::PSStatUp);
     double eVarPSDown     = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::PSStatDown);

     double eVarLowUp      = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::LowPtUp);
     double eVarLowDown    = ers.applyEnergyCorrection(eta, e, ptype, egRescaler::EnergyRescalerUpgrade::LowPtDown);

     if(i==0) {
       gZee->SetPointError(i,0,alphaUp-alpha);
       if( !subtractAverage )
	 gZee->SetPoint(i,40.*cosh(eta),alpha);
       else
	 gZee->SetPoint(i,40.*cosh(eta),0);
     }
     
     if( !subtractAverage ) {
       
       gNom->SetPoint(i, e/GeV, (e-eNominal)/eNominal); 
       
       gZeeAllUp->SetPoint(i, e/GeV, (e-eVarZeeAllUp)/eVarZeeAllUp); 
       gZeeAllDown->SetPoint(i, e/GeV, (e-eVarZeeAllDown)/eVarZeeAllDown); 
       
       gZeeUp->SetPoint(i, e/GeV, (e-eVarZeeUp)/eVarZeeUp); 
       gZeeDown->SetPoint(i, e/GeV, (e-eVarZeeDown)/eVarZeeDown); 
       
       gZeeGenUp->SetPoint(i, e/GeV, (e-eVarGenUp)/eVarGenUp); 
       gZeeGenDown->SetPoint(i, e/GeV, (e-eVarGenDown)/eVarGenDown); 
       
       gMATUp->SetPoint(i, e/GeV, (e-eVarMATUp)/eVarMATUp); 
       gMATDown->SetPoint(i, e/GeV, (e-eVarMATDown)/eVarMATDown); 
       
       gPSUp->SetPoint(i, e/GeV, (e-eVarPSUp)/eVarPSUp); 
       gPSDown->SetPoint(i, e/GeV, (e-eVarPSDown)/eVarPSDown); 
       
       gLowUp->SetPoint(i, e/GeV, (e-eVarLowUp)/eVarLowUp); 
       gLowDown->SetPoint(i, e/GeV, (e-eVarLowDown)/eVarLowDown); 
       
     } else {
       
       gNom->SetPoint(i, e/GeV, (e-eNominal)/eNominal - alpha); 
     
       gZeeAllUp->SetPoint(i, e/GeV, (e-eVarZeeAllUp)/eVarZeeAllUp - alpha); 
       gZeeAllDown->SetPoint(i, e/GeV, (e-eVarZeeAllDown)/eVarZeeAllDown - alpha); 
       
       gZeeUp->SetPoint(i, e/GeV, (e-eVarZeeUp)/eVarZeeUp - alpha); 
       gZeeDown->SetPoint(i, e/GeV, (e-eVarZeeDown)/eVarZeeDown - alpha); 
       
       gZeeGenUp->SetPoint(i, e/GeV, (e-eVarGenUp)/eVarGenUp - alpha); 
       gZeeGenDown->SetPoint(i, e/GeV, (e-eVarGenDown)/eVarGenDown - alpha); 
       
       gMATUp->SetPoint(i, e/GeV, (e-eVarMATUp)/eVarMATUp - alpha); 
       gMATDown->SetPoint(i, e/GeV, (e-eVarMATDown)/eVarMATDown - alpha); 
       
       gPSUp->SetPoint(i, e/GeV, (e-eVarPSUp)/eVarPSUp - alpha); 
       gPSDown->SetPoint(i, e/GeV, (e-eVarPSDown)/eVarPSDown - alpha); 
       
       gLowUp->SetPoint(i, e/GeV, (e-eVarLowUp)/eVarLowUp - alpha); 
       gLowDown->SetPoint(i, e/GeV, (e-eVarLowDown)/eVarLowDown - alpha); 
       
     }

     i++;

     if( logx )
       e *= exp(eStep);
     else
       e += eStep;

   }
 
   TCanvas* canvas = new TCanvas();
   if( logx )
     canvas->SetLogx();

   if( !subtractAverage ) {
     gZeeAllUp->SetMaximum(alpha+.01);
     gZeeAllUp->SetMinimum(alpha+-.015);
   } else {
     gZeeAllUp->SetMaximum(.01);
     gZeeAllUp->SetMinimum(-.015);
   }

   char grafname[99];
   if( ptype==0 )
     sprintf(grafname,"ES  Linearity (%s), #eta = %4.2f, Electrons", year.c_str(), eta);
   else if( ptype==1 )
     sprintf(grafname,"ES  Linearity (%s), #eta = %4.2f, Unconverted photons", year.c_str(), eta);
   else if( ptype==2 )
     sprintf(grafname,"ES  Linearity (%s), #eta = %4.2f, Converted photons", year.c_str(), eta);
   gZeeAllUp->SetTitle(grafname);
   
   gZeeAllUp->GetXaxis()->SetTitleOffset(1.2*gNom->GetXaxis()->GetTitleOffset());
   gZeeAllUp->GetYaxis()->SetTitleOffset(1.2*gNom->GetYaxis()->GetTitleOffset());
   gZeeAllUp->GetXaxis()->SetTitle("E [GeV]");
   gZeeAllUp->GetYaxis()->SetTitle("#alpha");
   gZeeAllUp->SetLineWidth(2);
   gZeeAllUp->SetLineStyle(1);
   gZeeAllUp->SetLineColor(1);
   gZeeAllUp->Draw("AL");

   gZee->SetMarkerStyle(20);
   gZee->Draw("P");

   gZeeAllDown->SetLineWidth(2);
   gZeeAllDown->SetLineStyle(2);
   gZeeAllDown->SetLineColor(1);
   gZeeAllDown->Draw("L");

   gZeeUp->SetLineWidth(2);
   gZeeUp->SetLineStyle(1);
   gZeeUp->SetLineColor(4);
   gZeeUp->Draw("L");
       
   gZeeDown->SetLineWidth(2);
   gZeeDown->SetLineStyle(2);
   gZeeDown->SetLineColor(4);
   gZeeDown->Draw("L");

   gZeeGenUp->SetLineWidth(2);
   gZeeGenUp->SetLineStyle(1);
   gZeeGenUp->SetLineColor(8);
   gZeeGenUp->Draw("L");

   gZeeGenDown->SetLineWidth(2);
   gZeeGenDown->SetLineStyle(2);
   gZeeGenDown->SetLineColor(8);
   gZeeGenDown->Draw("L");

   gMATUp->SetLineWidth(2);
   gMATUp->SetLineStyle(1);
   gMATUp->SetLineColor(2);
   gMATUp->Draw("L");

   gMATDown->SetLineWidth(2);
   gMATDown->SetLineStyle(2);
   gMATDown->SetLineColor(2);
   gMATDown->Draw("L");

   gPSUp->SetLineWidth(2);
   gPSUp->SetLineStyle(1);
   gPSUp->SetLineColor(6);
   gPSUp->Draw("L");

   gPSDown->SetLineWidth(2);
   gPSDown->SetLineStyle(2);
   gPSDown->SetLineColor(6);
   gPSDown->Draw("L");

   gLowUp->SetLineWidth(2);
   gLowUp->SetLineStyle(1);
   gLowUp->SetLineColor(1);
   gLowUp->Draw("L");

   gLowDown->SetLineWidth(2);
   gLowDown->SetLineStyle(2);
   gLowDown->SetLineColor(1);
   gLowDown->Draw("L");

   gNom->Draw("L");

   gZee->SetMarkerStyle(20);
   gZee->Draw("P");

   TLegend* leg0 = new TLegend(.45,.12,.65,.27);
   leg0->SetBorderSize(0);
   leg0->SetFillColor(0);
   leg0->AddEntry(gZeeUp,"#alpha, Method up", "l");
   leg0->AddEntry(gZeeDown,"#alpha, Method down", "l");

   leg0->AddEntry(gZeeGenUp,"#alpha, Model up", "l");
   leg0->AddEntry(gZeeGenDown,"#alpha, Model down", "l");

   TLegend* leg1 = new TLegend(.65,.12,.85,.27);
   leg1->SetBorderSize(0);
   leg1->SetFillColor(0);
   leg1->AddEntry(gMATUp,"MAT, up", "l");
   leg1->AddEntry(gMATDown,"MAT, down", "l");
   leg1->AddEntry(gPSUp,"PS, up", "l");
   leg1->AddEntry(gPSDown,"PS, down", "l");

   TLegend* leg2 = new TLegend(.25,.12,.45,.27);
   leg2->SetName("leg2");
   leg2->SetBorderSize(0);
   leg2->SetFillColor(0);

   leg2->AddEntry(gZee,"#alpha, Z peak", "p");
   leg2->AddEntry(gLowUp,"Low-pt, Up", "l");
   leg2->AddEntry(gLowDown,"Low-pt, Down", "l");
   if( subtractAverage ) {
     TText* txt = new TText();
     txt->SetNDC();
     txt->SetTextAlign(22);
     txt->SetTextSize(.03);
     txt->SetText(.5,.92,"Average scale subtracted");
     txt->Draw();
   }

   if( overlay ) {
     for(int i=0; i<overlay->GetN(); i++) {
       double x,y;
       overlay->GetPoint(i,x,y);
       overlay->SetPoint(i,x,y);
     }
     overlay->SetMarkerStyle(21);
     overlay->SetMarkerColor(1);
     overlay->SetMarkerSize(gZee->GetMarkerSize());
     overlay->Draw("psame");
     leg2->AddEntry(overlay,"linearity, E/p", "p");
   }

   leg0->Draw("same");
   leg1->Draw("same");
   leg2->Draw("same");


   // char plotname[99];
   // sprintf(plotname,"linearity_%d_%2.1f.png",ptype,eta);
   //   canvas->Print("overlay.ps");
  
   return canvas;
}


void plotloop(int ptype, string year) {
  
  char psnameFill[99];
  sprintf(psnameFill,"scaleLinearity_%d_%s.ps",ptype,year.c_str());
  
  char psnameOpen[99];
  sprintf(psnameOpen,"scaleLinearity_%d_%s.ps[",ptype,year.c_str());
  
  char psnameClose[99];
  sprintf(psnameClose,"scaleLinearity_%d_%s.ps]",ptype,year.c_str());
  
 TCanvas c1;
  c1.Print(psnameOpen);
  
  for(int i=0; i<48; i++) {

    double eta = (double)i/10. - 2.4 + 0.05;
    TCanvas* cnv = plotLinearity_ERUpgrade(eta,10.,100.,ptype, year);

    cnv->Print(psnameFill);

    delete cnv;
  }

  c1.Print(psnameClose);
}
