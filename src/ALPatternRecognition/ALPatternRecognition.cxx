//*************************************************************************
//* ===================
//*  PatterRecognition
//* ===================
//*
//* (Description)
//* 16/03/2017  S. Mechbal		
//*
//* This is the main program that implements the pattern recognition routine. Hit
//* Hits are divided between bending and non-bending plane, and we implement two algorithms 
//* to choose the proper hits in the two fields of view.
//***************************************************************************

#include "ALPatternRecognition.h"

 using namespace std;

void ALPatternRecognition::FindPattern() 

{
  //Read ROOT MC output file
  
  bool flagMC = true;
  TFile*filein=new TFile("/home/sarah/AESOPLITE/ALanalysis-master/ALanalysis-master/RawEventMC.root","READ");
  TTree *tree;
  if(flagMC)tree = (TTree*)filein->Get("MC");
  
   //Define variables to read event
   ALEvent *e = new ALEvent;      			
 //Set address to access event data
  tree->SetBranchAddress("event",&e); 
 // Get number of event in Tree
 int nentries=tree->GetEntries();

 for (int i=0;i<300;i++) {
   tree->GetEntry(i);
   TObjArray *hits = new TObjArray;
	 
   //loop over hits in event
      int nnhits = (int)e->get_Nhits();
   if (nnhits > 6) {
  // cout << "Event " << i << "  Number of hits= " << nnhits << endl;
    
     for(int j=0;j<nnhits;j++) {
      	Float_t X=((e->get_hits().at(j))->get_xin()+(e->get_hits().at(j))->get_xout())/2;   //Mean of xin and xout
        Float_t Z=((e->get_hits().at(j))->get_zin()+(e->get_hits().at(j))->get_zout())/2;
        Float_t Y=((e->get_hits().at(j))->get_yin()+(e->get_hits().at(j))->get_yout())/2;
        TVector3 *xx = new TVector3;                      					//measurement vector
        xx->SetXYZ(X,Y,Z);	
		    hits->Add(xx);	
	  }
  
   
   TVector3 *xv;
   vector<double> xnonbend;
   vector<double> xbend;
   vector<double> ybend;
   vector<double> ynonbend;
   vector<double> zbend;
   vector<double> znonbend;
   vector<double> xtoplayer;
   vector<double> xbottomlayer;
   vector<double> ztoplayer;
   vector<double> zbottomlayer;
   Double_t xmid, zmid;
     
     
   for (int k=0; k<hits->GetEntries(); k++) {
        xv = (TVector3*)hits->At(k);
        Double_t xx = xv->X();
        Double_t yy = xv->Y();
        Double_t zz = xv->Z();
     
     //non-bending plane
     if ((zz > -2.0 && zz < -1.0) || (zz > -18.0 && zz < -17.0) || (zz > -22.0 && zz < -21.0)) {
        xnonbend.push_back(xx);
        znonbend.push_back(zz);
        if (zz > -2.0 && zz < -1.0) {
          xtoplayer.push_back(xx);
          ztoplayer.push_back(zz);
        }
       if (zz > -22.0 && zz < -21.0) {
         xbottomlayer.push_back(xx);
         zbottomlayer.push_back(zz);
       }
       if (zz > -18.0 && zz < -17.0) {
         xmid = xx;
         zmid = zz;
     }
     }
     else {
       ybend.push_back(yy);
       zbend.push_back(zz);       
     }
        //cout << "Original MC coordinates, event " << i << "  , hit #" << k << "    x=" << xx << "  y=" << yy << "  z=" << zz << endl;
   } 

     
     ///////////////////////////////////////////////////////////////////////////FIT AND GRAPH////////////////////////////////////////////////////////////////////////////////////////////
     
     //Make graphs in the non-bending plane of top/bottom layer
      Int_t ntop = ztoplayer.size();
      Int_t nbottom = zbottomlayer.size();
      Int_t xsize = xnonbend.size();                                                                                             //number of points in NB plane before Pattern Recognition
      Int_t NConf = nbottom*ntop; 
      vector<double> chisquare;        
         cout << "------------------------------------NEXT EVENT  " << i << "-----------------------------------" << endl;
      cout << "Number of possible configuration is " << NConf << endl;
      gROOT->SetStyle("Plain");
      gStyle->SetOptFit(0001);
			TCanvas *c1 = new TCanvas("fitted event", Form("Event %d", i),1200,1200);
      TGraph**gNonBending= new TGraph*[NConf];
      TF1**func = new TF1*[NConf];                  //arrays of fit functions
     	c1->Divide(NConf+1,2);		
	   
	   
	   
     //////////////////////////////////////////////////////BENDING PLANE//////////////////////////////////////
	   
	   
	   
      //draw event in bending plane, invert Y and Z axis to have abscissa go in increasing order for cubic spline interpolation
      TMultiGraph *mg1 = new TMultiGraph();
      mg1->SetTitle(Form("Bending plane, all points"));
      TGraph *gr1 = new TGraph(ybend.size(), &zbend[0], &ybend[0]);
   		gr1->SetTitle(Form("Event %d: Bending plane, inverted view", i));
   		gr1->GetXaxis()->SetTitle("Z (in cm)");
   		gr1->GetYaxis()->SetTitle("Y (in cm)");
			gr1->GetXaxis()->SetRangeUser(-30, 5);
			gr1->GetYaxis()->SetLimits(-20,20);
      gr1->SetMarkerStyle(kCircle);
      gr1->SetMarkerColor(kBlue);
			c1->cd(1);                            
      TSpline3 *spline3 = new TSpline3("Test", &zbend[0], &ybend[0], 4);
      spline3->SetLineColor(kGreen);
      gr1->Draw("AP");
      spline3->Draw("same");
     
      TGraph *gr2 = new TGraph(ybend.size(), &ybend[0], &zbend[0]);
   		gr2->SetMarkerStyle(kCircle);	
      gr2->SetMarkerColor(kBlue);	
      mg1->Add(gr2);
     
     //second graph on multigraph
     TGraph *gadded = new TGraph();
     //evaluate y position of hits in non-bending plane
                for (int l=0; l<znonbend.size();l++) {
             Double_t z = znonbend[l];
             double y = spline3->Eval(z);
             cout << "In non-bending plane, z = " << z << "   evaluated at z in bending plane y = " << y << endl;
             gadded->SetPoint(l, y, z);
           }

     
    gadded->SetMarkerStyle(kStar);
    gadded->SetMarkerColor(kRed);
    mg1->Add(gadded);
    c1->cd(2);   
    mg1->Draw("AP");
    c1->Modified(); 
    mg1->GetXaxis()->SetLimits(-20,20);
    mg1->GetYaxis()->SetRangeUser(-30, 5);
    mg1->GetXaxis()->SetTitle("Y (in cm)");
   	mg1->GetYaxis()->SetTitle("Z (in cm)");
    TLegend *leg = new TLegend(0.1,0.8,0.48,0.9);
    leg->AddEntry(gr2,"MC points", "p");
    leg->AddEntry(gadded, "Interpolated points", "p");
    leg->Draw("same");
    c1->Update();

     
     ////////////NON-BENDING PLANE/////////////////
      
     //draw all hits in non-bend
      TGraph *gr3 = new TGraph(xnonbend.size(), &xnonbend[0], &znonbend[0]);
      gr3->SetTitle(Form("Event %d: non-bending plane", i));
      gr3->GetXaxis()->SetTitle("X (in cm)");
   		gr3->GetYaxis()->SetTitle("Z (in cm)");
		  gr3->GetYaxis()->SetRangeUser(-30,5);
			gr3->GetXaxis()->SetLimits(-25,10);
   		gr3->SetMarkerStyle(kCircle);
      c1->cd(3);
   		gr3->Draw("AP");
     
     //fit all possible lines 
      int NConfig = 0;
        for(int m=0; m<ntop; m++) {
          for(int n=0; n<nbottom; n++) {     
          gNonBending[NConfig] = new TGraph(1);
          cout << "xtop = " << xtoplayer[m] << "   ztop = " << ztoplayer[m] << endl;
          cout << "xbottom = " << xbottomlayer[m] << "   zbottom = " << zbottomlayer[m] << endl;
          gNonBending[NConfig]->SetPoint(0, xtoplayer[m], ztoplayer[m]);
          gNonBending[NConfig]->SetPoint(1, xbottomlayer[n], zbottomlayer[n]);
         	gNonBending[NConfig]->SetTitle(Form("Non-bending config %d", NConfig));
          gNonBending[NConfig]->Fit("pol1");                                         //fit line between top and bottom layer
          func[NConfig] =  (TF1*)gNonBending[NConfig]->GetFunction("pol1");          //retrieve the fit function
          gNonBending[NConfig]->SetPoint(2, xmid, zmid);                             //add point from mid layer(have to loop over possible points if more than one)
          TFormula *formula = func[NConfig]->GetFormula();          
          Double_t p0 = formula->GetParameter(0);                                     
          Double_t p1 = formula->GetParameter(1);
          func[NConfig]->FixParameter(0, p0);                                        //fix fit parameters  
          func[NConfig]->FixParameter(1, p1);                                        //fix fit parameters
          gNonBending[NConfig]->Fit(func[NConfig]);                                  //fit fixed function with midpoint added and retrieve chisquare value of fit 
          func[NConfig]->SetLineColor(kBlue);
          Double_t chi2 = func[NConfig]->GetChisquare();
          chisquare.push_back(chi2);
          gNonBending[NConfig]->GetXaxis()->SetTitle("X (in cm)");
          gNonBending[NConfig]->GetYaxis()->SetTitle("Z (in cm)");
          gNonBending[NConfig]->GetYaxis()->SetRangeUser(-30,5);
          gNonBending[NConfig]->GetXaxis()->SetLimits(-25,10);
          gNonBending[NConfig]->SetMarkerStyle(kCircle);
          c1->cd(4+NConfig);
          gNonBending[NConfig]->Draw("AP");  
          NConfig++;
          }
      }

     
     
     //now retrieve good points 
          TMultiGraph *mg2 = new TMultiGraph();
          TGraph *gr4 = new TGraph(); 
          TGraph *gr5 = new TGraph();
          Int_t index = min_element(chisquare.begin(), chisquare.end()) - chisquare.begin();      //index of min. chisquare fit
          cout << " Index of minimum chi2 is i=" << index << endl;
          Int_t Npoints = gNonBending[index]->GetN();
          cout << "number of points in graph " << Npoints << endl;
          xnonbend.clear();                                                                          //clear arrays with old hits
          znonbend.clear();                                                                  //clear arrays with old hits
          for (int l=0; l<3; l++) {
          gNonBending[index]->GetPoint(l, xnonbend[l], znonbend[l]);
          Double_t x = xnonbend[l];                                                          //choose "good" hits
          Double_t z = znonbend[l];
          gr4->SetPoint(l, x, z);
          }
      gr4->SetMarkerStyle(kCircle);	
      gr4->SetMarkerColor(kBlue);	
      mg2->Add(gr4);
     
      //and do final fit

      gr4->Fit("pol1");
      TF1 *finalfit = gr4->GetFunction("pol1");                                               //save the final fit function
      finalfit->SetLineColor(kRed); 
      finalfit->SetLineStyle(3);
      finalfit->SetLineWidth(1);
        
     //evaluate x coordinate for bending plane point given z
           TFormula *fitformula = finalfit->GetFormula();
           Double_t p0 = fitformula->GetParameter(0);                                     
           Double_t p1 = fitformula->GetParameter(1);
           TFormula *inverseformula = new TFormula("inverse" , "(x-[0])/([1])", true);
           inverseformula->SetParameter(0,p0);
           inverseformula->SetParameter(1,p1);
           for (int l=0; l<zbend.size();l++) {
             Double_t z = zbend[l];
             double x = inverseformula->Eval(z);
             cout << "In bending plane, z = " << z << "   evaluated x = " << x << endl;
             gr5->SetPoint(l, x, z);                                                //add points from bending plane in the graph
           }
      gr5->SetMarkerStyle(kStar);
      gr5->SetMarkerColor(kRed);
      mg2->Add(gr5);
      c1->cd(NConfig+4);
      mg2->Draw("AP");
      c1->Modified();
      mg2->SetTitle(Form("Non-bending plane, all points"));
      mg2->GetXaxis()->SetTitle("X (in cm)");
   		mg2->GetYaxis()->SetTitle("Z (in cm)");
		  mg2->GetYaxis()->SetRangeUser(-30,5);
			mg2->GetXaxis()->SetLimits(-25,10);
      TLegend *leg1 = new TLegend(0.1,0.8,0.48,0.9);
      leg1->AddEntry(gr4,"MC points", "p");
      leg1->AddEntry(gr5, "Interpolated points", "p");
      leg1->Draw("same");
      c1->Update();
     
     
      c1->SaveAs(Form("./Event_%d.jpg", i));
      c1->Update();
     
   		}
	}
}
