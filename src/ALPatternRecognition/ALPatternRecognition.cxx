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

void ALPatternRecognition::FindPattern(ALEvent *re, int* TckReg) {
	

   int nnhits = re->get_Nhits(); 
   int i = re->get_eventnumber();
   double E0 = re->get_EkMC();
   vector<double> xnonbend, ynonbend, znonbend;					//NB
   vector<double> xbend, ybend, zbend;							//B
   vector<double> xtoplayer, ytoplayer, ztoplayer;				//NB
   vector<double> xbottomlayer, ybottomlayer, zbottomlayer;		//NB
   Double_t xmid, ymid, zmid;									//NB
   vector<double> xbend1, ybend1, zbend1;						//B
   vector<double> xbend2, ybend2, zbend2;						//B
   vector<double> xbend3, ybend3, zbend3;						//B
   vector<double> xbend4, ybend4, zbend4;						//B

   bool testlayer[7]={false,false,false,false,false,false,false};
     
cout << "Event " << i << "  Number of hits= " << nnhits << endl;
	
//FLUKA coordinates used in PatternRecognition
	for(int j=0;j<nnhits;j++) { 
		
      	Double_t xx=((re->get_hits().at(j))->get_xin()+(re->get_hits().at(j))->get_xout())/2; 
        Double_t yy=((re->get_hits().at(j))->get_yin()+(re->get_hits().at(j))->get_yout())/2;
        Double_t zz=((re->get_hits().at(j))->get_zin()+(re->get_hits().at(j))->get_zout())/2;

     //non-bending plane
     if ((zz > -2.0 && zz < -1.0) || (zz > -18.0 && zz < -17.0) || (zz > -22.0 && zz < -21.0)) {
        xnonbend.push_back(xx);
        znonbend.push_back(zz);
	//record the extra MC coordinate in both views
		ynonbend.push_back(yy);
        if (zz > -2.0 && zz < -1.0) {
          xtoplayer.push_back(xx);
          ztoplayer.push_back(zz);
		  ytoplayer.push_back(yy);
        }
       if (zz > -22.0 && zz < -21.0) {
         xbottomlayer.push_back(xx);
         zbottomlayer.push_back(zz);
		 ybottomlayer.push_back(yy);
       }
       if (zz > -18.0 && zz < -17.0) {
         xmid = xx;
         zmid = zz;
		 ymid = yy;
    	 } 			//end if
     }				//end if NB plane
     else {
       ybend.push_back(yy);
       zbend.push_back(zz);    
	//record the extra MC coordinate
	   xbend.push_back(xx);
        if (zz > -4.0 && zz < -3.0) {
			ybend1.push_back(yy);
			zbend1.push_back(zz);
			xbend1.push_back(xx);
    	 }
        if (zz > -10.0 && zz < -9.0) {
			ybend2.push_back(yy);
			zbend2.push_back(zz);
			xbend2.push_back(xx);    	 
		}
        if (zz > -16.0 && zz < -15.0) {
			ybend3.push_back(yy);
			zbend3.push_back(zz);
			xbend3.push_back(xx);
    	 }
        if (zz > -20.0 && zz < -19.0) {
			ybend4.push_back(yy);
			zbend4.push_back(zz);
			xbend4.push_back(xx);
		  }					//end if
  	 } 						//end else B plane
  }  							//end j  
 
///IF EVENT DOES NOT HAVE 7 HITS, JUST PLOT IT
	
	      gROOT->SetStyle("Plain");
          gStyle->SetOptFit(0001);
     	  TCanvas *c1 = new TCanvas("fitted event", Form("Event %d", i),1200,1200);

    if (re->get_Ti()==127) {
           
	///////////////////////////////////////////////////TEST ALL CONFIGURATIONS BENDING PLANE/////////////////////////////////////////
	 Int_t npointsB = ybend.size();
	 Int_t nlayer1 = ybend1.size();
	 Int_t nlayer2 = ybend2.size();
	 Int_t nlayer3 = ybend3.size();
	 Int_t nlayer4 = ybend4.size();
//cout << "There are " << npointsB << " points in the bending plane " << endl;
//cout << "Layer 1 = " << nlayer1 << " points " << endl;
//cout << "Layer 2 = " << nlayer2 << " points " << endl;
//cout << "Layer 3 = " << nlayer3 << " points " << endl;
//cout << "Layer 4 = " << nlayer4 << " points " << endl;
	 Int_t  ncomb   = nlayer1*nlayer2*nlayer3*nlayer4;
	 vector<double> chisquareB;
     TObjArray**xvalues = new TObjArray*[ncomb];		 //array of TVector3 corresponding to "good" hits in B plane
	 TGraph**gBending= new TGraph*[ncomb];
	 TGraph**gBendingInverted = new TGraph*[ncomb];
     TF1**parabolas = new TF1*[ncomb];                  //arrays of fit functions for B plane 
	 TF1**invertedparabolas = new TF1*[ncomb];

	   
 //fit all possible parabolas 
      int ncombination = 0;
        for(int m=0; m<nlayer1; m++) {
          for(int n=0; n<nlayer2; n++) {   
			  for(int o=0; o<nlayer3; o++) {
				  for (int p=0; p<nlayer4; p++) {
		  vector<double> y;
		  TVector3 *xx1 = new TVector3(xbend1[m], ybend1[m], zbend1[m]);
		  TVector3 *xx2 = new TVector3(xbend2[n], ybend2[n], zbend2[n]);	
		  TVector3 *xx3 = new TVector3(xbend3[o], ybend3[o], zbend3[o]);
		  TVector3 *xx4 = new TVector3(xbend4[p], ybend4[p], zbend4[p]);
		  y.push_back(ybend1[m]);
          y.push_back(ybend2[n]);
		  y.push_back(ybend3[o]);
		  y.push_back(ybend4[p]);
		  xvalues[ncombination] = new TObjArray();
		  xvalues[ncombination]->Add(xx1);
	      xvalues[ncombination]->Add(xx2);
          xvalues[ncombination]->Add(xx3);
		  xvalues[ncombination]->Add(xx4);
		  Int_t imin = min_element(y.begin(), y.end()) - y.begin();
	      Int_t imax = max_element(y.begin(), y.end()) - y.begin();
		  Double_t ymax = y[imax];
	      Double_t ymin = y[imin];
          gBending[ncombination] = new TGraph();
          gBending[ncombination]->SetPoint(0, ybend1[m], zbend1[m]);							//first layer
          gBending[ncombination]->SetPoint(1, ybend2[n], zbend2[n]);							//second layer
          gBending[ncombination]->SetPoint(2, ybend3[o], zbend3[o]);							//third layer
          gBending[ncombination]->SetPoint(3, ybend4[p], zbend4[p]);					 		//fourth layer		
          gBending[ncombination]->SetTitle(Form("Bending config %d", ncombination));
		  gBending[ncombination]->SetMarkerStyle(kCircle);
		  gBending[ncombination]->SetMarkerColor(kBlue);
		  parabolas[ncombination] = new TF1("parabola", "pol2", ymin - 3 , ymax +3);
		  parabolas[ncombination]->SetLineColor(kBlue);
		  parabolas[ncombination]->SetLineWidth(1);
		  parabolas[ncombination]->SetLineStyle(3);
		  gBending[ncombination]->Fit("parabola", "R");
          gBendingInverted[ncombination] = new TGraph();
          gBendingInverted[ncombination]->SetPoint(0, zbend1[m], ybend1[m]);							//first layer
          gBendingInverted[ncombination]->SetPoint(1, zbend2[n], ybend2[n]);							//second layer
          gBendingInverted[ncombination]->SetPoint(2, zbend3[o], ybend3[o]);							//third layer
          gBendingInverted[ncombination]->SetPoint(3, zbend4[p], ybend4[p]);
	      invertedparabolas[ncombination] = new TF1("invertedparabola", "pol2");
		  gBendingInverted[ncombination]->Fit("invertedparabola");
		  Double_t chi2 = parabolas[ncombination]->GetChisquare();
		  chisquareB.push_back(chi2);
          //c1->cd(5+ncombination);
          //gBending[ncombination]->Draw("AP");  
          ncombination++;
          }	//end p
      }	// end o
	} //end n
} //end m
	
			  
///////////////////////////////////////////////////TEST ALL CONFIGURATIONS NB PLANE//////////////////////////////////////////////
	  Int_t npointsNB = xnonbend.size();
      Int_t ntop = ztoplayer.size();
      Int_t nbottom = zbottomlayer.size();
      Int_t NConf = nbottom*ntop; 
//cout << "There are " << npointsNB << " points in the non-bending plane " << endl;
//cout << "Top Layer = " << ntop << " points " << endl;
//cout << "Bottom Layer = " << nbottom << " points " << endl;
		
      vector<double> chisquare;        
	  TObjArray**yvalues = new TObjArray*[NConf];					//array of TVector3 corresponding to "good" hits in NB plane
      TGraph**gNonBending= new TGraph*[NConf];
      TF1**line = new TF1*[NConf];                  //arrays of fit functions for NB plane     
	
 //fit all possible lines 
	
   Int_t NConfig = 0;
   for(int m=0; m<ntop; m++) {
          for(int n=0; n<nbottom; n++) {     
		  TVector3 *xx1 = new TVector3(xtoplayer[m], ytoplayer[m], ztoplayer[m]);
		  TVector3 *xx2 = new TVector3(xmid, ymid, zmid);	
		  TVector3 *xx3 = new TVector3(xbottomlayer[n], ybottomlayer[n], zbottomlayer[n]);
		  yvalues[NConfig] = new TObjArray();
		  yvalues[NConfig]->Add(xx1);
	      yvalues[NConfig]->Add(xx2);
		  yvalues[NConfig]->Add(xx3);
          gNonBending[NConfig] = new TGraph(1);
          gNonBending[NConfig]->SetPoint(0, xtoplayer[m], ztoplayer[m]);
          gNonBending[NConfig]->SetPoint(1, xbottomlayer[n], zbottomlayer[n]);
          gNonBending[NConfig]->SetTitle(Form("Non-bending config %d", NConfig));
          gNonBending[NConfig]->Fit("pol1");                                         //fit line between top and bottom layer
          line[NConfig] =  (TF1*)gNonBending[NConfig]->GetFunction("pol1");          //retrieve the fit function
          gNonBending[NConfig]->SetPoint(2, xmid, zmid);                             //add point from mid layer(have to loop over possible points if more than one)
          TFormula *formula = line[NConfig]->GetFormula();          
          Double_t p0 = formula->GetParameter(0);                                     
          Double_t p1 = formula->GetParameter(1);
          line[NConfig]->FixParameter(0, p0);                                        //fix fit parameters  
          line[NConfig]->FixParameter(1, p1);                                        //fix fit parameters
          gNonBending[NConfig]->Fit(line[NConfig]);                                  //fit fixed function with midpoint added and retrieve chisquare value of fit 
          line[NConfig]->SetLineColor(kBlue);
          line[NConfig]->SetLineWidth(1);
          line[NConfig]->SetLineStyle(3);
          Double_t chi2 = line[NConfig]->GetChisquare();
          chisquare.push_back(chi2);
          gNonBending[NConfig]->GetXaxis()->SetTitle("X (in cm)");
          gNonBending[NConfig]->GetYaxis()->SetTitle("Z (in cm)");
          gNonBending[NConfig]->GetYaxis()->SetRangeUser(-30,5);
          gNonBending[NConfig]->GetXaxis()->SetLimits(-25,10);
          gNonBending[NConfig]->SetMarkerStyle(kCircle);
          //c1->cd(4+NConfig+ncombination);
          //gNonBending[NConfig]->Draw("AP");  
          NConfig++;
          }
      }	

	
///////////////////////////////////////////////////////////GRAPH//////////////////////////////////////////////////////
	      c1->Divide(2,2);		

	      TMultiGraph *mg1 = new TMultiGraph();
          TMultiGraph *mg2 = new TMultiGraph();
	      TGraph *gr1 = new TGraph(ybend.size(), &ybend[0], &zbend[0]);						//all MC hits, bending plane
          TGraph *gr3 = new TGraph(xnonbend.size(), &xnonbend[0], &znonbend[0]);			//all MC hits, NB plane
          TGraph *gInterpolateB = new TGraph();												//interpolated B hits in NB view
	      TGraph *gMCbending = new TGraph(xbend.size(), &xbend[0], &zbend[0]);				//MC B hits NB view 
	      TGraph *gInterpolateNB = new TGraph();											//interpolated NB hits in NB view
	      TGraph *gfinalNB = new TGraph();													//final chosen NB points
		  TGraph *gfinalB = new TGraph();													//final chosen B points 
	      TGraph *gMCnonbending = new TGraph(xnonbend.size(), &ynonbend[0], &znonbend[0]);	//all MC NB points in B view 		

//GET BEST FIT PARAMETERS

          int index = min_element(chisquare.begin(), chisquare.end()) - chisquare.begin();      //index of min. chisquare fit
	      //cout << " index of min chisquare in NB plane i = " << index << endl; 	  
           Double_t p0 = line[index]->GetParameter(0);                                     
           Double_t p1 = line[index]->GetParameter(1);						//p1 = 1/tanl 
		   Double_t tanl = 1/p1;
		   Double_t lambda = atan(tanl);
	       Double_t cosl = cos(lambda);
//invert formula to interpolate missing coordinate 		 
           TFormula *inverseformula = new TFormula("inverse" , "(x-[0])/([1])", true);
           inverseformula->SetParameter(0,p0);
           inverseformula->SetParameter(1,p1);
	     
		 int indexB = min_element(chisquareB.begin(), chisquareB.end()) - chisquareB.begin();      //index of min. chisquare fit         
	      //cout << " index of min chisquare in BENDING plane i = " << indexB << endl; 
         Double_t a = parabolas[indexB]->GetParameter(0);
		 Double_t b = parabolas[indexB]->GetParameter(1);
		 Double_t c = parabolas[indexB]->GetParameter(2);

	
//INTERPOLATE NB POINTS IN THE BENDING PLANE (FIND Y GIVEN Z)
	      for(int l=0; l<3; l++) {
	  	  TVector3 *xNonBending =(TVector3*)yvalues[index]->At(l);
		  double yMC = xNonBending->Y();
          double zMC = xNonBending->Z();
		  gfinalNB->SetPoint(l, yMC, zMC);									
		  double yInter = invertedparabolas[indexB]->Eval(zMC);
		  gInterpolateNB->SetPoint(l, yInter, zMC);								 
		  }
//INTERPOLATE B POINTS IN THE NB PLANE (FIND X GIVEN Z)
	     for (int l=0; l<4;l++) {
		    TVector3 *xx = (TVector3*)xvalues[indexB]->At(l);
			double xMC = xx->X();
		  	double zMC = xx->Z();
		  	gfinalB->SetPoint(l, xMC, zMC);
			Double_t xInter = inverseformula->Eval(zMC);
            gInterpolateB->SetPoint(l, xInter, zMC);                                               		
		   }
	
	
      gr1->SetTitle(Form("Event %d: bending plane", i));
      gr1->GetXaxis()->SetTitle("Y (in cm)");
   	  gr1->GetYaxis()->SetTitle("Z (in cm)");
	  gr1->GetYaxis()->SetRangeUser(-30,5);
	  gr1->GetXaxis()->SetLimits(-25,10);
   	  gr1->SetMarkerStyle(kCircle);
      c1->cd(1);
   	  gr1->Draw("AP");
	  c1->Update();
	
	  gfinalB->SetMarkerStyle(kCircle);
	  gfinalB->SetMarkerColor(kBlack);
	  gInterpolateNB->SetMarkerStyle(kStar);
	  gInterpolateNB->SetMarkerColor(kRed);
	  gMCnonbending->SetMarkerStyle(kCircle);
	  gMCnonbending->SetMarkerColor(kGreen);
	  gfinalNB->SetMarkerStyle(kCircle);
	  gfinalNB->SetMarkerColor(kBlack);
//add graphs to multigraph
	  mg1->Add(gBending[indexB]);
	  mg1->Add(gInterpolateNB);
	  mg1->Add(gMCnonbending);
	  mg1->Add(gfinalNB);
	  c1->cd(2);	
	  mg1->Draw("AP");
      mg1->SetTitle(Form("Bending plane, all points"));
	  mg1->GetXaxis()->SetLimits(-20,20);
      mg1->GetYaxis()->SetRangeUser(-30, 5);
	  mg1->GetXaxis()->SetTitle("Y (in cm)");
      mg1->GetYaxis()->SetTitle("Z (in cm)");
      TLegend *leg = new TLegend(0.1,0.8,0.48,0.9);
	  leg->AddEntry(gBending[indexB], "Chosen B points", "p");
      leg->AddEntry(gMCnonbending,"All NB MC points", "p");
      leg->AddEntry(gInterpolateNB, "NB Interpolated points", "p");
	  leg->AddEntry(gfinalNB, "Chosen NB points", "p");
      leg->Draw("same");
      c1->Update();
	
//NB PLANE	
	
	
      gr3->SetTitle(Form("Event %d: non-bending plane", i));
      gr3->GetXaxis()->SetTitle("X (in cm)");
   	  gr3->GetYaxis()->SetTitle("Z (in cm)");
	  gr3->GetYaxis()->SetRangeUser(-30,5);
   	  gr3->GetXaxis()->SetLimits(-25,10);
   	  gr3->SetMarkerStyle(kCircle);
      c1->cd(3);
   	  gr3->Draw("AP");
	

      gInterpolateB->SetMarkerStyle(kStar);
      gInterpolateB->SetMarkerColor(kRed);
	  gMCbending->SetMarkerStyle(kCircle);
	  gMCbending->SetMarkerColor(kGreen);
	  mg2->Add(gNonBending[index]);
      mg2->Add(gInterpolateB);
	  mg2->Add(gMCbending);
      mg2->Add(gfinalB);
      c1->cd(4);
      mg2->Draw("AP");
      c1->Modified();
      mg2->SetTitle(Form("Non-bending plane, all points"));
      mg2->GetXaxis()->SetTitle("X (in cm)");
   	  mg2->GetYaxis()->SetTitle("Z (in cm)");
	  mg2->GetYaxis()->SetRangeUser(-30,5);
	  mg2->GetXaxis()->SetLimits(-25,10);
      TLegend *leg1 = new TLegend(0.1,0.8,0.48,0.9);
      leg1->AddEntry(gNonBending[index],"Chosen NB points", "p");
      leg1->AddEntry(gInterpolateB, " B Interpolated  points", "p");
	  leg1->AddEntry(gMCbending, "B plane MC points" , "p");
	  leg1->AddEntry(gfinalB, "Chosen B points" , "p");
      leg1->Draw("same");
      c1->Update();
	
// QUICK AND DIRTY LEAST SQUARES FIT 
     Double_t R = 1/(2*c);							  //radius of curvature in [m]
	 Double_t B = 0.3;								     //magnetic field in [T]
     Double_t pt = fabs(0.3*B*R);						//transverse momentum in [GeV]
	 Double_t p = pt*sqrt(1+(tanl*tanl));				//momentum in [GeV]
	 cout << "Incoming energy of the particle is " << E0*1000 << " MeV " << endl;
     cout << "Reconstructed momentum = " << p*1000 << " MeV " << endl;
	}
	
///PLOT EVENTS WITH LESS THAN ONE HIT PER LAYER	
	else {
	
	  TGraph *gr1 = new TGraph(ybend.size(), &ybend[0], &zbend[0]);						//all MC hits, bending plane
      TGraph *gr3 = new TGraph(xnonbend.size(), &xnonbend[0], &znonbend[0]);			//all MC hits, NB plane
	  c1->Divide(1,2);
      gr1->SetTitle(Form("Event %d: bending plane", i));
      gr1->GetXaxis()->SetTitle("Y (in cm)");
   	  gr1->GetYaxis()->SetTitle("Z (in cm)");
	  gr1->GetYaxis()->SetRangeUser(-30,5);
   	  gr1->GetXaxis()->SetLimits(-25,10);
   	  gr1->SetMarkerStyle(kCircle);
      c1->cd(1);
   	  gr1->Draw("AP");
		
      gr3->SetTitle(Form("Event %d: non-bending plane", i));
      gr3->GetXaxis()->SetTitle("X (in cm)");
   	  gr3->GetYaxis()->SetTitle("Z (in cm)");
	  gr3->GetYaxis()->SetRangeUser(-30,5);
   	  gr3->GetXaxis()->SetLimits(-25,10);
   	  gr3->SetMarkerStyle(kCircle);
      c1->cd(2);
   	  gr3->Draw("AP");
	}
     

      c1->SaveAs(Form("./PatternReco_Event_%d.jpg", i));
      c1->Update();
}

