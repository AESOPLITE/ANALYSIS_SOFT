
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
#include "TCollection.h"
#include "TBox.h"


 using namespace std;
 //Maximum of configuration to try fitting
 
 int MAXB=1000;
 int MAXNB=1000;
 //Strip pitch in cm
 double  strippitch=0.0228;
//Detector resolution in cm
 double sigma_reso = strippitch/TMath::Sqrt(12);			
static const double kMelectron = 0.5109989461e-3;		//electron mass in GeV
static const double kMuon = 0.105658371;				//muon mass in GeV
 float zz0=0;//in cm
 
 TCanvas*can;
 TMultiGraph*multi;
 TMultiGraph*multiB;
 TF1*tmpf;
 TF1*tmpfB;
 TF1*tmpfB2;

int ALPatternRecognition::FindPattern(ALEvent *re, int DataType) {
cout << "Calling function FindPattern() " << endl;	
	
int nnhits = re->get_Nhits(); 
int i = re->get_eventnumber();
double E0 = re->get_EkMC();
TObjArray *xnonbend = new TObjArray();
TObjArray *xbend = new TObjArray();
TObjArray *allhits = new TObjArray(nnhits+1);
TObjArray *xtoplayer = new TObjArray();
TObjArray *xbottomlayer = new TObjArray();
TObjArray *xmid = new TObjArray();
TObjArray *xbend1 = new TObjArray();
TObjArray *xbend2 = new TObjArray();
TObjArray *xbend3 = new TObjArray();
TObjArray *xbend4 = new TObjArray();
vector<int> ibend1, ibend2, ibend3, ibend4;//keep track of index of hit
vector<int> itop, imid, ibottom;	//keep track of index of hit

 //Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
 //Load detector geometry parameters
 float*zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 int*valthres=new int[4];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 for(int i=0;i<4;i++)valthres[i]=0;
	
  

//MC or data file
 if(DataType==0) {//MC file	
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh);
 	}
 else {
 string DataFile="../src/ALSim/Dataparameters.dat";	 
 LoadDataparameters(DataFile,zL,OffsetLL,OffsetRL,valthres);
 }
 
 int*LayerWithHit= new int[7];
 for(int i=0;i<7;i++)LayerWithHit[i]=0;
 re->get_Layers(LayerWithHit);
 uint8_t Ti=(uint8_t)re->get_Ti();
 //Number of layers with hit(s)
 int NL= re->get_NLayers();
 //Number of layers with hit(s) in bending/non-bending plane
 int NLB = re->get_Layer(1) + re->get_Layer(2) + re->get_Layer(2) + re->get_Layer(5);
 int NLNB =  re->get_Layer(0) + re->get_Layer(4) + re->get_Layer(6);

  //if not less than 5 layers were touched then don't try pattern recognition
  if(NL<5) return 0;  


	for(int j=0;j<nnhits;j++) { 
		
      	Double_t xx=(re->get_hits().at(j))->get_x();
        Double_t yy=(re->get_hits().at(j))->get_y();
        Double_t zz=(re->get_hits().at(j))->get_z();
//Smear the coordinate to simulate the finite resolution of the detector, before discretization
	xx  += gRandom->Gaus(0., sigma_reso);   // smearing x
        yy += gRandom->Gaus(0., sigma_reso);   // smearing y
	int kk = re->get_hits().at(j)->get_k();
	TVector3 *X = new TVector3(xx,yy,zz);
	allhits->AddAt(X,kk);
 //add region of hit
	int L = re->get_hits().at(j)->get_L();
        cout << "Hit " << j << "  , layer" << L << ", index " << kk << "    x=" << xx << "  y=" << yy << "  z=" << zz << endl;

	if(L==0||L==4||L==6)//non-bending plane
         {
          xnonbend->Add(X);
          switch(L)
           {
            case 0: xtoplayer->Add(X);    itop.push_back(kk);    break;
            case 4: xmid->Add(X);         imid.push_back(kk);    break;
            case 6: xbottomlayer->Add(X); ibottom.push_back(kk); break;   
           }  
         }//if NB
       else//bending plane
         {
          xbend->Add(X);
          switch(L)
           {
            case 1: xbend1->Add(X);ibend1.push_back(kk);break;
            case 2: xbend2->Add(X);ibend2.push_back(kk);break;
            case 3: xbend3->Add(X);ibend3.push_back(kk);break;   
            case 5: xbend4->Add(X);ibend4.push_back(kk);break;   
           }  
         }//else       
      }  //j 
	
	//GraphAllPoints(re);
	
    //TEST ALL CONFIGURATIONS BENDING PLANE
    int npointsB, layer1, layer2, layer3, layer4, ncomb;
	layer1=layer2=layer3=layer4=1;
    npointsB = xbend->GetEntries();
    if(xbend1->GetEntries()!=0) layer1 = xbend1->GetEntries() ; 
    if(xbend2->GetEntries()!=0) layer2 = xbend2->GetEntries() ; 
    if(xbend3->GetEntries()!=0) layer3 = xbend3->GetEntries() ; 
    if(xbend4->GetEntries()!=0) layer4 = xbend4->GetEntries() ; 
   // cout << "entries in layer1 = " << xbend1->GetEntries() << "  , layer2 = " << xbend2->GetEntries() << "  , layer3 = " << xbend3->GetEntries() << "  , layer4= " << xbend4->GetEntries() << endl;
     ncomb   = layer1*layer2*layer3*layer4; 
	
// if less then 3 layers with hits in bending plane, can't reconstruct
	if(NLB < 3) {
	//cout << "NLB < 3, quitting function" << endl;
	return 0;
	}
	    if(ncomb<=0||ncomb>MAXB)
      {
         //  cout << " ncomb<=0||ncomb>MAXB " << endl;
	   return 0;
      }
	
    vector<double> chisquareB;
    TArrayI**indicesB = new TArrayI*[ncomb]; //array of TVector3 corresponding to "good" hits in B plane
    int*nhitB = new int[ncomb]; 
    for(int ijk=0;ijk<ncomb;ijk++)nhitB[ijk]=0;
    TGraphErrors**gBending= new TGraphErrors*[ncomb];
    TGraphErrors**gBendingInverted = new TGraphErrors*[ncomb];
    TF1**parabolas = new TF1*[ncomb];  //arrays of fit functions for B plane 
    TF1**invertedparabolas = new TF1*[ncomb];
	
    //fit all possible parabolas 
    int ncombination = 0;
	//cout << "BENDING PLANE " << endl;
    for(int m=0; m<layer1; m++)
      {   
       for(int n=0; n<layer2; n++)
         {  
          for(int o=0; o<layer3; o++) 
            {
             for(int p=0;p<layer4;p++) 
               {
                indicesB[ncombination] = new TArrayI(4);
                gBending[ncombination] = new TGraphErrors();
                gBendingInverted[ncombination] = new TGraphErrors();
                if(xbend1->GetEntries()!=0)
                 {
                  TVector3 *xx1 =(TVector3*)xbend1->At(m);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx1->Y(), xx1->Z());//first layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx1->Z(), xx1->Y());//first layer
	          nhitB[ncombination]++;
                 }
                if(xbend2->GetEntries()!=0)
                 {
                  TVector3 *xx2 =(TVector3*)xbend2->At(n);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx2->Y(), xx2->Z());//second layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx2->Z(), xx2->Y());//second layer
	          nhitB[ncombination]++;
                 }
                if(xbend3->GetEntries()!=0)
                 {
                  TVector3 *xx3 =(TVector3*)xbend3->At(o);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx3->Y(), xx3->Z());//third layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx3->Z(), xx3->Y());//third layer
	          nhitB[ncombination]++;
                 }
                if(xbend4->GetEntries()!=0)
                 {
                  TVector3 *xx4 =(TVector3*)xbend4->At(p);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx4->Y(), xx4->Z());//fourth layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx4->Z(), xx4->Y());//fourth layer
		  nhitB[ncombination]++;
                 }
 
                
                 for(int ijk=0;ijk<gBending[ncombination]->GetN();ijk++)gBending[ncombination]->SetPointError(ijk,strippitch/TMath::Sqrt(12.),0.);
                gBending[ncombination]->SetTitle(Form("Bending config %d", ncombination));
                gBending[ncombination]->SetMarkerStyle(kCircle);
                gBending[ncombination]->SetMarkerColor(kBlue);

                for(int ijk=0;ijk<gBendingInverted[ncombination]->GetN();ijk++)gBendingInverted[ncombination]->SetPointError(ijk,0.,strippitch/TMath::Sqrt(12.));
                gBendingInverted[ncombination]->SetMarkerStyle(kCircle);

                //fit function
                invertedparabolas[ncombination] = new TF1(Form("BConfig%d",ncombination),"[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,100);
                invertedparabolas[ncombination]->FixParameter(3,zz0); //Position of the center of the magnet in z coordinates, it is arbitrary as long as the same zz0 is used later in the analyis code 
                invertedparabolas[ncombination]->SetLineColor(kBlue);
                invertedparabolas[ncombination]->SetLineWidth(2);
                invertedparabolas[ncombination]->SetLineStyle(1);
                
                //FIT
                gBendingInverted[ncombination]->Fit(invertedparabolas[ncombination],"QSN");
                double chi2=invertedparabolas[ncombination]->GetChisquare();
                chisquareB.push_back(chi2);
		
		        //Order hits
				int point=0;
				if(xbend4->GetEntries()!=0) {			   
				  int k4 = ibend4.at(p);
				  indicesB[ncombination]->AddAt(k4,point);
				  //cout << "PR, bending plane, point added at index " << k4 << endl;
				  point++;
						 }
				if(xbend3->GetEntries()!=0) {
					int k3 = ibend3.at(o);
				 	indicesB[ncombination]->AddAt(k3,point);
					//cout << "PR, bending plane, point added at index " << k3 << endl;
					point++;
				    }
				if(xbend2->GetEntries()!=0) {
				 	int k2 = ibend2.at(n);
			     	        indicesB[ncombination]->AddAt(k2,point);
					//cout << "PR, bending plane, point added at index " << k2 << endl;
                                        point++;
					 }
				if(xbend1->GetEntries()!=0) {
				 	int k1 = ibend1.at(m);
			 	 	indicesB[ncombination]->AddAt(k1,point);
					//cout << "PR, bending plane, point added at index " << k1 << endl;
					  }
		      ncombination++;

               }//end p	 
            }// end o
		 }//end n
      } //end m
	
			  
    //TEST ALL CONFIGURATIONS NON BENDING PLANE
	int npointsNB, ntop, nmid, nbottom, NConf;
	ntop=nmid=nbottom=1;
    npointsNB = xnonbend->GetEntries();
    if(xtoplayer->GetEntries()!=0) ntop = xtoplayer->GetEntries(); //else cout << "non-bending plane top layer no hits " << endl;
    if(xmid->GetEntries()!=0) nmid = xmid->GetEntries(); //else cout << "non-bending plane middle layer no hits " << endl;
    if(xbottomlayer->GetEntries()!=0) nbottom = xbottomlayer->GetEntries(); //else cout << "non-bending plane botttom layer no hits " << endl;
   //cout << "entries in toplayer = " << ntop << "  , midlayer = " << nmid << "  , nbottom = " << nbottom << endl;
    //cout << "npointsNB = " << npointsNB << endl;
	NConf = nbottom*ntop*nmid; 
	
// if less then 2 hits in non-bending plane, can't reconstruct
	if(NLNB < 2)return 0;
	
    if(NConf<=0||NConf>MAXNB) return 0;
      
    
  
    //cout << "There are " << NConf << " possible combinations in the non-bending plane " << endl;
    vector<double> chisquare;      
    TArrayI**indicesNB = new TArrayI*[NConf];//array of TVector3 corresponding to "good" hits in NB plane
    int*nhitNB=new int[NConf];
    for(int ijk=0;ijk<NConf;ijk++)nhitNB[ijk]=0;
    TGraphErrors**gNonBending= new TGraphErrors*[NConf];
    TGraphErrors**gNonBending2= new TGraphErrors*[NConf];
    TF1**line = new TF1*[NConf];                  //arrays of fit functions for NB plane     

 //fit all possible lines 
    Int_t NConfig = 0;
    for(int m=0; m<ntop; m++) 
      {
       for(int n=0; n<nbottom; n++)
         { 
          for(int o=0; o<nmid; o++)
            {
             indicesNB[NConfig] = new TArrayI(3);
             gNonBending[NConfig] = new TGraphErrors();
             gNonBending2[NConfig] = new TGraphErrors();

             if(xtoplayer->GetEntries()!=0)
              {
               TVector3 *xx1 = (TVector3*)xtoplayer->At(m);
               gNonBending2[NConfig]->SetPoint(gNonBending2[NConfig]->GetN(), xx1->X(), xx1->Z());
               gNonBending[NConfig]->SetPoint(gNonBending[NConfig]->GetN(), xx1->Z(), xx1->X());
               nhitNB[NConfig]++;
              }
             if(xbottomlayer->GetEntries()!=0) 
              {
               TVector3 *xx2 = (TVector3*)xbottomlayer->At(n);
               gNonBending2[NConfig]->SetPoint(gNonBending2[NConfig]->GetN(), xx2->X(), xx2->Z());
               gNonBending[NConfig]->SetPoint(gNonBending[NConfig]->GetN(), xx2->Z(), xx2->X());
               nhitNB[NConfig]++;
              }
             if(xmid->GetEntries()!=0) 
              {
               TVector3 *xx3 = (TVector3*)xmid->At(o);
               gNonBending2[NConfig]->SetPoint(gNonBending2[NConfig]->GetN(), xx3->X(), xx3->Z());
               gNonBending[NConfig]->SetPoint(gNonBending[NConfig]->GetN(), xx3->Z(), xx3->X());
               nhitNB[NConfig]++;
              }
                
             //Same error on the measure X for every hits (may change in the future)   
             for(int ijk=0;ijk<gNonBending2[NConfig]->GetN();ijk++)gNonBending2[NConfig]->SetPointError(ijk,strippitch/TMath::Sqrt(12.),0);  
                
             for(int ijk=0;ijk<gNonBending[NConfig]->GetN();ijk++)gNonBending[NConfig]->SetPointError(ijk,0,strippitch/TMath::Sqrt(12.));
             
             //Display options
             gNonBending[NConfig]->SetTitle(Form("Non-bending config %d",  NConfig));
             gNonBending[NConfig]->SetMarkerStyle(kCircle);
             gNonBending[NConfig]->SetMarkerColor(kBlue);
             gNonBending2[NConfig]->SetTitle(Form("Non-bending config %d",  NConfig));
             gNonBending2[NConfig]->SetMarkerStyle(kCircle);
             gNonBending2[NConfig]->SetMarkerColor(kBlue);
             
             //fit function
             line[NConfig] =  new TF1(Form("NBConfig%d",NConfig),"pol1",-40,40);
             line[NConfig]->SetLineColor(kBlue);
             line[NConfig]->SetLineWidth(1);
             line[NConfig]->SetLineStyle(3);
             
             //FIT             
             gNonBending[NConfig]->Fit(line[NConfig],"QSN");  
             double chi2=line[NConfig]->GetChisquare();
             chisquare.push_back(chi2);

             //Order hits 
             int index=0;
             if(xmid->GetEntries()!=0)
              {
               int k3 = imid.at(o);
               indicesNB[NConfig]->AddAt(k3,index);
	       //cout << "PR, non-bending plane, point added at index " << k3 << endl;
               index++;
              }
             if(xbottomlayer->GetEntries()!=0)
              {
               int k2 = ibottom.at(n);
               indicesNB[NConfig]->AddAt(k2,index);
	       //cout << "PR, non-bending plane, point added at index " << k2 << endl;
               index++;
              }
             if(xtoplayer->GetEntries()!=0)
              {
               int k1 = itop.at(m);
               indicesNB[NConfig]->AddAt(k1,index);
	       //cout << "PR, non-bending plane, point added at index " << k1 << endl;
               index++;
              }
                        
             NConfig++;
            }//end o
         }//end n
      }//end m


    //GET BEST FIT PARAMETERS
    //BENDING
    int indexB = min_element(chisquareB.begin(), chisquareB.end()) - chisquareB.begin();      
    //NON-BENDING
    int index = min_element(chisquare.begin(), chisquare.end()) - chisquare.begin();
	
    Double_t a = invertedparabolas[indexB]->GetParameter(2);
    Double_t b = invertedparabolas[indexB]->GetParameter(1);
    Double_t c = invertedparabolas[indexB]->GetParameter(0);
 
    re->set_chi2BPR(chisquareB[indexB]);
    re->set_aPR(a);
    re->set_bPR(b);
    re->set_cPR(c);

    
    //Incoming Straight particle    
    float lim=TckZPos[1];//z position of 2nd layer
    float diff=2*a*lim+2*a*zz0+b;
    
    //Outcoming Straight particle    
    float limo=TckZPos[5];//z position of 6th layer
    float diffout=2*a*limo+2*a*zz0+b;
    
    //Set Deflection
    double deflection=diffout-diff;    
    re->set_deflecPR(deflection);

    //NON-BENDING
    //  cout << " index of min chisquare in NB plane i = " << index << endl; 	  
    Double_t p0 = line[index]->GetParameter(0);                                     
    Double_t p1 = line[index]->GetParameter(1);						//p1 = 1/tanl 
    TF1*inverseline=new TF1("inverseline","(x-[0])/([1])",-10,10);
    inverseline->SetParameters(p0,p1);
    
    //Fill PR coordinates
    re->set_interPR(p0);
    re->set_slopePR(p1);
    re->set_chi2NBPR(chisquare[index]);

//Interpolate bending plane hits in the non-bending plane

    for (int l=0; l<nhitB[indexB];l++)
      {
       int kk = indicesB[indexB]->At(l);
       re->get_hits().at(kk)->set_flagPR(true);
       //cout << "In PR, hit index " << kk << " set_flagPR " << endl;
       re->get_hits().at(kk)->set_yPR(invertedparabolas[indexB]->Eval(re->get_hits().at(kk)->get_z()));   
       re->get_hits().at(kk)->set_zPR(re->get_hits().at(kk)->get_z());
       //extrapolate the xPR from the fit in the non bending plane
       if(NLNB>=2)re->get_hits().at(kk)->set_xPR(line[index]->Eval(re->get_hits().at(kk)->get_z()));       
      }//l
	

    for(int l=0; l<nhitNB[index]; l++)
      {
       //index of chosen points given by array indicesNB
       int kk = indicesNB[index]->At(l);
       re->get_hits().at(kk)->set_flagPR(true);
      // cout << "In PR, hit index " << kk << " set_flagPR " << endl;
       if(p1!=0) 
	re->get_hits().at(kk)->set_zPR(re->get_hits().at(kk)->get_z());
        re->get_hits().at(kk)->set_xPR(p1*re->get_hits().at(kk)->get_z()+p0);
       //extrapolate the yPR from the fit in the bending plane
        if(NLB>=3)re->get_hits().at(kk)->set_yPR(invertedparabolas[indexB]->Eval(re->get_hits().at(kk)->get_z()));

     }
		
 //For events with <1 hit per layer: create ghost hit
    int kold = re->get_Nhits();
    vector<ALTckhit*> Hh;
	ALTckhit *h = new ALTckhit();
	for(int l=0;l<7;l++) {
	 	if(LayerWithHit[l]==0) {		//if there is a layer without a hit
			double z = TckZPos[l];
			double x = line[index]->Eval(z);
			double y = invertedparabolas[indexB]->Eval(z);
			h = new ALTckhit();
			h->set_xPR(x);
			h->set_yPR(y);
			h->set_zPR(z);
			//cout << "created ghost track for layer " << l << " zPR = " << z << endl;
			h->set_fGhost(true);
			h->set_flagPR(true);
			h->set_L(l); 
			h->set_k(kold);
			Hh.push_back(h);
			kold++;
		} //end if
	} //end for
//add missing hits to event
    for(int ij=0;ij<(int)Hh.size();ij++) re->add_hit(Hh.at(ij));


//Calculate and set momentum and total energy of particle from PR fit
	   double slopeNB = re->get_slopePR();
	   double thetaNB = TMath::ATan(slopeNB);
           double aPR = re->get_aPR(); 
	   double bPR = re->get_bPR();
	   double cPR = re->get_cPR();
	   //Signed curvature at the 4 points of the bending plane
	   double zzz[4]={TckZPos[5],TckZPos[3],TckZPos[2],TckZPos[1]};
	   double curv[4]={0,0,0,0};
	   TF1* fcurv=new TF1("fcurv","2*[0]/TMath::Power(1+TMath::Power(2*[0]*x+[1],2),3./2.)",-20,20);
	   fcurv->SetParameter(0,aPR);
	   fcurv->SetParameter(1,2*a*zz0+bPR);
	   double CurvMean=0; 	    
	   for(int ij=0;ij<4;ij++)
		{
         curv[ij]= fcurv->Eval(zzz[ij]);   
	    CurvMean+=	curv[ij]/4;
		}   
	   double Rmean=1./curv[2];
	   if(CurvMean!=0)	 Rmean=1./CurvMean;			 
	   //Extract a simple estimation of energy from the average of the 4 curvature a radius	 
	   double B=0.3; //in T	
	   double Pt=0.3 * B *	0.01*Rmean; //in GeV
	   double p0PR= Pt / TMath::Cos(fabs(thetaNB));   //in GeV
	   double mass;
	   int type;
	   double deflec = re->get_deflecPR();
	   double Q =  TMath::Sign(1,deflec);
	   if (DataType==1) { 
	   if (re->get_T2()) {			//if CK fired, electrons or positrons
	      if (Q>0) type=4;
	      else type=3;
	      }
	    else   {			//if CK not fired
	      if (Q>0) type=10;
	      else type=11;	       			
	    }
	  } //if data
	  //if MC	
	   else type = re->get_typeMC();	//if MC	
	   if (type==3 || type == 4) mass = kMelectron;	//for electrons in GeV
	   else if(type==10 || type ==11)  mass = kMuon;	//for muons in GeV
	   double EkPR=TMath::Sqrt(p0PR*p0PR+mass*mass);
	   re->set_EkPR(EkPR);
	   re->set_p0PR(p0PR);
	
//Calculate and set directional cosines for each hit	
	for(int j=0;j<re->get_Nhits();j++) { 		
		bool flagPR = re->get_hits().at(j)->get_flagPR();
		int k = re->get_hits().at(j)->get_k();
	        if(flagPR) {
		    float zPR =  re->get_hits().at(k)->get_zPR();
                    double slopeB = 2*aPR*zPR;
		    double thetaB = TMath::ATan(slopeB);
		    float cxPR = TMath::Cos(TMath::Pi()*0.5+thetaNB);
		    float cyPR = TMath::Cos(TMath::Pi()*0.5-thetaB);
		    float czPR = TMath::Cos(TMath::Pi()-thetaNB);
		    re->get_hits().at(k)->set_cxPR(cxPR);
		    re->get_hits().at(k)->set_cyPR(cyPR);
		    re->get_hits().at(k)->set_czPR(czPR);
		} //if flag
	}	//end for

return 1;
 
}
