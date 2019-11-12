
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
static const double kProton = 0.93827231;				//proton mass in GeV
static const double kAlpha = 3.727379378;				//alpha mass in GeV
 float zz0=0;//in cm

 TCanvas*can;
 TMultiGraph*multi;
 TMultiGraph*multiB;
 TF1*tmpf;
 TF1*tmpfB;
 TF1*tmpfB2;

int ALPatternRecognition::FindPattern(ALEvent *re, int DataType, float* zL,float* OffsetLL,float* OffsetRL,float* TrigThresh) {

//cout << "Calling function FindPattern() " << endl;
gROOT->Reset();
int nnhits = re->get_Nhits();
//int i = re->get_eventnumber();
double E0 = re->get_EkMC();
TObjArray *xnonbend = new TObjArray();
TObjArray *xbend = new TObjArray();
//TObjArray *allhits = new TObjArray(nnhits+1);
TObjArray *xtoplayer = new TObjArray();
TObjArray *xbottomlayer = new TObjArray();
TObjArray *xmid = new TObjArray();
TObjArray *xbend1 = new TObjArray();
TObjArray *xbend2 = new TObjArray();
TObjArray *xbend3 = new TObjArray();
TObjArray *xbend4 = new TObjArray();
vector<int> ibend1, ibend2, ibend3, ibend4;//keep track of index of hit
vector<int> itop, imid, ibottom;	//keep track of index of hit
vector<int> nstrip1, nstrip2, nstrip3,nstrip4;   //keep track of number of strips that make up a hit
vector<int> nstriptop, nstripmid,nstripbottom;	 //keep  track of number of strips of make up a hit
vector<int> L6S575;

//Number of layers with hit(s) with max 3 strips
int NL= re->get_NLayersc();
//Number of layers with hit(s) in bending/non-bending plane
//int NLB = re->get_Layer(1) + re->get_Layer(2) + re->get_Layer(3) + re->get_Layer(5);
//int NLNB =  re->get_Layer(0) + re->get_Layer(4) + re->get_Layer(6);
bool NL0 = false;
bool NL1 = false;
bool NL2 = false;
bool NL3 = false;
bool NL4 = false;
bool NL5 = false;
bool NL6 = false;

int* NLPRtmp= new int[7];
for(int i=0;i<7;i++)NLPRtmp[i]=0;
//cout << "NL = " << NL << "  NLB = " << NLB << "  NLNB = " << NLNB << endl;
//if not less than 5 layers were touched then don't try pattern recognition
if(NL<5)
 {
  return 0;
 }

//Check if the strip 63, chip 8, L6 has signal
//Calculate thenumber of hits in the chip



//cout << "nnhits = " << nnhits << endl;
for(int j=0;j<nnhits;j++)
  {
   int L = re->get_hits().at(j)->get_L();
   int kk = re->get_hits().at(j)->get_k();
   int nstrip = re->get_hits().at(j)->get_nstrips()+DataType; //DataType==0 for MC and 1 for Data
   int fstripID=re->get_hits().at(j)->get_fstripID();

   Double_t xx=(re->get_hits().at(j))->get_x();
   Double_t yy=(re->get_hits().at(j))->get_y();
   Double_t zz=(re->get_hits().at(j))->get_z();


   if (L==6 && fstripID==575 &&nstrip==1)
	  {
	   xx = xx - 32.5*strippitch; //For flight 2018, put the hit in middle of chip 8
	  }

   TVector3 *X = new TVector3(xx,yy,zz);

   if(nstrip>3) continue;

	 if(L==0||L==4||L==6)//non-bending plane
    {
     xnonbend->Add(X);
     switch(L)
           {
            case 0: if (abs(re->get_hits().at(j)->get_x())<4.) {xtoplayer->Add(X);itop.push_back(kk); nstriptop.push_back(nstrip);NL0=true;NLPRtmp[L]++;}break;
            case 4: if (abs(re->get_hits().at(j)->get_x())<10.) {xmid->Add(X);imid.push_back(kk);  nstripmid.push_back(nstrip);NL4=true;NLPRtmp[L]++;}break;
            case 6: if (abs(re->get_hits().at(j)->get_x())<10.){xbottomlayer->Add(X);ibottom.push_back(kk); nstripbottom.push_back(nstrip);NL6=true;NLPRtmp[L]++;}
                    if (fstripID==575 && nstrip==1) L6S575.push_back(1.);//Hits with single strip L6S575
                    else L6S575.push_back(0.);
                    break;
           }//end switch
    }//if NB
   else//bending plane
       {
        xbend->Add(X);
        switch(L)
           {
            case 1: if (abs(re->get_hits().at(j)->get_y())<6.) {xbend1->Add(X);ibend1.push_back(kk);nstrip1.push_back(nstrip);NL1=true;NLPRtmp[L]++;}break;
            case 2: if (abs(re->get_hits().at(j)->get_y())<6.25) {xbend2->Add(X);ibend2.push_back(kk);nstrip2.push_back(nstrip);NL2=true;NLPRtmp[L]++;}break;
            case 3: if (abs(re->get_hits().at(j)->get_y())<10) {xbend3->Add(X);ibend3.push_back(kk);nstrip3.push_back(nstrip);NL3=true;NLPRtmp[L]++;}break;
            case 5: if (abs(re->get_hits().at(j)->get_y())<10) {xbend4->Add(X);ibend4.push_back(kk);nstrip4.push_back(nstrip);NL5=true;NLPRtmp[L]++;}break;
           }
        }//else
    // cout << "Hit " << j << "  , L=" << L << ", fstripID " << fstripID << "  nstrip =  " << nstrip << "    x=" << xx << "  y=" << yy << "  z=" << zz << endl;

 }  //j

 int TiPRtmp=0;
 for(int i=0;i<7;i++)
   {
    re->set_NLPR(i,NLPRtmp[i]);
    if(NLPRtmp[i]>0)TiPRtmp+=(int)TMath::Power(2,i);
   }
 re->set_TiPR(TiPRtmp);

	//GraphAllPoints(re);

    //TEST ALL CONFIGURATIONS BENDING PLANE
 int npointsB, layer1, layer2, layer3, layer4, ncomb;
 layer1=layer2=layer3=layer4=1;
 npointsB = xbend->GetEntries();
 int NLB = (int)NL1 + (int)NL2 + (int)NL3 +(int)NL5;
 if((xbend1->GetEntries())!=0) layer1 = xbend1->GetEntries() ;
 if((xbend2->GetEntries())!=0) layer2 = xbend2->GetEntries() ;
 if((xbend3->GetEntries())!=0) layer3 = xbend3->GetEntries() ;
 if((xbend4->GetEntries())!=0) layer4 = xbend4->GetEntries() ;
 ncomb = layer1*layer2*layer3*layer4;
// cout << " npointsB = " << npointsB << ", NLB = " << NLB << endl;
 // if less then 3 layers with hits in bending plane, can't reconstruct
	if( NLB <3)
   {
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
  double sigma1,sigma2,sigma3,sigma4;

  //fit all possible parabolas
  int ncombination = 0;
  //	cout << "BENDING PLANE " << endl;
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
		              //  cout << " xbend1 " << endl;
                  TVector3 *xx1 =(TVector3*)xbend1->At(m);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx1->Y(), xx1->Z());//first layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx1->Z(), xx1->Y());//first layer
                  if(nstrip1.at(m)==3) sigma1= 2*sigma_reso;
                  else if(nstrip1.at(m)==2) sigma1= sigma_reso/(TMath::Sqrt(2.));
		              else sigma1=sigma_reso;
		              //  cout << "nstrip1 = " << nstrip1.at(m) << " sigma1 = " << sigma1 << endl;
		              gBending[ncombination]->SetPointError(gBending[ncombination]->GetN()-1,sigma1,0.);
		              gBendingInverted[ncombination]->SetPointError(gBendingInverted[ncombination]->GetN()-1,0.,sigma1);
		              nhitB[ncombination]++;
                 }
                if(xbend2->GetEntries()!=0)
                 {
		              //cout << "xbend2 " << endl;
                  TVector3 *xx2 =(TVector3*)xbend2->At(n);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx2->Y(), xx2->Z());//second layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx2->Z(), xx2->Y());//second layer
                  if(nstrip2.at(n)==3) sigma2= 2*sigma_reso;
                  else if(nstrip2.at(n)==2) sigma2= sigma_reso/(TMath::Sqrt(2.));
		              else sigma2=sigma_reso;
		              //  cout << "nstrip2 = " << nstrip2.at(n) << " sigma2 = " << sigma2 << endl;
		              gBending[ncombination]->SetPointError(gBending[ncombination]->GetN()-1,sigma2,0.);
		              gBendingInverted[ncombination]->SetPointError(gBendingInverted[ncombination]->GetN()-1,0.,sigma2);
		              nhitB[ncombination]++;
                 }

                if(xbend3->GetEntries()!=0)
                 {
		              // cout << "xbend3" << endl;
                  TVector3 *xx3 =(TVector3*)xbend3->At(o);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx3->Y(), xx3->Z());//third layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx3->Z(), xx3->Y());//third layer
                  if(nstrip3.at(o)==3) sigma3= 2*sigma_reso;
                  else if(nstrip3.at(o)==2) sigma3= sigma_reso/(TMath::Sqrt(2.));
		              else sigma3=sigma_reso;
		              //cout << "nstrip3 = " << nstrip3.at(o) << " sigma3 = " << sigma3 << endl;
		              gBending[ncombination]->SetPointError(gBending[ncombination]->GetN()-1,sigma3,0.);
		              gBendingInverted[ncombination]->SetPointError(gBendingInverted[ncombination]->GetN()-1,0.,sigma3);
		              nhitB[ncombination]++;
                 }
                if(xbend4->GetEntries()!=0)
                 {
		              // cout << " xbend4 " << endl;
                  TVector3 *xx4 =(TVector3*)xbend4->At(p);
                  gBending[ncombination]->SetPoint(gBending[ncombination]->GetN(), xx4->Y(), xx4->Z());//fourth layer
                  gBendingInverted[ncombination]->SetPoint(gBendingInverted[ncombination]->GetN(), xx4->Z(), xx4->Y());//fourth layer
                  if(nstrip4.at(p)==3) sigma4= 2*sigma_reso;
                  else if(nstrip4.at(p)==2) sigma4= sigma_reso/(TMath::Sqrt(2.));
		              else sigma4=sigma_reso;
		              // cout << "nstrip4 = " << nstrip4.at(p) << " sigma4 = " << sigma4 << endl;
		              gBending[ncombination]->SetPointError(gBending[ncombination]->GetN()-1,sigma4,0.);
		              gBendingInverted[ncombination]->SetPointError(gBendingInverted[ncombination]->GetN()-1,0.,sigma4);
		              nhitB[ncombination]++;
                 }

                gBending[ncombination]->SetTitle(Form("Bending config %d", ncombination));
                gBending[ncombination]->SetMarkerStyle(kCircle);
                gBending[ncombination]->SetMarkerColor(kBlue);
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
		            //cout << "chi2 bending = " << chi2 << endl;

		            //Order hits
			        	int point=0;
				        if(xbend4->GetEntries()!=0)
                 {
				          int k4 = ibend4.at(p);
				          indicesB[ncombination]->AddAt(k4,point);
			          	//  cout << "PR, bending plane, point added at index " << k4 << endl;
				         point++;
						    }
				       if(xbend3->GetEntries()!=0)
                {
				       	 int k3 = ibend3.at(o);
				 	       indicesB[ncombination]->AddAt(k3,point);
			           //	cout << "PR, bending plane, point added at index " << k3 << endl;
					       point++;
				        }
				       if(xbend2->GetEntries()!=0)
                {
				 	       int k2 = ibend2.at(n);
			     	     indicesB[ncombination]->AddAt(k2,point);
					       //cout << "PR, bending plane, point added at index " << k2 << endl;
                 point++;
					      }
				       if(xbend1->GetEntries()!=0)
                {
				 	       int k1 = ibend1.at(m);
			 	 	       indicesB[ncombination]->AddAt(k1,point);
					       //cout << "PR, bending plane, point added at index " << k1 << endl;
			          }
	             //delete gBending[ncombination];
	             //delete gBendingInverted[ncombination];
		           ncombination++;

              }//end p
           }// end o
         }//end n
       }//end m


    //TEST ALL CONFIGURATIONS NON BENDING PLANE
  int npointsNB, ntop, nmid, nbottom, NConf;
  ntop=nmid=nbottom=1;
  npointsNB = xtoplayer->GetEntries() + xmid->GetEntries() + xbottomlayer->GetEntries();
  int NLNB = (int)NL0 + (int)NL4 + (int)NL6;
  if(xtoplayer->GetEntries()!=0) ntop = xtoplayer->GetEntries(); //else cout << "non-bending plane top layer no hits " << endl;
  if(xmid->GetEntries()!=0) nmid = xmid->GetEntries(); //else cout << "non-bending plane middle layer no hits " << endl;
  if(xbottomlayer->GetEntries()!=0) nbottom = xbottomlayer->GetEntries(); //else cout << "non-bending plane botttom layer no hits " << endl;
  NConf = nbottom*ntop*nmid;
  double sigmatop,sigmamid,sigmabottom;
 // cout << " npointsNB = " << npointsNB << ", NLNB = " << NLNB << endl;

// if less then 2 hits in non-bending plane, can't reconstruct
	if (NLNB <2)
   {
	//  cout << "here we are" << endl;
	  return 0;
	 }

  if(NConf<=0 || NConf>MAXNB)
   {
	  return 0;
   }


   // cout << "There are " << NConf << " possible combinations in the non-bending plane " << endl;
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
               if(nstriptop.at(m)==3) sigmatop=2*sigma_reso;
               else if(nstriptop.at(m)==2) sigmatop=sigma_reso/(TMath::Sqrt(2.));
	             else sigmatop=sigma_reso;
	             //cout << " nstriptop = " << nstriptop.at(m) << " sigmatop = " << sigmatop << endl;
	            // cout << "x = " << xx1->X() << " z = " << xx1->Z() << endl;
	             gNonBending2[NConfig]->SetPointError(gNonBending2[NConfig]->GetN()-1,sigmatop,0);
	             gNonBending[NConfig]->SetPointError(gNonBending[NConfig]->GetN()-1,0,sigmatop);
               nhitNB[NConfig]++;
              }
             if(xbottomlayer->GetEntries()!=0)
              {
               TVector3 *xx2 = (TVector3*)xbottomlayer->At(n);
               gNonBending2[NConfig]->SetPoint(gNonBending2[NConfig]->GetN(), xx2->X(), xx2->Z());
               gNonBending[NConfig]->SetPoint(gNonBending[NConfig]->GetN(), xx2->Z(), xx2->X());
               if(nstripbottom.at(n)==3) sigmabottom=2*sigma_reso;
               else if(nstripbottom.at(n)==2) sigmabottom=sigma_reso/(TMath::Sqrt(2.));
	             else sigmabottom=sigma_reso;
               if(L6S575.at(n)==1) sigmabottom=64*sigma_reso;

	             // cout << " nstripbottom = " << nstripbottom.at(n) << " sigmabottom = " << sigmabottom << endl;
	            // cout << "x = " << xx2->X() << " z = " << xx2->Z() << endl;
	             gNonBending2[NConfig]->SetPointError(gNonBending2[NConfig]->GetN()-1,sigmabottom,0);
	             gNonBending[NConfig]->SetPointError(gNonBending[NConfig]->GetN()-1,0,sigmabottom);
               nhitNB[NConfig]++;
              }
             if(xmid->GetEntries()!=0)
              {
               TVector3 *xx3 = (TVector3*)xmid->At(o);
               gNonBending2[NConfig]->SetPoint(gNonBending2[NConfig]->GetN(), xx3->X(), xx3->Z());
               gNonBending[NConfig]->SetPoint(gNonBending[NConfig]->GetN(), xx3->Z(), xx3->X());
               if(nstripmid.at(o)==3) sigmamid=2*sigma_reso;
               else if(nstripmid.at(o)==2) sigmamid=sigma_reso/(TMath::Sqrt(2.));
	             else sigmamid=sigma_reso;
	             //  cout << " nstripmid = " << nstripmid.at(o) << " sigmamid = " << sigmamid << endl;
	            // cout << "x = " << xx3->X() << " z = " << xx3->Z() << endl;
	             gNonBending2[NConfig]->SetPointError(gNonBending2[NConfig]->GetN()-1,sigmamid,0);
	             gNonBending[NConfig]->SetPointError(gNonBending[NConfig]->GetN()-1,0,sigmamid);
               nhitNB[NConfig]++;
              }


             //Display options
             //gNonBending[NConfig]->SetTitle(Form("Non-bending config %d",  NConfig));
             //gNonBending[NConfig]->SetMarkerStyle(kCircle);
            // gNonBending[NConfig]->SetMarkerColor(kBlue);
             //gNonBending2[NConfig]->SetTitle(Form("Non-bending config %d",  NConfig));
            // gNonBending2[NConfig]->SetMarkerStyle(kCircle);
             //gNonBending2[NConfig]->SetMarkerColor(kBlue);

             //fit function
             line[NConfig] =  new TF1(Form("NBConfig%d",NConfig),"pol1",-20,20);
             //line[NConfig]->SetLineColor(kBlue);
             //line[NConfig]->SetLineWidth(1);
             //line[NConfig]->SetLineStyle(3);
             //FIT
             gNonBending[NConfig]->Fit(line[NConfig],"QSN");
             double chi2=line[NConfig]->GetChisquare();
             chisquare.push_back(chi2);
             //cout << " chi2 = " << chi2 << endl;
             //Order hits
             int index=0;
             if(xmid->GetEntries()!=0)
              {
               int k3 = imid.at(o);
               indicesNB[NConfig]->AddAt(k3,index);
	             //  cout << "PR, non-bending plane, point added at index " << k3 << endl;
               index++;
              }
             if(xbottomlayer->GetEntries()!=0)
              {
               int k2 = ibottom.at(n);
               indicesNB[NConfig]->AddAt(k2,index);
	             //   cout << "PR, non-bending plane, point added at index " << k2 << endl;
               index++;
              }
             if(xtoplayer->GetEntries()!=0)
              {
               int k1 = itop.at(m);
               indicesNB[NConfig]->AddAt(k1,index);
	             //    cout << "PR, non-bending plane, point added at index " << k1 << endl;
               index++;
              }
	           //   delete gNonBending[NConfig];
	           // delete gNonBending2[NConfig];
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
    float lim=zL[1];//z position of 2nd layer
    float zT1=33.49968+0.25; //z position of middle of T1
    float zT3=0.+0.25; //z position of middle of T3
    float zT4=-25.59012+0.5; //z position of middle of T4
    float diff=2*a*lim+2*a*zz0+b;
    float aa=invertedparabolas[indexB]->Eval(lim);

    float YT1PR=(zT1-lim)*diff+aa; //PR Y-position in T1
    float YT3PR=(zT3-lim)*diff+aa; //PR Y-position in T3
    re->set_thBiPR(TMath::ATan(diff));
    re->set_YT1PR(YT1PR);
    re->set_YT3PR(YT3PR);

    //Outcoming Straight particle
    float limo=zL[5];//z position of 6th layer
    float diffout=2*a*limo+2*a*zz0+b;
    aa=invertedparabolas[indexB]->Eval(limo);
    re->set_thBoPR(TMath::ATan(diffout));
    float YT4PR=(zT4-limo)*diffout+aa; //PR Y-position in T4
    re->set_YT4PR(YT4PR);

    //Set Deflection
    double deflection=TMath::ATan(diffout)-TMath::ATan(diff);
    re->set_deflecPR(deflection);

    //NON-BENDING
   //   cout << " index of min chisquare in NB plane i = " << index << endl;
    Double_t p0 = line[index]->GetParameter(0);
    Double_t p1 = line[index]->GetParameter(1);						//p1 = 1/tanl
    float XT1PR=line[index]->Eval(zT1);
    float XT3PR=line[index]->Eval(zT3);
    float XT4PR=line[index]->Eval(zT4);
    re->set_XT1PR(XT1PR);
    re->set_XT3PR(XT3PR);
    re->set_XT4PR(XT4PR);


    //Fill PR coordinates
    re->set_interPR(p0);
    re->set_slopePR(p1);
    re->set_chi2NBPR(chisquare[index]);
    float thetaNB = TMath::ATan(p1);
    re->set_thNBPR(thetaNB);

//Interpolate bending plane hits in the non-bending plane

    for (int l=0; l<nhitB[indexB];l++)
      {
       int kk = indicesB[indexB]->At(l);
       re->get_hits().at(kk)->set_flagPR(true);
       // cout << "In PR, hit index " << kk << " set_flagPR " << endl;
       re->get_hits().at(kk)->set_yPR(invertedparabolas[indexB]->Eval(re->get_hits().at(kk)->get_z()));
       re->get_hits().at(kk)->set_zPR(re->get_hits().at(kk)->get_z());
       //extrapolate the xPR from the fit in the non bending plane
       if(NLNB>=2)re->get_hits().at(kk)->set_xPR(line[index]->Eval(re->get_hits().at(kk)->get_z()));
      }//l

//cout <<"Interpolation in NB plane done " << endl;
    for(int l=0; l<nhitNB[index]; l++)
      {
       //index of chosen points given by array indicesNB
       int kk = indicesNB[index]->At(l);
       //cout << "kk = " << kk << endl;
       re->get_hits().at(kk)->set_flagPR(true);
       //cout << "In PR, hit index " << kk << " set_flagPR " << endl;
       if(p1!=0)
	       re->get_hits().at(kk)->set_zPR(re->get_hits().at(kk)->get_z());
       re->get_hits().at(kk)->set_xPR(p1*re->get_hits().at(kk)->get_z()+p0);
       //extrapolate the yPR from the fit in the bending plane
       if(NLB>=3)re->get_hits().at(kk)->set_yPR(invertedparabolas[indexB]->Eval(re->get_hits().at(kk)->get_z()));
       //cout << "next one " << endl;
      }
//cout << "Interpolation in B plane done" << endl;

//Calculate and set directional cosines for each hit
    for(int j=0;j<re->get_Nhits();j++)
      {
		  bool flagPR = re->get_hits().at(j)->get_flagPR();
		  int k = re->get_hits().at(j)->get_k();
	    if(flagPR)
       {
		    float zPR =  re->get_hits().at(k)->get_zPR();
        double slopeB = 2*a*zPR;
		    double thetaB = TMath::ATan(slopeB);
		    float cxPR = TMath::Cos(TMath::Pi()*0.5+thetaNB);
		    float cyPR = TMath::Cos(TMath::Pi()*0.5-thetaB);
		    float czPR = TMath::Cos(TMath::Pi()-thetaNB);
		    // cout << "cxPR = " << cxPR << " cyPR = " << cyPR << " czPR = " << czPR << endl;
		    re->get_hits().at(k)->set_cxPR(cxPR);
		    re->get_hits().at(k)->set_cyPR(cyPR);
		    re->get_hits().at(k)->set_czPR(czPR);
		   } //if flag
	    }	//end for

//Calculate and set momentum and total energy of particle from PR fit

	//Signed curvature at the 4 points of the bending plane
	double zzz[4]={zL[5],zL[3],zL[2],zL[1]};
	double curv[4]={0,0,0,0};
	TF1* fcurv=new TF1("fcurv","2*[0]/TMath::Power(1+TMath::Power(2*[0]*x+[1],2),3./2.)",-20,20);
	fcurv->SetParameter(0,a);
	fcurv->SetParameter(1,2*a*zz0+b);
	double CurvMean=0;
	for(int ij=0;ij<4;ij++)
	  {
     curv[ij]= fcurv->Eval(zzz[ij]);
	   CurvMean+=	curv[ij]/4;
		}
	double Rmean=fabs(1./curv[2]);
	if(CurvMean!=0)	 Rmean=fabs(1./CurvMean);
	//Extract a simple estimation of energy from the average of the 4 curvature a radius
	double B=0.3; //in T
	double Pt=0.3 * B *	0.01*Rmean; //in GeV
	double p0PR= Pt / TMath::Cos(fabs(thetaNB));   //in GeV
	double mass=0;
	int type;
	double Q =  TMath::Sign(1,deflection);
	if (DataType==1)
   {
	  if (re->get_T2())
     {			//if CK fired, electrons or positrons
	    if (Q>0) type=4;
	    else type=3;
	   }
	  else
     {			//if CK not fired
	    if (Q>0) type=10;
	    else type=11;
	   }
	 } //if data
  else type = re->get_typeMC();	//if MC
  if (type==3 || type == 4) mass = kMelectron;	//for electrons in GeV
  else if(type==10 || type ==11)  mass = kMuon;	//for muons in GeV
  else if(type==1)  mass = kProton;	//for protons in GeV
  else if(type==-6)  mass = kAlpha;	//for alphas in GeV
  double EkPR=TMath::Sqrt(p0PR*p0PR+mass*mass);			//total energy
  re->set_EkPR(EkPR);
  re->set_p0PR(p0PR);

//cout << "end of PR function " << endl;

//delete pointers
//cout << "bending? " << endl;

for(int z =0; z<ncomb;z++)
  {
   //	cout << "bending, z = " << z << endl;
	 delete indicesB[z];
  }

for(int z =0; z<NConf;z++)
  {
   //cout << "nonbending z = " << z << endl;
   delete indicesNB[z];
	}

  delete[] parabolas;
  delete[] invertedparabolas;
  delete[] line;
  delete[] nhitB;
  delete[] nhitNB;
  delete[] gBending;
  delete[] gBendingInverted;
  delete[] gNonBending;
  delete[] gNonBending2;
  delete[] indicesB;
  delete[] indicesNB;
//   gROOT->gObjectTable->Print()

  return 1;

}
