////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 11, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "MakeEventData.h"
#include "ALPatternRecognition.h"
#include "TBox.h"
int MakeEventData(string filename)
{

 
 //Load configuration parameter
/* int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg);
*/

 //Strip pitch in cm
 float  strippitch=0.0228;

 float zz0=-9.6;//in cm
 
 TCanvas*can;
 TMultiGraph*multi;
 TMultiGraph*multiB;
 TF1*tmpf;
 TF1*tmpfB;
 TF1*tmpfB2;


 //Input root file 
 TFile*file=new TFile(Form("%s.root",filename.c_str()),"READ");
 cout << "Input file is open" <<endl;

 //Input root file 
 TFile*fileout=new TFile(Form("%s.EVENT.root",filename.c_str()),"RECREATE");
 cout << "Output file is created" <<endl;
 
 
 //Get Tree from the input file
 TTree *tree = (TTree*)file->Get("BPD");
 //Define variables to read event
 ALEvent *e = new ALEvent();
 //Set address to access event data
 tree->SetBranchAddress("event",&e);  
 // Create a TTree
 TTree *DEtree = new TTree("Data"," Data Event");
 //Define variables to make the Reco event
 ALEvent *de=new ALEvent();
 // Create a branch with event
 DEtree->Branch("event",&de); 
 // Get number of events in Tree
 int nentries=tree->GetEntries();
 cout << "Number  of events: " << nentries << endl;
 //Loop over the events
 // This is the main loop
 
 for (int k=0;k<nentries;k++)
   {
    tree->GetEntry(k); //Load the entry k in the variable e  
    //cout << "Make Event Data: "<<k <<endl;
    //Copy the raw event into a Data event with same structure 
    //if(k!=3) continue;
    de=new ALEvent();
    de->Copy(e);  
     
    //loop over the number of clusters
    int nnhits = (int)de->get_Nhits();

    for(int i=0;i<nnhits;i++)
      {
       float x=-999.;
       float y=-999.;
       float z=-999.;
       float offsetLL=0;//Left ladder offset
       float offsetRL=0;//Right ladder offset
       float center=0;
       int L=((de->get_hits()).at(i))->get_L();//0 to 6
       int nstrips=((de->get_hits()).at(i))->get_nstrips(); //0 if only one strip
       int fstripID=((de->get_hits()).at(i))->get_fstripID();
 
       //Determine the coordinate x,y,z parameters
       //This should be in a configuration files
       switch(L)
         {
          case 0: z=-1.46 ;offsetLL=0.8889;offsetRL=9.8472;break;//Board A
          case 1: z=-3.46 ;offsetLL=0.8866;offsetRL=9.8447;break;//Board B
          case 2: z=-9.46 ;offsetLL=0.8904;offsetRL=9.8478;break;//Board C
          case 3: z=-15.46;offsetLL=0.8871;offsetRL=9.8465;break;//Board D
          case 4: z=-17.46;offsetLL=0.8886;offsetRL=9.8458;break;//Board E
          case 5: z=-19.46;offsetLL=0.8873;offsetRL=9.8460;break;//Board H
          case 6: z=-21.46;offsetLL=0.8884;offsetRL=9.8466;break;//Board G
         }//switch
         
       //Determine the coordinates from the strips
       //Equations from Sarah's email of Spetember 4 2017.
       //Just calculate the mean value of the strip positions
       //Assume that the center of the coodinates is 0 is equidistant from the 
       //two ladders
       center=(offsetRL+(offsetLL+384.*strippitch))/2.;
       //In non-bending plane: Coordinate X in cm
       if(L==0||L==4||L==6)
        {
         //First 6 chips: 0 to 5; strip number 1 to 384
         if(fstripID>0 &&fstripID<=384)
          {
           x=(fstripID-0.5)*strippitch+offsetLL-center; 
          }
         //Last 6 chips: 6 to 11; strip number 385 to 768
         if(fstripID>384)
          {
           x=(fstripID-384.5)*strippitch+offsetRL-center;
          }
          
         //Calculate the center of the cluster depending on the number of strips 
    
         x+=nstrips*strippitch/2.;    
              
          
         y=-999.; 
         
        }
       
       //In bending plane: Coordinate Y in cm
       if(L==1||L==2||L==3||L==5)
        {
         //First 6 chips: 0 to 5; strip number 1 to 384
         if(fstripID>0 &&fstripID<=384)
          {
           y=(fstripID-0.5)*strippitch+offsetLL-center; 
       
          }
         //Last 6 chips: 6 to 11; strip number 385 to 768
         if(fstripID>384)
          {
           y=(fstripID-384.5)*strippitch+offsetRL-center;
       
          }
         y+=nstrips*strippitch/2.;    
 
         x=-999.; 
        }

        
       //Fill the coordinates    
       ((de->get_hits()).at(i))->set_x(x);
       ((de->get_hits()).at(i))->set_y(y);
       ((de->get_hits()).at(i))->set_z(z);          
      }  //i

    ////////////////////////////////////
    //TRIGGER
    ////////////////////////////////////  
      
    //if not all layer touched then don't try pattern recognition
    //cout << "Internal trigger: " << de->get_Ti() <<endl;
    if(de->get_Ti()!=127)
     {   
      /////////////////////    
      //Fill the output file 
      /////////////////////
      DEtree->Fill();
      //Free memory    
      delete de;
      continue;
    }
    //T1&T3&T4
    if(de->get_EneT1().at(0)<0 || de->get_EneT3().at(0)<0 ||de->get_EneT4().at(0)<0)
     {   
      /////////////////////    
      //Fill the output file 
      /////////////////////
      DEtree->Fill();
      //Free memory    
      delete de;
      continue;
    }
    
      
    ////////////////////////////////////
    //Pattern Recognition     Sarah's Code
    ////////////////////////////////////
    
    TObjArray *xnonbend = new TObjArray();
    TObjArray *xbend = new TObjArray();
    TObjArray *allhits = new TObjArray();
    TObjArray *xtoplayer = new TObjArray();
    TObjArray *xbottomlayer = new TObjArray();
    TObjArray *xmid = new TObjArray();
    TObjArray *xbend1 = new TObjArray();
    TObjArray *xbend2 = new TObjArray();
    TObjArray *xbend3 = new TObjArray();
    TObjArray *xbend4 = new TObjArray();
    vector<int> ibend1, ibend2, ibend3, ibend4;//keep track of index of hit
    vector<int> itop, imid, ibottom;//keep track of index of hit

    for(int i=0;i<nnhits;i++)
      {
       float x=(de->get_hits().at(i))->get_x();
       float y=(de->get_hits().at(i))->get_y();
       float z=(de->get_hits().at(i))->get_z();
       int L=((de->get_hits()).at(i))->get_L();//0 to 6
       int fstripID=((de->get_hits()).at(i))->get_fstripID();//1 to 768
       int noisy=((de->get_hits()).at(i))->get_noisy();//0 or 1
      
       if(noisy) continue; //remove cluster with noisy strips

       de->get_hits().at(i)->set_k(i);
       int kk = de->get_hits().at(i)->get_k();
       cout << "add hit index kk = " << kk << "   x = " << x << "  y = " << y << "   z = " << z << endl; 
       TVector3 *X = new TVector3(x,y,z);
       allhits->AddAt(X,kk);
       
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
      }  //i  
    
    //TEST ALL CONFIGURATIONS BENDING PLANE
    
    Int_t npointsB = xbend->GetEntries();
    Int_t layer1 = xbend1->GetEntries();
    Int_t layer2 = xbend2->GetEntries();
    Int_t layer3 = xbend3->GetEntries();
    Int_t layer4 = xbend4->GetEntries();
    Int_t  ncomb   = layer1*layer2*layer3*layer4;
    
    cout << "There are " << ncomb << " possible combinations in the bending plane " << endl;
    if(ncomb<=0||ncomb>100)
      {
       DEtree->Fill();
       //Free memory    
       delete de;
       continue;
      }
    vector<double> chisquareB;
    TArrayI**indicesB = new TArrayI*[ncomb]; //array of TVector3 corresponding to "good" hits in B plane
    TGraph**gBending= new TGraph*[ncomb];
    TGraph**gBendingInverted = new TGraph*[ncomb];
    TF1**parabolas = new TF1*[ncomb];  //arrays of fit functions for B plane 
    TF1**invertedparabolas = new TF1*[ncomb];
    

    //fit all possible parabolas 
    int ncombination = 0;
    for(int m=0; m<layer1; m++)
      {  
       for(int n=0; n<layer2; n++)
         {  
          for(int o=0; o<layer3; o++) 
            {
             for(int p=0;p<layer4;p++) 
               {
                TVector3 *xx1 =(TVector3*)xbend1->At(m);
                TVector3 *xx2 =(TVector3*)xbend2->At(n);
                TVector3 *xx3 =(TVector3*)xbend3->At(o);
                TVector3 *xx4 =(TVector3*)xbend4->At(p);
                indicesB[ncombination] = new TArrayI(4);
                //cout << "BENDING PLANE COMBINATION " << ncombination << endl;
                //cout << "xx1 vector x = " << xx1->X() << "  y = " << xx1->Y() << " z = " << xx1->Z() << endl;
                int k1 = ibend1.at(m);
                indicesB[ncombination]->AddAt(k1,0);
                //cout << "index of vector xx1 is k1 = " << k1 << endl;
                //cout << "xx2 vector x = " << xx2->X() << "  y = " << xx2->Y() << " z = " << xx2->Z() << endl;
                int k2 = ibend2.at(n);
                indicesB[ncombination]->AddAt(k2,1);
                //cout << "index of vector xx2 is k2 = " << k2 << endl;
                //cout << "xx3 vector x = " << xx3->X() << "  y = " << xx3->Y() << " z = " << xx3->Z() << endl;
                int k3 = ibend3.at(o);
                //cout << "index of vector xx3 is k3 = " << k3 << endl;
                indicesB[ncombination]->AddAt(k3,2);
                //cout << "xx4 vector x = " << xx4->X() << "  y = " << xx4->Y() << " z = " << xx4->Z() << endl;
                int k4 = ibend4.at(p);
                //   cout << "index of vector xx4 is k4 = " << k4 << endl;
                indicesB[ncombination]->AddAt(k4,3);
 
                gBending[ncombination] = new TGraph();
                gBending[ncombination]->SetPoint(0, xx1->Y(), xx1->Z());//first layer
                gBending[ncombination]->SetPoint(1, xx2->Y(), xx2->Z());//first layer
                gBending[ncombination]->SetPoint(2, xx3->Y(), xx3->Z());//first layer
                gBending[ncombination]->SetPoint(3, xx4->Y(), xx4->Z());//first layer
                
                gBending[ncombination]->SetTitle(Form("Bending config %d", ncombination));
                gBending[ncombination]->SetMarkerStyle(kCircle);
                gBending[ncombination]->SetMarkerColor(kBlue);
                //parabolas[ncombination] = new TF1("parabola", "pol2");
                //parabolas[ncombination]->SetLineColor(kBlue);
                //parabolas[ncombination]->SetLineWidth(1);
                //parabolas[ncombination]->SetLineStyle(3);
               // gBending[ncombination]->Fit("parabola");
                gBendingInverted[ncombination] = new TGraph();
                gBendingInverted[ncombination]->SetPoint(0, xx1->Z(), xx1->Y());//first layer//first layer
                gBendingInverted[ncombination]->SetPoint(1, xx2->Z(), xx2->Y());//first layer//first layer
                gBendingInverted[ncombination]->SetPoint(2, xx3->Z(), xx3->Y());//first layer//first layer
                gBendingInverted[ncombination]->SetPoint(3, xx4->Z(), xx4->Y());//first layer//first layer
                invertedparabolas[ncombination] = new TF1(Form("BConfig%d",ncombination),"[2]*(x+[3])*(x+[3])+[1]*(x+[3])+[0]",-100,40);
                invertedparabolas[ncombination]->FixParameter(3,zz0); //Postion of the center of the magnet in z coordinates
                gBendingInverted[ncombination]->Fit(invertedparabolas[ncombination],"WQSN");
                invertedparabolas[ncombination]->SetLineColor(kBlue);
                gBendingInverted[ncombination]->SetMarkerStyle(kCircle);
                invertedparabolas[ncombination]->SetLineWidth(2);
                invertedparabolas[ncombination]->SetLineStyle(1);
                Double_t chi2 = invertedparabolas[ncombination]->GetChisquare();
                chisquareB.push_back(chi2);
                ncombination++;
               }//end p
            }// end o
         }//end n
      } //end m
    
    //TEST ALL CONFIGURATIONS NON BENDING PLANE
    Int_t npointsNB = xnonbend->GetEntries();
    Int_t ntop = xtoplayer->GetEntries();
    Int_t nbottom = xbottomlayer->GetEntries();
    Int_t nmid = xmid->GetEntries();
    Int_t NConf = nbottom*ntop*nmid; 
    if(NConf<=0||NConf>100)
      {
       DEtree->Fill();
       //Free memory    
       delete de;
       continue;
      }
    
    
    cout << "There are " << NConf << " possible combinations in the non-bending plane " << endl;
    vector<double> chisquare;      
    TArrayI**indicesNB = new TArrayI*[NConf];//array of TVector3 corresponding to "good" hits in NB plane
    TGraph**gNonBending= new TGraph*[NConf];
    TGraph**gNonBending2= new TGraph*[NConf];
    TF1**line = new TF1*[NConf];                  //arrays of fit functions for NB plane     
   
    //fit all possible lines 
    Int_t NConfig = 0;
    for(int m=0; m<ntop; m++) 
      {
       for(int n=0; n<nbottom; n++)
         { 
          for(int o=0; o<nmid; o++)
            {
             cout << "NON BENDING PLANE COMBINATION " << NConfig << endl;
             TVector3 *xx1 = (TVector3*)xtoplayer->At(m);
             TVector3 *xx2 = (TVector3*)xbottomlayer->At(n);
             TVector3 *xx3 = (TVector3*)xmid->At(o);
             indicesNB[NConfig] = new TArrayI(3);
             cout << "xx1 vector x = " << xx1->X() << "  y = " << xx1->Y() << " z = " << xx1->Z() << endl;
             int k1 = itop.at(m);
             indicesNB[NConfig]->AddAt(k1,0);
             cout << "index of vector xx1 is k1 = " << k1 << endl;
             cout << "xx2 vector x = " << xx2->X() << "  y = " << xx2->Y() << " z = " << xx2->Z() << endl;
             int k2 = ibottom.at(n);
             cout << "index of vector xx2 is k2 = " << k2 << endl;			
             indicesNB[NConfig]->AddAt(k2,1);
             cout << "xx3 vector x = " << xx3->X() << "  y = " << xx3->Y() << " z = " << xx3->Z() << endl;
             int k3 = imid.at(o);
             cout << "index of vector xx3 is k3 = " << k3 << endl;			
             indicesNB[NConfig]->AddAt(k3,2);
            
             gNonBending[NConfig] = new TGraph();
             gNonBending2[NConfig] = new TGraph();
             gNonBending2[NConfig]->SetPoint(0, xx1->X(), xx1->Z());
             gNonBending2[NConfig]->SetPoint(1, xx2->X(), xx2->Z());
             gNonBending2[NConfig]->SetPoint(2, xx3->X(), xx3->Z());
             gNonBending[NConfig]->SetPoint(0, xx1->Z(), xx1->X());
             gNonBending[NConfig]->SetPoint(1, xx2->Z(), xx2->X());
             gNonBending[NConfig]->SetPoint(2, xx3->Z(), xx3->X());
             gNonBending[NConfig]->SetTitle(Form("Non-bending config %d",  NConfig));
             gNonBending[NConfig]->SetMarkerStyle(kCircle);
             gNonBending[NConfig]->SetMarkerColor(kBlue);
             gNonBending2[NConfig]->SetTitle(Form("Non-bending config %d",  NConfig));
             gNonBending2[NConfig]->SetMarkerStyle(kCircle);
             gNonBending2[NConfig]->SetMarkerColor(kBlue);
             line[NConfig] =  new TF1(Form("NBConfig%d",NConfig),"pol1",-40,40);
 //            line[NConfig] =  new TF1(Form("NBConfig%d",NConfig),"pol1",-10,10);
             gNonBending[NConfig]->Fit(line[NConfig],"WQSN");                                         
             //line[NConfig] =  (TF1*)gNonBending[NConfig]->GetFunction("pol1");          
             line[NConfig]->SetLineColor(kBlue);
             line[NConfig]->SetLineWidth(1);
             line[NConfig]->SetLineStyle(3);
             Double_t chi2 = line[NConfig]->GetChisquare();
             chisquare.push_back(line[NConfig]->GetChisquare());
             cout << "chi2: " << chi2 <<endl;
             NConfig++;
            }//end o
         }//end n
      }//end m
    //GRAPH
    TMultiGraph *mg1 = new TMultiGraph();
    TMultiGraph *mg2 = new TMultiGraph();
    TGraph *gr1 = new TGraph();//all MC hits, bending plane
    TGraph *invertedgr1 = new TGraph();//all MC hits, bending plane
    TGraph *gr3 = new TGraph();//all MC hits, NB plane
    TGraph *gInterpolateB = new TGraph();//interpolated B hits in NB view
    TGraph *gMCbending = new TGraph();//MC B hits NB view 
    TGraph *gInterpolateNB = new TGraph();//interpolated NB hits in NB view
    TGraph *gfinalNB = new TGraph();//final chosen NB points
    TGraph *gfinalB = new TGraph();//final chosen B points 
    TGraph *gMCnonbending = new TGraph();//all MC NB points in B view  
   
    for(int i=0; i<xbend->GetEntries(); i++)
      {
       TVector3 *x = (TVector3*)xbend->At(i);
       gr1->SetPoint(i, x->Y(), x->Z());
       invertedgr1->SetPoint(i, x->Z(), x->Y());
       gMCbending->SetPoint(i, x->X(), x->Z());
      }//i
    invertedgr1->SetMarkerStyle(5);
    for(int i=0; i<xnonbend->GetEntries(); i++) 
      {
       TVector3 *x = (TVector3*)xnonbend->At(i);
       gr3->SetPoint(i, x->X(), x->Z());
       gMCnonbending->SetPoint(i, x->Y(), x->Z());
      }//i
  
    //GET BEST FIT PARAMETERS
    //BENDING
    //Get the reciproque  of the parabola ax^2+bx+c
    //Need 2 functions to plot
    int indexB = min_element(chisquareB.begin(), chisquareB.end()) - chisquareB.begin();      //index of min. chisquare fit         
    Double_t a = invertedparabolas[indexB]->GetParameter(2);
    Double_t b = invertedparabolas[indexB]->GetParameter(1);
    Double_t c = invertedparabolas[indexB]->GetParameter(0);
    
    de->set_chi2B(chisquareB[indexB]);
    de->set_a(a);
    de->set_b(b);
    de->set_c(c);

    
    Double_t b2a=0;  
    Double_t b4ac=0;
    if(a!=0)
     {
      b2a= b/(2*a);  
      b4ac=(b*b-4*a*c)/(4*a);
     }   
    else 
     {
       DEtree->Fill();
       //Free memory    
       delete de;
       continue;
     }   
    if(a>=0)tmpfB=new TF1("tmpfB","TMath::Sqrt((x+[0])/[1])-[2]",gBendingInverted[indexB]->Eval(-b/(2*a)),1000);
    if(a<0)tmpfB=new TF1("tmpfB","TMath::Sqrt((x+[0])/[1])-[2]",-1000,gBendingInverted[indexB]->Eval(-b/(2*a)));
    tmpfB->FixParameter(0,b4ac);
    tmpfB->FixParameter(1,a);
    tmpfB->FixParameter(2,b2a);
    if(a>=0)tmpfB2=new TF1("tmpfB2","-TMath::Sqrt((x+[0])/[1])-[2]",gBendingInverted[indexB]->Eval(-b/(2*a)),1000);
    if(a<0)tmpfB2=new TF1("tmpfB2","-TMath::Sqrt((x+[0])/[1])-[2]",-1000,gBendingInverted[indexB]->Eval(-b/(2*a)));
    tmpfB2->FixParameter(0,b4ac);
    tmpfB2->FixParameter(1,a);
    tmpfB2->FixParameter(2,b2a);

    
    tmpfB->SetLineColor(kBlue);
    tmpfB->SetLineWidth(2);
    tmpfB->SetLineStyle(1);
    tmpfB2->SetLineColor(kOrange);
    tmpfB2->SetLineWidth(2);
    tmpfB2->SetLineStyle(1);   
    
    for (int l=0; l<4;l++)
      {
       int kk = indicesB[indexB]->At(l);
       de->get_hits().at(kk)->set_flagPR(true);
       //float yPR  = (float) yMC;
       //float zPR = (float) zMC;
       de->get_hits().at(kk)->set_yPR(invertedparabolas[indexB]->Eval(de->get_hits().at(kk)->get_z()));       
       de->get_hits().at(kk)->set_zPR(invertedparabolas[indexB]->GetX(de->get_hits().at(kk)->get_y()));          
      }//l
    
    //Incoming Straight particle    
    float lim=-3.46;//z position of 2nd layer
    TF1*incomingB=new TF1("incomingB","pol1",lim,40);
    incomingB->SetLineStyle(2);
    incomingB->SetLineWidth(2);
    incomingB->SetLineColor(kBlue);
    float aa=invertedparabolas[indexB]->Eval(lim);
    float diff=2*a*lim+2*a*zz0+b;
    incomingB->FixParameter(0,aa-diff*lim);
    incomingB->FixParameter(1,diff);
    //Outcoming Straight particle    
    float limo=-19.46;//z position of 6th layer
    TF1*outcomingB=new TF1("outcomingB","pol1",-30,limo);
    float aaout=invertedparabolas[indexB]->Eval(limo);
    float diffout=2*a*limo+2*a*zz0+b;
    outcomingB->FixParameter(0,aaout-diffout*limo);
    outcomingB->FixParameter(1,diffout);
    outcomingB->SetLineWidth(2);
    outcomingB->SetLineColor(kBlue);
     
    double deflection=diffout-diff;    
    de->set_deflec(deflection);
    
    //Signed curvature at the 3 points in and around magnets
    float zzz[3]={-15.46,-9.46,-3.46};
    float curv[3]={3,3,3};
    TF1* fcurv=new TF1("fcurv","2*[0]/TMath::Power(1+TMath::Power(2*[0]*x+[1],2),3./2.)",-20,20);
    fcurv->SetParameter(0,a);
    fcurv->SetParameter(1,2*a*zz0+b);
    for(int ij=0;ij<3;ij++)
      {
       curv[ij]= fcurv->Eval(zzz[ij]);   
      }  
    
    //NON-BENDING
    int index = min_element(chisquare.begin(), chisquare.end()) - chisquare.begin();      //index of min. chisquare fit
    //   cout << " index of min chisquare in NB plane i = " << index << endl; 	  
    Double_t p0 = line[index]->GetParameter(0);                                     
    Double_t p1 = line[index]->GetParameter(1);						//p1 = 1/tanl 
    Double_t tanl = 1/p1;
    Double_t lambda = atan(tanl);
    Double_t cosl = cos(lambda);
    
    
    //Fill Fit parameter
    //Fill PR coordinates
    de->set_inter(p0);
    de->set_slope(p1);
    de->set_chi2NB(chisquare.at(0));

    cout << "chi2: " << chisquare.at(index) <<endl;
    cout << "p0: " << p0 <<endl;
    cout << "p1: " << p1 <<endl;

    for(int l=0; l<3; l++)
      {
       //index of chosen points given by array indicesNB
       int kk = indicesNB[index]->At(l);    
       de->get_hits().at(kk)->set_flagPR(true);
 //      if(p1!=0)
 //       de->get_hits().at(kk)->set_xPR((de->get_hits().at(kk)->get_z()-p0)/p1);
 //       de->get_hits().at(kk)->set_zPR(p1*de->get_hits().at(kk)->get_x()+p0);
       if(p1!=0)
        de->get_hits().at(kk)->set_zPR((de->get_hits().at(kk)->get_x()-p0)/p1);
        de->get_hits().at(kk)->set_xPR(p1*de->get_hits().at(kk)->get_z()+p0);
      }
      
      
    //PLOT  
    if(de->get_eventnumber()<100)
     {  
      can=new TCanvas(Form("Event %d",de->get_eventnumber()),Form("Event %d",de->get_eventnumber()),200,10,1600,800);
      can->Divide(2,1);
      //Detector layout
      //Layers
      TLine**Line=new TLine*[7];
      Line[0]=new TLine(-9,-1.46 ,9,-1.46 );
      Line[1]=new TLine(-9,-3.46 ,9,-3.46 );
      Line[2]=new TLine(-9,-9.46 ,9,-9.46 );
      Line[3]=new TLine(-9,-15.46,9,-15.46);
      Line[4]=new TLine(-9,-17.46,9,-17.46);
      Line[5]=new TLine(-9,-19.46,9,-19.46);
      Line[6]=new TLine(-9,-21.46,9,-21.46);
      
      for(int ijk=0;ijk<7;ijk++)
        {
         Line[ijk]->SetLineColor(kBlack);  
         Line[ijk]->SetLineWidth(1);  
        }
       
      //Magnets 
      TBox*boxM1=new TBox(-15,-8.6,-6.703,-4.6);
      TBox*boxM2=new TBox(-15,-14.6,-6.703,-10.6);
      TBox*boxM3=new TBox(6.703,-8.6,15,-4.6);
      TBox*boxM4=new TBox(6.703,-14.6,15,-10.6);
      boxM1->SetFillColor(kGray);
      boxM2->SetFillColor(kGray);
      boxM3->SetFillColor(kGray);
      boxM4->SetFillColor(kGray);
      
      //T1
      TBox*boxT1=new TBox(-13.,33.,13.,33.5);
      boxT1->SetFillColor(kRed);
      if(de->get_EneT1().at(0)>0)boxT1->SetFillColor(kGreen);
      TPaveText*PHT1=new TPaveText(15,31,25,35);
      PHT1->AddText(Form("PHT1=%d",(int)de->get_EneT1().at(0)));
      PHT1->SetFillStyle(0);PHT1->SetBorderSize(0);

      //T2
      TGraph* grT2 = new TGraph();
      grT2->SetPoint(0,-6.5,2.05588);//Bottom
      grT2->SetPoint(1,6.5,2.05588);//Botton
      grT2->SetPoint(2,13.5,27.8824);//Top
      grT2->SetPoint(3,-13.5,27.8824);//Top
      grT2->SetPoint(4,-6.5,2.05588);//Bottom again
      grT2->SetFillStyle(0);
      grT2->SetLineWidth(3);
      grT2->SetLineColor(kRed);
      if(de->get_EneT2().at(0)>0)grT2->SetLineColor(kGreen); 
      TPaveText*PHT2=new TPaveText(15,10,25,14);
      PHT2->AddText(Form("PHT2=%d",(int)de->get_EneT2().at(0)));
      PHT2->SetFillStyle(0);PHT2->SetBorderSize(0);     
      //T3
      TBox*boxT3=new TBox(-3.5,0.,3.5,0.5);
      boxT3->SetFillColor(kRed);
      if(de->get_EneT3().at(0)>0)boxT3->SetFillColor(kGreen);
      TPaveText*PHT3=new TPaveText(15,1,25,5);
      PHT3->AddText(Form("PHT3=%d",(int)de->get_EneT3().at(0)));
      PHT3->SetFillStyle(0);PHT3->SetBorderSize(0);

      //Guard
      TBox*boxG1=new TBox(-13.5,0.,-3.5,0.5);
      boxG1->SetFillColor(kRed);
      if(de->get_Eneg().at(0)>0)boxG1->SetFillColor(kGreen);
      TBox*boxG2=new TBox(+3.5,0.,13.5,0.5);
      boxG2->SetFillColor(kRed);
      if(de->get_Eneg().at(0)>0)boxG2->SetFillColor(kGreen);     
      TPaveText*PHG=new TPaveText(15,-3,25,1);
      PHG->AddText(Form("PHG=%d",(int)de->get_Eneg().at(0)));
      PHG->SetFillStyle(0);PHG->SetBorderSize(0);
      
      //T4
      TBox*boxT4=new TBox(-18.,-24.,18.,-23.);
      boxT4->SetFillColor(kRed);
      if(de->get_EneT4().at(0)>0)boxT4->SetFillColor(kGreen);
      TPaveText*PHT4=new TPaveText(15,-22,25,-18);
      PHT4->AddText(Form("PHT4=%d",(int)de->get_EneT4().at(0)));
      PHT4->SetFillStyle(0);PHT4->SetBorderSize(0);
      
      //Bending plot   
      can->cd(1);
      multi=new TMultiGraph();multi->SetTitle(Form("Event: %d, Non bending plane",(int)de->get_eventnumber()));
      gr3->SetMarkerStyle(5);
      multi->Add(gr3,"p");
      multi->Add(gNonBending2[index],"p");
      multi->Add(grT2,"");
      multi->Draw("a");

      multi->GetXaxis()->SetTitle("X (cm)");
      multi->GetYaxis()->SetTitle("Z (cm)");
      multi->GetXaxis()->SetLimits(-25,25);
      multi->SetMaximum(40);
      multi->SetMinimum(-30);
      gPad->Update();
      //Magnet
      boxM1->Draw("l same");
      boxM2->Draw("l same");
      boxM3->Draw("l same");
      boxM4->Draw("l same");       

      //T1
      boxT1->Draw("l same");
      PHT1->Draw();
      //T2
      PHT2->Draw();
      //T3
      boxT3->SetLineWidth(1);
      boxT3->SetLineColor(1);
      boxT3->Draw("l same");
      PHT3->Draw();
      //Guard
      boxG1->SetLineWidth(1);
      boxG2->SetLineWidth(1);
      boxG1->SetLineColor(1);
      boxG2->SetLineColor(1);

      boxG1->Draw("l same");
      boxG2->Draw("l same");
      PHG->Draw();
      //T4
      boxT4->Draw("l same");
      PHT4->Draw();
      //Layers
      for(int ijk=0;ijk<7;ijk++) Line[ijk]->Draw("same");  
      
      tmpf=new TF1("tmpf","pol1",-20,20);
      if(p1!=0)tmpf->FixParameter(0,-p0/p1);
      if(p1!=0)tmpf->FixParameter(1,1./p1);
      tmpf->SetLineColor(kBlue);
      tmpf->SetLineWidth(2);
      tmpf->SetLineStyle(1);
     // line[index]->Draw("same");
      if(p1!=0)tmpf->Draw("same");
      //Bending plot   
      //Detector layout
      //Layers
      TLine**LineB=new TLine*[7];
      LineB[0]=new TLine(-1.46 ,-9,-1.46 ,9);
      LineB[1]=new TLine(-3.46 ,-9,-3.46 ,9);
      LineB[2]=new TLine(-9.46 ,-9,-9.46 ,9);
      LineB[3]=new TLine(-15.46,-9,-15.46,9);
      LineB[4]=new TLine(-17.46,-9,-17.46,9);
      LineB[5]=new TLine(-19.46,-9,-19.46,9);
      LineB[6]=new TLine(-21.46,-9,-21.46,9);
      
      for(int ijk=0;ijk<7;ijk++)
        {
         LineB[ijk]->SetLineColor(kBlack);  
         LineB[ijk]->SetLineWidth(1);  
        }
      //Magnets 
      TBox*boxM1B=new TBox(-8.6,-15,-4.6, -6.703);
      TBox*boxM2B=new TBox(-14.6,-15,-10.6, -6.703);
      TBox*boxM3B=new TBox(-8.6,6.703,-4.6,15);
      TBox*boxM4B=new TBox(-14.6,6.703,-10.6,15);
      boxM1B->SetFillColor(kGray);
      boxM2B->SetFillColor(kGray);
      boxM3B->SetFillColor(kGray);
      boxM4B->SetFillColor(kGray);       
      //T1
      TBox*boxT1B=new TBox(33.,-13.,33.5,13.);
      boxT1B->SetFillColor(kRed);
      if(de->get_EneT1().at(0)>0)boxT1B->SetFillColor(kGreen);

      //T2
      TGraph* grT2B = new TGraph();
      grT2B->SetPoint(0,2.05588,-6.5);//Bottom
      grT2B->SetPoint(1,2.05588,6.5);//Botton
      grT2B->SetPoint(2,27.8824,13.5);//Top
      grT2B->SetPoint(3,27.8824,-13.5);//Top
      grT2B->SetPoint(4,2.05588,-6.5);//Bottom again
      grT2B->SetFillStyle(0);
      grT2B->SetLineWidth(3);
      grT2B->SetLineColor(kRed);
      if(de->get_EneT2().at(0)>0)grT2B->SetLineColor(kGreen); 
      
      //T3
      TBox*boxT3B=new TBox(0.,-3.5,0.5,3.5);
      boxT3B->SetFillColor(kRed);
      if(de->get_EneT3().at(0)>0)boxT3B->SetFillColor(kGreen);
      //Guard
      TBox*boxG1B=new TBox(0.,-13.5,0.5,-3.5);
      boxG1B->SetFillColor(kRed);
      if(de->get_Eneg().at(0)>0)boxG1B->SetFillColor(kGreen);
      TBox*boxG2B=new TBox(0.,3.5,0.5,13.5);
      boxG2B->SetFillColor(kRed);
      if(de->get_Eneg().at(0)>0)boxG2B->SetFillColor(kGreen);     
      
      //T4
      TBox*boxT4B=new TBox(-24.,-18.,-23.,18.);
      boxT4B->SetFillColor(kRed);
      if(de->get_EneT4().at(0)>0)boxT4B->SetFillColor(kGreen);


      can->cd(2);
      multiB=new TMultiGraph(); multiB->SetTitle("Bending plane");
      invertedgr1->SetMarkerStyle(5);
      multiB->Add(invertedgr1,"p");
      multiB->Add(gBendingInverted[indexB],"p");
      multiB->Add(grT2B,"");
      multiB->Draw("a");
      multiB->GetXaxis()->SetTitle("Z (cm)");
      multiB->GetYaxis()->SetTitle("Y (cm)");
      multiB->GetXaxis()->SetLimits(-30,40);
      multiB->SetMaximum(25);
      multiB->SetMinimum(-25);
      gPad->Update();
      //Magnet
      boxM1B->Draw("l same");
      boxM2B->Draw("l same");
      boxM3B->Draw("l same");
      boxM4B->Draw("l same");       
 
      //T1
      boxT1B->Draw("l same");
      //T2
      //T3
      boxT3B->SetLineWidth(1);
      boxT3B->SetLineColor(1);
      boxT3B->Draw("l same");
      //Guard
      boxG1B->SetLineWidth(1);
      boxG1B->SetLineColor(1);
      boxG2B->SetLineWidth(1);
      boxG2B->SetLineColor(1);
      boxG1B->Draw("l same");
      boxG2B->Draw("l same");
      //T4
      boxT4B->Draw("l same");
      //Layers
      for(int ijk=0;ijk<7;ijk++) LineB[ijk]->Draw("same");      
      
      
      outcomingB->SetLineStyle(2);
      outcomingB->Draw("same");
      incomingB->Draw("same");
      invertedparabolas[indexB]->DrawF1(limo,lim,"same");

      TPaveText*Def=new TPaveText(10,20,30,24);
      Def->AddText(Form("Deflection: %6.4f",(float)deflection));
      Def->SetFillStyle(0);Def->SetBorderSize(0);
      Def->Draw();
      TPaveText*Cur=new TPaveText(-20,-23,20,-16);
      Cur->AddText(Form("Curvature at L3=%5.4f",curv[0]));
      Cur->AddText(Form("Curvature at L2=%5.4f",curv[1]));
      Cur->AddText(Form("Curvature at L1=%5.4f",curv[2]));
      Cur->SetFillStyle(0);Cur->SetBorderSize(0);
      Cur->Draw();
      
      fileout->cd();
      can->Write();
      if(abs(deflection)>=0.1&&de->get_EneT2().at(0)<=0&&de->get_Eneg().at(0)<=0) can->SaveAs(Form("%s_Event%06d_deflec.pdf",filename.c_str(),(int)de->get_eventnumber()));
     }
    
    /////////////////////    
    //Fill the output file 
    /////////////////////
    DEtree->Fill();
    //Free memory    
    delete de;
   }//k End loop on events
 delete e;
 //Write tree in output file
 fileout->cd();
 DEtree->Write();
 //Close files
 fileout->Close();
 file->Close();
      
 return 0;
}




