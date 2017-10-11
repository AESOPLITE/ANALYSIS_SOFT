////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, September 11, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "MakeEventData.h"
#include "ALPatternRecognition.h"
#include "LoadDataparameters.h"
#include "TBox.h"
int MakeEventData(string filename)
{

 
 //Load configuration parameter
 float* zL=new float[7];
 float*OffsetLL=new float[7];
 float*OffsetRL=new float[7];
 for(int i=0;i<7;i++)zL[i]=OffsetLL[i]=OffsetRL[i]=0;
 string paramfile="../src/ALSim/Dataparameters.dat"; 
 
 LoadDataparameters(paramfile,zL,OffsetLL,OffsetRL);

  for(int i=0;i<7;i++)
   {
    cout << "L"<<i <<", zL:" << zL[i] ;
    cout << ", OffsetLL:" << OffsetLL[i] ;
    cout << ", OffsetRL:" << OffsetRL[i] << endl;
   }  
 
 //Maximum of configuration to try fitting
 
 int MAXB=1000;
 int MAXNB=1000;

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
    if(k%10000==0) cout << "Event " << k << endl;
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
       //Get the values from configuration file read at the beginning of th efunction
       z=zL[L];
       offsetLL=OffsetLL[L];
       offsetRL=OffsetRL[L];
         
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
   /* if(de->get_EneT1().at(0)<0 || de->get_EneT3().at(0)<0 ||de->get_EneT4().at(0)<0)
     {   
      /////////////////////    
      //Fill the output file 
      /////////////////////
      DEtree->Fill();
      //Free memory    
      delete de;
      continue;
    }
    */
      
    ////////////////////////////////////
    //Pattern Recognition     Sarah's Code
    ////////////////////////////////////
    
    TObjArray *xnonbend = new TObjArray();
    TObjArray *xbend = new TObjArray();
    TObjArray *allhits = new TObjArray(nnhits);
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
       //int noisy=((de->get_hits()).at(i))->get_noisy();//0 or 1
      
       //if(noisy) continue; //remove cluster with noisy strips

       de->get_hits().at(i)->set_k(i);
       int kk = de->get_hits().at(i)->get_k();
       //cout << "add hit index kk = " << kk << "   x = " << x << "  y = " << y << "   z = " << z << endl; 
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
    
    //cout << "There are " << ncomb << " possible combinations in the bending plane " << endl;
    if(ncomb<=0||ncomb>MAXB)
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
    if(NConf<=0||NConf>MAXNB)
      {
       DEtree->Fill();
       //Free memory    
       delete de;
       continue;
      }
    
    
    //cout << "There are " << NConf << " possible combinations in the non-bending plane " << endl;
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
             // cout << "NON BENDING PLANE COMBINATION " << NConfig << endl;
             TVector3 *xx1 = (TVector3*)xtoplayer->At(m);
             TVector3 *xx2 = (TVector3*)xbottomlayer->At(n);
             TVector3 *xx3 = (TVector3*)xmid->At(o);
             indicesNB[NConfig] = new TArrayI(3);
             // cout << "xx1 vector x = " << xx1->X() << "  y = " << xx1->Y() << " z = " << xx1->Z() << endl;
             int k1 = itop.at(m);
             indicesNB[NConfig]->AddAt(k1,0);
             // cout << "index of vector xx1 is k1 = " << k1 << endl;
             //cout << "xx2 vector x = " << xx2->X() << "  y = " << xx2->Y() << " z = " << xx2->Z() << endl;
             int k2 = ibottom.at(n);
             //cout << "index of vector xx2 is k2 = " << k2 << endl;			
             indicesNB[NConfig]->AddAt(k2,1);
             // cout << "xx3 vector x = " << xx3->X() << "  y = " << xx3->Y() << " z = " << xx3->Z() << endl;
             int k3 = imid.at(o);
             // cout << "index of vector xx3 is k3 = " << k3 << endl;			
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
             //cout << "chi2: " << chi2 <<endl;
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
      
    for (int l=0; l<4;l++)
      {
       int kk = indicesB[indexB]->At(l);
       de->get_hits().at(kk)->set_flagPR(true);
       de->get_hits().at(kk)->set_yPR(invertedparabolas[indexB]->Eval(de->get_hits().at(kk)->get_z()));       
       de->get_hits().at(kk)->set_zPR(invertedparabolas[indexB]->GetX(de->get_hits().at(kk)->get_y()));          
      }//l
    
    //Incoming Straight particle    
    float lim=zL[1];//z position of 2nd layer
    float diff=2*a*lim+2*a*zz0+b;
    
    //Outcoming Straight particle    
    float limo=zL[5];//z position of 6th layer
    float diffout=2*a*limo+2*a*zz0+b;
    
    //Set Deflection
    double deflection=diffout-diff;    
    de->set_deflec(deflection);
    
    //NON-BENDING
    int index = min_element(chisquare.begin(), chisquare.end()) - chisquare.begin();      //index of min. chisquare fit
    //   cout << " index of min chisquare in NB plane i = " << index << endl; 	  
    Double_t p0 = line[index]->GetParameter(0);                                     
    Double_t p1 = line[index]->GetParameter(1);//p1 = 1/tanl 
    Double_t tanl = 1/p1;
    Double_t lambda = atan(tanl);
    Double_t cosl = cos(lambda);
    
    
    //Fill PR coordinates
    de->set_inter(p0);
    de->set_slope(p1);
    de->set_chi2NB(chisquare.at(0));

    //cout << "chi2: " << chisquare.at(index) <<endl;
    // cout << "p0: " << p0 <<endl;
    // cout << "p1: " << p1 <<endl;
 
    for(int l=0; l<3; l++)
      {
       //index of chosen points given by array indicesNB
       int kk = indicesNB[index]->At(l);    
       de->get_hits().at(kk)->set_flagPR(true);

       if(p1!=0)
        de->get_hits().at(kk)->set_zPR((de->get_hits().at(kk)->get_x()-p0)/p1);
        de->get_hits().at(kk)->set_xPR(p1*de->get_hits().at(kk)->get_z()+p0);
      }
    
    /////////////////////    
    //RECONSTRUCTION
    /////////////////////    
    
    //Include Kalman Filter
    
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
 DEtree->Write(0,TObject::kOverwrite);
 //Close files
 fileout->Close();
 file->Close();
      
 return 0;
}




