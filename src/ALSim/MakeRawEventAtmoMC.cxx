////////////////////////////////////////////////////////////////////////////////////////////////////
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
///			   Sarah Mechbal, smechbal@ucsc.edu
///	   Santa Cruz Institute for Particle Physics, University of California, Santa Cruz, March 2019 
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "MakeRawEventAtmoMC.h"


int MakeRawEventAtmoMC(int typeT,int layer, int cycle,string Inppath,string Inppath2,string Outpath,string startfile,string endfile)
{
 cout << "Calling MakeRawEventAtmoMC" << endl;   
 float OffsetLL=0;
 float OffsetRL=0;
    
//Load region numbers used in the MC geometry 
 int*TckReg=new int[7];
 int*TrigReg=new int[4];
 int*GReg=new int[1];
 int*ShellReg=new int[2];
 float*TckZPos=new float[7];
 float*TrigThresh=new float[4];
 float*GuardThresh=new float[1];
 for(int i=0;i<7;i++)TckReg[i]=0;
 for(int i=0;i<4;i++)TrigReg[i]=0;
 for(int i=0;i<2;i++)ShellReg[i]=0;
 for(int i=0;i<1;i++)GReg[i]=0;
 for(int i=0;i<7;i++)TckZPos[i]=0;
 for(int i=0;i<4;i++)TrigThresh[i]=0;
 for(int i=0;i<1;i++)GuardThresh[i]=0;
 int nL0 = 0;
 int nL1 = 0;
 int nL2 = 0;
 int nL3 = 0;
 int nL4 = 0;
 int nL5 = 0;
 int nL6 = 0;
 int mreg11=0;
 int mreg12=0;
 int mreg13=0;
 int mreg14=0;
 int mreg15=0;
 int mreg16=0;
 int mreg17=0;
 int rL0 = 0;
 int rL1 = 0;
 int rL2 = 0;
 int rL3 = 0;
 int rL4 = 0;
 int rL5 = 0;
 int rL6 = 0;
	
 string MCparamfile="../src/ALSim/MCparameters.dat"; 
 
 LoadMCparameters(MCparamfile,TckReg,TrigReg,GReg,TckZPos,TrigThresh,GuardThresh,ShellReg);
      
 
 float ladderOffsetLeft[7] = { 8.889, 8.866, 8.904, 8.871, 8.886, 8.873, 8.884 };
 float ladderOffsetRight[7] = { 98.472, 98.447, 98.478, 98.465, 98.460, 98.466 };
 
 
 
 TFile *file;
 TFile *fileout;
 //Input file 

 file=new TFile(Form("%s/%d/%s/%s_%d_%d_%03d%s.root",Inppath.c_str(),typeT,Inppath2.c_str(),startfile.c_str(),typeT,layer,cycle,endfile.c_str()),"READ");
	cout << "Input file is open" <<endl;
   //Output file 
   fileout=new TFile(Form("%s/%d/RawEvent_%s_%d_%d_%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,layer,cycle,endfile.c_str()),"RECREATE");
   cout << "Output file " << Form("%s/%d/RawEvent_%s_%d_%d_%03d%s.root",Outpath.c_str(),typeT,startfile.c_str(),typeT,layer,cycle,endfile.c_str())<< " is created" <<endl;
  
  
 //Get ntuple from the input file
 TNtuple*ntuple=(TNtuple*)file->Get("Track");
 cout << "Got the ntuple from the input file" <<endl;

 //Define variables to read ntuple
 float col1=0;
 float col2=0;
 float col3=0;
 float col4=0;
 float col5=0;
 float col6=0;
 float col7=0;
 float col8=0;
 float col9=0;
 float col10=0;
 float col11=0;
 float col12=0;
 float col13=0;
 float col14=0;
 float col15=0;
 float col16=0;
 float col17=0;
		
 //Variables for PP
 float typePP=0;
 float EkPP=0;
 float ZenPP=0;
 float AziPP=0;
 float CoLatSP=0;
 float CoLonSP=0;
 float flagSP=0;
 //Variables for SP
 float ncase=0; 
 float mreg=0;
 float mtrack=0;
 float type=0;
 float EkMC=0;
 float eMC=0; 
 float pMC=0;
 float x=0;
 float y=0;
 float z=0;
 float age=0; 
 float cx=0;
 float cy=0;
 float cz=0; 
 float Edep=0; 
 float flag=0; 
	
 //Set addresses to access ntuple data
 ntuple->SetBranchAddress("col1",&col1); 
 ntuple->SetBranchAddress("col2",&col2); 
 ntuple->SetBranchAddress("col3",&col3); 
 ntuple->SetBranchAddress("col4",&col4); 
 ntuple->SetBranchAddress("col5",&col5); 
 ntuple->SetBranchAddress("col6",&col6); 
 ntuple->SetBranchAddress("col7",&col7); 
 ntuple->SetBranchAddress("col8",&col8); 
 ntuple->SetBranchAddress("col9",&col9); 
 ntuple->SetBranchAddress("col10",&col10); 
 ntuple->SetBranchAddress("col11",&col11); 
 ntuple->SetBranchAddress("col12",&col12); 
 ntuple->SetBranchAddress("col13",&col13); 
 ntuple->SetBranchAddress("col14",&col14); 
 ntuple->SetBranchAddress("col15",&col15); 
 ntuple->SetBranchAddress("col16",&col16); 
 ntuple->SetBranchAddress("col17",&col17); 
 
 cout << "Set Branch Address done" <<endl;

 
 // Get number of entries in ntuple (number of lines in the output fluka file) 
 int nentries=ntuple->GetEntries();
 cout << "Number  of entries: " << nentries << endl;
 int ievt=-1; // Identify the event
 // Create a TTree
 TTree *tree = new TTree("MC","Raw event MC");
 ALEvent *e = new ALEvent();
 // Create a branch with event
 tree->Branch("event",&e);  

 int j=0;
 int nh=0;
 double EneFoam=0;
 double EneShell=0;
 double EneT1=0;
 double EneT2=0; 
 double EneT3=0; 
 double EneT4=0; 
 double Eneg=0; 
 double timeFoam=0;
 double timeShell=0;
 double timeT1=0;
 double timeT2=0; 
 double timeT3=0; 
 double timeT4=0;
 int nOptPhCK=0;  //Number of Optical phtoon produced in T2
 int prevreg=0; 
 double timeg=0; 
 double nFoam=0;				//number of segments of track in insulation foam
 double nShell=0;				//number of segments of track in shell
 double nT1=0; 						//number of segments of track in T1
 double nT2=0;
 double nT3=0; 						
 double nT4=0; 
 double ng=0; 
 int iFoam=0;
 int iShell=0;
 int iT1=0; 
 int iT2=0;
 int iT3=0; 
 int iT4=0; 
 int ig=0;
 
 float oldx=0;
 float oldy=0;
 float oldz=0;
 
//Temporary variable for internal triggers 
 int*Titmp=new int[7];
 for (int i=0;i<7;i++)Titmp[i]=0;
 
 int nL=0;
 int iL=0;
 ALTckhit* h=new ALTckhit();

 //Vector of segments to discretize
 //Clear every time there is a new layer crossed
 vector<float> X;
 vector<float> Y;
 vector<float> Z;
 vector<float> CZ;
 vector<float> T;
 vector<float> EDEP;
 float CoordDisc=-999;
 int fstrip=-1;
 int fstripID=-1;
 int nstrip=0;
 int chip=-1;
 //Loop over the entries and create the events
 // This is the main loop
 for (int i=0;i<nentries;i++)
   // for (int i=0;i<30;i++)
	{
    ntuple->GetEntry(i); //Load the entry i in the variables  
    if(i%100000==0) cout << "entry: " << i <<endl;
   // cout << "i = " << i << ", ncase " << ncase << ", mreg " << mreg << " type = " << type << " p = " << pMC*1000 <<" MeV,  x=" << x << ", y=" << y << ", z = " << z << endl; 

	 if(col17==-1) {
		 ncase=col1;
		 type=col2;
		 EkMC=col3;
		 age=col4;
		 x=col5;
		 y=col6;
		 z=col7;
		 cx=col8;
		 cy=col9;
		 cz=col10;
		 typePP=col11;
		 EkPP=col12;
		 ZenPP=col13;
		 AziPP=col14;
		 CoLatSP=col15;
		 CoLonSP=col16;
		 flagSP=col17;}

		else {
		 ncase=col1;
		 mreg=col2;
		 mtrack=col3;
		 type=col4;
		 age=col5;
		 eMC=col6;
		 x=col7;
		 y=col8;
		 z=col9;
		 cx=col10;
		 cy=col11;
		 cz=col12;
		 Edep=col13;
		 flag=col14;
			}
/*
   pMC = TMath::Sqrt(EkMC*EkMC + 2*mass);    
   e->add_posX(x);
   e->add_posY(y);
   e->add_posZ(z);
   e->add_posType(type);
   e->add_posAge(age);
   e->add_posP(pMC);
   pMC = TMath::Sqrt(EkMC*EkMC + 2*mass); 
*/
   if(mreg==11) mreg11++;
   if(mreg==12) mreg12++;
   if(mreg==13) mreg13++;
   if(mreg==14) mreg14++;
   if(mreg==15) mreg15++;
   if(mreg==16) mreg16++;
   if(mreg==17) mreg17++;

    ////////////////////////////////// 
    // Triggers
    ////////////////////////////////// 
    if(prevreg==ShellReg[0] && prevreg!=mreg)//Foam insulation
     {
      if(nFoam>0)
       {
        e->add_EneIsofoam(EneFoam); e->add_timeIsofoam(timeFoam);
		EneFoam=0;timeFoam=0;nFoam=0;iFoam=0;
       }
     }	 
    if(prevreg==ShellReg[1] && prevreg!=mreg)//Aluminium Shell
     {
      if(nShell>0)
       {
        e->add_EneShell(EneShell); e->add_timeShell(timeShell);
	EneShell=0;timeShell=0;nShell=0;iShell=0;
       }
     }	
    if(prevreg==TrigReg[0] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T1
     {
      if(nT1>0)
       {
        e->add_EneT1(EneT1); e->add_timeT1(timeT1);
        if (EneT1 > TrigThresh[0]) e->set_T1(true);
        EneT1=0;timeT1=0;nT1=0;iT1=0;
       }
     }
   if(prevreg==TrigReg[1] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T2
     {
      if(nT2>0)
       {
        e->add_EneT2(EneT2);
        e->add_timeT2(timeT2);
        if (EneT2 > TrigThresh[1]) e->set_T2(true);
        e->set_NphCK(nOptPhCK);
        EneT2=0;timeT2=0;nT2=0;iT2=0;nOptPhCK=0;
       }
     }
    if(prevreg==TrigReg[2] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T3
     {
      if(nT3>0)
       {
        e->add_EneT3(EneT3); e->add_timeT3(timeT3);
        if(EneT3 > TrigThresh[2]) e->set_T3(true);
        EneT3=0;timeT3=0;nT3=0;iT3=0;
       }
     }
    if(prevreg==TrigReg[3] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//T4
     {
      if(nT4>0)
       {
        e->add_EneT4(EneT4); e->add_timeT4(timeT4);
        if (EneT4 > TrigThresh[3]) e->set_T4(true);
        EneT4=0;timeT4=0;nT4=0;iT4=0;
       }
     }
    if(prevreg==GReg[0] && (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))//guard
     {
      if(ng>0)
       {
        e->add_Eneg(Eneg); e->add_timeg(timeg);
        if (Eneg > TrigThresh[4]) e->set_guard(true);
        Eneg=0;timeg=0;ng=0;ig=0;
       }
     }

    ////////////////////////////////// 
    // Tracker
    ////////////////////////////////// 
     
//    if(prevreg>=TckReg[0] && prevreg<=TckReg[6]&& (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))
    if(prevreg>=TckReg[0] && prevreg<=TckReg[6]&& (prevreg!=mreg || abs(oldx-x)>0.1 || abs(oldy-y)>0.1 || abs(oldz-z)>0.1))
     {       
      
      //APPLY DISCRETISATION
      CoordDisc=-999.;
      fstrip=-1;
      nstrip=0;
      chip=-1;
//      cout << "size of vector x = " << (int)X.size() << endl;      
      if((int)X.size()!=0) 
       {
//	cout << "applying discretization" << endl;
        if(h->get_L()==0) nL0++;
        if(h->get_L()==1) nL1++;
        if(h->get_L()==2) nL2++;
        if(h->get_L()==3) nL3++;
        if(h->get_L()==4) nL4++;
        if(h->get_L()==5) nL5++;
        if(h->get_L()==6) nL6++;

        //OffsetLL=ladderOffsetLeft[h->get_L()];
       // OffsetRL=ladderOffsetRight[h->get_L()];
        OffsetLL=OffsetRL=0;
	CoordDisc=Discretize(h->get_L(),X,Y,Z,CZ,T,EDEP,&chip,&fstrip,&fstripID,&nstrip,OffsetLL,OffsetRL,true);
//        cout << "CoordDisc=" <<  CoordDisc <<endl; 
        if(CoordDisc!=-999.)
         {   
          //Update parameters of the hit after discretisation
          if(h->get_L()==0||h->get_L()==4||h->get_L()==6) h->set_x(CoordDisc); //Non bending plane
          if(h->get_L()==1||h->get_L()==2||h->get_L()==3||h->get_L()==5) h->set_y(CoordDisc);//Bending plane
          h->set_fstripID(fstripID);//from 0 to 767
          h->set_fstrip(fstrip);//from 0 to 63
          h->set_nstrips(nstrip);//number of strips in cluster
          h->set_chip(chip);//chip 0 to 11 of the first strip (strip lowest index)
          if(nL>0){h->set_cx(h->get_cx()/nL);h->set_cy(h->get_cy()/nL);h->set_cz(h->get_cz()/nL);}
          //End of the hit.
          //Add the hit to the event
	  //index k of hit starts at 0
          h->set_k(nh);
	  nh++;
          if(h->get_x()==-999) cout << "X coord is -999. Event " << i <<endl;
          if(h->get_y()==-999) cout << "Y coord is -999. Event " << i <<endl;
          if(h->get_z()==-999) cout << "Z coord is -999. Event " << i <<endl;
          e->add_hit(h);
          Titmp[(int)h->get_L()]=1;
         } 
	else
         {
	  //  if(h->get_L()==0||h->get_L()==4||h->get_L()==6)  cout << " L = " << h->get_L() << ", x = " << x << ", z = "  << z << endl;
          // if(h->get_L()==1||h->get_L()==2||h->get_L()==3||h->get_L()==5)  cout << " L = " << h->get_L() << ", y = " << y << ", z = "  << z << endl;

	  if(h->get_L()==0) rL0++;
          if(h->get_L()==1) rL1++;
          if(h->get_L()==2) rL2++;
          if(h->get_L()==3) rL3++;
          if(h->get_L()==4) rL4++;
          if(h->get_L()==5) rL5++;
          if(h->get_L()==6) rL6++;
         }
		
        X.clear();
        Y.clear();
        Z.clear();
        CZ.clear();
        T.clear();
        EDEP.clear();
       //Check the layer for internal trigger       
       }
 

      //Reset the hit structure
      h=new ALTckhit();
      nL=0;
      iL=0;
     }  
    if(ncase!=ievt)//Start a new event
     {
      j++;
      nh=0;

      if(i==0)//if first event
       {
        ievt=ncase; 
        delete e;
        e=new ALEvent();
        e->set_eventnumber(j);
        e->set_ncase(ncase);
        //First line  with a given ncase contains the source particle info
        e->set_typeMC(type);
        //Calculate kinetic energy from momentum at point of injection
        float mass=0.000511;						//electron mass in GeV
        if(type==11 || type ==10)	mass=0.10566;			//muon mass in GeV
	if(type==1)			mass=0.93827;			//proton mass in GeV
	if(type==-6)			mass=3.72739;			//alpha mass in GeV 
        pMC = TMath::Sqrt(EkMC*EkMC+2*mass);   
        e->set_pMC(pMC);   
        e->set_EkMC(EkMC);
        e->set_X0MC(x);
        e->set_Y0MC(y);
        e->set_Z0MC(z);
        e->set_CX0MC(cx);
        e->set_CY0MC(cy);
        e->set_CZ0MC(cz);
	e->set_typePP(typePP);
	e->set_EkPP(EkPP);
	e->set_ZenPP(ZenPP);
	e->set_AziPP(AziPP);
	e->set_CoLatSP(CoLatSP);
	e->set_CoLonSP(CoLonSP);
        e->add_posX(x);
        e->add_posY(y);
        e->add_posZ(z);
        e->add_posCX(cx);
        e->add_posCY(cy);
        e->add_posCZ(cz);
        e->add_posType(type);
        e->add_posAge(age);
        e->add_posP(pMC);
		 
	for(int ij=0;ij<7;ij++)Titmp[ij]=0;
       }//if
      else
       {
        //Fill internal trigger information of the previous event
        int t=0;
        for (int ij=0;ij<7;ij++) t+=Titmp[ij]*(int)TMath::Power(2,ij);
        e->set_Ti(t);  
        //Fill the previous event in the output file  
        tree->Fill(); 
        delete e;
        e=new ALEvent();
        h=new ALTckhit();
        nh=0;
        ievt=ncase; 

        e->set_eventnumber(j);
        e->set_ncase(ncase);
        //First line with a given ncase contains the source particle info
	e->set_typeMC(type);
        e->set_pMC(pMC);
        //Calculate momentum from kinetic energy at point of injection
        float mass=0.000511;						//electron mass in GeV
        if(type==11 || type ==10)	mass=0.10566;			//muon mass in GeV
	if(type==1)			mass=0.93827;			//proton mass in GeV
	if(type==-6)			mass=3.72739;			//alpha mass in GeV 
        pMC = TMath::Sqrt(EkMC*EkMC + 2*mass);
        e->set_pMC(pMC);       
        e->set_EkMC(EkMC);
        e->set_X0MC(x);
        e->set_Y0MC(y);
        e->set_Z0MC(z);
        //cout << " x0 = " << x << " y0 = " << y << "  z0 = " << z << endl;
        e->set_CX0MC(cx);
        e->set_CY0MC(cy);
        e->set_CZ0MC(cz);
	//  Fill in info for PP
	e->set_typePP(typePP);
	e->set_EkPP(EkPP);
	e->set_ZenPP(ZenPP);
	e->set_AziPP(AziPP);
	e->set_CoLatSP(CoLatSP);
	e->set_CoLonSP(CoLonSP);
        e->add_posX(x);
        e->add_posY(y);
        e->add_posZ(z);
        e->add_posCX(cx);
        e->add_posCY(cy);
	e->add_posCZ(cz);        
	e->add_posType(type);
        e->add_posAge(age);
        e->add_posP(pMC);

        //Reset internal trigger
        for (int ij=0;ij<7;ij++)Titmp[ij]=0;
       }     //else
     }  //if
    
    ////////////////////////////////// 
    //Check to see if the particle has crossed one of the scintallator 
    //(r19=air, r1=T1, r6=T3, r7=guard, r11=Tracker1 ,..., r17=Tracker7, r18=T4)
    ////////////////////////////////// 
	//Shell & Insulation
    if (mreg == ShellReg[0]&& type!=7 && type !=0 && type!=8 && type!=211) //Insulating foam
      {
       if(i==iFoam+1){EneFoam+=Edep;nFoam++;}
       else{nFoam=1;EneFoam=Edep;timeFoam=age;}
       iFoam=i;
      }
    if (mreg == ShellReg[1]&& type!=7 && type !=0 && type!=8 && type!=211)//Aluminium Shell
      {
       if(i==iShell+1){EneShell+=Edep;nShell++;}
       else{nShell=1;EneShell=Edep;timeShell=age;}
       iShell=i;
      }	     
    //Triggers
    if (mreg == TrigReg[0]&& type!=7 && type !=0 && type!=8 && type!=211)//T1
      {
       if(i==iT1+1){EneT1+=Edep;nT1++;}
       else{nT1=1;EneT1=Edep;timeT1=age;}
       iT1=i;
      }
    if (mreg == TrigReg[1]&& type!=7 && type !=0 && type!=8 && type!=211)//T2
      {
       if(i==iT2+1)
        {
         EneT2+=Edep;
         nT2++;
        }
       else
        {
         nT2=1;
         EneT2=Edep;
         timeT2=age;
        }
       if(type==-1 && mtrack==1) //Production of Opticl Photon in CK, if mtrack ==0, this is the line that corresponds to the absorption of the  optical photon
        {
         nOptPhCK+=1;   
        }  
       iT2=i;
      }

    if (mreg == TrigReg[2]&& type!=7 && type !=0 && type!=8 && type!=211)//T3
      {
       if(i==iT3+1){EneT3+=Edep;nT3++;}
       else{nT3=1;EneT3=Edep;timeT3=age;}
       iT3=i;
      }
    if (mreg == TrigReg[3]&& type!=7 && type !=0 && type!=8 && type!=211)//T4
      {
       if(i==iT4+1){EneT4+=Edep;nT4++;}
       else{nT4=1;EneT4=Edep;timeT4=age;}
       iT4=i;
      }
    if (mreg == GReg[0]&& type!=7 && type !=0 && type!=8 && type!=211)//guard
      {
       if(i==ig+1){Eneg+=Edep;ng++;}
       else{ng=1;Eneg=Edep;timeg=age;}
       ig=i;
      }
    
    for (int ii=0;ii<7;ii++)
      {
//       if(mreg == TckReg[ii])
       if(mreg == TckReg[ii] && type!=7 && type !=0 && type!=8 && type!=211) //no photon,no neutron, no zero, no EM-deposit
        {
         if(i==iL+1 && abs(oldx-x)<=0.1 && abs(oldy-y)<=0.1 && abs(oldz-z)<=0.1 ) //if the segment of track is not the first of the hit 
          {
           h->set_DeltaE(h->get_DeltaE()+Edep);// Add the energy deposited in the segment
           h->set_cx(h->get_cx()+cx);
           h->set_cy(h->get_cy()+cy); 
           h->set_cz(h->get_cz()+cz);
           h->set_xout(x);
           h->set_yout(y);   
           h->set_zout(z);
           h->set_x((h->get_xin() + x)/2.);
           h->set_y((h->get_yin() + y)/2.);
           h->set_z((h->get_zin() + z)/2.);	
           X.push_back(x);
           Y.push_back(y);
           Z.push_back(z);
           T.push_back(type);
           CZ.push_back(cz);
           EDEP.push_back(Edep);
           nL++;//Increment the number of segments that compose the hit
          }//if
         else //if this segment of track is the first part of the hit
          {
           //Fill the data at the entrance of the particle in the tracker layer
           h->set_mregMC(mreg); 
           //set layer number (for MC and data)
           int L = (int)mreg%11;
	   h->set_L(L);
           h->set_mtrackMC(mtrack);
           h->set_typeMC((int)type);
           h->set_eMC(eMC); 				//total energy of particle crossing region
           h->set_flag(flag);
           h->set_k(nh);
           h->set_age(age);
           h->set_cx(cx);
           h->set_cy(cy); 
           h->set_cz(cz);
           h->set_xin(x);
           h->set_yin(y);   
           h->set_zin(z);
           h->set_xout(x);
           h->set_yout(y);   
           h->set_zout(z);
	   // if(h->get_mregMC()!=mreg)   cout << "Region mreg  = " << mreg << " , mregMC " << h->get_mregMC() <<  endl;
           h->set_x(x);
           h->set_y(y);
           h->set_z(z);
           X.push_back(x);
           Y.push_back(y);
           Z.push_back(z);
           CZ.push_back(cz);
           T.push_back(type);
           EDEP.push_back(Edep);
           h->set_DeltaE(Edep);
           nL=1;
          } //else
         iL=i; 
         ii=7;
        }  //if
      }//ii
    prevreg=mreg;
    oldx=x;
    oldy=y;
    oldz=z;
   }//i

  //Fill the last event
  if(nFoam>0)
   {
    e->add_EneIsofoam(EneFoam); e->add_timeIsofoam(timeFoam); 
   }  	 
  if(nShell>0)
   {
    e->add_EneShell(EneShell); e->add_timeShell(timeShell);
   }
  if(nT1>0)
   {
    e->add_EneT1(EneT1);e->add_timeT1(timeT1); 
    if (EneT1 > TrigThresh[0]) e->set_T1(true);
   }
  if(nT2>0)
   {
    e->add_EneT2(EneT2);
    //e->add_timeT2(timeT2);
    if (EneT2 > TrigThresh[1]) e->set_T2(true);
   }
  if(nT3>0)
   {
    e->add_EneT3(EneT3);e->add_timeT3(timeT3);
    if (EneT3 > TrigThresh[2]) e->set_T3(true);
   }
  if(nT4>0)
   {
    e->add_EneT4(EneT4);e->add_timeT4(timeT4);
    if(EneT4 > TrigThresh[3]) e->set_T4(true);
   }
  if(ng>0)
   {
    e->add_Eneg(Eneg);e->add_timeg(timeg);
    if (Eneg > TrigThresh[4]) e->set_guard(true);
   }
 
 
 //APPLY DISCRETISATION
 CoordDisc=-999.;
 fstrip=-1;
 fstripID=-1;
 nstrip=0;
 chip=-1;
 
 if((int)X.size()!=0) 
  {
   //CoordDisc=Discretize(h->get_L(),X,Y,Z,CZ,T,EDEP,&chip,&fstrip,&fstripID,&nstrip,OffsetLL,OffsetRL,true);
  
   OffsetLL=ladderOffsetLeft[h->get_L()];
   OffsetRL=ladderOffsetRight[h->get_L()];
   OffsetLL=OffsetRL=0;

   CoordDisc=Discretize(h->get_L(),X,Y,Z,CZ,T,EDEP,&chip,&fstrip,&fstripID,&nstrip,OffsetLL,OffsetRL,true);   
      
   
   cout << "CoordDisc2=" <<  CoordDisc <<endl; 
   
   if(CoordDisc!=-999.)
    {   
     //Update parameters of the hit after discretisation
     if(h->get_L()==0||h->get_L()==4||h->get_L()==6) h->set_x(CoordDisc); //Non bending plane
     if(h->get_L()==1||h->get_L()==2||h->get_L()==3||h->get_L()==5) h->set_y(CoordDisc);//Bending plane
     h->set_fstripID(fstripID);//from 0 to 767
     h->set_fstrip(fstrip);//from 0 to 62
     h->set_nstrips(nstrip);//number of strips in cluster
     h->set_chip(chip);//chip 0 to 11
     if(nL>0){h->set_cx(h->get_cx()/nL);h->set_cy(h->get_cy()/nL);h->set_cz(h->get_cz()/nL);}//add the last hit 
     //End of the hit.
     //Add the hit to the event
     h->set_k(nh);
     nh++;
     e->add_hit(h);
    //Check the layer for internal trigger
     Titmp[(int)h->get_L()]=1;
    } 
   X.clear();
   Y.clear();
   Z.clear();
   CZ.clear();
   T.clear();
   EDEP.clear();
//   delete h;
      
  }

 int tt=0;
 for (int ij=0;ij<7;ij++) tt+=Titmp[ij]*(int)TMath::Power(2,ij);
 //cout << "before setTi" <<endl;

 e->set_Ti(tt);
 //cout << "before tree fill" <<endl;

 tree->Fill();
// cout << "before delete e" <<endl;

 delete e;
 //cout << "after delete e" <<endl;

 //Write tree in output file 
 tree->Write();

 //Close files
 fileout->Close();
 file->Close();
 //Delete TTree
 // tree->Delete(); int*TckReg=new int[7];
 delete[] TckReg;
 delete[] TrigReg;
 delete[] GReg;
 delete[] ShellReg;
 delete[] TckZPos;
 delete[] TrigThresh;
 delete[] GuardThresh;
 delete[] Titmp;


 return 1;
}

