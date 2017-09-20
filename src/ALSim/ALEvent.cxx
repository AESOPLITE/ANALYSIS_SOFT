////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 

#include "ALEvent.h"
ClassImp(ALEvent)
ClassImp(ALTckhit)

//Constructors
ALEvent::ALEvent()// Default
{
 eventnumber=0; //Event number
 
 yPHA=-1;//Year from PHA line linked to the event
 mPHA=-1;//Month from PHA line linked to the event
 dPHA=-1;//Day from PHA line linked to the event
 hPHA=-1;//Hour from PHA line linked to the event
 miPHA=-1;//Minute from PHA line linked to the event
 sPHA=-1;//Second from PHA line linked to the event
 GoPHA=-1;//Go counter from PHA line linked to the event
 tPHA=-1;//timer from from PHA line linked to the event
 yEVT=-1;//Year from EVT line linked to the event
 mEVT=-1;//Month from EVT line linked to the event
 dEVT=-1;//Day from EVT line linked to the event
 hEVT=-1;//Hour from EVT line linked to the event
 miEVT=-1;//Minute from EVT line linked to the event
 sEVT=-1;//Second from EVT line linked to the event  
 EVT="";//Data from EVT line linked to the event
 for(int i=0;i<7;i++) L[i]=string();//Data from  ASI lines of the event
   
 ncase=0; 
 typeMC=-99; //type of particle
 EkMC=0;   //kinetic energy of the particle
 X0MC=Y0MC=Z0MC=0;//Coordinates of the partcle at the injection point 
 CX0MC=CY0MC=CZ0MC=0; //Incidence cosines of the partcle at the injection point 
 Nhits=0; //Number of hits in the event
 typereco=-999; //type of particle
 Ekreco=-999;   //kinetic energy of the particle
 p0reco=-999;  //momentum of the particle
 X0reco=Y0reco=Z0reco=0;//Coordinates of the partcle at the injection point 
 CX0reco=CY0reco=CZ0reco=0; //Incidence cosines of the partcle at the injection point 
 ndf=0;
 chi2=cl=0;
 d0=phi0=cpa=dz=tanl=0;
 d0err2=phi0err2=cpaerr2=dzerr2=tanlerr2=0;
	
 EkPR=-999; 
 p0PR =-999;
 a=b=c=0;
 inter=slope=0;
 chi2B=chi2NB=clB=clNB=0;
 deflec=0;
 //Triggers: default is false
 T1=false;
 T2=false;
 T3=false;
 T4=false;
 guard=false;
 Ti=0;
}


void ALEvent::Copy(ALEvent* e)
{
  //Single variables 
  
  eventnumber =e->get_eventnumber();
  yPHA=e->get_yPHA();
  mPHA=e->get_mPHA();
  dPHA=e->get_dPHA();
  hPHA=e->get_hPHA();
  miPHA=e->get_miPHA();
  sPHA=e->get_sPHA();
  GoPHA=e->get_GoPHA();
  tPHA=e->get_tPHA();
  yEVT=e->get_yEVT();
  mEVT=e->get_mEVT();
  dEVT=e->get_dEVT();
  hEVT=e->get_hEVT();
  miEVT=e->get_miEVT();
  sEVT=e->get_sEVT();
  EVT=e->get_EVT();
  for(int i=0;i<7;i++)L[i]=e->get_L(i);



   ncase =e->get_ncase();
   typeMC =e->get_typeMC();
   EkMC =e->get_EkMC();
   X0MC =e->get_X0MC();
   Y0MC =e->get_Y0MC();
   Z0MC =e->get_Z0MC();
   CX0MC =e->get_CX0MC();
   CY0MC =e->get_CY0MC();
   CZ0MC =e->get_CZ0MC();
   Nhits =e->get_Nhits();
   typereco =e->get_typereco();
   Ekreco =e->get_Ekreco();
   p0reco =e->get_p0reco();
   X0reco =e->get_X0reco();
   Y0reco =e->get_Y0reco();
   Z0reco =e->get_Z0reco();
   CX0reco =e->get_CX0reco();
   CY0reco =e->get_CY0reco();
   CZ0reco =e->get_CZ0reco();
   ndf =e->get_ndf();
   chi2 =e->get_chi2();
   cl =e->get_cl();
   d0 =e->get_d0();
   phi0 =e->get_phi0();
   cpa =e->get_cpa();
   dz =e->get_dz();
   tanl =e->get_tanl();
   d0err2 =e->get_d0err2();
   phi0err2 =e->get_phi0err2();
   cpaerr2 =e->get_cpaerr2();
   dzerr2 =e->get_dzerr2();
   tanlerr2 =e->get_tanlerr2();

   EkPR = e->get_EkPR();
   p0PR = e->get_p0PR();
   a = e->get_a();
   b = e->get_b();
   c = e->get_c();
   inter = e->get_inter();
   slope = e->get_slope();
   chi2B = e->get_chi2B();
   chi2NB = e->get_chi2NB();
   clB = e->get_clB();
   clNB = e->get_clNB();
   deflec = e->get_deflec();
   
   T1 =e->get_T1();
   T2 =e->get_T2();
   T3 =e->get_T3();
   T4 =e->get_T4();
   guard =e->get_guard();
   Ti =e->get_Ti();

   
   //Vectors of double
   for(int i=0;i<(int)(e->get_EneT1()).size();i++) EneT1.push_back((e->get_EneT1()).at(i));
   for(int i=0;i<(int)(e->get_EneT2()).size();i++) EneT2.push_back((e->get_EneT2()).at(i));
   for(int i=0;i<(int)(e->get_EneT3()).size();i++) EneT3.push_back((e->get_EneT3()).at(i));
   for(int i=0;i<(int)(e->get_EneT4()).size();i++) EneT4.push_back((e->get_EneT4()).at(i));
   for(int i=0;i<(int)(e->get_Eneg()).size();i++) Eneg.push_back((e->get_Eneg()).at(i));
   for(int i=0;i<(int)(e->get_timeT1()).size();i++) timeT1.push_back((e->get_timeT1()).at(i));
   for(int i=0;i<(int)(e->get_timeT2()).size();i++) timeT2.push_back((e->get_timeT2()).at(i));
   for(int i=0;i<(int)(e->get_timeT3()).size();i++) timeT3.push_back((e->get_timeT3()).at(i));
   for(int i=0;i<(int)(e->get_timeT4()).size();i++) timeT4.push_back((e->get_timeT4()).at(i));
   for(int i=0;i<(int)(e->get_timeg()).size();i++) timeg.push_back((e->get_timeg()).at(i));
 
   //Vectors of ALTckhit
   for(int i=0;i<(int)(e->get_hits()).size();i++) hits.push_back((e->get_hits()).at(i));


 
}
