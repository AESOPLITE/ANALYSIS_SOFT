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
 ncase=0; 
 typeMC=-99; //type of particle
 EkMC=0;   //kinetic energy of the particle
 X0MC=Y0MC=Z0MC=0;//Coordinates of the partcle at the injection point 
 CX0MC=CY0MC=CZ0MC=0; //Incidence cosines of the partcle at the injection point 
 Nhits=0; //Number of hits in the event
 typereco=-999; //type of particle
 Ekreco=-999;   //kinetic energy of the particle
 X0reco=Y0reco=Z0reco=0;//Coordinates of the partcle at the injection point 
 CX0reco=CY0reco=CZ0reco=0; //Incidence cosines of the partcle at the injection point 
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
   X0reco =e->get_X0reco();
   Y0reco =e->get_Y0reco();
   Z0reco =e->get_Z0reco();
   CX0reco =e->get_CX0reco();
   CY0reco =e->get_CY0reco();
   CZ0reco =e->get_CZ0reco();
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

