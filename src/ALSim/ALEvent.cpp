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
 
}


