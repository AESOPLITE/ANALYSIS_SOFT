////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, January 17, 2017
////////////////////////////////////////////////////////////////////////////////////////// 

#include "ALTckhit.h"

ClassImp(ALTckhit)


//Constructors
ALTckhit::ALTckhit()// Default
{
   mregMC=0; // region if MC
   mtrackMC=0; // mtrack of track
   typeMC=0;   //type of particle if MC
   eMC=0;     //total energy if MC (entrance of the tracker layer)
   xin=0;      //x coordinate of hit MC (entrance of the tracker layer)
   yin=0;      //y coordinate of hit MC (entrance of the tracker layer)
   zin=0;      //z coordinate of hit MC (entrance of the tracker layer)
   xout=0;     //x coordinate of hit MC (exit the tracker layer)
   yout=0;     //y coordinate of hit MC (exit the tracker layer)
   zout=0;     //z coordinate of hit MC (exit the tracker layer)
   age=0;    //Age of particle if MC  (entrance of the tracker layer)
   cx=0;     //cosineX of momentum MC (Average in the tracker layer)
   cy=0;     //cosineY of momentum MC (Average in the tracker layer)
   cz=0;     //cosineZ of momentum MC (Average in the tracker layer)
   flag=0;     //flag value: Track=1, Hit=0
   DeltaE=0;     //Energy along the track or at hit location Hit=0  

  //Reconstructed information
   xreco=0;      //x coordinate of hit 
   yreco=0;      //y coordinate of hit 
   zreco=0;      //z coordinate of hit 
   agereco=0;    //Age of particle  
   cxreco=0;     //cosineX of momentum 
   cyreco=0;     //cosineY of momentum 
   czreco=0;     //cosineZ of momentum    
   ereco=0;     //kinetic energy 
   k=0;     //kth hit in event
}