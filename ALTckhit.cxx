
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

  
   // Time data from the corressponding "ASI" LINE
   y=-1;//Year from first ASI line of the event
   m=-1;//Month from first ASI line of the event
   d=-1;//Day from first ASI line of the event
   hour=-1;//Hour from first ASI line of the event
   mi=-1;//Minute from first ASI line of the event
   s=-1;//Second from first ASI line of the event
  
  
   //Raw Data from cluster information
   L=-1;         //Layer from 0 to 6 top to bottom
   chip=-1;      //Chip ID: 0 to 11
   nstrips=-1;   //Number of strips in the cluster
   nstripsNC=-1;   //Number of strips in the next chip if it is a boundary cluster (except boundary chips 5-6)
   fstrip=-1;    //First strip ID from 0 to 63 
   fstripID=-1;  //First strip on the layer from 0 767
   noisy=0;
   parityerr[0]=parityerr[1]=-1;//Parity error of the clusters that make the hit
   chiperr[0]=chiperr[1]=-1;//Chip error of the clusters that make the hit
   overflow[0]=overflow[1]=-1;//overflow of the clusters that make the hit

   //Coordinates of the cluster in cm determined from the raw data
   x=-999.;
   y=-999.;
   z=-999.;
 
  //Pattern Reco info
   xPR=0;
   yPR=0;
   zPR=0;
   cxPR=0;
   cyPR=0;
   czPR=0;
   fGhost=false;
   flagPR=false;

  //Reconstructed information
   xreco=-999;      //x coordinate of hit 
   yreco=-999;      //y coordinate of hit 
   zreco=-999;      //z coordinate of hit 
   agereco=-99;    //Age of particle  
   cxreco=-999;     //cosineX of momentum 
   cyreco=-999;     //cosineY of momentum 
   czreco=-999;     //cosineZ of momentum    
   ereco=-999;     //kinetic energy 
   fUsed=false;    //if site accepted by filter
   k=-999;     //kth hit in event

}


//Copy function
void ALTckhit::Copy(ALTckhit* h)
{
  //Single variables

   mregMC=h->get_mregMC();
   mtrackMC=h->get_mtrackMC();
   typeMC=h->get_typeMC();
   eMC=h->get_eMC();     
   xin=h->get_xin();     
   yin=h->get_yin();      
   zin=h->get_zin();      
   xout=h->get_xout();     
   yout=h->get_yout();     
   zout=h->get_zout();     
   age=h->get_age();    
   cx=h->get_cx();     
   cy=h->get_cy();     
   cz=h->get_cz();     
   flag=h->get_flag();    
   DeltaE=h->get_DeltaE();    

  
   // Time data from the corressponding "ASI" LINE
   y=h->get_year();
   m=h->get_m();
   d=h->get_d();
   hour=h->get_hour();
   mi=h->get_mi();
   s=h->get_s();
  
  
   //Raw Data from cluster information
   L=h->get_L();         
   chip=h->get_chip();   
   nstrips=h->get_nstrips();  
   nstripsNC=h->get_nstripsNC();  
   fstrip= h->get_fstrip();     
   fstripID=h->get_fstripID();    
   noisy=h->get_noisy();  
   parityerr[0]=parityerr[1]=h->get_parityerr(0);  
   chiperr[0]=chiperr[1]=h->get_chiperr(0);  
   overflow[0]=overflow[1]=h->get_overflow(0);  

   //Coordinates of the cluster in cm determined from the raw data
   x=h->get_x();
   y=h->get_y();
   z=h->get_z();
 
  //Pattern Reco info
   xPR=h->get_xPR();
   yPR=h->get_yPR();
   zPR=h->get_zPR();
   cxPR=h->get_cxPR();
   cyPR=h->get_cyPR();
   czPR=h->get_czPR();
   fGhost=h->get_fGhost();
   flagPR=h->get_flagPR();

  //Reconstructed information
   xreco=h->get_xreco();
   yreco=h->get_yreco();
   zreco=h->get_zreco();     
   agereco=h->get_agereco();    
   cxreco=h->get_cxreco();    
   cyreco=h->get_cyreco();    
   czreco=h->get_czreco();    
   ereco=h->get_ereco();     
   fUsed=h->get_fUsed(); 
   k=h->get_k();     


}
