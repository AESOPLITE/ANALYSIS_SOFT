////////////////////////////////////////////////////////////////////////////////////////// 
///    Author: Pierre-Simon Mangeard, psmangeard@gmail.com
///    Department of Physics and Astronomy, University of Delaware, October 28, 2016
////////////////////////////////////////////////////////////////////////////////////////// 
#ifndef __ALEVENT__
#define __ALEVENT__
#include "headers.h"
#include "ALTckhit.h"

/*
struct ALHit //Track or hits
{
 float mregMC; // region if MC
 float mtrackMC; // mtrack of track
 int typeMC;   //type of particle if MC
 float eMC;     //total energy if MC (entrance of the tracker layer)
 float xin;      //x coordinate of hit MC (entrance of the tracker layer)
 float yin;      //y coordinate of hit MC (entrance of the tracker layer)
 float zin;      //z coordinate of hit MC (entrance of the tracker layer)
 float xout;     //x coordinate of hit MC (exit the tracker layer)
 float yout;     //y coordinate of hit MC (exit the tracker layer)
 float zout;     //z coordinate of hit MC (exit the tracker layer)
 float age;    //Age of particle if MC  (entrance of the tracker layer)
 float cx;     //cosineX of momentum MC (Average in the tracker layer)
 float cy;     //cosineY of momentum MC (Average in the tracker layer)
 float cz;     //cosineZ of momentum MC (Average in the tracker layer)
 int flag;     //flag value: Track=1, Hit=0
 float DeltaE;     //Energy along the track or at hit location Hit=0
 //Reconstructed information
 float xreco;      //x coordinate of hit 
 float yreco;      //y coordinate of hit 
 float zreco;      //z coordinate of hit 
 float agereco;    //Age of particle  
 float cxreco;     //cosineX of momentum 
 float cyreco;     //cosineY of momentum 
 float czreco;     //cosineZ of momentum    
 float ereco;     //kinetic energy 
 int k;     //kth hit in event 
};
*/

class ALEvent:public TObject
{
 private:
  
   int eventnumber; //Event number
   //Monte Carlo information: Truth variable names finish with MC 
   int ncase; 
   int typeMC; //type of particle
   double EkMC;   //kinetic energy of the particle
   double X0MC,Y0MC,Z0MC;//Coordinates of the partcle at the injection point 
   double CX0MC,CY0MC,CZ0MC; //Incidence cosines of the partcle at the injection point 
   // Hits information
   int Nhits; //Number of hits in the event
   
   
   //Pattern Recognition info
   
   double EkPR, p0PR;			//kinetic energy and momentum of particle from least squares fit
   double chi2NB, chi2B, clNB, clB;	//chi2 of parabolic/linear fit in bending/nonbending plane
   double a, b, c;			//parameters of parabolic fit ( y(x) = a + bx + cx*x)
   double slope, inter;			//parameters of linear fit 
  
  
  //Reconstruction information: variables finish with reco
   int typereco;                    //type of particle
   double Ekreco, p0reco;           //kinetic energy and momentum of the particle
   double X0reco,Y0reco,Z0reco;     //Coordinates of the partcle at the injection point 
   double CX0reco,CY0reco,CZ0reco;  //Incidence cosines of the partcle at the injection point 
   int ndf;                         // number of degrees of freedom
   double chi2, cl;                 //chi2 of fit and confidence level (cl = Prob(chi2, ndf)
   double d0, phi0, cpa, dz, tanl;  //reconstructed helical track parameters
   double d0err2, phi0err2, cpaerr2, dzerr2, tanlerr2; //err^2 of track parameters
   
   //Hits information
   std::vector<ALTckhit*> hits;  
   
   //External Triggers
   bool T1;
   bool T2;
   bool T3;
   bool T4;
   bool guard;
   std::vector<double> EneT1;  
   std::vector<double> EneT2;  
   std::vector<double> EneT3;
   std::vector<double> EneT4;
   std::vector<double> Eneg;   
   std::vector<double> timeT1;  
   std::vector<double> timeT2;  
   std::vector<double> timeT3;
   std::vector<double> timeT4;
   std::vector<double> timeg;     
   //Internal Triggers 
   //Tracker layer for 0 to 6. top layer is 0
   //integer value coded with 7 bits values
   //\Sum\Limits_{k=0}^{6} 2^k
   int Ti=0;
   
 public: 
   //Constructors
   ALEvent();// Default
   //Destructor
   ~ALEvent(){};   
   ///////////////////////////////
   // Methods
   ///////////////////////////////
   void Copy(ALEvent*);
   ////////////////////////////////
   //"setting" member methods
   ////////////////////////////////
   void set_eventnumber(int a){eventnumber=a;}
   void set_ncase(int a){ncase=a;}
   void set_typeMC(int a){typeMC=a;}
   void set_EkMC(double a){EkMC=a;}
   void set_X0MC(double a){X0MC=a;}
   void set_Y0MC(double a){Y0MC=a;}
   void set_Z0MC(double a){Z0MC=a;}
   void set_CX0MC(double a){CX0MC=a;}
   void set_CY0MC(double a){CY0MC=a;}
   void set_CZ0MC(double a){CZ0MC=a;}
   ////////////////////////////////
   void set_Nhits(int a){Nhits=a;}
   void add_Nhits(){Nhits++;}
   ////////////////////////////////  
   
   void set_EkPR(double a){EkPR=a;}
   void set_p0PR(double a){p0PR=a;}
   void set_a(double b){a=b;}
   void set_b(double a){b=a;}
   void set_c(double a){c=a;}
   void set_inter(double a){inter=a;}
   void set_slope(double a){slope=a;}
   void set_chi2B(double a){chi2B=a;}
   void set_chi2NB(double a){chi2NB=a;}
   void set_clB(double a){clB=a;}
   void set_clNB(double a){clNB=a;}
   
   
   ////////////////////////////////
   void set_typereco(int a){typereco=a;}
   void set_Ekreco(double a){Ekreco=a;}
   void set_p0reco(double a){p0reco=a;}
   void set_X0reco(double a){X0reco=a;}
   void set_Y0reco(double a){Y0reco=a;}
   void set_Z0reco(double a){Z0reco=a;}
   void set_CX0reco(double a){CX0reco=a;}
   void set_CY0reco(double a){CY0reco=a;}
   void set_CZ0reco(double a){CZ0reco=a;}
   void set_ndf(int a){ndf=a;}
   void set_chi2(double a){chi2=a;}
   void set_cl(double a){cl=a;}
   void set_d0(double a){d0=a;}
   void set_phi0(double a){phi0=a;}
   void set_cpa(double a){cpa=a;}
   void set_dz(double a){dz=a;}
   void set_tanl(double a){tanl=a;}
   void set_d0err2(double a){d0err2=a;}
   void set_phi0err2(double a){phi0err2=a;}
   void set_cpaerr2(double a){cpaerr2=a;}
   void set_dzerr2(double a){dzerr2=a;}
   void set_tanlerr2(double a){tanlerr2=a;}
   
   void add_hit(ALTckhit* h){hits.push_back(h);Nhits++;}
   
   //Set Pattern Reco variable at the position of the hit of index k
   void set_hxPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_xPR(a);}
   void set_hyPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_yPR(a);}
   void set_hzPR(int k, float a){if(k<(int)hits.size())(hits.at(k))->set_zPR(a);}

   //Set reconstruted variable at the position of the hit of index k
   void set_hxreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_xreco(a);}
   void set_hyreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_yreco(a);}
   void set_hzreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_zreco(a);}
   void set_hagereco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_agereco(a);}
   void set_hcxreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_cxreco(a);}
   void set_hcyreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_cyreco(a);}
   void set_hczreco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_czreco(a);}
   void set_hereco(int k,float a){if(k<(int)hits.size())(hits.at(k))->set_ereco(a);}
   
   
   void set_T1(bool a){T1=a;}
   void set_T2(bool a){T2=a;}
   void set_T3(bool a){T3=a;}
   void set_T4(bool a){T4=a;}
   void set_guard(bool a){guard=a;}
   void add_EneT1(double a){EneT1.push_back(a);}
   void add_EneT2(double a){EneT2.push_back(a);}
   void add_EneT3(double a){EneT3.push_back(a);}
   void add_EneT4(double a){EneT4.push_back(a);}
   void add_Eneg(double a){Eneg.push_back(a);}
   void add_timeT1(double a){timeT1.push_back(a);}
   void add_timeT3(double a){timeT3.push_back(a);}
   void add_timeT4(double a){timeT4.push_back(a);}
   void add_timeg(double a){timeg.push_back(a);}   
   void set_Ti(int a){Ti=a;}
   
   ////////////////////////////////
   //"Getting" member methods
   ////////////////////////////////
   int get_eventnumber(){return eventnumber;}
   int get_ncase(){return ncase;}
   int get_typeMC(){return typeMC;}
   double get_EkMC(){return EkMC;}
   double get_X0MC(){return X0MC;}
   double get_Y0MC(){return Y0MC;}
   double get_Z0MC(){return Z0MC;}
   double get_CX0MC(){return CX0MC;}
   double get_CY0MC(){return CY0MC;}
   double get_CZ0MC(){return CZ0MC;}


   ////////////////////////////////
   int get_Nhits(){return Nhits;}
   ////////////////////////////////
   
   double get_EkPR(){return EkPR;}
   double get_p0PR(){return p0PR;}
   double get_a(){return a;}
   double get_b(){return b;}
   double get_c(){return c;}
   double get_inter(){return inter;}
   double get_slope(){return slope;}
   double get_chi2B(){return chi2B;}
   double get_chi2NB(){return chi2NB;}
   double get_clB(){return clB;}
   double get_clNB(){return clNB;}
 
   ////////////////////////////////
   int get_typereco(){return typereco;}
   double get_Ekreco(){return Ekreco;}
   double get_p0reco(){return p0reco;}
   double get_X0reco(){return X0reco;}
   double get_Y0reco(){return Y0reco;}
   double get_Z0reco(){return Z0reco;}
   double get_CX0reco(){return CX0reco;}
   double get_CY0reco(){return CY0reco;}
   double get_CZ0reco(){return CZ0reco;}
   int get_ndf(){return ndf;}
   double get_chi2(){return chi2;}
   double get_cl(){return cl;}
   double get_d0(){return d0;}
   double get_phi0(){return phi0;}
   double get_cpa(){return cpa;}
   double get_dz(){return dz;}
   double get_tanl(){return tanl;}
   double get_d0err2(){return d0err2;}
   double get_phi0err2(){return phi0err2;}
   double get_cpaerr2(){return cpaerr2;}
   double get_dzerr2(){return dzerr2;}
   double get_tanlerr2(){return tanlerr2;}
   std::vector<ALTckhit*> get_hits(){return hits;}
   bool get_T1(){return T1;}
   bool get_T2(){return T2;}
   bool get_T3(){return T3;}
   bool get_T4(){return T4;}
   bool get_guard(){return guard;}
   std::vector<double>  get_EneT1(){return EneT1;}
   std::vector<double>  get_EneT2(){return EneT2;}
   std::vector<double>  get_EneT3(){return EneT3;}
   std::vector<double>  get_EneT4(){return EneT4;}
   std::vector<double>  get_Eneg(){return Eneg;}
   std::vector<double>  get_timeT1(){return timeT1;}
   std::vector<double>  get_timeT2(){return timeT2;}
   std::vector<double>  get_timeT3(){return timeT3;}
   std::vector<double>  get_timeT4(){return timeT4;}
   std::vector<double>  get_timeg(){return timeg;}
   int get_Ti(){return Ti;}

   ////////////////////////////////
   ClassDef(ALEvent,1)

};

#endif
