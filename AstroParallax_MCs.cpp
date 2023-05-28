#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
using namespace std;

const int Num=10000;
const double MaxD=65.0;///kpc
const double step=MaxD/(double)Num/1.0;///step in kpc



const double pi= M_PI;
const double RA=180.0/M_PI;
const double Hp=6.62607004*pow(10.0,-34);
const double KP=3.08568025*pow(10.,19); // in meter.
const double G=6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double M_sun=1.98892*pow(10.,30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=364.0;
const double binary_fraction=double(2.0/3.0);
const double Avks=double(8.20922);
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
const double Dsun=8.0; 

///============================ Besancon constant ==============///
const double R_sun=8.5;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]=  {3.1,2.5,3.1,3.1};

///============================ WFIRST & OGLE & KMTNet ===================
const int M=5;///number of filters  VIKH,W149
const double satu[M]= {12.0, 12.0, 13.0, 13.0, 14.8}; //it should be changed
const double thre[M]= {20.0, 21.0, 21.0, 21.0, 26.0};
const double FWHM[M]= {0.33, 0.33, 0.33, 0.33 , 0.33};//3*pixel_size (0.11") of WFIRST, VIKH W149 Z087
const double AlAv[M]= {1.009,0.600,0.118,0.184,0.225};///From besancon model[VIKH W149]
const double sigma[M]={0.022,0.022,0.02,0.025,0.025};//MOAاستفاده از مقاله کاردلی
const double Akv=0.118;
const double cade1=double(15.16/60.0/24.0);///W149_cadence in  day
const double dt=cade1;//double(7.58/60.0/24.0);//[days]
const int coun=int(5.0*year/dt) + 2;

////========================== LMC =============================
const double bLMC= -32.9; 
const double lLMC= 280.5; 
const double RaLMC= 80.89375;
const double DecLMC=-68.2438888888889;  
const double sLMC= 5.0;///degree (half of size of LMC)
const double DLMC= 49.97;///KPC 

const int nn=int(455);
const double dec1= -72.5; 
const double dec2= -64.8; 
const double rig1=float((4.0+27.0/60.0)*15.0 );
const double rig2=float((5.0+50.0/60.0)*15.0 );
const double ddec=5.0*fabs(dec2-dec1)/nn/1.0;
const double drig=5.0*fabs(rig2-rig1)/nn/1.0;

const int nex=23104; 
const int nrd=10000; 

///**************  Number constant  ********************///
const int Na=int(38);
const int Nw=int(123);///sigma_WFIRST.dat
const int YZ=3615;///number of rows in file yzma
const int N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
/// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///======================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI, lat, lon;
    double right, decli;
    double od_disk,od_ThD,od_bulge,od_halo,opd;///optical  depth
    double od_dlmc, od_blmc, od_hlmc; 
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double rho_hlmc[Num], rho_dlmc[Num],rho_blmc[Num]; 
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs, Romins,nstart, nstarti;
    double Nblend[M], blend[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M];
    double type, Tstar, logl, col, Rstar, mass,vs;
    double Romax,ro_star, deltao;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double xv, yv, zv; 
    double tetv; 
};
struct lens{
    int numl,struc;
    double Ml, Dl, vl , Vt, xls;
    double rhomaxl,tE,RE;
    double mul, u0;
    double t0;
    double tetE;
    double piE; 
};
struct astromet{
   double tetp;
   double omegae;
   double vearth;
   double fact;
   double Ve_n1, Ve_n2;
   double sx0, sy0, lx0, ly0;
   double sx1, sy1, lx1, ly1;
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[5][N1],Rs_d[N1],mass_d[N1],type_d[N1]; int cl_d[N1];  ///thin disk
    double Teff_b[N2],logl_b[N2],Mab_b[5][N2],Rs_b[N2],mass_b[N2],type_b[N2]; int cl_b[N2];  /// bulge
    double Teff_t[N3],logl_t[N3],Mab_t[5][N3],Rs_t[N3],mass_t[N3],type_t[N3]; int cl_t[N3];  ///thick disk
    double Teff_h[N4],logl_h[N4],Mab_h[5][N4],Rs_h[N4],mass_h[N4],type_h[N4]; int cl_h[N4];  /// halo
};
/*struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
};*/
struct roman{
     double magw[Nw], errw[Nw];
     double maga[Na], erra[Na]; 
};
struct galactic{
   double l[nrd], b[nrd], RA[nrd], Dec[nrd]; 
};
///===================== FUNCTION ==============================
void read_cmd(CMD & cm);
void func_source(source & s, CMD & cm);
void func_lens(lens & l, source & s);
void vrel(source & s , lens & l);
void Disk_model(source & s, int );
void optical_depth(source & s);
double RandN(double sigma, double Nn);
double RandR(double down, double up);
void astrometric(source & s, astromet & as, double tim);
double errwfirst(roman & ro, double ghadr, int flag);
//////////////////////////////////////////////////////
    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;disk_
    FILE * _randStream;
///////////////////////////////////////////////////////

///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//

int main()
{
///****************************************************************************
//gettimeofday(&newTime,NULL);
//ftime(&tp);
	time(&_timeNow);
	_randStream = fopen("/dev/urandom", "r");
	_dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
///****************************************************************************
     source s;
     lens l;
     astromet as;
     CMD cm;
     galactic ga;
     roman ro; 
     read_cmd(cm);
     
     FILE* magg;  
     FILE* data3;    

///===========================================================

    FILE * wfirstf;
    wfirstf=fopen("./files/sigma_WFIRST.txt","r");
    for(int i=0; i<Nw;++i){
    fscanf(wfirstf,"%lf     %lf\n",&ro.magw[i],&ro.errw[i]);}
    fclose(wfirstf);
    
    FILE * roman;
    roman=fopen("./files/roman_astro2.txt","r");
    for(int i=0; i<Na;++i){
    fscanf(roman,"%lf     %lf\n",&ro.maga[i],&ro.erra[i]);
    ro.erra[i]= pow(10.0,ro.erra[i])*1000.0;///miliarcs
    }
    fclose(roman);
    //cout<<"*****  sigma_wfirst was read *****  "<<endl;
///===========================================================

   FILE* convert; 
   convert=fopen("./files/convert_coordinate_2.dat","r");
   int coo=0; 
   for(int i=0; i<100; ++i){
   for(int j=0; j<100; ++j){
   fscanf(convert,"%lf  %lf   %lf  %lf \n",&ga.RA[coo], &ga.Dec[coo], &ga.l[coo], &ga.b[coo]);  coo++; }}   
   fclose(convert); 
   


    FILE* fil1;
    fil1=fopen("./files/MONT/IBH_ROMANMC.txt","w");
    fclose(fil1);



    double obs1[6][2]={0.0};
    double obs2[4][2]={0.0};
    obs1[0][0]=0.0*year+0.0;    obs1[0][1]=0.0*year+62.0;
    obs1[1][0]=0.0*year+182.0;  obs1[1][1]=0.0*year+182.0+62.0;
    obs1[2][0]=1.0*year+0.0;    obs1[2][1]=1.0*year+62.0;

    obs2[0][0]=1.0*year+182.0;  obs2[0][1]=1.0*year+182.0+62.0;
    obs2[1][0]=2.0*year+0.0;    obs2[1][1]=2.0*year+62.0;
    obs2[2][0]=2.0*year+182.0;  obs2[2][1]=2.0*year+182.0+62.0;
    obs2[3][0]=3.0*year+0.0;    obs2[3][1]=3.0*year+62.0;

    obs1[3][0]=3.0*year+182.0;  obs1[3][1]=3.0*year+182.0+62.0;
    obs1[4][0]=4.0*year+0.0;    obs1[4][1]=4.0*year+62.0;
    obs1[5][0]=4.0*year+182.0;  obs1[5][1]=4.0*year+182.0+62.0;
    
   
    double cade2[2000];
    for(int i=0; i<2000; i+=4){
    cade2[i]=  cade1; 
    cade2[i+1]=cade1; 
    cade2[i+2]=cade1; 
    cade2[i+3]= 10.0;}


    as.tetp=double(M_PI/3.0);
    as.omegae=double(2.0*M_PI/year); ///radian per day
    as.vearth= as.omegae;//*0.001/(24.0*3600.0); ///km/s
    as.fact=double(1000.0*3600.0*24.0/AU);
  
    
    int flagf, ii;
    int ndw, w1;
    char   filnam1[40], filnam2[40];
    double magnio, erro, test, dist;
    double magni[M];
    double timp1,lonn, chi, chi1, tim;
    double flag0, flag1, flag2;
    double t_0, t_1;
    double dchi, mh;
    double Nmlw=0.0;
    int    flag_peak=0, flag_det=0, flag_mlw=0, hh;
    double lonn0, lat0, siga, snr, errx, erry;
    double u0, u1, Vt, Astar0, Astar1, Dx0, Dy0, Dx1, Dy1, Dx2, Dy2; 
    int    counter=0;
    double sxi, syi, mag_w0,flago, mag_w1, mindd; 
    double mus1, mus2, mul1, mul2;
    double errm, erra; 

    s.right= RaLMC;  
    s.decli= DecLMC;   
   
    hh=-1; mindd=10000.0; 
    for(int i=0; i<10000; ++i){
    dist= double(s.right- ga.RA[i])*(s.right- ga.RA[i]) + (s.decli- ga.Dec[i])*(s.decli- ga.Dec[i]); 
    dist=sqrt(dist); 
    if(dist<=mindd){ mindd=dist; hh= i;   }}
    if(hh<0){cout<<"ERROR:right:    "<<s.right<<"\t s.dec:  "<<s.decli<<endl;   int uue;  cin>>uue; }
    s.lon=double(ga.l[hh]);   
    s.lat=double(ga.b[hh]);
    s.TET=(360.0-s.lon)/RA;///radian
    s.FI=  s.lat/RA; 
    Disk_model(s, 0); 
   


     for(int icon=0; icon<10000; ++icon) {
     cout<<"***** Step:  "<<icon<<endl;
     //for (int li=1; li<=7; li++){
    // if(li==1) {lonn0=1.3;  lat0=-0.875; }
    // if(li==2) {lonn0=0.9;  lat0=-0.875; }
    // if(li==3) {lonn0=1.3;  lat0=-1.625; }
    // if(li==4) {lonn0=0.9;  lat0=-1.675; }
    // if(li==5) {lonn0=0.5;  lat0=-1.675; }
    // if(li==6) {lonn0=0.1;  lat0=-1.675; }
     //if(li==7) {lonn0=-0.3; lat0=-1.675; }

     //s.lat= lat0+ double(RandR(-0.375 , 0.375));
     //lonn =lonn0+ double(RandR(-0.200 , 0.200));
     //cout<<"lat:  "<<s.lat<<"\t lonn:  "<<lonn<<endl;
     //cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;

    // if(lonn<=0.0)   s.lon=360.0+lonn;
    // else            s.lon=lonn;
    // s.TET=double(360.0-s.lon)/RA;///radian
    // s.FI=s.lat/RA;///radian
    // Disk_model(s, icon);
     
    // if(int(Extinction(ex,s))==1){
     do{
     func_source(s, cm);
     func_lens(l, s);
     }while(l.tE<=0.1 or l.tE>3000.0);
     optical_depth(s);
    
    
     mus1= double(s.SV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
     mus2= double(s.SV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
     mul1= double(s.LV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
     mul2= double(s.LV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days 
     l.piE=double(1.0/l.Dl -1.0/s.Ds)/l.tetE; 
    
     test=(double)rand()/((double)(RAND_MAX)+(double)(1.0));
     if(test<=s.blend[4] and s.magb[4]<=thre[4]){
     counter+=1; 

     flagf=0;
     if(counter%1==0){
     flagf=1;
     sprintf(filnam1,"./files/MONT/%c%c%c%c%d.dat",'m','a','g','_', counter);
     sprintf(filnam2,"./files/MONT/%c%c%c%c%d.dat",'d','a','t','_', counter);
     magg=fopen(filnam1,"w");    data3=fopen(filnam2,"w");}
     
     
     double *magn=new double[coun];
     double *timn=new double[coun];
     double *soux=new double[coun];
     double *souy=new double[coun];
///*********************************************************************************
  
     int numn=0;   
     for(int i=0; i<coun; ++i){magn[i]=timn[i]=-1.0; soux[i]=souy[i]=0.0;}
     sxi= -l.u0*l.tetE*sin(s.tetv);
     syi= +l.u0*l.tetE*cos(s.tetv);
     as.lx1=as.ly1=as.lx0=as.ly0=0.0;//Lens is at the origin 
     as.sx1=as.sx0=sxi;//n1
     as.sy1=as.sy0=syi;//n2


     for(tim=l.t0;  tim>0.0;  tim-=dt){
     numn+=1;
     astrometric(s,as, tim);
     as.sx1 -=  double(s.SV_n1- s.VSun_n1)*as.fact*dt/s.Ds - as.Ve_n1*dt/s.Ds;///mili arcs
     as.sy1 -=  double(s.SV_n2- s.VSun_n2)*as.fact*dt/s.Ds - as.Ve_n2*dt/s.Ds;///mili arcs
     as.lx1 -=  double(s.LV_n1- s.VSun_n1)*as.fact*dt/l.Dl - as.Ve_n1*dt/l.Dl;///mili arcs
     as.ly1 -=  double(s.LV_n2- s.VSun_n2)*as.fact*dt/l.Dl - as.Ve_n2*dt/l.Dl;///mili arcs
     
     as.sx0 -=  double(s.SV_n1- s.VSun_n1)*as.fact*dt/s.Ds;///mili arcs
     as.sy0 -=  double(s.SV_n2- s.VSun_n2)*as.fact*dt/s.Ds;///mili arcs
     as.lx0 -=  double(s.LV_n1- s.VSun_n1)*as.fact*dt/l.Dl;///mili arcs
     as.ly0 -=  double(s.LV_n2- s.VSun_n2)*as.fact*dt/l.Dl;///mili arcs
     
     u0= sqrt((as.sx0-as.lx0)*(as.sx0-as.lx0)+(as.sy0-as.ly0)*(as.sy0-as.ly0))/l.tetE;//without dimention 
     u1= sqrt((as.sx1-as.lx1)*(as.sx1-as.lx1)+(as.sy1-as.ly1)*(as.sy1-as.ly1))/l.tetE;//without dimention
     
     Astar0= (u0*u0+2.0)/sqrt(u0*u0*(u0*u0+4.0));
     Astar1= (u1*u1+2.0)/sqrt(u1*u1*(u1*u1+4.0));
     
     mag_w0= s.magb[4]- 2.5*log10(Astar0*s.blend[4] + 1.0 - s.blend[4]);//W149
     mag_w1= s.magb[4]- 2.5*log10(Astar1*s.blend[4] + 1.0 - s.blend[4]);//W149
     
     Dx0=double(as.sx0-as.lx0)/(u0*u0+2.0); 
     Dy0=double(as.sy0-as.ly0)/(u0*u0+2.0);
     Dx1=double(as.sx1-as.lx1)/(u1*u1+2.0); 
     Dy1=double(as.sy1-as.ly1)/(u1*u1+2.0);
     
     Dx2=double(as.sx0-as.lx0)/(u1*u1+2.0); 
     Dy2=double(as.sy0-as.ly0)/(u1*u1+2.0);
     
     
     ii=int(tim/dt)+1;
     if(ii>=0 and ii<coun){
     magn[ii]= Astar1;  timn[ii]=tim; soux[ii]=Dx1;  souy[ii]=Dy1;}//  soux[ii]= as.sx1+Dx0;  souy[ii]=as.sy1+Dy0;}
     if(flagf>0 and numn%5==0)  
     
     fprintf(magg,"%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf  %.4lf %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf  %.4lf %.4lf %d\n",
     tim,mag_w1,mag_w0,Astar1,Astar0,as.sx1,as.sy1,as.lx1,as.ly1,as.sx0,as.sy0,as.lx0,as.ly0, Dx0, Dy0, Dx1, Dy1,Dx2,Dy2,sxi,syi,0);}//22
     
  

///*********************************************************************************

     as.lx1=as.ly1=as.lx0=as.ly0=0.0;///Lens is at the origin 
     as.sx1=as.sx0= sxi;///n1
     as.sy1=as.sy0= syi;///n2
     for(tim=l.t0; tim<(5.0*year); tim+=dt){
     numn+=1;
    
     astrometric(s,as, tim);
      
     as.sx1 +=  double(s.SV_n1- s.VSun_n1)*as.fact*dt/s.Ds - as.Ve_n1*dt/s.Ds;///mili arcs
     as.sy1 +=  double(s.SV_n2- s.VSun_n2)*as.fact*dt/s.Ds - as.Ve_n2*dt/s.Ds;///mili arcs
     as.lx1 +=  double(s.LV_n1- s.VSun_n1)*as.fact*dt/l.Dl - as.Ve_n1*dt/l.Dl;///mili arcs
     as.ly1 +=  double(s.LV_n2- s.VSun_n2)*as.fact*dt/l.Dl - as.Ve_n2*dt/l.Dl;///mili arcs
     
     as.sx0 +=  double(s.SV_n1- s.VSun_n1)*as.fact*dt/s.Ds;///mili arcs
     as.sy0 +=  double(s.SV_n2- s.VSun_n2)*as.fact*dt/s.Ds;///mili arcs
     as.lx0 +=  double(s.LV_n1- s.VSun_n1)*as.fact*dt/l.Dl;///mili arcs
     as.ly0 +=  double(s.LV_n2- s.VSun_n2)*as.fact*dt/l.Dl;///mili arcs
        
     u0= sqrt((as.sx0-as.lx0)*(as.sx0-as.lx0)+(as.sy0-as.ly0)*(as.sy0-as.ly0))/l.tetE;//without dimention 
     u1= sqrt((as.sx1-as.lx1)*(as.sx1-as.lx1)+(as.sy1-as.ly1)*(as.sy1-as.ly1))/l.tetE;//without dimention
     

     Astar0= (u0*u0+2.0)/sqrt(u0*u0*(u0*u0+4.0));
     Astar1= (u1*u1+2.0)/sqrt(u1*u1*(u1*u1+4.0));
     mag_w0= s.magb[4]- 2.5*log10(Astar0*s.blend[4] + 1.0 - s.blend[4]);//W149
     mag_w1= s.magb[4]- 2.5*log10(Astar1*s.blend[4] + 1.0 - s.blend[4]);//W149
     
     
     Dx0= double(as.sx0-as.lx0)/(u0*u0+2.0); 
     Dy0= double(as.sy0-as.ly0)/(u0*u0+2.0);
     Dx1= double(as.sx1-as.lx1)/(u1*u1+2.0); 
     Dy1= double(as.sy1-as.ly1)/(u1*u1+2.0);
     
     Dx2=double(as.sx0-as.lx0)/(u1*u1+2.0); 
     Dy2=double(as.sy0-as.ly0)/(u1*u1+2.0);
     
    
     ii=int(tim/dt)+1;
     if(ii>=0 and ii<coun){
     magn[ii]= Astar1;  timn[ii]=tim; soux[ii]=Dx1;  souy[ii]=Dy1;}// soux[ii]= as.sx1+Dx0;  souy[ii]=as.sy1+Dy0;}
     if(flagf>0 and numn%5==0)  
     fprintf(magg,"%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf  %.4lf %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf %.4lf %.4lf %.4lf %.4lf  %.4lf %.4lf  %.4lf  %d\n",
     tim,mag_w1,mag_w0,Astar1,Astar0,as.sx1,as.sy1,as.lx1,as.ly1,as.sx0,as.sy0,as.lx0,as.ly0, Dx0, Dy0,Dx1, Dy1,Dx2,Dy2,sxi,syi,1);}//22
    //0     1      2     3       4      5     6      7       8     9     10       11     12    13   14   15  16  17
  
     
///************************** WFIRST  *********************************************
     flag0=flag1=flag2=0.0;
     ndw=flag_peak=flag_det=0;
     t_1=t_0=0.0;
     chi=chi1=0.0;
     timp1=0.0;
     int vx=0, datf; 
     
     
     
     for(int j=0; j<coun; ++j){
     if(magn[j]>0.0){
     tim=timn[j];  

    
     
     for(int i=0;i<M; ++i){
     magni[i] =s.magb[i]-2.5*log10(magn[j]*s.blend[i] + 1.0 - s.blend[i]);}
     
     
     
     flago=0.0;
     for(int kl=0; kl<6; ++kl){
     if((obs1[kl][0]-tim)*(obs1[kl][1]-tim)<=0.0){flago=-1.0; break;}}
     if(flago>-0.5){
     for(int kl=0; kl<4; ++kl){
     if((obs2[kl][0]-tim)*(obs2[kl][1]-tim)<=0.0){flago=1.0; break;}} }
     //test=(double)rand()*100.0/((double)(RAND_MAX)+(double)(1.0));
      
     
     
     
     datf=0; 
     if(flago<-0.5){
     timp1 +=dt;
     if(timp1>=cade1){    datf=1;  timp1-=cade1;} }
     
     else if(flago>+0.5){
     timp1 +=dt;
     if(timp1>=cade2[vx]){datf=0;  timp1-=cade2[vx]; vx+=1;      if(vx>1998) vx=0; } }
     
     
     
     if(magni[4]>=satu[4] and magni[4]<=thre[4] and datf>0){ 
     errm= errwfirst(ro, magni[4], 0); 
     erra= errwfirst(ro, magni[4], 1); 
     magnio=magni[4] + RandN(errm ,3.0);
     chi  +=fabs( (magnio- s.magb[4])*(magnio- s.magb[4])/(errm*errm));
     chi1 +=fabs( (magnio-  magni[4])*(magnio-  magni[4])/(errm*errm));
     errx= RandN(erra ,3.0);
     erry= RandN(erra ,3.0);
     fprintf(data3,"%.8lf  %.8lf  %.8lf    %.8lf    %.8lf    %.8lf   %.1lf\n",tim, magnio, errm, soux[j]+errx, souy[j]+erry, erra, flago);
     flag2=0.0;
     if(fabs(magnio-s.magb[4])>fabs(4.0*errm) )    flag2=1.0;
     if(ndw>1 and float(flag0+flag1+flag2)>2.0)    flag_det=1; 
     flag0=flag1;
     flag1=flag2;
     ndw+=1;}

     }}

    if(flagf>0){fclose(magg); fclose(data3);}
    cout<<"End of loop time ******************"<<endl;
    ///************ WFIRST *********************************************************************
    dchi=fabs(chi-chi1);
    if(dchi>800.0 and flag_det>0 and ndw>5){
   
    
    Nmlw+=1.0;
    fil1=fopen("./files/MONT/IBH_ROMANMC.txt","a+");
    fprintf(fil1,
   "%d  %.5lf   %.5lf  "///3
   "%d  %.5lf   %.5lf   %.5lf  "///7
   "%d   %d  %.5lf   %.5lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf   %.4lf   %.4lf  %.4lf  "   //21
   "%.5lf  %.9lf  %.5lf  %.9lf   %.2lf  %.2lf   %.4lf  %.4lf   "  //29
   "%.8lf  %.6lf  %.6lf  %.8lf    %.7lf  %.8lf  %.4lf  %.8lf  %.9lf " ///38
   "%d  %d  %.1lf  %d   %.9lf  %.9lf   %.9lf  %.9lf  %.9lf    %.9lf\n", //49
   counter,s.lat, s.lon, //3
   l.struc, l.Ml, l.Dl, l.vl, //7
   s.struc, s.cl, s.mass, s.Ds, s.Tstar, s.Rstar, s.logl, s.type, s.col, s.vs, s.Mab[1], s.Mab[4], s.Map[1], s.Map[4],//21
   s.magb[1], s.magb[4], s.blend[1], s.blend[4], s.Nblend[1], s.Nblend[4], s.Ai[1], s.Ai[4], //29
   l.tE, l.RE/AU, l.t0, l.mul, l.Vt, l.u0, s.opd*1.0e6, s.ro_star, l.tetE,//38
   flagf, flag_det, dchi, ndw, mus1, mus2, s.tetv,mul1, mul2, l.piE);//49
   fclose(fil1);
   cout<<"************* End of saving in the file ***********************"<<endl;}
   if(int(counter)%1==0){
   cout<<"============================================================="<<endl;
   cout<<"counter:     "<<counter<<endl;
   cout<<"lat:  "<<s.lat<<"\t lon:  "<<s.lon<<endl;
   cout<<"********************** SOURCE **************************"<<endl;
   cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<"\t strucs:  "<<s.struc<<endl;
   cout<<"cl:  "<<s.cl<<"\t mass:  "<<s.mass<<"\t Tstar:  "<<s.Tstar<<endl;
   cout<<"Rstar:  "<<s.Rstar<<"\t logl:  "<<s.logl<<"\t type:  "<<s.type<<endl;
   cout<<"col:  "<<s.col<<"\t vs:  "<<s.vs<<"\t Mag_I:  "<<s.Mab[1]<<"\t Mag_W149:  "<<s.Mab[4]<<endl;
   cout<<"********************** LENS **************************"<<endl;
   cout<<"Dl:  "<<l.Dl<<"\t numl:  "<<l.numl<<"\t strucl:  "<<l.struc<<endl;
   cout<<"mass:  "<<l.Ml<<"vl:  "<<l.vl<<endl;
   cout<<"*********************** LENSING ***********************"<<endl;
   cout<<"tE:  "<<l.tE<<"\t RE(AU):  "<<l.RE/AU<<"\t t0:  "<<l.t0<<endl;
   cout<<"Vt:  "<<l.Vt<<"\t mul(mas/days):  "<<l.mul<<"\t u0:  "<<l.u0<<endl;
   cout<<"************ WFIRST & OGLE ***************************"<<endl;
   cout<<"flag_det:  "<<flag_det<<"\t dchi:  "<<dchi<<"\t ndw:  "<<ndw<<endl;
   cout<<"Nmlw:    "<<Nmlw<<"\t flag_peak:  "<<flag_peak<<endl;
   cout<<"==============================================================="<<endl;}
   
   delete [] magn, timn, soux, souy;
   
   }//if wfirst>0
   //}//end of EXtinction
   //}//loop icon
   } ///loop il
   return(0);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double errwfirst(roman & ro, double ghadr, int flag){

     double error=0.0, shib=0.0; 
     if(flag==0){///calculating photometry errors ROMAN
     if(ghadr<ro.magw[0] or  ghadr==ro.magw[0])          error= ro.errw[0];
     
     else if(ghadr>ro.magw[Nw-1] or ghadr==ro.magw[Nw-1]) {
     shib=(ro.errw[Nw-1]-ro.errw[Nw-2])/(ro.magw[Nw-1]-ro.magw[Nw-2]); 
     error=ro.errw[Nw-1]+shib*(ghadr-ro.magw[Nw-1]);    }
     
     else{
     for(int i=1; i<Nw; ++i){
     if(double((ghadr-ro.magw[i])*(ghadr-ro.magw[i-1]))<0.0 or  ghadr==ro.magw[i-1]){
     shib=(ro.errw[i]-ro.errw[i-1])/(ro.magw[i]-ro.magw[i-1]); 
     error=ro.errw[i-1]+shib*(ghadr-ro.magw[i-1]); 
     break;}}}
     }
     
     
     if(flag==1){///calculating astrometry error ROMAN
     if(ghadr<ro.maga[0] or  ghadr==ro.maga[0])            error= ro.erra[0]; 
     else if(ghadr>ro.maga[Na-1] or ghadr==ro.maga[Na-1]){
     shib=(ro.erra[Na-1]-ro.erra[Na-2])/(ro.maga[Na-1]-ro.maga[Na-2]); 
     error=ro.erra[Na-1]+shib*(ghadr-ro.maga[Na-1]);   }   
       

     else{
     for(int i=1; i<Na; ++i){
     if(double((ghadr-ro.maga[i])*(ghadr-ro.maga[i-1]))<0.0 or ghadr==ro.maga[i-1]){
     shib=(ro.erra[i]-ro.erra[i-1])/(ro.maga[i]-ro.maga[i-1]);
     error=ro.erra[i-1] +  shib*(ghadr-ro.maga[i-1]); 
     break;}}}
     }
 
     return(error);     
}     
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         AStrometric calculations                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void astrometric(source & s, astromet & as , double tim)
{
   double vex, vey, VSE_R, VSE_T, VSE_Z, Ve_x;  
   vex= as.vearth * cos(as.omegae* tim + M_PI/2.0);
   vey= as.vearth * sin(as.omegae* tim + M_PI/2.0);
   VSE_R = double( cos(as.tetp)*vex );
   VSE_T = double( vey );
   VSE_Z = double( sin(as.tetp)*vex );
 
   as.Ve_n1= VSE_R*sin(s.deltao)  - VSE_T*cos(s.deltao);
   Ve_x=    -VSE_R*cos(s.deltao)  - VSE_T*sin(s.deltao);
   as.Ve_n2=-sin(s.FI)*Ve_x + cos(s.FI)*VSE_Z;
    
   
    if(fabs(as.Ve_n1)>fabs(2.0*as.vearth) or fabs(as.Ve_n2)>fabs(2.0*as.vearth)){
    cout<<"ERROR Vsun_as_l:  "<<as.Ve_n1<<"\t VSun_s:  "<<s.VSun_n1<<"\t Earth_V:  "<<as.vearth<<endl;
    cout<<"Vsun_as_b:  "<<as.Ve_n2<<"\t VSun_s:  "<<s.VSun_n2<<endl;
    cout<<"TET:  "<<s.TET<<"\t  FI:  "<<s.FI<<endl;
    cout<<"tetp:  "<<as.tetp<<endl;
    cout<<"s.deltao:  "<<s.deltao<<endl;
    int uue ; cin>>uue;}
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*M_sun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD  &  cm)
{
    int num,struc,nums,yye;
    double rho,rf;
    double Ds,Ai[M],Av;
    double Mab[M],Map[M];
    double maxnb=0.0;


    for(int i=0; i<M; ++i){
    s.Fluxb[i]=s.Nblend[i]=0.0;
    s.Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    s.Nblend[i]=s.Nblend[i]+RandN(sqrt(s.Nblend[i]),1.5);
    if(s.Nblend[i]<=1.0) s.Nblend[i]=1.0;
    if(s.Nblend[i]>maxnb) maxnb=s.Nblend[i];}
    for(int i=0; i<M; ++i){s.magb[i]=s.Ai[i]=s.Map[i]=s.Mab[i]=0.0;}
    //cout<<"Nblend[4]:  "<<s.Nblend[4]<<"\t maxnb: "<<maxnb<<endl;


    for(int k=1; k<=int(maxnb); ++k){

    do{
    num=int((Num-15.00)*(double)rand()/(double)(RAND_MAX+1.) +10.00);
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    Ds=(double)(num*step);
    }while(rho>s.Rostari[num] or Ds<0.1 or Ds>MaxD);///distance larger than 50.
    if(Ds>MaxD or Ds<0.1){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<num<<endl; cin>>yye;}
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}
   // cout<<"Ds:  "<<Ds<<endl;


     
    rf= RandR(0.0, s.Rostar0[nums]); 
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]))struc=4;///halo_lmc
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums]+s.rho_hlmc[nums]+s.rho_dlmc[nums]))struc=5;//disk_lmc
else struc=6;///bulge_lmc


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(struc==0 or struc==5){//thin disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0));
    for(int i=0; i<5; ++i){Mab[i]=cm.Mab_d[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_d[num];
    s.type= cm.type_d[num];
    s.mass= cm.mass_d[num];
    s.Tstar=cm.Teff_d[num];
    s.logl= cm.logl_d[num];
    s.col=Mab[0]-Mab[1];
    s.cl= cm.cl_d[num];}
    if(cm.mass_d[num]<0.0 or int(cm.cl_d[num])>5 or cm.Teff_d[num]<0.0 or cm.type_d[num]>8.0 or cm.mass_d[num]>1000.0){
    cout<<"Error(thin disk) temp: "<<cm.Teff_d[num]<<"\t mass: "<<cm.mass_d[num]<<"\t counter: "<<num<<endl; cin>>yye; }}


    if(struc==1 or struc==6){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_b[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_b[num];
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Tstar=cm.Teff_b[num];
    s.col= Mab[0]-Mab[1];
    s.logl=cm.logl_b[num];
    s.cl= cm.cl_b[num];}
    if(cm.mass_b[num]<0.0 or int(cm.cl_b[num])>5 or cm.Teff_b[num]<0.0 or cm.type_b[num]>8.0 or cm.mass_b[num]>10000.0){
    cout<<"Error(bulge) temp: "<<cm.Teff_b[num]<<"\t mass: "<<cm.mass_b[num]<<"\t counter: "<<num<<endl;   cin>>yye; }}




    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_t[i][num]; }
    if(k==1){
    s.Rstar= cm.Rs_t[num];
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Tstar=cm.Teff_t[num];
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_t[num];
    s.cl= cm.cl_t[num];}
    if( cm.mass_t[num]<0.0 or int(cm.cl_t[num])>5 or cm.Teff_t[num]<0.0 or cm.type_t[num]>8.0 or cm.mass_t[num]>100000.0){
    cout<<"Error(thick disk) temp: "<<cm.Teff_t[num]<<"\t mass: "<<cm.mass_t[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}



    if(struc==3 or struc==4){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_h[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_h[num];
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Tstar=cm.Teff_h[num];
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_h[num];
    s.cl= cm.cl_h[num];}
    if(cm.mass_h[num]<0.0 or int(cm.cl_h[num])>5 or cm.Teff_h[num]<0.0 or cm.type_h[num]>8.0 or cm.mass_h[num]>10000.0){
    cout<<"Error(Galactic halo) temp: "<<cm.Teff_h[num]<<"\t mass: "<<cm.mass_h[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}



    Mab[4]=(Mab[2]+Mab[3]+Mab[4])/3.0;//absolute magnitude in W149: (K+H+J)/3
    //ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    //Av=ex.Aks*Avks;
 
Av=0.1;    
    if(Av<0.0)    Av=0.0;

    if(Av>20.0 or Av<-0.00365 or Ds>MaxD or Ds<0.0){
    cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl;
    cout<<"Nblend[4]:  "<<s.Nblend[4]<<"\t Nblend[1]:  "<<s.Nblend[1]<<endl;
    cout<<"\t latitude:   "<<s.lat<<"\t longtide:  "<<s.lon<<endl; }


    for(int i=0;  i<M; ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i] + 5.0*log10(Ds*100.0) + Ai[i];
    if(s.Nblend[i]>=k){ s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));} }
    if(k==1){
    for(int i=0; i<M;  ++i){s.Ai[i]=Ai[i]; s.Map[i]=Map[i]; s.Mab[i]= Mab[i];}
    s.col=s.col+s.Ai[0]-s.Ai[1];}
    }///loop over the stars


    for(int i=0; i<M; ++i){
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.magb[i]= -2.5*log10(fabs(s.Fluxb[i]));    
    if(int(s.Nblend[i])<1.0 or (s.Nblend[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or
    s.Fluxb[i]<-0.009543 or s.Ai[i]<0.0){
    cout<<"BIGG ERRROR nsbl: "<<s.Nblend[i]<<"\t Nblend[i]: "<<s.Nblend[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<"\t ext:  "<<s.Ai[i]<<endl; cin>>yye;} }
    if(s.type>8.0 or s.type<1.0 or s.mass<=0.0 or s.Tstar<0.0 or s.Rstar<0.0 or s.mass>10000.0 or
    s.nums>Num or s.nums<=0 or Av<0.0 or s.cl<0 or s.cl>=6){
    cout<<"ERROR(source):  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<"\t cl:  "<<s.cl<<endl;   cin>>yye;}
   // cout<<"************** End of func_source  ****************"<<endl;
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s)
{
    double f,test;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    int il, yye;
    double mmin=2.0;
    double mmax=600.0;

    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}


    do{
    l.numl=(int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test  =((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){
    cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);



double randflag=RandR(0.0,s.Rostar0[l.numl]);
if (randflag<=fabs(s.rho_disk[l.numl]) ) l.struc=0;//disk
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1;//bulge
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2;//thick
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]) ) l.struc=3;//halo
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]+s.rho_hlmc[l.numl]))l.struc=4;//halo_lmc
else if (randflag<=fabs(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl]+s.rho_halo[l.numl]+s.rho_hlmc[l.numl]+s.rho_dlmc[l.numl])) l.struc=5;//disk_lmc
else if (randflag<=fabs( s.Rostar0[l.numl])) l.struc=6;//bar_lmc
else {  cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}

    //do{
    //l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(mmax - mmin)) + mmin;
    //test=fabs((double)rand()/(double)(RAND_MAX+1.)*57.0);
    //if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
    //if(l.Ml> 1.0) f=pow(l.Ml,-3.0);
    //}while(test>f);
    l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(mmax - mmin)) + mmin;
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*M_sun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    s.ro_star=s.Rstar*Rsun*l.xls/l.RE;
    l.mul=l.Vt*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
    l.u0=RandR(0.0,1.0);
    l.tetE= double(l.RE/AU/l.Dl);///marcs

    l.t0=RandR(2.0*year,3.0*year);
    //l.pt1=0.0;
    //l.pt2=5.0*year; //+2.5*l.tE;
  
    if(s.ro_star<=0.0 or l.tE<=0.0 or l.Dl>s.Ds or l.Vt<=0.0 or l.Ml<0.0){
    cout<<"ERROR ro_star:  "<<s.ro_star<<endl;
    cout<<"Vt:  "<<l.Vt<<endl;
    cout<<"RE: "<<l.RE/AU<<"\t xls:  "<<l.xls<<"\t tE: "<<l.tE<<endl;
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<endl;
    cout<<"numl:  "<<l.numl<<"\t Vt:  "<<l.Vt<<"\t mul:  "<<l.mul<<endl;  cin>>yye;}
   // cout<<"ksi:   "<<l.ksi<<"\t ro_star:  "<<s.ro_star<<endl;
    //cout<<"************** End of func_Lens  ****************"<<endl;
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye, uui, k1, k2, h, g;
    double metal, age,gravity, MB;
    char filename[40];
    FILE *fp2;


    double Age1[YZ]={0.0}; double B1[YZ]={0.0};  double M1[YZ]={0.0};   double mm1[YZ]={0.0};
    double Age2[YZ]={0.0}; double B2[YZ]={0.0};  double M2[YZ]={0.0};   double mm2[YZ]={0.0};
    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0};
    FILE *meta;
    meta=fopen("./files/CMD_WFIRST/metal.txt","r");
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf   %d  %d\n",&Metal[i],&count[i],&number[i]);
    if((Metal[i]<Metal[i-1] and i>0) or float(Metal[i])<-0.001 or number[i]==0 or count[i]>YZ or (abs(count[i-1]+number[i-1]-count[i])>2 and i>0)){
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta);
    FILE *hks;
    hks=fopen("./files/CMD_WFIRST/HKS.txt", "r");
    for(int i=0; i<YZ; ++i){
    fscanf(hks,"%lf  %lf  %lf  %lf\n",&Age1[i],&mm1[i],&B1[i],&M1[i]);
    if(Age1[i]<0.0 or mm1[i]<0.0 or fabs(B1[i])>0.3 or M1[i]<0.5 or Age1[i]>18.0){
    cout<<"ERROR Age(HKS): "<<Age1[i]<<"\t metal: "<<mm1[i]<<"\t B[i]"<<B1[i]<<"\t M[i]: "<<M1[i]<<"\t i: "<<i<<endl; cin>>uui; }}
    fclose(hks);
    FILE *ji;
    ji=fopen("./files/CMD_WFIRST/JI.txt", "r");
    for(int i=0; i<YZ; ++i){
    fscanf(ji,"%lf   %lf   %lf  %lf\n",&Age2[i],&mm2[i],&B2[i],&M2[i]);
    if(Age2[i]<0.0 or mm2[i]<0.0 or fabs(B2[i])>1.7 or M2[i]<0.5 or Age2[i]>18.0  or Age1[i]!=Age2[i] or mm1[i]!=mm2[i]){
    cout<<"ERROR Age(JI): "<<Age2[i]<<"\t metal: "<<mm2[i]<<"\t B[i]"<<B2[i]<<"\t M[i]: "<<M2[i]<<"\t i: "<<i<<endl;
    cout<<"Age1[i]:  "<<Age1[i]<<"\t mm1[i]:  "<<mm1[i]<<endl;  cin>>uui;}}
    fclose(ji);




////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','i','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Rs_d[j],&MB,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.cl_d[j],&cm.type_d[j]);
    ///*******************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0])         h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_d[3][j]= double(B1[g]+M1[g]*cm.Mab_d[2][j]*1.05263157894737); ///H-band  Mks/Mk=1.05263157894737
    cm.Mab_d[4][j]= double(B2[g]+M2[g]*cm.Mab_d[1][j]);   ///J-band
    if(fabs(cm.Mab_d[3][j]-cm.Mab_d[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_d[4][j]-cm.Mab_d[1][j])>1.5){
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[4][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
    ///********************************************************
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.Teff_d[j]<0.0 or metal>0.1 or age>10 or
    int(cm.cl_d[j])==6 or float(cm.type_d[j])>8.0 or cm.type_d[j]<1.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;






////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','b','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Rs_b[j],&MB,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.cl_b[j],&cm.type_b[j]);
    ///*****************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_b[3][j]= double(B1[g]+M1[g]*cm.Mab_b[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_b[4][j]= double(B2[g]+M2[g]*cm.Mab_b[1][j]);   ///J-band
    if(fabs(cm.Mab_b[3][j]-cm.Mab_b[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_b[4][j]-cm.Mab_b[1][j])>1.5){
    cout<<"ERROR: Mab_b(y-band): "<<cm.Mab_b[4][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///*****************************************************
    if(cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.Teff_b[j]<0.0 or age>10 or metal>0.9 or cm.cl_b[j]==6 or
    cm.type_b[j]>=8.0 or (cm.cl_b[j]==5 and int(cm.type_b[j])>8.0) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) ){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl;  cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;






////=================================== THICK DISK =============================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','k','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Rs_t[j],&MB,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.cl_t[j],&cm.type_t[j]);
    ///*********************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_t[3][j]= double(B1[g]+M1[g]*cm.Mab_t[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_t[4][j]= double(B2[g]+M2[g]*cm.Mab_t[1][j]);   ///J-band
    if(fabs(cm.Mab_t[3][j]-cm.Mab_t[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_t[4][j]-cm.Mab_t[1][j])>1.5){
    cout<<"ERROR: Mab_t(y-band): "<<cm.Mab_t[4][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///********************************************************
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or metal>0.2||cm.cl_t[j]==6|| cm.type_t[j]>=8.0 or
    (cm.cl_t[j]==5 and int(cm.type_t[j])>8) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;








////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','h','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Rs_h[j],&MB,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.cl_h[j],&cm.type_h[j]);
    ///************************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_h[3][j]= double(B1[g]+M1[g]*cm.Mab_h[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_h[4][j]= double(B2[g]+M2[g]*cm.Mab_h[1][j]);   ///J-band
    if(fabs(cm.Mab_h[3][j]-cm.Mab_h[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_h[4][j]-cm.Mab_h[1][j])>1.5){
    cout<<"ERROR: Mab_h(y-band): "<<cm.Mab_h[4][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///***********************************************************
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<0 or cm.cl_h[j]<0  or cm.cl_h[j]==6  or  cm.Teff_h[j]<0.0 or
    metal>0.1 || cm.cl_h[j]>7|| cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>8) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or
    (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;


}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s, int yy)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nnf=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   s.Romins=10000000000.0;
   double fd=1.0;///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;

   double alm= 2.0; ////2KPC,  Mancini 2004
   double rlm_0= 1.76*0.01;/// solar mass/ Pc^-3
   double M_dlm=2.6;//[Msun]mass of disk LMC  from Mancini 2004 paper
   double M_blm=0.15*M_dlm;///[Msun] mass of bar LMC
   double xb0, yb0, zb0;///KPC
   double frac=0.05; // fraction of halo in the form of compact objects
   double qq=0.688; 
   double Rd0=1.54;
   double xol, yol, zol, x0, y0, z0; 
   double r0,  pos1, inc1, Rlm; 

  char filename[40];
  FILE *fill;     
  FILE *fil1;
  fil1=fopen("./files/density/right.txt", "a+"); 
  int flagf=0;
  if(fabs(s.decli+68.5)<0.8){
  flagf=1; 
  sprintf(filename,"./files/density/%c%d%c%d.dat",'d',0,'_',yy);
  fill=fopen(filename,"w");
   //if(!fill){cout<<"cannot open file longtitude : "<<s.lg<<"\t latitude: "<<s.bg<<endl;  exit(0);}
  }


for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;
   s.rho_hlmc[i]=s.rho_dlmc[i]=s.rho_blmc[i]=0.0; 


   x=double(i*step);
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = Dsun- x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);
   double rsun= sqrt(zb*zb+ yb*yb+ xb*xb); 


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///Msun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nnf)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nnf)*exp(-fabs(zb)/0.8)/(1.0+0.5*nnf);///Msun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=frac*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///Msun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///Msun/pc^3
///==================================================================



///=====================  LMC density profiles ====================== 
   x0 = -x*cos(s.decli/RA)*sin(s.right/RA-RaLMC/RA);
   y0 =  x*sin(s.decli/RA)*cos(DecLMC/RA)-x*cos(s.decli/RA)*sin(DecLMC/RA)*cos(s.right/RA-RaLMC/RA);
   z0 = DLMC- x*cos(s.decli/RA)*cos(DecLMC/RA)*cos(s.right/RA-RaLMC/RA)-x*sin(s.decli/RA)*sin(DecLMC/RA);
   if(fabs(x-DLMC)<double(step*1.5)){s.xv=x0;  s.yv=y0;   s.zv=z0;}
   r0= sqrt(x0*x0+y0*y0+z0*z0); 
 
  
///=====================  LMC density HALO ======================
   if(r0<15.0) s.rho_hlmc[i]=fabs(rlm_0*frac/(1.0+r0*r0/alm/alm) );  ///Msun/PC^3
   else        s.rho_hlmc[i]=0.0;   
  
  
///=====================  LMC density Disk ======================
   pos1=(170.0-90.0)*M_PI/180.0; 
   inc1=34.7*M_PI/180.0;  
   xol= x0*cos(pos1)+ y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   Rlm= sqrt(xol*xol + yol*yol ); 
   double zd0=0.3;//KP Kim's paper
   double Rd1=1.8;///KPc  Kim's paper
   s.rho_dlmc[i]=fabs(M_dlm/(4.0*M_PI*zd0*Rd1*Rd1)*exp(-Rlm/Rd1)*exp(-fabs(zol/zd0))); ////kim [Msun/Pc^3]  Kim 2000
   

///=====================  LMC density Bulge ======================
   xb0= 1.2;  yb0=zb0= 0.44;
   pos1=(110.0-90.0)*M_PI/180.0; 
   inc1=0.0*M_PI/180.0;  
   xol= x0* cos(pos1) +  y0* sin(pos1);
   yol=-x0*sin(pos1)*cos(inc1) + y0*cos(pos1)*cos(inc1) -z0*sin(inc1);
   zol=-x0*sin(pos1)*sin(inc1) + y0*cos(pos1)*sin(inc1) +z0*cos(inc1);    
   s.rho_blmc[i]=M_blm*pow(2.0*M_PI,-1.5)/(xb0*yb0*zb0)*exp(-0.5*(xol*xol/xb0/xb0+yol*yol/yb0/yb0+zol*zol/zb0/zb0));///[Msun/Pc^3]  
///===============================================================================



s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]) +fabs(s.rho_hlmc[i])+fabs(s.rho_dlmc[i])+ fabs(s.rho_blmc[i]);///[Msun/pc^3]

s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Msun/deg^2]

s.Nstari[i]= binary_fraction*((s.rho_disk[i]+s.rho_dlmc[i])*fd/0.403445+s.rho_ThD[i]*fh/0.4542+(s.rho_halo[i]+s.rho_hlmc[i])*fh/0.4542+(s.rho_bulge[i]+s.rho_blmc[i])*fb/0.308571);////[Nt/pc^3] 

s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=  s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]

if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
if(s.Rostari[i]<s.Romins) s.Romins=s.Rostari[i];///source selection

//cout<<"Rostari[i]:  "<<s.Rostari[i]<<"\t Ds:  "<<x<<endl;

if(flagf>0)
fprintf(fill,"%e   %e   %e   %e   %e  %e  %e   %e   %e   %e\n",
  x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.rho_dlmc[i],s.rho_blmc[i],s.rho_hlmc[i], s.Rostar0[i],s.Nstari[i]); 
 // cout<<"rho_disk(LMC):  "<<s.rho_dlmc[i]<<"\t rho_bar(LMC): "<<s.rho_blmc[i]<<endl;
  //cout<<"rho_halo(LMC):  "<<s.rho_hlmc[i]<<endl;
}

if(flagf>0) fprintf(fil1,"%.5lf  %.5lf  %.5lf  %.5lf\n",s.right, s.decli, s.TET*RA, s.FI*RA);

//cout<<"Romins:  "<<s.Romins<<"\t Romaxs:  "<<s.Romaxs<<endl;

if(flagf>0)   {fclose(fill); fclose(fil1);}
//int rre; cin>>rre;

   //cout<<"xLMC:  "<<xLMC<<"\t yLMC:  "<<yLMC<<"\t zLMC:  "<<zLMC<<endl;
   //cout<<"LMC_distance:  "<<sqrt(xLMC*xLMC + yLMC*yLMC +  zLMC*zLMC)<<endl;
   //cout<<"xol:  "<<xol<<"\t yol:  "<<yol<<"\t  zol:  "<<zol<<endl;
   //cout<<"distance_LMC_Center:  "<<rlm<<"\t projected_distance:  "<<Rlm<<endl;
  // cout<<"Nstart [Nt/deg^2]: "<<s.Nstart<<"\t Ro_star [Mass/deg^2]: "<<s.Rostart<<endl;
   //cout<<">>>>>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<"<<endl;
   //exit(0);
   
}
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nnd; ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return rr;
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
/*
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0 or Lon<0.0 or Lon<0.1  or s.lon<0.1) Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 || Lon<0.25 || fabs(Lat)>10.0 || (Lon>100 && Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     //cout<<"Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;

     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 && fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  || ex.dis[i]>50.0 || ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }}
     fclose(fpd);}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<flag<<endl;
     return(flag);
}
*/
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;



  double NN=2.5;
  double sigma_R_Disk, sigma_T_Disk, sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){ Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}

 
for (int i=0;i<2; ++i){
 test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])                       {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))              {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))       {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])){sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))
                                       {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5]))
                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                    {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}
    
    
    double v_R_lmc=-57.0;
    double v_T_lmc=-226.0; 
    double v_Z_lmc= 221.0;
    double sigma_LMC=20.2; 
    double err_rlmc= 13.0; ///error of global velocity
    double err_tlmc= 15.0; 
    double err_zlmc= 19.0; 


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    
    
    else if(s.struc>3){
    SVR= RandN(sigma_LMC, NN); 
    SVT= RandN(sigma_LMC, NN); 
    SVZ= RandN(sigma_LMC, NN); 
    SVZ +=  v_Z_lmc +RandR(-err_zlmc, err_zlmc);
    SVR +=  v_R_lmc +RandR(-err_rlmc, err_rlmc);
    SVT +=  v_T_lmc +RandR(-err_tlmc, err_tlmc);}
    //s.vs=sqrt(SVT*SVT+SVZ*SVZ+SVR*SVR);
        
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/R_sun,0.0394) + 0.00712);
    
    s.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   
   
   else if(l.struc>3){
   LVR= RandN(sigma_LMC, NN); 
   LVT= RandN(sigma_LMC, NN); 
   LVZ= RandN(sigma_LMC, NN); 
   LVZ +=  v_Z_lmc +RandR(-err_zlmc, err_zlmc);
   LVR +=  v_R_lmc +RandR(-err_rlmc, err_rlmc);
   LVT +=  v_T_lmc +RandR(-err_tlmc, err_tlmc); }
   
   
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/R_sun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(R_sun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(R_sun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= pi-fabs(tetd);
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    

    if(vls1>0.0 and vls2>=0.0)         s.tetv=atan(fabs(vls2)/fabs(vls1));
    else if(vls1<=0.0 and vls2>0.0)    s.tetv=atan(fabs(vls1)/fabs(vls2)) + M_PI/2.0;
    else if(vls1<0.0 and vls2<=0.0)    s.tetv=atan(fabs(vls2)/fabs(vls1)) + M_PI;
    else if(vls1>=0.0 and vls2<0.0)    s.tetv=atan(fabs(vls1)/fabs(vls2)) + 3.0*M_PI/2.0;
    else if(s.tetv>2.0*M_PI)           s.tetv-= 2.0*M_PI; 
    
    if(vls1==0.0 and vls2==0.0){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    int iee;  cin>>iee;}

    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<s.vs<<endl;
    cout<<"source type:  "<<s.type<<"\t s.struc:  "<<s.struc<<endl;
    int yee; cin>>yee;}
    cout<<"(vrel_func):   s.deltao:  "<<s.deltao<<"FI: "<<s.FI<<endl;
    
///*****************************************************
}



   




///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
