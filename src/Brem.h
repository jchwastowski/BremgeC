/*
The generator of the bremsstrahlung events for: e-p, e-gas and e-p beam size effect modified
 - BGREMGE Fortran generator of L. Suszycki 
 - The beam size effect corrections by K. Piotrzkowski
references:
   Proc. of the Physics at HERA Workshop, eds. W. Buchmueller and G. Ingelman
   Hamburg, Oct. 29-30, 1991, p.1463
 - C++ implentation: 
       authors: J. J. Chwastowski, IFJ PAN, Krakow, Poland
                Janusz.Chwastowski@ifj.edu.pl
		B. Pawlik, AGH University of Krakow, Krakow, Poland
                pawlik.bogdan@gmail.com
*/
#include <ctime>
#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <random>

#include <math.h>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "TROOT.h"
#include "TMath.h"
#include "TRint.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"


//extern unsigned long int gsl_rng_default_seed;

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "BremConst.h"

using namespace HepMC3;
using namespace std;


class Brem {

private:

  double   E_p0;           //  incident proton energy (GeV)
  double   E_e0;           //  incident electron energy (GeV)
  double   Eg_min;         //  minimum energy of the radiated photon (GeV)
  double   Eg_max;         //  maximum energy of the radiated photon (GeV)
  int      Seed_Mode;      //  <0/0/>0  set seed from previousrun/gsl default/seed given by user
  int      ModeBr;         //  generator mode 0/1/2 ep bremsstrahlung/ ep+beam size effect/e-gas
  std::string outname;     //  name of the hepmc3 output file
  double   Tg_max;         //  radiated photon polar angle theta max 
  double   Z_gas;          //  e-gas
  double   xVTX, yVTX, zVTX;

  
  unsigned long long Seed; //  seed feeded to gsl_rng  
  int      ntry;           //  number of iterations
  double   Q2;             //  Q^2 of the bremsstrahlung process
  
  // Lorentz gammas of the proton, electron and the scattered electron
  double gamma_p0;
  double gamma_e0;
  double gamma_e;
  double SigTot;     // Bremge calculated cross-section
  // internal
  double delta, delta_max, delta_pr, d2_min, d2_max, du3;
  double Eg, tgx, tgy, gth, Ee, tex, tey, eth;
  double mDET_POSITION;                        // distance to the face of epic LumiDirectPC detector [mm]
  double xpos, ypos, zpos;                     // radiated photon position (x, y) at mDET_POSITION
  double x_ep, y_ep, x2_ep, y2_ep, xy2_ep;     // beam size
  double beam_sigmaX, beam_sigmaY;             // beam angular spread
  double Z13, BE12, BE22, BE32, SCRMAX;
  double fGeGp;                                // factor containing the Lorentz gammas
  // auxiliaries
  double Ehl, SCR,  EBound[2], DsDuma;
  TLorentzVector gamma_vec;
  gsl_rng *r;
  GenEvent hepevt;

  double unif(){
    return gsl_rng_uniform (r);
  }
  
  double gaus(){
    return gsl_ran_gaussian_ziggurat(r,beam_sigmaX);    
  }

public:
  double proton_in[4];
  double lepton_in[4];
  double proton_out[4],lepton_out[4],gamma_out[4];
  TString outname_root;
  TString outname_hepmc3;
  TFile* rootfile;
  TTree *BremTree;
  TH1D  *hphot_en, *hphot_theta,*hele_en, *hele_theta;
  TH2D  *hphot_edepxy, *hphot_thxthy;
  TString RunNumberStr;
  int     RunNumber;
  // seters
  void set_Seed(unsigned long int s) {gsl_rng_set( r, s) ; Seed = s; }
  void set_Modebr(int m) { ModeBr = m; }
  void set_Eg_min(double e) { Eg_min = e;}
  void set_Eg_max(double e) { Eg_max = e;}
  void set_E_p0(double e)   { E_p0 = e;  }
  void set_E_e0(double e)   { E_e0= e; }
  void set_detector_position( double pos) { mDET_POSITION = pos; }
  void set_VTX( double x, double y, double z) { xVTX=x; yVTX=y; zVTX=z; }
  void set_beam_spread( double sigx, double sigy ){ beam_sigmaX=sigx; beam_sigmaY=sigy; }
  // geters
  long unsigned int  get_Seed() { return Seed;}
  int    get_Modebr(){ return ModeBr; }
  double get_Z_gas() { return Z_gas;}
  double get_Eg_min(){ return Eg_min;}
  double get_Eg_max(){ return Eg_max;}
  double get_E_p0()  {return E_p0;}
  double get_E_e0()  {return E_e0;}
  double get_Tg_max(){return Tg_max;}
  int    get_Ntry()  {return ntry;}
  double get_SigTot(){return SigTot;}
  void   get_VTX( double &x, double &y, double &z) { x=xVTX; y=yVTX; z=zVTX; } 
  double get_Q2()    {return Q2;}         // to be called after the call to bremge
  gsl_rng* get_gsl_r(){ return r; } 
  //
  
  Brem(double Ee0 = 10., double Ep0 = 275.,
       double Egmin = 0.1, double Egmax = 9.9,
       std::string output_name = "bremgeout.hepmc3", int seed_mode = -1,
       double Tgmax = 0.001, double Zgas = 4.2 ) :E_e0(Ee0),E_p0(Ep0),
						  Eg_min(Egmin),Eg_max(Egmax),
						  outname(output_name),
						  Seed_Mode(seed_mode),
						  Tg_max(Tgmax),
						  Z_gas(Zgas)
  {

    // Standard constructor setting the default values

    mDET_POSITION = -65000;
    
    ModeBr = 0;
#ifdef BREM_MODE_1
    ModeBr = 1;
#endif
#ifdef BREM_MODE_2
    modeBr = 2;
#endif

    xVTX=0.;
    yVTX=0.;
    zVTX=0.;
#ifdef LUMI_DIRECT_TEST
    zVTX = mDET_POSITION;
#endif

    beam_sigmaX = 0.0002;
    beam_sigmaY = 0.0002;
    
    x_ep = 0.00001;
    y_ep = 0.00001;
    x2_ep = x_ep*x_ep;
    y2_ep = y_ep*y_ep;
    xy2_ep = x2_ep*y2_ep*CONV;

    gamma_e0 = E_e0/e_mass;
    gamma_p0 = E_p0/p_mass;

    delta_max = Tg_max*gamma_e0+delta_min;
    d2_min = delta_min*delta_min;
    d2_max = delta_max*delta_max;
    du3 = (1+d2_max)/(d2_max-d2_min);

    if(ModeBr==2) {
      fGeGp = 4.*gamma_e0;
      Z13 = pow(Z_gas,1./3.);
      BE12 = (Z13/BCON*B1)*(Z13/BCON*B1);
      BE22 = (Z13/BCON*B2)*(Z13/BCON*B2);
      BE32 = (Z13/BCON*B3)*(Z13/BCON*B3);
    }
    else {
      fGeGp = 4.0*gamma_e0*gamma_p0;
    }
    Ehl = log(Eg_max/Eg_min);
    SCR = 1.0;
    
    rootfile=nullptr;
    BremTree=nullptr;
    
    // gsl random generator setup and Run Number
    //
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    T = gsl_rng_ranlux; 
    r = gsl_rng_alloc (T);
    
    // get today in form last dof the year and the number this might be eg 23001 (first day of 2023)   
    int mDate;
    int runNum;
    std::ifstream instrm;
    system("date +%g%j 1> tmp.txt");    
    instrm.open("tmp.txt",ifstream::in);
    if ( instrm.is_open() ) { instrm >> mDate; }
    system("rm tmp.txt");
    instrm.close();
    // check if LastRunIfo exists
    instrm.open("LastRunInfo.txt",ifstream::in);
    if( instrm.is_open() ){
      instrm >> RunNumber >> Seed ;
      instrm.close();
      if ( RunNumber/1000 == mDate ) { RunNumber++; RunNumberStr = TString( to_string(RunNumber)); RunNumberStr.Prepend("run");}     
    } else {
	RunNumber  = mDate*1000;
	RunNumberStr = TString( to_string(RunNumber) );
	set_Seed(0);
	cout<<"\t Brem::Brem: new LastRunInfo.txt will be created on exit. \n";
    }

    if(Seed_Mode < 0) {
      set_Seed(Seed);
      cout<<"\n\t Brem::Brem: ranlux seed " << gsl_rng_get(r)<<" from LastRunInfo.txt \n";
    }else if( Seed_Mode == 0 ){
      set_Seed(0); 
      cout<<"\n\t Brem::Brem: Seed_Mode "<<Seed_Mode<<" ranlux seeded with:\t"<< gsl_rng_get(r) <<"\t user given seed, not saved "<<endl;      
    }else if( Seed_Mode > 0) {
      set_Seed( Seed_Mode );
      cout<<"n\t Brem::Brem: Seed_Mode "<<Seed_Mode<<" ranlux seeded with:\t"<< gsl_rng_get(r) <<"\t user given seed, not saved "<<endl;  
    }
       
    //
    // output file name preparation
    
    outname_root=TString(outname);
    int pos1 = outname_root.Last('/')+1;
    int pos2 = outname_root.Last('.');
    TString dir = TString( outname_root(0, pos1));
    outname_root = TString( outname_root(pos1,pos2-pos1));      				     
    outname_root.Prepend(RunNumberStr + "_");
    outname_root.Prepend(dir);
    outname_hepmc3 = outname_root;
    //
    outname_root.Append(".root");
    outname_hepmc3.Append(".hepmc3");
   if ( dir != "" ){
      TString cmd = "mkdir -p ";
      cmd  += dir;
      system(cmd.Data());
    }

    Brem_Crossection();
    Brem_Root_Init();
  }

  void Brem_Root_Init(){
    
      cout<< "\t Brem::Brem_root_init: ROOT output file: "<<outname_root<<endl;
    
      rootfile = new TFile(outname_root,"RECREATE");
      
      BremTree = new TTree("BremTree", "bremstrahlung ep events e_beam=10GeV pbeam=275GeV photon_theta_max=0.001[rad]");      
      BremTree->Branch("el_en",&Ee);
      BremTree->Branch("el_thx",&tex);
      BremTree->Branch("el_thy",&tey);
      BremTree->Branch("el_theta",&eth);
      
      BremTree->Branch("phot_en",&Eg);
      BremTree->Branch("phot_thx",&tgx);
      BremTree->Branch("phot_thy",&tgy);
      BremTree->Branch("phot_theta",&gth);
      
      BremTree->Branch("phot_posx",&xVTX);
      BremTree->Branch("phot_posy",&yVTX);
      
      double thmax = Tg_max*1000. ; //radians to miliradians
      hphot_en      = new TH1D ("hphot_en", "brem photon energy; E_{photon}/E_{electron beam}; N_{evt}" , 105, 0., 1.05);
      hphot_theta   = new TH1D ("hphot_theta", "brem photon polar angle; #theta_{#gamma} [mrad]; N_{evt}" , 100, 0., thmax);
      hphot_edepxy  = new TH2D ("hphot_edepxy","EdepXY;X [mm];Y [mm]",240,-120.,120.,240,-120.,120.);
      hphot_thxthy = new TH2D ("hphot_thxthy","EdepThxThy;#theta_{X} [mrad];#theta_{Y} [mrad]",200,-thmax,thmax,200,-thmax,thmax);
      hphot_en->Sumw2();
      hphot_theta->Sumw2();
   
      hele_en      = new TH1D ("hele_en", "out electron energy; E_{electron_out/E_{electron beam}; N_{evt}" , 105, 0., 1.05);
      hele_theta   = new TH1D ("hele_theta", "out electron polar angle; #theta [mrad]; N_{evt} " ,100,0.,thmax ); 
      hele_en->Sumw2();
      hele_theta->Sumw2();
  }
  
  void Brem_Crossection(){
    /*
      get the totalbremsstrahlung cross-section in [Eg_min, Eg_max]
      Trapezoidal integration in log(E_gamma)
    */

    EBound[0] = Eg_min;
    EBound[1] = Eg_max;
    SigTot = 0.0;
    for(int i=0;i<2;i++) {
      double Eg = EBound[i];
      double Ee = E_e0-Eg;
      double EE0EE = E_e0/Ee;
      double WE;
      if(ModeBr==2) {
	WE = (E_e0/(E_e0-Eg)+(E_e0-Eg)/E_e0-2./3.)*
	  (Z_gas*Z_gas*log(184./pow(Z_gas,1./3.))
           +Z_gas*log(1194./pow(Z_gas,2./3.)))
	  +Z_gas*(Z_gas+1.0)/9.0;
      }
      else {
	WE = (E_e0/(E_e0-Eg)+(E_e0-Eg)/E_e0-2./3.)*(log(fGeGp*Ee/Eg)-1./2.);
      }
      double ds = 4*ALR2MB*WE*(E_e0-Eg)/E_e0;
      SigTot += ds;
    }
    SigTot = SigTot*log(Eg_max/Eg_min)/2.0;

    /*
    set the density cut to 10.*SigTot: approx. 10% efficienncy of the generator
    */

    DsDuma = 10.0*SigTot;


  }
  //
  // generator
  int Bremge(double &Eg, double &tgx, double &tgy, double &Ee, double &tex, double &tey) {
    /*
      bremge output:
      returns gives the  number of the iterations
      Eg  - the photon energy
      tgx - arctan(pgx/pgz), wher pgz, pgy, pgz are the components of the photon momentum
      tgy - arctan(pgy/pgz)
      Ee  - the scattered electron energy
      tex - arctan(pex/pez)
      tey - arctan(pey/pez)
    */
    ntry = 0;
    while(1) {
      ntry++;

      // get the photon and electron energies

      double u1 = unif();
      Eg = Eg_min*exp(u1*log(Eg_max/Eg_min));
      Ee = E_e0-Eg;

      double EgEe = Eg/(Ee);
      double XI2 = Eg*EgEe/E_e0/2.0;
      gamma_e =(Ee)/e_mass;

      /* 
	 get delta = Th_g*gamma_e0
	 Th_g is the photon polar angle w.r.t. the initial electron direction
      */
      double u2 = unif();
      delta = sqrt((u2+d2_min*du3)/(du3-u2));
      double d21 = delta*delta+1.0;
      /*
	now the angle between the photon and the scattered electron planes
      */
      double Dfim = d21*EgEe/fGeGp;
      double Fim = Dfim/delta;
      double Pifim = log(1.+M_PI/Fim);
      double u3 = unif();
      double Fi = Fim*(exp(u3*Pifim)-1.);
      double Dfifim = (Fi+Fim)*delta;
      double AU1 = (delta_max-delta)/Dfifim;
      AU1 = atan(AU1);
      double BU1 = (delta-delta_min)/Dfifim;
      BU1 = atan(BU1);
      /*
	get delta-delta' (DDPR)
	delta' = Th_e*gamma_e 
      */
      double u4 = unif();
      double DDPR= Dfifim*tan(u4*(BU1+AU1)-AU1);
      double DDPR2 = DDPR*DDPR;
      double YAFIU = log(Eg_max/Eg_min)*d21*d21*Pifim*(AU1+BU1)*(Dfifim*Dfifim+DDPR2)/(2.*(1.+d2_min)*delta*delta*du3);
      /*
	now calculate dSigma/(dFi dTh_e dTh_g dE_g)
      */
      double delta_pr = delta-DDPR;
      double DMDPR = delta*delta_pr;
      double dpr21 = delta_pr*delta_pr+1.;
      double DDPR21 = d21*dpr21;
      double WQ = (-d21+E_e0/(Ee)*dpr21)/fGeGp;
      double Q2M2 = WQ*WQ+DDPR2;
      double WW0 = DDPR*(1.-DMDPR)/DDPR21;
      WW0 = WW0*WW0+XI2*DDPR2/DDPR21;
      double WW1 = (1+XI2)/DDPR21;
      double WFI = DMDPR*2.*(1.-cos(Fi));
      if(Fi<=0.01) {
	double WFI = DMDPR*Fi*Fi;
      }
      Q2M2 = Q2M2+WFI;
      /*
       process Q2
      */
      Q2 = Q2M2*e_mass*e_mass;
      double Q2TR = WFI+DDPR2;
      if(ModeBr==2) {
	double SCR = Q2M2* (AL1/(BE12+Q2M2)+AL2/(BE22+Q2M2)+AL3/(BE32+Q2M2));
	SCR = SCR*SCR;
	u1 = unif();
	if(u1*SCRMAX>SCR) continue;
      }
      double DSDFI = AR28PI*DMDPR*(WW0+WW1*WFI)/Q2M2/Q2M2;
      double DSDU = DSDFI*2.*Ee/E_e0*YAFIU;
      u1 = unif();
      if(DSDU<DsDuma*u1) continue;
      double TG = delta/gamma_e0;
      double TE = delta_pr/gamma_e;
      double TX = TE*sin(Fi);
      if(unif()>0.5) TX = -TX;
      double TY =-TE*cos(Fi)+TG;

      double PHI = 2.*M_PI*unif();
      double CPHI = cos(PHI);
      double SPHI = sin(PHI);
      tgx = -TG*SPHI;
      tgy = TG*CPHI;
      tex = TX*CPHI-TY*SPHI;
      tey = TX*SPHI+TY*CPHI;
      /* beam-size effect */
      if(ModeBr==1) {
        double FEXP=(x2_ep*CPHI*CPHI+y2_ep*SPHI*SPHI)/(xy2_ep*Q2TR);
	if(FEXP>30) {
	  Ee = 0.;
	  return(0);
	}else {
	  if(unif()<exp(-FEXP)) return(0);
	}
      }
      return(0);
    }
  }
  //
  //  converter
  void Brem_Event_HepMC3() {
    
    int j = Bremge(Eg, tgx, tgy, Ee, tex, tey);

    hepevt = GenEvent(Units::GEV,Units::MM);
    // proton in
    double p_pz = sqrt(get_E_p0()*get_E_p0()-p_mass*p_mass);
    proton_in[0] = 0.0;
    proton_in[1] = 0.0;
    proton_in[2] = p_pz;
    proton_in[3] = get_E_p0();
    // electron in
    double e_pz = -sqrt(get_E_e0()*get_E_e0()-e_mass*e_mass);
    lepton_in[0] = 0.0;
    lepton_in[1] = 0.0;
    lepton_in[2] = e_pz;
    lepton_in[3] = get_E_e0();
    // bremsstrahlung photon
#ifdef BEAM_SPREAD_THETA
    tgx += gaus();
    tgy += gaus();
    //    double COSTH_PHOT =1./sqrt(1.+tan(tgx)*tan(tgx)+tan(tgy)*tan(tgy));
    //#else
#endif
    double COSTH_PHOT =1./sqrt(1.+tan(tgx)*tan(tgx)+tan(tgy)*tan(tgy));
    //#endif
    double g_pz = -Eg*COSTH_PHOT;
    double g_px =  g_pz*tan(tgx);
    double g_py =  g_pz*tan(tgy);
    gamma_out[0] = g_px;
    gamma_out[1] = g_py;
    gamma_out[2] = g_pz;
    gamma_out[3] = Eg;    
    gth = acos( COSTH_PHOT );
    xVTX = zVTX*tan(tgx);
    yVTX = zVTX*tan(tgy);
    
    hphot_en->Fill(Eg/get_E_e0());
    hphot_theta->Fill(gth*1000.);
    hphot_edepxy->Fill(xVTX, yVTX);
    hphot_thxthy->Fill(tgx*1000.,tgy*1000.);
    //
    // electron out
    double e_pp = sqrt(Ee*Ee-e_mass*e_mass);
    double CTH_ELE = 1./sqrt(1.+tan(tex)*tan(tex)+tan(tey)*tan(tey));
    eth = acos(1./sqrt(1.+tan(tex)*tan(tex)+tan(tey)*tan(tey)));
    e_pz = -e_pp*CTH_ELE;
    double e_px = e_pz*tan(tex);
    double e_py = e_pz*tan(tey);
    lepton_out[0] = e_px;
    lepton_out[1] = e_py;
    lepton_out[2] = e_pz;
    lepton_out[3] = Ee;

    hele_en->Fill(Ee/get_E_e0());
    hele_theta->Fill(eth*1000.);
    
    /*
    Virtual photon using HEPMC3 convention
    */

    double q[4];
    q[0] = lepton_in[0]-e_px-g_px;
    q[1] = lepton_in[1]-e_py-g_py;
    q[2] = lepton_in[2]-e_pz-g_pz;
    q[3] = lepton_in[3]-Ee-Eg;
    /*
    double sumx = (g_px+e_px+q[0])*1.e9; // eV
    double sumy = (g_py+e_py+q[1])*1.e9; // eV
    double sumz = (g_pz+e_pz+q[2]); // GeV
    double sumE = (Eg+Ee); // GeV
    printf("\t  test: momenta sum ele_out+photon+virt_gamma (px, py, pz, E) = ( %6.2f eV, %6.2f eV, %6.1f GeV ,%6.1f GeV )\n",sumx,sumy,sumz, sumE);
    */
 #ifdef DEBUG   
    if(q[3]!=0.0) {
      cout<<"\n\n 0  q   "<<q[3]<<'\t'<<q[2]<<'\t'<<q[1]<<'\t'<<q[0]<<endl;
      cout<<" 0b en "<<lepton_in[3]<<'\t'<<Ee<<'\t'<<Eg<<'\t'<<(Ee+Eg)/get_E_e0().<<endl;
      double qq = q[3]*q[3]-q[0]*q[0]-q[1]*q[1]-q[2]*q[2];
      cout<<" 1  qq    "<<-qq<<" g Q2 "<<get_Q2()<<" ratio "<<-qq/get_Q2()<<endl;
      q[3] = 0.0;
      qq = q[3]*q[3]-q[0]*q[0]-q[1]*q[1]-q[2]*q[2];
      cout<<" 2 qq(4) "<<-qq<<" g Q2 "<<get_Q2()<<" ratio "<<-qq/get_Q2()<<endl; 
    }
#endif
    proton_out[0] = proton_in[0]+q[0];
    proton_out[1] = proton_in[1]+q[1];
    proton_out[2] = proton_in[2]+q[2];
    proton_out[3] = proton_in[3]+q[3];

    GenParticlePtr p1 = make_shared<GenParticle>( FourVector(proton_in[0],proton_in[1],proton_in[2],proton_in[3]),2212,  4 );
    GenParticlePtr p2 = make_shared<GenParticle>( FourVector(lepton_in[0],lepton_in[1],lepton_in[2],lepton_in[3]),  11,  4 );
    GenParticlePtr p3 = make_shared<GenParticle>( FourVector(q[0],   q[1], q[2], q[3]),  22,  3 ); //<-- virtual photon

    GenParticlePtr p4 = make_shared<GenParticle>( FourVector( gamma_out[0], gamma_out[1], gamma_out[2], gamma_out[3]),  22,  1 );
    GenParticlePtr p5 = make_shared<GenParticle>( FourVector(lepton_out[0],lepton_out[1],lepton_out[2],lepton_out[3]),  11,  1 );
    GenParticlePtr p6 = make_shared<GenParticle>( FourVector(proton_out[0],proton_out[1],proton_out[2],proton_out[3]),2212,  1 );


#ifdef LUMI_DIRECT_TEST
    GenVertexPtr v1 = make_shared<GenVertex>(FourVector(xVTX, yVTX, zVTX,0.0));
    v1->add_particle_in(p2);   // _in(p1); electron in
    v1->add_particle_out(p4);  // _in(p3); photon out
#else
    GenVertexPtr v1 = make_shared<GenVertex>(FourVector(xVTX, yVTX, zVTX,0.0));
    v1->add_particle_in(p2);   // _in(p1); electron in
    v1->add_particle_out(p5);  // _in(p2); electron out  
    v1->add_particle_out(p4);  // _in(p3); photon out
#endif
    hepevt.add_vertex(v1);
#ifndef LUMI_DIRECT_TEST
    GenVertexPtr v2 = make_shared<GenVertex>(FourVector(xVTX, yVTX, zVTX,0.0));
    v2->add_particle_in(p1);        //_out(p4); proton in
    v2->add_particle_out(p6);       //          proton out
    /*    v2->add_particle_out(p5);  */
    hepevt.add_vertex(v2);
#endif

#ifdef SAVE_P7
    GenParticlePtr p7 = make_shared<GenParticle>(FourVector(get_Q2(),0.,(double)get_Ntry(),0.),22,3);
    GenVertexPtr v3 = make_shared<GenVertex>();
    v3->add_particle_out (p7);
    hepevt.add_vertex(v3);
#endif
    BremTree->Fill();
  }
  
  void run_Brem(int Nevents ) {
        
    WriterAscii output_hepmc3(outname_hepmc3.Data());
  
    for(int i=0;i<Nevents;i++){
    
      Brem_Event_HepMC3();
    
      output_hepmc3.write_event(hepevt);      // Save event to output file
    }
  
    output_hepmc3.close();
    end_Brem();
  }
  
  void end_Brem(){    
    rootfile->cd();
    BremTree->Write();
    hele_en->Write();
    hele_theta->Write();
    hphot_en->Write();
    hphot_theta->Write();
    hphot_edepxy->Write();
    hphot_thxthy->Write();
    rootfile->Close();
    //
    ulong long sseed;  
    sseed = ( Seed_Mode < 0)? gsl_rng_get(r) : (Seed_Mode == 0 )? 0 :(Seed_Mode > 0 )? get_Seed(): 0; 
    ofstream ostrm("LastRunInfo.txt",ofstream::out);
    if ( ostrm.is_open() ){
      ostrm<< '\t' << RunNumber << '\t' << sseed << '\n';
      ostrm.close();
    }else{
      cout<< " can not open file LastRunInfo.txt for writing" << endl;
    }
    
      Print();
      gsl_rng_free (r);
  }
  void Print(){
  // print parameters used in this run
    
    printf("\n\n\t end of run number %10d  \n", RunNumber);
    cout<<"--------------------------------------------------------------------\n\n";
    printf(" - - electron beam energy        (      e_beam = %4.1f [GeV]) \n", E_e0);
    printf(" - - proton beam energy          (      p_beam = %4.1f [GeV]) \n",E_p0);     
    printf(" - - photon minimum energy       (e_photon_min = %4.1f [GeV]) \n",Eg_min); 
    printf(" - - photon maximum energy       (e_photon_max = %4.1f [GeV]) \n",Eg_max);
    printf(" - - photon maximum polar angle  (       tgmax = %4.3f [radians]) \n",Tg_max);
    //
    printf("\n--------------------------------------------------------------------\n");
    printf("\t generator reported cross section: %6.3f [ub] \n", SigTot);
    printf("\t - HepMC3 output file: %16s \n",outname_hepmc3.Data());
    printf("\t - ROOT   output file: %16s \n",outname_root.Data());
    cout<<"--------------------------------------------------------------------\n\n";

  }
  
};
