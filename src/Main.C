#include "Brem.h"

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

//extern unsigned long int gsl_rng_default_seed;

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include "TString.h"

using namespace HepMC3;
using namespace std;


int main(int argc, char **argv) {

  std::string Me( argv[0] );
  int p = Me.rfind('/');
  Me.erase(0, p+1);
  // defaults
  std::string LastSeedFile("LastRunInfo.txt");
  double e_beam = 10.;
  double p_beam = 275.;
  double e_photon_min = 0.1 ;
  double e_photon_max = 9.99;   
  int seed_mode = -1;
  int Nevents = 100;
  double tgmax  = 0.001;       // default theta_max of the bremstrahlung photon
  //
  if( argc<3 ) {
    printf("\t Usage: %s :  <number_of_events> <HepMC3_output_file> <photon_theta_max> <generator_seed>\n",Me.c_str());
    printf("\t\t\t       -  <number_of_events>       obligatory \n"); 
    printf("\t\t\t       -  <HepMC3_output_file>     obligatory \n");
    printf("\t\t\t       -  <photon_theta_max>       optional photon theta_max, if ommitted default = 0.001 [radian] \n");
    printf("\t\t\t       -  <generator_seed>         optional may have three values (default is -1) : \n");
    printf("\t\t\t                             if = 0 -> Start ranlux with the default seed \n"); 
    printf("\t\t\t                             if > 0 -> Seed ranlux with default given seed \n");
    printf("\t\t\t                             if < 0 -> or ommited use seed saved in previous run in file %s \n",LastSeedFile.c_str());
    printf("                        ===========================================================================================\n");
    exit(-1);
  }else if(argc==3) {
    cout<<"\n\n\t"<<argv[0]<<" called for:"
	<<"\n\t\t   -events: "<<argv[1]
	<<"\n\t\t   -output file: "<<argv[2]
	<<"\n\t\t   -max. photon theta:"<<tgmax
    	<<"\n\t\t   -using seed from file : "<<LastSeedFile
	<<endl;
    Nevents = atoi(argv[1]);
  }else if(argc==4) {
    cout<<"\n\n\t"<<argv[0]<<" called for\t"<<argv[1]<<" events\n "
	<<"\t\t output file: "<<argv[2]
	<<"\n\t\t max. photon theta:"<<argv[3]
    	<<"\n\t\t using seed from file "<< LastSeedFile
	<< endl;
    Nevents = atoi(argv[1]);
    tgmax = atof(argv[3]);
  }else if(argc==5){    
    cout<<"\n\n\t"<<argv[0]<<" called for\t"<<argv[1]<<" events\n "
	<<"\t\t output file: "<<argv[2]
	<<"\n\t\t max. photon theta:"<<argv[3]
	<<"\n\t\t         seed mode:"<<argv[4]
	<< endl;
    Nevents = atoi(argv[1]);
    tgmax = atof(argv[3]);
    seed_mode = atoi(argv[4]);
  }else if(argc>5) {
    cout<<"\n\n\t Error!!! "<<argv[0]<<" called with more than four parameters. argc = "<<argc-1<<"\n ";
    exit(1);
  }

  /* 
     Brem * b = new Brem(e_beam, p_beam,
     e_photon_min, e_photon_max,
     seed, output_name,
     tgmax);
  */
  std::string output_name( argv[2] );
  

  Brem * brem = new Brem(e_beam, p_beam, e_photon_min, e_photon_max, output_name, seed_mode, tgmax);
  
  brem->run_Brem( Nevents  );
  return 0;
}
