BremgeC is a generator of bremsstrahlung events in ep (electron-proton) and eA (electron-Ion) scattering, 
ported to C++.
It is based on the original Bremge written in Fortran in 1991 for simulation of luminosity detector at HERA/DESY.
The detail informationon the Bremge and physics ground can be found in sub-directory doc. 
Recently BremgeC is used for development of the luminosity detector (LumiDirectPCal) part of the EIC project detector.

*** 1. Overview. ***

........ to be filled

** 2. Compilation. **

The ** HepMC3, GSL ** and **ROOT** libraries are required to compile the code. 

Assuming you have those, execute bash script make_bremge.sh provided in the scripts directory:

 bash: ./scripts/make_bremge.sh
 
 
This, if done successfully, creates executable ** bremgenC ** , and puts it into created bin directory. 
The script takes one command line argument 0/1/2/3. 
For general purpose flag 0 or equivalently none can be provided, flags 1,2 and 3 are for special use.  
 
** 3. Runnig bremgeC. **

The executable bremgenC takes up to four command line arguments:

	 Usage: bremgenC :  <number_of_events> <HepMC3_output_file> <photon_theta_max> <generator_seed>
			         -  <number_of_events>       obligatory 
			         -  <HepMC3_output_file>     obligatory 
			         -  <photon_theta_max>       optional photon theta_max,
				                                 if ommitted default = 0.001 [radian] 
			         -  <generator_seed>         optional may have three values (default is -1): 
			                                     if = 0 -> Start ranlux with the default seed 
			                                     if > 0 -> Seed ranlux with default user given seed 
			                                     if < 0 -> or ommited use seed saved in previous run
						                                    in the file LastRunInfo.txt

First two, number of events and name of the HepMC3 output file, are obligatory. 
The default values of parameters needed to instantiate class Brem and can not be changed 
from the command line are set in the src/Main.C code, these are following:

  - - electron beam energy        (e_beam =  10. GeV)
  - - proton beam energy          (p_beam = 275. GeV)      
  - - photon minimum energy       (e_photon_min = 0.1 [GeV])  
  - - photon maximum energy       (e_photon_max = 9.9 [GeV])
  - - photon maximum polar angle  (tgmax = 0.001 [radians])

These values, except the latter one, can be changed only by editing Main.C, and recompiling code. 
On the output, after succesfull run (aproximately 17 usec/event), are two files:

  - - runYYDDDrr_your_file_name.hepmc3            ------ Bremge output converted to HepMC3 format
  - - runYYDDDrrr_your_file_name.root             ------ contains BremTree and some diagnostic histograms

To avoid accidental overwriting files generated earlier, your_file_name is prepended with the strin "runYYDDDrrr_", 
here YY last 2 digits of the year, DDD the day of the year and rrr is run number done that day.

