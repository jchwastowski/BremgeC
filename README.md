BremgeC is a generator of bremsstrahlung events in ep (electron-proton) and eA (electron-Ion) scattering, 
ported to C++.
It is based on the original Bremge written in Fortran in 1991 to simulate the luminosity detector at ZEUS/HERA/DESY.
Detailed information on the Bremge and physics ground can be found in the sub-directory doc. 
Presently, BremgeC is used to develop the luminosity detector (LumiDirectPCal) part of the ePIC/EIC detector.

# 1. Overview

Bremsstrahlung events in the electron-proton collisions are generated according to the ultra-relativistic differential cross-section
in Born and small-angle approximation [1]. Higher-order diagrams and inelastic contributions are neglected since they are
well below 1% [2, 3]. Also, the energy transfer to the target of the order of m_e/M_target is neglected leading to the relation

		E_{\gamma} + E'_{e} = E_{e,beam}
In the case of the electron-nucleus bremsstrahlung the screening of the nucleus field by the atomic electrons is taken into account 
by applying the Thomas-Fermi-Moliere form factor [4].
The beam-size effect [5] is related to interactions with the large impact parameter value comparable to the transverse dimensions
of the beams. It is simulated by suppressing events with small momentum transfer.

1. V. B. Berestetskii, E. M. Lifshitz and L. P. Pitaevskii, Quantum Electrodynamics, Pergamon Press, 1971.
2. M. van der Horst, Theoretical Calculation for Electron-Proton Scattering, Ph. D. thesis, U. Amsterdam, 1990;
   L.J.G. Gaemers and M. van der Hosrt, Nucl. Phys. B316 (1989) 269;
   M. van der Horst, Nucl. Phys. B347 (1990) 149 and Phys. Lett. B244 (1990) 07.
3. A. A. Akhundov et al., Electron-State Bremsstranlung Process ep->γep(X) at HERA, DESY-90-130.
4. Y-S. Tsai, Pair Production and Bremsstrahlung of Charged Leptons, Rev. of Modern Phys. 46 (1974) 815.
5. G.L. Kotkin et al., Influence of the Transverse Beam Sizes on the ep→epγ Cross-section at the HERA 
   and a Future {CERN} Electron-Proton Collider, Z. Phys. C39 (1988) 61. 




# 2. Compilation.

The **HepMC3**, **GSL**, and **ROOT** libraries are required to compile the code. 

Assuming you have those, execute the bash script make_bremge.sh provided in the scripts directory:

 bash: `./scripts/make_bremge.sh`
 
 
If done successfully, this creates executable ** bremgenC ** and stores it in the created bin directory. 
The script takes one command line argument 0/1/2/3. 
For general purposes flag 0 or equivalently none can be provided, flags 1,2 and 3 are for special use.  

======
# 3. Running BremgeC

The executable bremgenC takes up to four command line arguments:

	 Usage: bremgenC :  <number_of_events> <HepMC3_output_file> <photon_theta_max> <generator_seed>
			         -  <number_of_events>       obligatory 
			         -  <HepMC3_output_file>     obligatory 
			         -  <photon_theta_max>       optional photon theta_max,
				                                 if omitted default = 0.001 [radian] 
			         -  <generator_seed>         optional may have three values (default is -1): 
			                                     if = 0 -> Start ranlux with the default seed 
			                                     if > 0 -> Seed ranlux with default user given seed 
			                                     if < 0 -> or omitted use the seed saved in the previous run
						                                    in the file LastRunInfo.txt

The first two, the number of events and the name of the HepMC3 output file are obligatory. 
The default values of parameters needed to instantiate class Brem and can not be changed 
from the command line are set in the src/Main.C code, these are the following:

  - - electron beam energy        (e_beam =  10. GeV)
  - - proton beam energy          (p_beam = 275. GeV)      
  - - photon minimum energy       (e_photon_min = 0.1 [GeV])  
  - - photon maximum energy       (e_photon_max = 9.9 [GeV])
  - - photon maximum polar angle  (tgmax = 0.001 [radians])

These values, except the latter one, can be changed only by editing Main.C, and recompiling code. 
After a successful run (approximately 17 usec/event) two output files can be found:

  - - runYYDDDrr_your_file_name.hepmc3            ------ Bremge output converted to HepMC3 format
  - - runYYDDDrrr_your_file_name.root             ------ contains BremTree and some diagnostic histograms

To avoid accidental overwriting files generated earlier, your_file_name is prepended with the string "runYYDDDrrr_", 
where YY denotes the last 2 digits of the year, DDD -- the day of the year and rrr is the run number done that day.

