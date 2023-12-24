/**
 *
 * @author  Naveen Kumar Kaliannan
 *
 * @Reach me via naveenkumar5892@gmail.com
 *
 * @ Upon request, a minimal assistance will
 * @ be provided for setting up this program.
 *
 */


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>
#include <vector> 
#include <numeric>
#include <sstream>  
#include <cmath> 
#include "../include/io.h"
#include "../include/numerics.h"
#include "../include/dipol.h"
#include "../include/velocity.h"
#include "../include/tinker.h"
#include "../include/molar.h"
#include "../include/rdf.h"
#include "../include/assign.h"
#include "../include/birefriengence.h"
#include "../include/dipolparameterization.h"
#include "../include/cp2k.h"

using namespace std;


int main ( int argc, char** argv )
{
  uint natoms = 0, nsteps = 0, id = 0, nmol = 0;
  string temp, xyzfilename, psffilename, fieldfilename;
  float dt = 0;
  vector<float> L (3,0.0);

  //input arguments
  xyzfilename = argv[1];
  psffilename = argv[2];
  nsteps = atoi(argv[3]);
  dt = atof(argv[4]);
  L[0] = atof(argv[5]);
  L[1] = atof(argv[5]);
  L[2] = atof(argv[5]);
  natoms = atoi(argv[6]);
  nmol = atoi(argv[7]);
  fieldfilename = argv[8];

  vector<Atom> r(natoms*nsteps);
  vector<Molecular> mol ;
  vector<Vector> E (nsteps);


  readtrajectory_gro(r, nsteps, natoms, xyzfilename, L); 
  BringintoBox(r, nsteps, natoms, L);
   
  /*Computer Polarizability anisotropy*/
  readExternalfield(E, nsteps, fieldfilename);
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[9]);
  Induced_dipole(mol, nsteps, nmol, L, 500, E);
  Induced_polarisability(mol, nsteps, nmol, L, 500, E);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[10]);
  Induced_polarisabilityduetofirsthyperpolarizability(mol, nsteps, nmol, L, 500, E);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[11]);
  Induced_polarisabilityduetosecondhyperpolarizability(mol, nsteps, nmol, L, 500, E);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[12]);

  /*Computer velocities*/
  computeatomicvelocity(r, nsteps, natoms, L, dt);
  computecomvelocity(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computepositionmolframe(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computeangularvelocity(r, nsteps, natoms, L, dt);

  /*CCF tensor components*/
  Print_intramolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[13]);
  Print_intermolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[14]);
  transform_angular_com_velocity_into_mol_frame(r, nsteps, natoms, L, dt);
  Print_intramolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[15]);
  Print_intermolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[16]);
 
  Print_power_Spectra(r, nsteps,natoms,L, dt, argv[17]);
  Print_power_Spectra_intermolecular(r, nsteps,natoms,L, dt, argv[18]);

  PrintKEnCosine(r, nsteps, natoms, L, dt, argv[19]);

  Assigncoordinationforpurewater(r, nsteps, natoms, L, dt);
  population_hbonds(r, nsteps, natoms, L, dt, argv[20]);

/*
  ------ Reading different format trajectory---------
  readtrajectory(r, nsteps, natoms, xyzfilename, L); 
  readtrajectory_gro(r, nsteps, natoms, xyzfilename, L); 
  readtrajectory_tinker(r, nsteps, natoms, xyzfilename, L); 
  readtrajectory_tinkerhp(r, nsteps, natoms, xyzfilename, L);
  readpsf(r, nsteps,  natoms, psffilename); 
  ------ Reading different format trajectory---------
  
  ------ Writing different format trajectory---------
  Print(r, nsteps, natoms, L, mol, nmol, dt, "PBC-trajectory.xyz", "ATM");
  Print(r, nsteps, natoms, L, mol, nmol, dt, "COM.xyz", "MOL");
  Print_tinker(r, nsteps, natoms, L, mol, nmol, dt, "water.xyz", "TINKER");
  Print_tinkerhp(r, nsteps, natoms, L, mol, nmol, dt, "../water2.arc", "TINKER");
  ------ Writing different format trajectory---------

  ------ Applying PBC ---------
  BringintoBox(r, nsteps, natoms, L);
  ------ Applying PBC ---------

  ---Estimate the simulation cell size and number of water and ions for different concentrations---
  EstimateCellSizenNumberofResidue();
  ---Estimate the simulation cell size and number of water and ions for different concentrations---

  -----Assign atomic mass to each atom----------------
  AssignAtomicMass(r, nsteps, natoms);
  -----Assign atomic mass to each atom----------------

  ----Transforming atomic to molecular or atomic to atomic------
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  TransformAtomictoAtomic(r, nsteps, natoms, L, mol, nmol);
  ----Transforming atomic to molecular or atomic to atomic------

  ----- Prints radial distribution function ----
  Printrdf(r, nsteps, natoms, L, dt, argv[9]);
  ----- Prints radial distribution function ----

  ---- Computing molecular dipole [debye], molecular polarisability [angstrom3], molecular polarizability anisotropy [angstrom3] ----     
  readExternalfield(E, nsteps, fieldfilename);
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[9]);
  Induced_dipole(mol, nsteps, nmol, L, 500, E);
  Induced_polarisability(mol, nsteps, nmol, L, 500, E);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[10]);
  Induced_polarisabilityduetofirsthyperpolarizability(mol, nsteps, nmol, L, 500, E);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[11]);
  Induced_polarisabilityduetosecondhyperpolarizability(mol, nsteps, nmol, L, 500, E);
  PrintOpticalBirefringence(mol, nsteps, nmol, L, dt, argv[12]);
  ---- Computing molecular dipole [debye], molecular polarisability [angstrom3], molecular polarizability anisotropy [angstrom3] ----     

  ----- Assinging coordination number and asymmetry gamma parameter for pure water -------------
  classifywater(r, nsteps, natoms, L, dt);
  PrintOOdistance1(r, nsteps, natoms, L, dt, "OOdistance");
  PrintOOdistance(r, nsteps, natoms, L, dt, "OOdistance");
  Assigncoordinationforpurewater(r, nsteps, natoms, L, dt); 
  Assigngammaforpurewater(r, nsteps, natoms, L, dt); 
  ----- Assinging coordination number and asymmetry gamma parameter for pure water -------------

   ---- Print cosine and Kinetic energy----------------
   PrintKEnCosine(r, nsteps, natoms, L, dt, argv[13])
   ---- Print cosine and Kinetic energy----------------


   ------ Print Rotational, translational and trans-rot correlation function --------------------
  computeatomicvelocity(r, nsteps, natoms, L, dt);
  computecomvelocity(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computepositionmolframe(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computeangularvelocity(r, nsteps, natoms, L, dt);

  // Local Frame
  Print_intramolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[9]);
  Print_intermolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[10]);

  // Molecular Frame
  transform_angular_com_velocity_into_mol_frame(r, nsteps, natoms, L, dt);
  Print_intramolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[11]);
  Print_intermolecular_coupling_btwn_translation_and_rotation(r, nsteps, natoms, L, dt,  argv[12]);
   ------ Print Rotational, translational and trans-rot correlation function --------------------

   ------ Print population of HBonds --------------------
  readtrajectory(r, nsteps, natoms, xyzfilename, L);
  BringintoBox(r, nsteps, natoms, L);
  Assigncoordinationforpurewater(r, nsteps, natoms, L, dt);
  population_hbonds(r, nsteps, natoms, L, dt, argv[9]);
   ------ Print population of HBonds --------------------

   ------ Print HB strength as a function of time --------------------
  readtrajectory(r, nsteps, natoms, xyzfilename, L);
  BringintoBox(r, nsteps, natoms, L);
  Print_ALMO_data(r, nsteps, natoms, L, dt, argv[9], argv[10], argv[11]);
   ------ Print HB strength as a function of time --------------------


   ------ Print power spectra, spectra of rotational and COM  --------------------
  readtrajectory(r, nsteps, natoms, xyzfilename, L); 
  BringintoBox(r, nsteps, natoms, L);
  computeatomicvelocity(r, nsteps, natoms, L, dt);
  computecomvelocity(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computepositionmolframe(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computeangularvelocity(r, nsteps, natoms, L, dt);  
  Print_power_Spectra(r, nsteps,natoms,L, dt, "RAMAN_IR_SPECTRA.dat");	   
   ------ Print power spectra, spectra of rotational and COM  --------------------

  */

  return 0;
}
