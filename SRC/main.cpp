/**
 *
 * @author  Naveen Kumar Kaliannan
 * @
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

  readtrajectory(r, nsteps, natoms, xyzfilename, L);
  readExternalfield(E, nsteps, fieldfilename);
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  Induced_dipole(mol, nsteps, nmol, L, 500, E);
  Induced_polarisability(mol, nsteps, nmol, L, 500, E);
  //Induced_polarisabilityduehyperpolarizability(mol, nsteps, nmol, L, 500, E);
  Print(r, nsteps, natoms, L, mol, nmol, dt, "Total.data", "DIP-T");

  /*printing birefriengence total, permanent and induced components*/
//  Print_birefriengenceT_purewater(mol, nsteps, nmol, L, dt, argv[9]);
//  Print_birefriengenceP_purewater(mol, nsteps, nmol, L, dt, argv[10]);
//  Print_birefriengenceI_purewater(mol, nsteps, nmol, L, dt, argv[11]);

  /*printing cosine and decomposed translation and roational KE*/
//  Print_KErot(r, nsteps, natoms, L, dt, argv[12]);
//  Print_KEtrans(r, nsteps, natoms, L, dt, argv[13]);
//  Print_Cosine2(r, nsteps, natoms, L, dt, argv[14]);
//  Print_Cosine(r, nsteps, natoms, L, dt, argv[15]);


  /*printing cosine and decomposed translation and roational KE*/
  //Print_KErot(r, nsteps, natoms, L, dt, argv[12]);
  //Print_KEtrans(r, nsteps, natoms, L, dt, argv[13]);
  //Print_Cosine2(r, nsteps, natoms, L, dt, argv[14]);
  //Print_Cosine(r, nsteps, natoms, L, dt, argv[15]);

  /* assinging coordination number and asymmetry gamma parameter for pure water
  Assigncoordinationforpurewater(r, nsteps, natoms, L, dt); 
  Assigngammaforpurewater(r, nsteps, natoms, L, dt); 
  */
  
  /*Reading Trajectories
  readtrajectory(r, nsteps, natoms, xyzfilename, L); 
  readtrajectory_tinker(r, nsteps, natoms, xyzfilename, L); 
  readtrajectory_tinkerhp(r, nsteps, natoms, xyzfilename, L);
  readpsf(r, nsteps,  natoms, psffilename); 
  */

  /*Applying PBC
  BringintoBox(r, nsteps, natoms, L);
  */

  /*Assign atomic mass to each atom
  AssignAtomicMass(r, nsteps, natoms);
  */

  /*Transforming atomic to molecular or atomic to atomic
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  TransformAtomictoAtomic(r, nsteps, natoms, L, mol, nmol);
  */

  /*Computing molecular dipole [debye] and polarisability [angstrom3]
  readExternalfield(E, nsteps, fieldfilename);
  Induced_dipole(mol, nsteps, nmol, L, 500, E);
  Induced_polarisability(mol, nsteps, nmol, L, 500, E);
  Print(r, nsteps, natoms, L, mol, nmol, dt, "Permanet.data", "DIP-P");
  Print(r, nsteps, natoms, L, mol, nmol, dt, "Induced.data", "DIP-I");
  Print(r, nsteps, natoms, L, mol, nmol, dt, "Total.data", "DIP-T");


  Print_birefriengenceT_purewater(mol, nsteps, nmol, L, dt, argv[10]);
  Print_birefriengenceP_purewater(mol, nsteps, nmol, L, dt, argv[10]);
  Print_birefriengenceI_purewater(mol, nsteps, nmol, L, dt, argv[10]);
  */

  /*Printing cosine, Kinetic energy, rdf
  Printrdf(r, nsteps, natoms, L, dt, argv[9]);
  PrintKEnCosine(r, nsteps, natoms, L, dt, argv[9]);

  Print_KErot(r, nsteps, natoms, L, dt, argv[9]);
  Print_KEtrans(r, nsteps, natoms, L, dt, argv[9]);
  Print_Cosine2(r, nsteps, natoms, L, dt, argv[9]);
  Print_Cosine(r, nsteps, natoms, L, dt, argv[9]);
  */


  /*Printing trajectories
  Print(r, nsteps, natoms, L, mol, nmol, dt, "PBC-trajectory.xyz", "ATM");
  Print(r, nsteps, natoms, L, mol, nmol, dt, "COM.xyz", "MOL");
  Print_tinker(r, nsteps, natoms, L, mol, nmol, dt, "water.xyz", "TINKER");
  Print_tinkerhp(r, nsteps, natoms, L, mol, nmol, dt, "../water2.arc", "TINKER");
  */

  /*Estimate the cell size and number of water and ions
  EstimateCellSizenNumberofResidue();
  */

  return 0;
}






