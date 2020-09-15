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

  //AssignAtomicMass(r, nsteps, natoms);
  //BringintoBox(r, nsteps, natoms, L);
  computeatomicvelocity(r, nsteps, natoms, L, dt);
  PrintKEs(r, nsteps, natoms, L, dt, argv[9]);
  Printcosine(r, nsteps, natoms, L, dt, argv[10]);

  //TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  //readExternalfield(E, nsteps, fieldfilename);
  //Induced_dipole_pol(mol, nsteps, nmol, L, 500, E);

 // TransformAtomictoAtomic(r, nsteps, natoms, L, mol, nmol);
  //TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  //readpsf(r, nsteps,  natoms, psffilename);
  //readmullikencharges(r, nsteps, natoms, "test" ) ;

  //Print(r, nsteps, natoms, L, mol, nmol, dt, "PBC-trajectory.xyz", "ATM");
  //Print(r, nsteps, natoms, L, mol, nmol, dt, "COM.xyz", "MOL");
  //Print(r, nsteps, natoms, L, mol, nmol, dt, "Permanet.data", "DIP-P");
  //Print(r, nsteps, natoms, L, mol, nmol, dt, "Induced.data", "DIP-I");

  //Print_tinker(r, nsteps, natoms, L, mol, nmol, dt, "water.xyz", "TINKER");
  //Print(r, nsteps, natoms, L, mol, nmol, dt, "PBC-trajectory.xyz", "ATM");
  //Print(r, nsteps, natoms, L, mol, nmol, dt, argv[9], "DIP-T");

  return 0;
}






