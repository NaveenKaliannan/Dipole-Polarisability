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

using namespace std;

#define PI 3.14159265


int main ( int argc, char** argv )
{
  uint natoms = 0, nsteps = 0, id = 0, nmol = 0;
  string temp, xyzfilename, psffilename;
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
  vector<Atom> r(natoms*nsteps);
  vector<Molecular> mol ;

  readtrajectory(r, nsteps, natoms, xyzfilename, L);
  AssignAtomicMass(r, nsteps, natoms);
  BringintoBox(r, nsteps, natoms, L);
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  //computevelocity(r, nsteps, natoms, L, dt);
  //readpsf(r, nsteps,  natoms, psffilename);
  Print(r, nsteps, natoms, L, mol, nmol,  "new-traj1.xyz", "ATM");
  Print(r, nsteps, natoms, L, mol, nmol,  "new-traj2.xyz", "MOL");
  Print(r, nsteps, natoms, L, mol, nmol, "new-traj3.xyz", "DIP");

  return 0;
}






