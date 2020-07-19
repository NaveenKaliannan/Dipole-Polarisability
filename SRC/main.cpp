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


using namespace std;

#define PI 3.14159265



int main ( int argc, char** argv )
{
  uint natoms = 0, nsteps = 0, id = 0, nmol = 0;
  string temp, xyzfilename, psffilename;
  float dt = 0, L = 0;

  //input arguments
  xyzfilename = argv[1];
  psffilename = argv[2];
  nsteps = atoi(argv[3]);
  dt = atof(argv[4]);
  L = atof(argv[5]);
  natoms = atoi(argv[6]);
  nmol = atoi(argv[7]);
  vector<Atom> r(natoms*nsteps);
  vector<Molecular> mol ;

  readtrajectory(r, nsteps, natoms, xyzfilename, L);
  //readpsf(r, nsteps,  natoms, psffilename);
  Print(r, nsteps, natoms, L, mol, nmol,  "new-traj1.xyz", "ATM");
  TransformAtomictoMolecular(r, nsteps, natoms, L, mol, nmol);
  Print(r, nsteps, natoms, L, mol, nmol,  "new-traj.xyz", "MOL");
  Print(r, nsteps, natoms, L, mol, nmol, "new-traj.xyz", "DIP");

  return 0;
}






