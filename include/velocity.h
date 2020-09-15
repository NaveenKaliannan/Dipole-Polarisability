#ifndef VELOCITY_H
#define VELOCITY_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"
#include "../include/io.h"

using namespace std;


#define vAperfmstoamu   0.0460908 // converts velocity [angstrom/fms] to atomic units  
#define jtohartree  2.2937104

void computeatomicvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);

void PrintKEs(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename) ; 

void Printcosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);

#endif
