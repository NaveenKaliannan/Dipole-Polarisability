#ifndef PBC_H
#define PBC_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"
#include "../include/io.h"

using namespace std;

void replica(const vector<float> & L, uint ncell, vector<float> & PB_L, vector<Vector_int> & imageno) ;
void dist(vector<Molecular> &mol, uint idi, uint idj, const vector<float> & L,  vector<float> & PB_L, vector<Vector_int> & imageno, uint index,  double &x, double & y, double &z );
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L) ; 

#endif
