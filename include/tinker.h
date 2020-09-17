#ifndef TINKER_H
#define TINKER_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"

using namespace std;

void Print_tinker(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE);

#endif
