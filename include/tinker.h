#ifndef TINKER_H
#define TINKER_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"

using namespace std;
void readtrajectory_tinker(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L);
void Print_tinker(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE);


void readtrajectory_tinkerhp(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L);
void Print_tinkerhp(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE); 
#endif
