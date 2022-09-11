#ifndef CP2K_H
#define CP2K_H

#include<vector>
#include<fstream>
#include "../include/io.h"

using namespace std;

void readmullikencharges(vector<Atom> &r, uint nsteps, uint natoms, string filename) ; 
void Print_ALMO_data(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string donorfile, string acceptorfile, string filename);

#endif
