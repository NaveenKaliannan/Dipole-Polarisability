#ifndef RDF_H
#define RDF_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"
#include "../include/io.h"

using namespace std;


void Printrdf(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);
void population_hbonds(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);

#endif
