#ifndef DIPOL_H
#define DIPOL_H

#include<vector>
#include<fstream>
#include "../include/io.h"
#include "../include/numerics.h" 

using namespace std;

void parameters(Molecular &mol);
void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol);


#endif
