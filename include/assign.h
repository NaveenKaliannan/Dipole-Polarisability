#ifndef ASSIGN_H
#define ASSIGN_H

#include<vector>
#include<fstream>
#include "../include/io.h"
#include "../include/numerics.h" 



using namespace std;

void AssignAtomicMass(vector<Atom> &r, uint nsteps, uint natoms);
void Assigncoordinationforpurewater(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt); 
void Assigngammaforpurewater(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt) ;

void parameters(Molecular &mol);
void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol);
void TransformAtomictoAtomic(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol);

#endif
