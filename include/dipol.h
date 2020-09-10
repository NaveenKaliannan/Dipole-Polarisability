#ifndef DIPOL_H
#define DIPOL_H

#include<vector>
#include<fstream>
#include "../include/io.h"
#include "../include/numerics.h" 



using namespace std;

void FieldduetoPermanentMultipoles(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field);
void Fieldduetodipole(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field);
void FieldduetoExternalField(vector<Molecular> &mol, uint t, uint nmol, const vector<Vector> &E,  vector <Vector> & Field);
void Induced_dipole_pol(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter,vector<Vector> &E);

void parameters(Molecular &mol);
void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol);
void TransformAtomictoAtomic(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol);


void dist(vector<Molecular> &mol, uint idi, uint idj, const vector<float> & L,  vector<float> & PB_L, vector<Vector_int> & imageno, uint index,  float &x, float & y, float &z );
void replica(const vector<float> & L, uint ncell, vector<float> & PB_L, vector<Vector_int> & imageno);


#endif
