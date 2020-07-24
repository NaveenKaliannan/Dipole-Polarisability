#ifndef DIPOL_H
#define DIPOL_H

#include<vector>
#include<fstream>
#include "../include/io.h"
#include "../include/numerics.h" 


#define sl 2.4380 // screeing length thole 

using namespace std;


void parameters(Molecular &mol);
void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol);
void Induced_dipole_pol(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter,vector<Vector> &E);
void dipoletensorfield(Matrix &Tij, float rij, float x, float y, float z) ;
void Tij_dipole(const Matrix & Tij, const Vector & dipole, Vector & dummy) ;
void Tij_Pol(const Matrix & Tij, const Matrix & Pol, Matrix & dummy);
void copydata(const vector<Molecular> &mol, uint nsteps, uint nmol,  vector<Matrix> & TP,  vector<Vector> & TD );
void Pol_Efield(const Matrix & A, const Vector & b, Vector & dummy);
void Pol_Ifield(const Matrix & A, const Vector & b, Vector & dummy) ;
#endif
