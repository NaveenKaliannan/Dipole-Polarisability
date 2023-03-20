#ifndef DIPOL_H
#define DIPOL_H

#include<vector>
#include<fstream>
#include "../include/io.h"
#include "../include/numerics.h" 



using namespace std;

// dipole
void FieldduetoPermanentMultipoles(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field);
void Fieldduetodipole(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field);
void FieldduetoExternalField(vector<Molecular> &mol, uint t, uint nmol, const vector<Vector> &E,  vector <Vector> & Field);
void Induced_dipole(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E); 

//polarizability
void TensorduetoExternalField(vector<Molecular> &mol, uint t, uint nmol, const vector<Vector> &E,  vector <Vector> & Field) ;
void TensorduetoPermanentMultipoles(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field) ; 
void Induced_polarisability(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E) ; 
void Induced_polarisabilityduetohyperpolarizability(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E) ;
void Induced_polarisabilityduetofirsthyperpolarizability(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E);
void Induced_polarisabilityduetosecondhyperpolarizability(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E);

void TensorduetoPolarisability(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Matrix> & C) ;

void PrintOpticalBirefringence(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename) ;
void Print_IR_RAMAN_Spectra(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename);
void Print_power_Spectra(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);

#endif
