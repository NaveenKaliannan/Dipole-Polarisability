#ifndef BIREFRIENGENCE_H
#define BIREFRIENGENCE_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"
#include "../include/io.h"

using namespace std;


void Print_birefriengenceI_purewater(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename);
void Print_birefriengenceP_purewater(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename);
void Print_birefriengenceT_purewater(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename);

void Print_Cosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);
void Print_Cosine2(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);
void Print_KEtrans(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename) ;
void Print_KErot(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename); 
#endif
