#ifndef VELOCITY_H
#define VELOCITY_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"
#include "../include/io.h"

using namespace std;


#define vAperfmstoamu   0.0460908 // converts velocity [angstrom/fms] to atomic units  
#define jtohartree  2.2937104

void computeatomicvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);
void computecomposition(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);
void computecomvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);
void computeangularvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);
void computepositionmolframe(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);
void Print_trans_rot_transrot_cf(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);


float angle_btwn_3points(const vector<Atom> &r, uint i, uint j1, uint j2, const vector<float> & L ) ; 
void PrintKEnCosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename) ; 
void classifywater(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt) ; 
void PrintOOdistance(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);
void Printwaterioncoordination(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename) ;
void Printtrans_rot_ccfn(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename); 
void PrintOOdistance1(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename);
#endif
