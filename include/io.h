#ifndef IO_H
#define IO_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"

using namespace std;


#define PI 3.14159265
#define jtohartree  2.2937104
#define amu  1.66053906
#define amutoangstrom  0.52917720859
#define vAperfmstoamu   0.0460908 // converts velocity [angstrom/fms] to atomic units  
#define Unitvectortodebye  1.93  // easily comparable to cp2k results. The dipole moment slightly changes due to the variation in H-bonding environments.
#define polfieldtodebye   17.14131  // E-30 (angstrom3 to m3) * 5.142 E11 (field a.u. to V/m or J C-1 m-1) * 8.85418 E-12 * 4 * 3.14 (permentivty of free space C2 m-1 J-1) / 3.336E-30 (coversion from C.m to debye)
#define polpointchargetodebye 4.8021  // E-30 (angstrom3 to m3) * 1.602E-19 (elementary charge C) / E-20 (angstrom2 to m2) * 3.336E-30 (coversion from C.m to debye)
#define amu3toangstrom3 0.148035889  // pow(0.52917720859,3)
#define angstrom3toamu3  6.75511  // pow(0.52917720859,-3)
#define pointchargedistancetodebye 4.8021  // E-10 (angstrom to m) * 1.602E-19 (elementary charge C) / 3.336E-30 (coversion from C.m to debye)  = 16.02E-30 / 3.336 E-30
#define amutodebye  	 2.541161873  //  0.52917720859 * E-10 (angstrom to m) * 1.602E-19 (elementary charge C) / 3.336E-30 (coversion from C.m to debye)  = 16.02E-30 / 3.336 E-30

#define amufieldtoangstrom  3.57637  //  5.14220642 E11 (field a.u. to V/m or J C-1 m-1) * 8.85418 E-12 * 4 * 3.14 (permentivty of free space C2 m-1 J-1) = 57.2935C m-2
                                              //  57.2935 E-20 C angstrom-2 / 16.02E-20 

struct Matrix
{
  double xx,xy,xz,
        yy,yx,yz,
        zz,zx,zy;  
};


struct Vector
{
  double x,y,z;
};

struct Vector_int
{
  uint x,y,z;
};



struct Atom
{
  string symbol;
  double x,y,z;
  double vx,vy,vz;
  uint index;
  string segname;
  uint resid;
  string resname;
  string atomname;
  string atomtype;
  float charge;
  float atomicmass;
};


struct Molecular
{
  string MOL;
  double x,y,z;  // Center of mass
  float q;      // Charge
  float m;      // mass
  float sl;     //screening length thole
  float vdwr;     //vanderwall radii

  Vector PD;    // dipole (permanent) 
  Matrix PPol;  // Polarisability matrix  (permanent)

  Vector ID;    // dipole (induced) 
  Matrix IPol;  // Polarisability matrix  (induced)

  Vector TD;    // dipole (total) 
  Matrix TPol;  // Polarisability matrix  (total)

};

void AssignAtomicMass(vector<Atom> &r, uint nsteps, uint natoms);
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L);
void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L);
void readmullikencharges(vector<Atom> &r, uint nsteps, uint natoms, string filename);
void readExternalfield(vector<Vector> &E, uint nsteps, string fieldfilename);
void computevelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt);
void readpsf(vector<Atom> &r, uint nsteps,  uint natoms, string psffilename);
void Print(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE);
void Print_tinker(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE);

void init_matrix_zero(Matrix &dummyM);
void init_vector_zero(Vector & dummyv);
void init_Matrix_zero(vector<Matrix> & Tij, uint nsteps, uint nmol);
void init_Vector_zero(vector<Vector> & dipole, uint nsteps, uint nmol);

void Mat_vec(const Matrix & A, const Vector & b, Vector & dummy);
void Mat_Mat(const Matrix & A, const Matrix & B, Matrix & C) ;

#endif
