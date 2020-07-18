#ifndef IO_H
#define IO_H

#include<vector>
#include<fstream>

using namespace std;


#define PI 3.14159265
#define jtohartree  2.2937104
#define amu  1.66053906
#define Unitvectortodebye   2.08  // easily comparable to cp2k results
#define polfieldtodebye   17.14131  // E-30 (angstrom3 to m3) * 5.142 E11 (field a.u. to V/m or J C-1 m-1) * 8.85418 E-12 * 4 * 3.14 (permentivty of free space C2 m-1 J-1) / 3.336E-30 (coversion from C.m to debye)
#define polpointchargetodebye 4.8021  // E-30 (angstrom3 to m3) * 1.602E-19 (elementary charge C) / E-20 (angstrom2 to m2) * 3.336E-30 (coversion from C.m to debye)


struct Atom
{
  string symbol;
  float x,y,z;
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
  float x,y,z;                    // Center of mass
  float q;                        // Charge
  float m;                        // mass

  float PD_x,PD_y,PD_z;           // dipole (permanent) 
  float Pol_xx,Pol_xy,Pol_xz,
        Pol_yy,Pol_yx,Pol_yz,
        Pol_zz,Pol_zx,Pol_zy;     // Polarisability matrix  (permanent)

};

void AssignAtomicMass(vector<Atom> &r, uint nsteps, uint natoms);
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, float L);
void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, float L);
void readpsf(vector<Atom> &r, uint nsteps,  uint natoms, string psffilename);
void Print(vector<Atom> &r, uint nsteps, uint natoms, float L, string filename);
void Printmol(vector<Molecular> &r, uint nsteps, uint nmol, float L, string filename);
void Printdipol_pol(vector<Molecular> &r, uint nsteps, uint nmol, float L, string filename);
void parameters(Molecular &mol);
void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, float L, vector<Molecular> &mol, uint nmol);

#endif
