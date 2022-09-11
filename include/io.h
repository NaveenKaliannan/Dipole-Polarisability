#ifndef IO_H
#define IO_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"

using namespace std;


#define deltat 1

#define vAperfmstoamu   0.0460908
#define PI 3.14159265
#define amu  1.66053906
#define amutoangstrom  0.52917720859
#define Unitvectortodebye  1.93  // easily comparable to cp2k results. The dipole moment slightly changes due to the variation in H-bonding environments.
#define polfieldtodebye   17.14131  // E-30 (angstrom3 to m3) * 5.142 E11 (field a.u. to V/m or J C-1 m-1) * 8.85418 E-12 * 4 * 3.14 (permentivty of free space C2 m-1 J-1) / 3.336E-30 (coversion from C.m to debye)
#define polpointchargetodebye 4.8021  // E-30 (angstrom3 to m3) * 1.602E-19 (elementary charge C) / E-20 (angstrom2 to m2) * 3.336E-30 (coversion from C.m to debye)
#define amu3toangstrom3 0.148035889  // pow(0.52917720859,3)
#define amu5toangstrom5 0.041495944  // pow(0.52917720859,5) 1 au = 0.320662E-52 C3m3 J-2 = (0.320662E-52 C3m3 J-2 * 1.602E-19 C * E50 (conversion from m5 to angstrom5) (1/(4*4*3.14*3.14*8.85418 E-12*8.85418 E-12))) angstrom**5
#define angstrom3toamu3  6.75511  // pow(0.52917720859,-3)
#define pointchargedistancetodebye 4.8021  // E-10 (angstrom to m) * 1.602E-19 (elementary charge C) / 3.336E-30 (coversion from C.m to debye)  = 16.02E-30 / 3.336 E-30
#define amutodebye  	 2.541161873  //  0.52917720859 * E-10 (angstrom to m) * 1.602E-19 (elementary charge C) / 3.336E-30 (coversion from C.m to debye)  = 16.02E-30 / 3.336 E-30
#define debyetoangstrom 0.20823

#define amufieldtoangstrom2  3.57637  //  5.14220652 E11 (field a.u. to V/m or J C-1 m-1) * 8.85418 E-12 * 4 * PI (permentivty of free space C2 m-1 J-1) = 57.2935 C m-2
                                              //  57.2935 E-20 C angstrom-2 / 16.02E-20 

#define amufieldtokcalpermol 1185.8592256  // 5.14220652 E11 (field a.u. to V/m or J C-1 m-1) *  1.602E-19 C 
                              // 5.14220652 E11 * E-3 (J to KJ) * E-10 (m-1 to angstrom-1) * 1.602E-19 * 6.023E23  KJ per angstrom per mol = 4961.635 KJ per mol per angstrom 
                              // 4961.635  KJ per mol per angstrom  * 0.239006 = 1185.8592256 kcal per mol per angstrom


#define Kfactor 332.0637 // 8.98 E9 (permetivity of free space Nm2C-2 or JmC-2) * (1.602E-19 C) * (1.602E-19 C) / (A2)  
                        //= 8.98 E9 * sqr(1.602E-19) * 6.023E23 (to per mol) * E-3 (J to KJ) * E10 (m to angstrom)  KJ/(mol. Angstrom) = 1389.552 KJ per mol per angstrom * 0.239006 =  332.0637 Kcal per mol per angstrom

struct rank3tensor
{
  double xxx,yyy,zzz,
         xxy,xxz,yyx,yyz,zzx,zzy,
         yxx,zxx,xyy,zyy,xzz,yzz,
         xyx,xzx,yxy,yzy,zxz,zyz,
         xyz,xzy,yxz,yzx,zxy,zyx;
};


struct rank4tensor
{
  double xxxx, xxyy, xxzz,
         yyxx, yyyy, yyzz,
         zzxx, zzyy, zzzz;
};

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
  string gro;
  string symbol;
  double x,y,z;
  double comx, comy, comz;
  double vx,vy,vz;
  double comvx, comvy, comvz;
  double angvx, angvy, angvz;

  uint index;
  string segname;
  uint resid;
  string resname;
  string atomname;
  string atomtype;
  float charge;
  float atomicmass;

  /*only for water*/
  uint totalhbonds;
  uint totaldonorhbonds;
  uint totalacceptorhbonds;
  float gamma_d, gamma_a;  
  uint D1, D2, A1, A2;
  float ED1, ED2, EA1, EA2;
  float CTD1, CTD2, CTA1, CTA2;

  /*classifywater*/
  string watertype;
};


struct Molecular
{
  string MOL;
  double x,y,z;  // Center of mass
  double vx,vy,vz; // center of mass velocity
  float q;      // Charge
  float m;      // mass
  double sl;     //screening length thole
  float vdwr;     //vanderwall radii

  Vector PD;    // dipole (permanent) 
  Matrix PPol;  // Polarisability matrix  (permanent)

  Vector ID;    // dipole (induced) 
  Matrix IPol;  // Polarisability matrix  (induced)

  Vector TD;    // dipole (total) 
  Matrix TPol;  // Polarisability matrix  (total)

  rank3tensor hyperpol;
  rank4tensor shyperpol;

  /*only for water*/
  uint totalhbonds;
  uint totaldonorhbonds;
  uint totalacceptorhbonds;
  float gamma_d, gamma_a;
};


void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L);
void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L);
void readExternalfield(vector<Vector> &E, uint nsteps, string fieldfilename);
void readpsf(vector<Atom> &r, uint nsteps,  uint natoms, string psffilename);
void Print(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE);
void readtrajectory_gro(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L) ; 

void init_matrix_zero(Matrix &dummyM);
void init_vector_zero(Vector & dummyv);
void init_Matrix_zero(vector<Matrix> & Tij, uint nsteps, uint nmol);
void init_Vector_zero(vector<Vector> & dipole, uint nsteps, uint nmol);

void Mat_vec(const Matrix & A, const Vector & b, Vector & dummy);
void Mat_Mat(const Matrix & A, const Matrix & B, Matrix & C) ;
void Thirdranktensor_vec(const rank3tensor & A, const Vector & b, Matrix & C) ;
void Fourthranktensor_vec(const rank4tensor & A, const Vector & b, Matrix & C);

#endif
