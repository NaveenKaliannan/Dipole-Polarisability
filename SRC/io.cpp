
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>
#include <vector> 
#include <numeric>
#include <sstream>  
#include <cmath>
#include "../include/numerics.h" 
#include "../include/io.h"


using namespace std;


//Reading the xyz trajectory
void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L)
{
  string temp;
  ifstream xyzfile(xyzfilename);
  for(uint t = 0; t < nsteps; ++t )
    { 
      xyzfile >> natoms; 
      getline(xyzfile, temp); 
      getline(xyzfile, temp); 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          xyzfile >> r[id].symbol >> r[id].x  >> r[id].y  >> r[id].z;           
          //cout << r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << endl;
        }
    }
  xyzfile.close();
  xyzfile.clear();
}

void readtrajectory_gro(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L)
{
    string temp;
    ifstream xyzfile(xyzfilename);
    for(uint t = 0; t < nsteps; ++t )
      {
        getline(xyzfile, temp);
        xyzfile >> natoms;
        getline(xyzfile, temp);
        for(uint i = 0;i < natoms;++i)
          {
            uint id = natoms*t+i;
            xyzfile >> temp >> r[id].symbol >> temp >> r[id].x  >> r[id].y  >> r[id].z;
            r[id].x *= 10;
            r[id].y *= 10;
            r[id].z *= 10;
            //cout << t << "  " << i << "  " <<  r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << endl;
          }
        xyzfile >> temp >>  temp >> temp;
        getline(xyzfile, temp);
     }
  xyzfile.close();
  xyzfile.clear();
}



//Reading the psf file
void readpsf(vector<Atom> &r, uint nsteps,  uint natoms, string psffilename)
{
  string input;
  ifstream is(psffilename);
  getline(is, input);
  if(input.substr(0,3) != "PSF")
     cout <<"PSF detected a non-PSF file" << endl;
  getline(is, input); 
  getline(is, input);
  getline(is, input);
  getline(is, input);
  getline(is, input);
  uint num_atoms;
  if(num_atoms == natoms)
        cout <<"PSF and xyz file doesnt have same number of atoms" << endl;
  stringstream(input) >> num_atoms;
  for(uint i = 0; i < num_atoms; i++) 
    {
      getline(is, input);
      stringstream ss(input);
      ss >> r[i].index ; 
      ss >> r[i].segname ; 
      ss >> r[i].resid ; 
      ss >> r[i].resname ; 
      ss >> r[i].atomname ; 
      ss >> r[i].atomtype ; 
      ss >> r[i].charge ; 
      ss >> r[i].atomicmass ;
      //cout <<  r[i].index << "  " << r[i].segname <<  "  " << r[i].resid << "  " << r[i].resname << "  " <<  r[i].atomname << "  " << r[i].atomtype << "  " << r[i].charge << " " << r[i].atomicmass <<  endl;
    }
  is.close();
  is.clear();
}


//Reading the externally applied field over time
void readExternalfield(vector<Vector> &E, uint nsteps, string fieldfilename)
{
  string temp;
  ifstream file(fieldfilename); 
  for(uint t = 0; t < nsteps; ++t )
    { 
      uint id = t; 
      file >> temp >> E[id].x  >> E[id].y  >> E[id].z;       
    }
  file.close();
  file.clear();
}

void Print(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE)
{
  ofstream outfile(filename);
  for(uint t = 0; t < nsteps; t += deltat)
    { 
      float a =0, b = 0, c = 0, a1 =0, b1 = 0, c1 = 0;
      float axy =0, axz = 0, ayz = 0, ayx =0, azx = 0, azy = 0;

      if(TYPE[0] == 'A' && TYPE[1] == 'T' && TYPE[2] == 'M')
      {
        uint n = 1;
        outfile << natoms * pow(n,3) << endl ;
        outfile << "BOX Length " << L[0]*n << "  " << L[1]*n << "  " << L[2]*n << "\n";
        for(uint x = 0; x < 1; ++x)
          {
            for(uint y = 0; y < n; ++y)
              {
                for(uint z = 0; z < n; ++z)
                  {
                    for(uint i = 0;i < natoms;++i)
                      {
                        uint id = natoms*t+i; 
                        outfile << r[id].symbol <<  "  " << r[id].x + L[0] * x << "  " << r[id].y + L[1] * y << "  " << r[id].z + L[2] * z << endl;
                      }
                  }
              }
          }
      }
      else if(TYPE[0] == 'M' && TYPE[1] == 'O' && TYPE[2] == 'L')
      {
        outfile << nmol << endl ;
        outfile << "BOX Length " << L[0] << "  " << L[1] << "  " << L[2] << "\n";
        for(uint i = 0;i < nmol;++i)
          {
            uint id = nmol*t+i; 
            outfile << mol[id].MOL << " " << mol[id].x << "  " << mol[id].y << "  " << mol[id].z << endl;
          }
      }
      else if(TYPE[0] == 'D' && TYPE[1] == 'I' && TYPE[2] == 'P' && TYPE[3] == '-')
      {
        for(uint i = 0;i < nmol;++i)
          {
            uint id = nmol*t+i;
             if(TYPE[4] == 'P')
               {
                 a  += mol[id].PD.x ; b  += mol[id].PD.y ; c  += mol[id].PD.z ;
                 a1 += mol[id].PPol.xx ; b1 += mol[id].PPol.yy ; c1 += mol[id].PPol.zz ;
                 axy += mol[id].PPol.xy ; axz += mol[id].PPol.xz ; ayz += mol[id].PPol.yz ; 
                 ayx += mol[id].PPol.yx ; azx += mol[id].PPol.zx ; azy += mol[id].PPol.zy ;
               }
             else if(TYPE[4] == 'I')
               {
                 a  += mol[id].ID.x ; b  += mol[id].ID.y ; c  += mol[id].ID.z ;
                 a1 += mol[id].IPol.xx ; b1 += mol[id].IPol.yy ; c1 += mol[id].IPol.zz ;
                 axy += mol[id].IPol.xy; axz +=  mol[id].IPol.xz; ayz +=  mol[id].IPol.yz ; 
                 ayx += mol[id].IPol.yx ; azx += mol[id].IPol.zx; azy +=  mol[id].IPol.zy;
               }
             else if(TYPE[4] == 'T' )
               {
                 a += mol[id].PD.x + mol[id].ID.x ;
                 b += mol[id].PD.y + mol[id].ID.y ; 
                 c += mol[id].PD.z + mol[id].ID.z ;

                 a1 += mol[id].PPol.xx + mol[id].IPol.xx; b1 += mol[id].PPol.yy + mol[id].IPol.yy; c1 += mol[id].PPol.zz + mol[id].IPol.zz;
                 axy += mol[id].PPol.xy + mol[id].IPol.xy; axz += mol[id].PPol.xz + mol[id].IPol.xz; ayz += mol[id].PPol.yz + mol[id].IPol.yz ; 
                 ayx += mol[id].PPol.yx + mol[id].IPol.yx ; azx += mol[id].PPol.zx + mol[id].IPol.zx; azy += mol[id].PPol.zy + mol[id].IPol.zy;
               }
          }

outfile <<  t * dt  << " " <<  a   << " " <<  b  << " "  << c << "  " //<<  (a-0.5*(b+c) )/nmol
                          << " " <<  a1  << " " <<  b1 << " "  << c1 //<< "   " << (a1-0.5*(b1+c1)) /nmol
                          << " " <<  axy << " " << axz << " "  << ayz
                          << " " <<  ayx << " " << azx << " "  << azy  << endl;

      /*  outfile <<  t * dt << " " <<  a   << " " <<  b  << " "  << c << "  " << a - 0.5 * (b+c)
                           << " " <<  a1  << " " <<  b1 << " "  << c1 << "   " << a1 - 0.5 * (b1+c1)
                           << " " <<  axy << " " << axz << " "  << ayz
                           << " " <<  ayx << " " << azx << " "  << azy  << endl; */
        
       // outfile << "#Time  " << t * dt << "\n" << "## Total Mol [string] " << nmol << " Total \u03BC  [Debye] in x y z " << a << "  " << b << "  " << c << "   "  <<  " Total \u03B1 [Angstrom] xx yy zz  " << a1 << "  " <<  b1  << "  " << c1  << endl;
        //outfile << "##Mol [string]  \u03BC  [Debye] in x y z \u03B1 [Angstrom] xx yy zz  \n";
        for(uint i = 0;i < nmol;++i)
          {
            uint id = nmol*t+i;
             if(TYPE[4] == 'P')
               {
                 outfile << id <<  "  " <<  mol[id].MOL <<  "  " << mol[id].PD.x << "  " << mol[id].PD.y << "  " << mol[id].PD.z << "  "<< mol[id].PPol.xx << "  " << mol[id].PPol.yy << "  " << mol[id].PPol.zz << endl;
outfile << mol[id].x << "  " << mol[id].y << endl;
               }
             else if(TYPE[4] == 'I')
               {
                 //outfile << mol[id].MOL <<  "  " <<  mol[id].ID.x << "  " <<  mol[id].ID.y << "  " <<  mol[id].ID.z << "  "<<  mol[id].IPol.xx << "  " <<  mol[id].IPol.yy << "  " <<  mol[id].IPol.zz << endl;
               }
             else if(TYPE[4] == 'T')
               {
             //    outfile << mol[id].MOL <<  "  " << mol[id].PD.x << "  " << mol[id].PD.y << "  " << mol[id].PD.z << "  "<< mol[id].PPol.xx << "  " << mol[id].PPol.yy << "  " << mol[id].PPol.zz << endl;
            //     outfile << mol[id].MOL <<  "  " << mol[id].PD.x + mol[id].ID.x << "  " << mol[id].PD.y + mol[id].ID.y << "  " << mol[id].PD.z + mol[id].ID.z << "  "<< mol[id].PPol.xx + mol[id].IPol.xx << "  " << mol[id].PPol.yy + mol[id].IPol.yy << "  " << mol[id].PPol.zz + mol[id].IPol.zz << endl;
               }
          }
      }
    }
  outfile.close();
  outfile.clear();
}


void init_matrix_zero(Matrix &dummyM)
{
  dummyM.xx= 0;dummyM.yy= 0;dummyM.zz= 0;
  dummyM.xy= 0;dummyM.xz= 0;dummyM.yz= 0;
  dummyM.yx= 0;dummyM.zx= 0;dummyM.zy= 0;
}


void init_vector_zero(Vector & dummyv)
{
  dummyv.x = 0; dummyv.y = 0; dummyv.z = 0; 
}

void init_Matrix_zero(vector<Matrix> & Tij, uint nsteps, uint nmol)
{
  for(uint t = 0; t < nsteps; ++t )
    {        
      for(uint i = 0;i < nmol;++i)
        {
          uint id = nmol*t+i; 
          Tij[id].xx= 0;Tij[id].yy= 0;Tij[id].zz= 0;
          Tij[id].xy= 0;Tij[id].xz= 0;Tij[id].yz= 0;
          Tij[id].yx= 0;Tij[id].zx= 0;Tij[id].zy= 0;
        }
    }
}


void init_Vector_zero(vector<Vector> & dipole, uint nsteps, uint nmol)
{
  for(uint t = 0; t < nsteps; ++t )
    {        
      for(uint i = 0;i < nmol;++i)
        {
          uint id = nmol*t+i; 
          dipole[id].x = 0 ;
          dipole[id].y = 0 ;
          dipole[id].z = 0 ;
        }
    }
}


void Mat_vec(const Matrix & A, const Vector & b, Vector & dummy)
{  
  dummy.x = A.xx * b.x + A.xy * b.y + A.xz * b.z;
  dummy.y = A.yx * b.x + A.yy * b.y + A.yz * b.z;
  dummy.z = A.zx * b.x + A.zy * b.y + A.zz * b.z;
}


void Mat_Mat(const Matrix & A, const Matrix & B, Matrix & C)
{
  C.xx = A.xx * B.xx + A.xy * B.yx + A.xz * B.zx;
  C.xy = A.xx * B.xy + A.xy * B.yy + A.xz * B.zy;
  C.xz = A.xx * B.xz + A.xy * B.yz + A.xz * B.zz;

  C.yx = A.yx * B.xx + A.yy * B.yx + A.yz * B.zx;
  C.yy = A.yx * B.xy + A.yy * B.yy + A.yz * B.zy;
  C.yz = A.yx * B.xz + A.yy * B.yz + A.yz * B.zz;

  C.zx = A.zx * B.xx + A.zy * B.yx + A.zz * B.zx;
  C.zy = A.zx * B.xy + A.zy * B.yy + A.zz * B.zy;
  C.zz = A.zx * B.xz + A.zy * B.yz + A.zz * B.zz;
}

//https://www.somewhereville.com/2006/04/
void Thirdranktensor_vec(const rank3tensor & A, const Vector & b, Matrix & C)
{
  C.xx += A.xxx * b.x + A.xxy * b.y + A.xxz * b.z;
  C.xy += A.xyx * b.x + A.xyy * b.y + A.xyz * b.z;
  C.xz += A.xzx * b.x + A.xzy * b.y + A.xzz * b.z;

  C.yx += A.yxx * b.x + A.yxy * b.y + A.yxz * b.z;
  C.yy += A.yyx * b.x + A.yyy * b.y + A.yyz * b.z;
  C.yz += A.yzx * b.x + A.yzy * b.y + A.yzz * b.z;

  C.zx += A.zxx * b.x + A.zxy * b.y + A.zxz * b.z;
  C.zy += A.zyx * b.x + A.zyy * b.y + A.zyz * b.z;
  C.zz += A.zzx * b.x + A.zzy * b.y + A.zzz * b.z;
}


//https://csmbrannon.net/2015/08/10/fourth-order-tensor-tutorial-excerpts-voigt-and-mandel-representations-as-well-as-isotropy-topics/
void Fourthranktensor_vec(const rank4tensor & A, const Vector & b, Matrix & C)
{
  double Exx, Eyy, Ezz;
  Exx = 0.5 * b.x * b.x ;
  Eyy = 0.5 * b.y * b.y ;
  Ezz = 0.5 * b.z * b.z ;

  C.xx += A.xxxx * Exx + A.xxyy * Eyy + A.xxzz * Ezz;
  C.yy += A.yyxx * Exx + A.yyyy * Eyy + A.yyzz * Ezz;
  C.zz += A.zzxx * Exx + A.zzyy * Eyy + A.zzzz * Ezz; 

}


