
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
#include "../include/dipol.h"
#include "../include/pbc.h"
#include "../include/assign.h"

#define ncell  1
#define rcut  7



using namespace std;

void Induced_dipole(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E)
{
  uint idi = 0, idj = 0;
  double a = 0.0, sum = 0.0, b = 0.0, eps = 0.0;
  vector<Vector> Field (nmol);
  vector<Vector> rsd (nmol), zrsd (nmol), conj (nmol), vec(nmol);

  for(uint t = 0; t < nsteps;t += deltat )
    {
      FieldduetoPermanentMultipoles(mol, t, nmol, L, Field);
      FieldduetoExternalField(mol, t, nmol, E,  Field);
      for(uint i = 0;i < nmol;++i)
        {
          idi = nmol*t+i;            
          Mat_vec(mol[idi].PPol, Field[i], mol[idi].ID);
        }

      init_Vector_zero(rsd, 1, nmol);
      init_Vector_zero(zrsd, 1, nmol);
      init_Vector_zero(conj, 1, nmol);
      init_Vector_zero(vec, 1, nmol);
      Fieldduetodipole(mol, t, nmol, L, Field);
      for(uint i = 0;i < nmol;++i)
        {
          idi = nmol*t+i;  
          rsd[i] =  Field[i] ;
          Mat_vec(mol[idi].PPol, rsd[i], zrsd[i]);
          conj[i] = zrsd[i];
        }

      // Induced dipole is computed via conjugate gradient
      // (CG) approach, exactly as in  Tinker Molecular
      // Modeling software package
      for(uint iter = 0;iter < niter;++iter)
        {
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              vec[i] = mol[idi].ID;
              mol[idi].ID = conj[i]; 
            }
          Fieldduetodipole(mol, t, nmol, L, Field);
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              mol[idi].ID = vec[i];
              vec[i].x = (conj[i].x/mol[idi].PPol.xx) - Field[i].x;
              vec[i].y = (conj[i].y/mol[idi].PPol.yy) - Field[i].y;
              vec[i].z = (conj[i].z/mol[idi].PPol.zz) - Field[i].z;
            }

          a = 0.0 ;
          sum = 0.0;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              a   = a   + conj[i].x * vec[i].x + conj[i].y * vec[i].y + conj[i].z * vec[i].z ;
              sum = sum + rsd[i].x * zrsd[i].x + rsd[i].y * zrsd[i].y + rsd[i].z * zrsd[i].z ;
            }
          a = (sum/a) ;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              mol[idi].ID.x = mol[idi].ID.x + a * conj[i].x ; 
              mol[idi].ID.y = mol[idi].ID.y + a * conj[i].y ; 
              mol[idi].ID.z = mol[idi].ID.z + a * conj[i].z ;

              rsd[i].x = rsd[i].x - a * vec[i].x;
              rsd[i].y = rsd[i].y - a * vec[i].y;
              rsd[i].z = rsd[i].z - a * vec[i].z;
            }

          b = 0.0;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              Mat_vec(mol[idi].PPol, rsd[i], zrsd[i]);
              b = b + rsd[i].x * zrsd[i].x + rsd[i].y * zrsd[i].y + rsd[i].z * zrsd[i].z ;
            }
          b = (b/sum) ; 
          eps = 0.0;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              conj[i].x = zrsd[i].x + b * conj[i].x ;
              conj[i].y = zrsd[i].y + b * conj[i].y ;
              conj[i].z = zrsd[i].z + b * conj[i].z ;
              eps = eps + rsd[i].x * rsd[i].x + rsd[i].y * rsd[i].y + rsd[i].z * rsd[i].z; 
            } 
          if (eps < 0.000001)  { cout << t <<  "  " << iter << endl; iter = niter;}   
          if (iter == niter - 1)  { cout << "Frame " << t <<  " is not converged (dipole) " << endl;  }                                 
        }
    }
      
}

void Induced_polarisability(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E)
{
  uint idi = 0, idj = 0;
  double eps = 0.0, b = 0.0, a = 0.0;
  vector<Matrix> tensor (nmol), dummy_ipolar (nmol);
  vector<Vector> Field (nmol);

  for(uint t = 0; t < nsteps;t += deltat )
    {
      // computing the induced polarisability via self-consistent field (scf)
      for(uint iter = 0;iter < niter;++iter)
        {
          TensorduetoPolarisability(mol, t, nmol, L, tensor);
          a = 0.0;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i;            
              a += mol[idi].PPol.xx + mol[idi].IPol.xx ; 
              a += mol[idi].PPol.yy + mol[idi].IPol.yy ; 
              a += mol[idi].PPol.zz + mol[idi].IPol.zz ; 
            }            
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i;  
              Mat_Mat(mol[idi].PPol, tensor[i], mol[idi].IPol);
            }
          b = 0.0;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i;            
              b += mol[idi].PPol.xx + mol[idi].IPol.xx ;   
              b += mol[idi].PPol.yy + mol[idi].IPol.yy ; 
              b += mol[idi].PPol.zz + mol[idi].IPol.zz ;  
            }

          eps = a - b  ; 
          if (abs(eps) < 0.000001)  { cout << t <<  "  " << iter << endl; iter = niter;}   
          if (iter == niter - 1)  { cout << "Frame " << t <<  " is not converged (polarisability) " << endl;  }                                 
        }
    }      
}

void Induced_polarisabilityduehyperpolarizability(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E)
{
  uint idi = 0, idj = 0;
  double eps = 0.0, b = 0.0, a = 0.0;
  vector<Vector> Field (nmol);

  for(uint t = 0; t < nsteps;t += deltat )
    {
      TensorduetoPermanentMultipoles(mol, t, nmol, L, Field);
      TensorduetoExternalField(mol, t, nmol, E,  Field);
      for(uint i = 0;i < nmol;++i)
        {
          idi = nmol*t+i;  
          Thirdranktensor_vec(mol[idi].hyperpol, Field[i], mol[idi].IPol );
        }
    }      
}



// Tensor change due to externally applied field
void TensorduetoExternalField(vector<Molecular> &mol, uint t, uint nmol, const vector<Vector> &E,  vector <Vector> & Field)
{
  uint idi = 0;
  for(uint i = 0;i < nmol;++i)
    {
      idi = nmol*t+i;  
      Field[i].x += amufieldtoangstrom2 * E[t].x ;
      Field[i].y += amufieldtoangstrom2 * E[t].y ;
      Field[i].z += amufieldtoangstrom2 * E[t].z ;
    }
}


// Tensor change due to permanent multipoles [charge, dipole, quadrupole]
void TensorduetoPermanentMultipoles(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field)
{
  uint idi = 0, idj = 0;
  double x = 0, y = 0, z = 0, rij = 0;
  double u = 0.0, ar = 0.0, st1 = 0.0, st2 = 0.0, r3 = 0.0, r5 = 0.0, pc = 0.0, pd = 0.0;
  float mean_alpha1 = 0, mean_alpha2 = 0;
  vector<float> PB_L(6,0.0);
  vector<Vector_int> imageno(pow(ncell,3));

  Vector dipole;
  Matrix Tij;
  init_Vector_zero(Field, 1, nmol); 

  replica(L, ncell, PB_L, imageno);
  for(uint i = 0;i < nmol;++i)
    {
      idi = nmol*t+i;  
      for(uint j = 0;j < nmol;++j)
        {
          idj = nmol*t+j;
          for(uint cell = 0; cell < pow(ncell,3); ++cell)
            {
              dist(mol, idi, idj, L, PB_L, imageno, cell, x, y, z );
              rij = mindis(x,y,z,PB_L);       
              if ((i != j && cell == 0 && rij < rcut  ) || (cell > 0 && rij < rcut))
              {               
                r3  = pow(rij, -3.0)  ; 
                r5  = 3.0 * pow(rij, -5.0) ;

                pc =  debyetoangstrom * mol[idj].q * pointchargedistancetodebye;
                pd =  debyetoangstrom * ((mol[idj].PD.x + mol[idj].ID.x ) * x + (mol[idj].PD.y + mol[idj].ID.y ) * y + (mol[idj].PD.z + mol[idj].ID.z ) * z ) ; 

                Field[i].x += x * ( r3 * pc + r5 * pd ) ;
                Field[i].y += y * ( r3 * pc + r5 * pd ) ;
                Field[i].z += z * ( r3 * pc + r5 * pd ) ;

                Tij.xx = -r3 + r5 * x * x ;
                Tij.yy = -r3 + r5 * y * y ; 
                Tij.zz = -r3 + r5 * z * z ;
                Tij.xy =  0  + r5 * x * y ;  
                Tij.xz =  0  + r5 * x * z ;
                Tij.yz =  0  + r5 * y * z ; 
                Tij.yx =  0  + r5 * y * x ;    
                Tij.zx =  0  + r5 * z * x ;
                Tij.zy =  0  + r5 * z * y ; 

                dipole.x = debyetoangstrom * (mol[idj].PD.x + mol[idj].ID.x ); 
                dipole.y = debyetoangstrom * (mol[idj].PD.y + mol[idj].ID.y ); 
                dipole.z = debyetoangstrom * (mol[idj].PD.z + mol[idj].ID.z ); 

                Field[i].x += Tij.xx * dipole.x + Tij.xy * dipole.y + Tij.xz * dipole.z;
                Field[i].y += Tij.yx * dipole.x + Tij.yy * dipole.y + Tij.yz * dipole.z;
                Field[i].z += Tij.zx * dipole.x + Tij.zy * dipole.y + Tij.zz * dipole.z;
              }                          
            }
        }
    }
}



//Field Tensor due to permanent and induced polarisability
void TensorduetoPolarisability(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Matrix> & C)
{
  uint idi = 0, idj = 0;
  double x = 0, y = 0, z = 0, rij = 0;
  double u = 0.0, ar = 0.0, st1 = 0.0, st2 = 0.0, r3 = 0.0, r5 = 0.0;
  float mean_alpha1 = 0, mean_alpha2 = 0;
  vector<float> PB_L(6,0.0);
  vector<Vector_int> imageno(pow(ncell,3));
  Matrix polar;
  Matrix Tij;

  init_Matrix_zero(C, 1, nmol);
  replica(L, ncell, PB_L, imageno);
  for(uint i = 0;i < nmol;++i)
    {
      idi = nmol*t+i;  
      for(uint j = 0;j < nmol;++j)
        {
          idj = nmol*t+j;
          for(uint cell = 0; cell < pow(ncell,3); ++cell)
            {
              dist(mol, idi, idj, L, PB_L, imageno, cell, x, y, z );
              rij = mindis(x,y,z,PB_L);       
              if ((i != j && cell == 0 && rij < rcut ) || (cell > 0 && rij < rcut))
              {
                r3  = pow(rij, -3.0)  ; 
                r5  = 3.0 * pow(rij, -5.0) ;
                /*
                  without damping function 
                  r3  = pow(rij, -3.0)  ; 
                  r5  = 3.0 * pow(rij, -5.0) ;
 
                  Thole damping function
                  mean_alpha1 = (mol[idi].PPol.xx + mol[idi].PPol.yy + mol[idi].PPol.zz ) / 3.0;
                  mean_alpha2 = (mol[idj].PPol.xx + mol[idj].PPol.yy + mol[idj].PPol.zz ) / 3.0;
                  u = rij / pow(mean_alpha1 * mean_alpha2, 0.16666);
                  ar  = mol[idj].sl * pow(u,3.0)  ;   
                  st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                  st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar);                   
                  r3  = pow(rij, -3.0)  * (st1 + 0.0); 
                  r5  = 3.0 * pow(rij, -5.0)  * (st2 + 0.0);

                  Wu and Yang damping function 
                  float cdamp = 3.54 ; 
                  float r_vdwr = mol[idi].vdwr + mol[idj].vdwr ;
                  float WY1 =  pow( 1.0  - exp( -1.0 * cdamp * pow( rij / r_vdwr ,3)) , 2);
                  r3  = pow(rij, -3.0)  * WY1 ; 
                  r5  = 3.0 * pow(rij, -5.0)  * WY1;
                */


                Tij.xx = -r3 + r5 * x * x ;
                Tij.yy = -r3 + r5 * y * y ; 
                Tij.zz = -r3 + r5 * z * z ;
                Tij.xy =  0  + r5 * x * y ;  
                Tij.xz =  0  + r5 * x * z ;
                Tij.yz =  0  + r5 * y * z ; 
                Tij.yx =  0  + r5 * y * x ;    
                Tij.zx =  0  + r5 * z * x ;
                Tij.zy =  0  + r5 * z * y ; 

                polar.xx = mol[idj].PPol.xx + mol[idj].IPol.xx ;
                polar.yy = mol[idj].PPol.yy + mol[idj].IPol.yy ; 
                polar.zz = mol[idj].PPol.zz + mol[idj].IPol.zz ;
                polar.xy = mol[idj].PPol.xy + mol[idj].IPol.xy ;  
                polar.xz = mol[idj].PPol.xz + mol[idj].IPol.xz ;
                polar.yx = mol[idj].PPol.yx + mol[idj].IPol.yx ;    
                polar.yz = mol[idj].PPol.yz + mol[idj].IPol.yz ; 
                polar.zx = mol[idj].PPol.zx + mol[idj].IPol.zx ;
                polar.zy = mol[idj].PPol.zy + mol[idj].IPol.zy ;

                C[i].xx += Tij.xx * polar.xx + Tij.xy * polar.yx + Tij.xz * polar.zx;
                C[i].xy += Tij.xx * polar.xy + Tij.xy * polar.yy + Tij.xz * polar.zy;
                C[i].xz += Tij.xx * polar.xz + Tij.xy * polar.yz + Tij.xz * polar.zz;

                C[i].yx += Tij.yx * polar.xx + Tij.yy * polar.yx + Tij.yz * polar.zx;
                C[i].yy += Tij.yx * polar.xy + Tij.yy * polar.yy + Tij.yz * polar.zy;
                C[i].yz += Tij.yx * polar.xz + Tij.yy * polar.yz + Tij.yz * polar.zz;

                C[i].zx += Tij.zx * polar.xx + Tij.zy * polar.yx + Tij.zz * polar.zx;
                C[i].zy += Tij.zx * polar.xy + Tij.zy * polar.yy + Tij.zz * polar.zy;
                C[i].zz += Tij.zx * polar.xz + Tij.zy * polar.yz + Tij.zz * polar.zz;
              }                          
            }
        }
    }
}



// Field due to externally applied field
void FieldduetoExternalField(vector<Molecular> &mol, uint t, uint nmol, const vector<Vector> &E,  vector <Vector> & Field)
{
  uint idi = 0;
  for(uint i = 0;i < nmol;++i)
    {
      idi = nmol*t+i; 
      Field[i].x += polfieldtodebye * E[t].x ;
      Field[i].y += polfieldtodebye * E[t].y ;
      Field[i].z += polfieldtodebye * E[t].z ;
    }
}



//Field due to induced dipoles
void Fieldduetodipole(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field)
{
  uint idi = 0, idj = 0;
  double x = 0, y = 0, z = 0, rij = 0;
  double u = 0.0, ar = 0.0, st1 = 0.0, st2 = 0.0, r3 = 0.0, r5 = 0.0, pc = 0.0, pd = 0.0;
  float mean_alpha1 = 0, mean_alpha2 = 0;
  vector<float> PB_L(6,0.0);
  vector<Vector_int> imageno(pow(ncell,3));
  Vector dipole;
  Matrix Tij;

  init_Vector_zero(Field, 1, nmol); 

  replica(L, ncell, PB_L, imageno);
  for(uint i = 0;i < nmol;++i)
    {
      idi = nmol*t+i; 
      for(uint j = 0;j < nmol;++j)
        {
          idj = nmol*t+j;
          for(uint cell = 0; cell < pow(ncell,3); ++cell)
            {
              dist(mol, idi, idj, L, PB_L, imageno, cell, x, y, z );
              rij = mindis(x,y,z,PB_L);       
              if ((i != j && cell == 0 && rij < rcut  ) || (cell > 0 && rij < rcut))
              {
                r3  = pow(rij, -3.0)  ; 
                r5  = 3.0 * pow(rij, -5.0) ;

                /*
                  without damping function 
                  r3  = pow(rij, -3.0)  ; 
                  r5  = 3.0 * pow(rij, -5.0) ;
 
                  Thole damping function
                  mean_alpha1 = (mol[idi].PPol.xx + mol[idi].PPol.yy + mol[idi].PPol.zz ) / 3.0;
                  mean_alpha2 = (mol[idj].PPol.xx + mol[idj].PPol.yy + mol[idj].PPol.zz ) / 3.0;
                  u = rij / pow(mean_alpha1 * mean_alpha2, 0.16666);
                  ar  = mol[idj].sl * pow(u,3.0)  ;   
                  st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                  st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar);                   
                  r3  = pow(rij, -3.0)  * (st1 + 0.0); 
                  r5  = 3.0 * pow(rij, -5.0)  * (st2 + 0.0);

                  Wu and Yang damping function 
                  float cdamp = 3.54 ; 
                  float r_vdwr = mol[idi].vdwr + mol[idj].vdwr ;
                  float WY1 =  pow( 1.0  - exp( -1.0 * cdamp * pow( rij / r_vdwr ,3)) , 2);
                  r3  = pow(rij, -3.0)  * WY1 ; 
                  r5  = 3.0 * pow(rij, -5.0)  * WY1;
                */

                Tij.xx = -r3 + r5 * x * x ;
                Tij.yy = -r3 + r5 * y * y ; 
                Tij.zz = -r3 + r5 * z * z ;
                Tij.xy =  0  + r5 * x * y ;  
                Tij.xz =  0  + r5 * x * z ;
                Tij.yz =  0  + r5 * y * z ; 
                Tij.yx =  0  + r5 * y * x ;    
                Tij.zx =  0  + r5 * z * x ;
                Tij.zy =  0  + r5 * z * y ; 

                dipole.x = mol[idj].ID.x; 
                dipole.y = mol[idj].ID.y; 
                dipole.z = mol[idj].ID.z; 

                Field[i].x += Tij.xx * dipole.x + Tij.xy * dipole.y + Tij.xz * dipole.z;
                Field[i].y += Tij.yx * dipole.x + Tij.yy * dipole.y + Tij.yz * dipole.z;
                Field[i].z += Tij.zx * dipole.x + Tij.zy * dipole.y + Tij.zz * dipole.z;        
              }                          
            }
        }
    }
}


// Field due to permanent multipoles [charge, dipole, quadrupole]
void FieldduetoPermanentMultipoles(vector<Molecular> &mol, uint t, uint nmol, const vector<float> & L, vector <Vector> & Field)
{
  uint idi = 0, idj = 0;
  double x = 0, y = 0, z = 0, rij = 0;
  double u = 0.0, ar = 0.0, st1 = 0.0, st2 = 0.0, r3 = 0.0, r5 = 0.0, pc = 0.0, pd = 0.0;
  float mean_alpha1 = 0, mean_alpha2 = 0;
  vector<float> PB_L(6,0.0);
  vector<Vector_int> imageno(pow(ncell,3));

  init_Vector_zero(Field, 1, nmol); 

  replica(L, ncell, PB_L, imageno);
  for(uint i = 0;i < nmol;++i)
    {
      idi = nmol*t+i;  
      for(uint j = 0;j < nmol;++j)
        {
          idj = nmol*t+j;
          for(uint cell = 0; cell < pow(ncell,3); ++cell)
            {
              dist(mol, idi, idj, L, PB_L, imageno, cell, x, y, z );
              rij = mindis(x,y,z,PB_L);       
              if ((i != j && cell == 0 && rij < rcut  ) || (cell > 0 && rij < rcut))
              {
                 /* vdw radii
                   1.  https://www.cgl.ucsf.edu/chimerax/docs/user/radii.html
                   2.  Consistent van der Waals Radii for the Whole Main Group and Van der Waals Radii of Elements
                   3. Consistent Treatment of Inter- and Intramolecular Polarization in Molecular Mechanics Calculations
                   4.  Set of van der Waals and Coulombic Radii of Protein  Atoms for Molecular and Solvent-Accessible Surface Calculation, Packing Evaluation, and Docking
                   5. https://bip.weizmann.ac.il/course/structbioinfo/databases/CCDC_Mercury/appa_glossary.4.73.html
                   6.  Bondi, J.Phys.Chem., 68, 441, 1964
                   7. Semiempirical GGA-Type Density Functional Constructed with a Long-Range Dispersion Correction
                 */
                 /* Damping functions
                   1. Wu and Yang damping function J. Chem. Phys. 2002, 116, 515-524 and Transferable ab initio intermolecular potentials. 1. Derivation from methanol dimer and trimer calculations
                   2. Effect of the Damping Function in Dispersion Corrected Density Functional Theory
                   3. A Universal Damping Function for Empirical Dispersion Correctionon Density Functional Theory
                   4.  https://www.pamoc.it/kw_dfd.html
                   5. Density Functional Theory Augmented with an Empirical Dispersion Term. Interaction Energies and Geometries of 80 Noncovalent Complexes Compared with Ab Initio Quantum Mechanics Calculations
                   6. Damping functions in the effective fragmentpotential method
                   7. On the performance of molecular polarization methods. II. Water and carbon tetrachloride close to a cation
                   8. Gaussian and Exponeital damping function:  The polarizable point dipoles method withelectrostatic damping: Implementation on amodel system
                   9.  Becke-Johnson  (BJ)  damping  function: A generally applicable atomic-chargedependent London dispersion correction
                   10. Thole damping:  Molecular and Atomic Polarizabilities: Thole’s Model Revisited
                   11. Thole damping: Molecular polarizabilities calculated with a modified dipole interaction
                   12. sum of vdwradii: Empirical D3 Dispersion as a Replacement for ab Initio Dispersion Terms in Density Functional Theory-Based Symmetry-Adapted Perturbation Theory
                   13. The polarizable point dipoles method withelectrostatic damping: Implementation on amodel system
                   14. Interatomic Methods for the Dispersion Energy Derived from the AdiabaticConnection Fluctuation-Dissipation Theorem
                   Note: Fermic  and Hesselmann daming function goes to zero quickly at smaller r, not suitable for here
                 */

                //Thole damping function
                r3  = pow(rij, -3.0)  ; 
                r5  = 3.0 * pow(rij, -5.0) ;

                /*
                  without damping function 
                  r3  = pow(rij, -3.0)  ; 
                  r5  = 3.0 * pow(rij, -5.0) ;
 
                  Thole damping function
                  mean_alpha1 = (mol[idi].PPol.xx + mol[idi].PPol.yy + mol[idi].PPol.zz ) / 3.0;
                  mean_alpha2 = (mol[idj].PPol.xx + mol[idj].PPol.yy + mol[idj].PPol.zz ) / 3.0;
                  u = rij / pow(mean_alpha1 * mean_alpha2, 0.16666);
                  ar  = mol[idj].sl * pow(u,3.0)  ;   
                  st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                  st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar);                   
                  r3  = pow(rij, -3.0)  * (st1 + 0.0); 
                  r5  = 3.0 * pow(rij, -5.0)  * (st2 + 0.0);

                  Wu and Yang damping function 
                  float cdamp = 3.54 ; 
                  float r_vdwr = mol[idi].vdwr + mol[idj].vdwr ;
                  float WY1 =  pow( 1.0  - exp( -1.0 * cdamp * pow( rij / r_vdwr ,3)) , 2);
                  r3  = pow(rij, -3.0)  * WY1 ; 
                  r5  = 3.0 * pow(rij, -5.0)  * WY1;
                */


                // Field due to point charge and dipole (quadrupole not considered)
                //point charge - http://www.physics.umd.edu/courses/Phys260/agashe/S10/notes/lecture18.pdf 
                //point dipole - http://www.physnet.org/modules/pdf_modules/m120.pdf 
                // Here we follow the standart physics convention to compare
                // with CP2K: dipole moment from - to +. hence i->j direction
                // In tinker, dipole moment from + to -. hence j->i direction, 
                // However, at the end in accelration computation (verlet.f), they multiply with -ve sign
                pc =  mol[idj].q * pointchargedistancetodebye;
                pd =  mol[idj].PD.x * x + mol[idj].PD.y * y + mol[idj].PD.z * z  ; 
                Field[i].x += x * ( r3 * pc + r5 * pd ) ;
                Field[i].y += y * ( r3 * pc + r5 * pd ) ;
                Field[i].z += z * ( r3 * pc + r5 * pd ) ;

                // Field due to dipole moment (quadrupole moment not considered)
                Field[i].x += -r3 * mol[idj].PD.x  ;
                Field[i].y += -r3 * mol[idj].PD.y  ;
                Field[i].z += -r3 * mol[idj].PD.z  ;
              }                          
            }
        }
    }
}


void PrintOpticalBirefringence(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_cation (nsteps, 0.0),
                total_OBF_anion  (nsteps, 0.0),
                total_OBF_both   (nsteps, 0.0),
                total_OBF_ions   (nsteps, 0.0), 
                total_OBF_rem    (nsteps, 0.0);

  float count = 0, count_cation = 0, count_anion = 0, count_both = 0, count_rem = 0, count_ions = 0 ;
  float rij = 0, temp = 0, x = 0, y = 0, z = 0 ;
  uint hbond_cation = 0, hbond_anion = 0;

  for(uint i = 0;i < nmol;++i)
    {
      if(1) 
        {
          hbond_cation = 0, hbond_anion = 0;
          for(uint j = 0;j < nmol;++j)
            { 
              uint idi = nmol*3500+i;  
              uint idj = nmol*3500+j;

              if((mol[j].MOL[0] == 'M' || mol[j].MOL[0] == 'N') && i != j )
                {
                  x = min_distance(mol[idj].x - mol[idi].x, L[0]);
                  y = min_distance(mol[idj].y - mol[idi].y, L[1]);
                  z = min_distance(mol[idj].z - mol[idi].z, L[2]); 
                  rij = mindis(x,y,z,L); 
                  if(rij < 3.0 && rij > 0)
                    {
                      hbond_cation += 1;  
                    }      
                }
              if((mol[j].MOL[0] == 'C' || mol[j].MOL[0] == 'F') && i != j ) 
                {
                  x = min_distance(mol[idj].x - mol[idi].x, L[0]);
                  y = min_distance(mol[idj].y - mol[idi].y, L[1]);
                  z = min_distance(mol[idj].z - mol[idi].z, L[2]); 
                  rij = mindis(x,y,z,L); 
                  if(rij < 3.0 && rij > 0)
                    {
                      hbond_anion += 1;
                    }      
                }
            }

          /*counting */
          count += 1;
          if(hbond_cation == 0 && hbond_anion == 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
            {
              count_rem   += 1; 
            }
          else if(hbond_cation > 0 && hbond_anion > 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
            {
              count_both  += 1; 
            }
          else if(hbond_cation == 0 && hbond_anion > 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
            {
              count_anion += 1; 
            }
          else if(hbond_cation > 0 && hbond_anion == 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
            {
              count_cation += 1; 
            }

          if((hbond_cation > 0 || hbond_anion > 0  || (hbond_cation + hbond_anion) > 0 ) && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O' )
            {
              count_ions += 1; 
            }

          for(uint t = 0; t < nsteps; t += deltat)
            {
              uint id = nmol*t+i;
              temp          =  (mol[id].PPol.xx + mol[id].IPol.xx) - 0.5 * (mol[id].PPol.yy + mol[id].IPol.yy + mol[id].PPol.zz + mol[id].IPol.zz )  ; 
              total_OBF[t] += temp ;

              if(hbond_cation == 0 && hbond_anion == 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
                {
                  total_OBF_rem[t] += temp ;
                }
              else if(hbond_cation > 0 && hbond_anion > 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
                {
                  total_OBF_both[t] += temp ;
                }
              else if(hbond_cation == 0 && hbond_anion > 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
                {
                  total_OBF_anion[t] += temp ;
                }
              else if(hbond_cation > 0 && hbond_anion == 0 && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O')
                {
                  total_OBF_cation[t] += temp ;
                }

               if((hbond_cation > 0 || hbond_anion > 0  || (hbond_cation + hbond_anion) > 0 ) && mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O' )
                {
                  total_OBF_ions[t] += temp ;
                }
            }
        } 
    }

  if(count == 0)        {count = 1;}
  if(count_rem == 0)    {count_rem = 1;}
  if(count_both == 0)   {count_both = 1;}
  if(count_anion == 0)  {count_anion = 1;}
  if(count_cation == 0) {count_cation = 1;}
  if(count_ions == 0)   {count_ions = 1;} 

  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 0; t < nsteps; t += deltat)
    {
      outfile << t*dt << "  " << total_OBF[t]        / count        << "  " << 
                                 total_OBF_cation[t] / count_cation << "  " << 
                                 total_OBF_anion[t]  / count_anion  << "  " << 
                                 total_OBF_both[t]   / count_both   << "  " << 
                                 total_OBF_ions[t]   / count_ions   << "  " << 
                                 total_OBF_rem[t]    / count_rem    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}




