
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


using namespace std;


void Induced_dipole_pol(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, uint niter, vector<Vector> &E)
{
  float x = 0, y = 0, z = 0, rij = 0;
  float rcut = 20;
  Matrix Tij;
  Vector Eion;
  vector<Vector> TD (nmol*nsteps), dummyTD (nmol) ; 
  vector<Matrix> TP (nmol*nsteps), dummyTP (nmol) ; 
  copydata(mol, nsteps, nmol, TP, TD );  
  for(uint t = 0; t < nsteps; ++t )
    {   
      //first induced dipole Mu_i due to External field E and ionic charges are computed: Mu_i = alpha_0 * E + alpha_0 * q * r /pow(rij,2)
      //rhat is the unit vector describing the line between i and j molecules
      for(uint i = 0;i < nmol;++i)
        {
          uint idi = nmol*t+i;  
          Pol_Efield(mol[idi].PPol, E[t], TD[idi]);
          for(uint j = 0;j < nmol;++j)
            {
              uint idj = nmol*t+j; 
              x = min_distance(mol[idi].x - mol[idj].x, L[0]);
              y = min_distance(mol[idi].y - mol[idj].y, L[1]);
              z = min_distance(mol[idi].z - mol[idj].z, L[2]);
              rij = mindis(x,y,z,L); 
              if(i == j){}
              else if (i != j && rij < rcut )
                {
                  //http://www.physics.umd.edu/courses/Phys260/agashe/S10/notes/lecture18.pdf
                  float b_x = x, b_y = y, b_z = z;          
                  float sum  =  pow(b_x * b_x  + b_y * b_y + b_z * b_z, 0.5);
                  b_x = b_x / sum ;
                  b_y = b_y / sum ;
                  b_z = b_z / sum ;

                  Eion.x = mol[idj].q * b_x / pow(rij,2.0) ;
                  Eion.y = mol[idj].q * b_y / pow(rij,2.0) ;
                  Eion.z = mol[idj].q * b_z / pow(rij,2.0) ;
                  Pol_Efield(mol[idi].PPol, Eion, TD[idi]);
                }
            } 
          mol[idi].ID.x += TD[idi].x ;
          mol[idi].ID.y += TD[idi].y ;
          mol[idi].ID.z += TD[idi].z ;        
        }

      for(uint iter = 0; iter < niter ; ++iter)
        {          
          init_Matrix_zero(dummyTP, 1, nmol);
          init_Vector_zero(dummyTD, 1, nmol);     
          for(uint i = 0;i < nmol;++i)
            {
              uint idi = nmol*t+i; 
              for(uint j = 0;j < nmol;++j)
                {
                  uint idj = nmol*t+j; 
                  x = min_distance(mol[idi].x - mol[idj].x, L[0]);
                  y = min_distance(mol[idi].y - mol[idj].y, L[1]);
                  z = min_distance(mol[idi].z - mol[idj].z, L[2]);
                  rij = mindis(x,y,z,L); 
                  if(i == j){}
                  else if (i != j && rij < rcut )
                  {
                    dipoletensorfield(Tij, rij, x, y, z); 
                    Tij_dipole(Tij, TD[idj], dummyTD[i]);
                    Tij_Pol(Tij, TP[idj], dummyTP[i]);   
                  }
                }
            }
          for(uint i = 0;i < nmol;++i)
            {
              uint idi = nmol*t+i; 
              Mat_vec(mol[idi].PPol, dummyTD[i], TD[idi]);
              Mat_Mat(mol[idi].PPol, dummyTP[i], TP[idi]);

              mol[idi].ID.x += TD[idi].x ;
              mol[idi].ID.y += TD[idi].y ;
              mol[idi].ID.z += TD[idi].z ;

              mol[idi].IPol.xx += TP[idi].xx;mol[idi].IPol.yy += TP[idi].yy;mol[idi].IPol.zz += TP[idi].zz;
              mol[idi].IPol.xy += TP[idi].xy;mol[idi].IPol.xz += TP[idi].xz;mol[idi].IPol.yz += TP[idi].yz;
              mol[idi].IPol.yx += TP[idi].yx;mol[idi].IPol.zx += TP[idi].zx;mol[idi].IPol.zy += TP[idi].zy;
            }
        } 
    }
}



void parameters(Molecular &mol)
{  
/* The parameter for the polarizability are generally set to axx = 1.626 angstrom,
   ayy = 1.495 angstrom,and azz = 1.286 angstrom. 
   While the Xaxis corresponds to the vector that connects H-H atoms, 
   the Yaxis is the water bisector, and
   the Zaxis is a vector perpendicular to both H-H vector and water bisector
   cite the following papers [Struct. Dyn. 2, 054102 (2015); https://doi.org/10.1063/1.4932597,  
   J. Chem. Phys. 141, 184201 (2014); https://doi.org/10.1063/1.4901216]
   Morita Hynes J. Phys. Chem. B 2002, 106, 673-685,
   Ni and Skinner  J. Chem. Phys. 143, 014502 (2015); https://doi.org/10.1063/1.4923462 
   Remi Khatib, J.Phys.Chem.C2016, 120, 18665−18673
   Yuki nagata papers
   Kenneth J Miller, Tholes, Siberstein and van Duijnen https://doi.org/10.1021/jp980221f research works
   Munir S. Skaf and Sérgio M. Vechi J. Chem. Phys 119, 2181 (2003); https://doi.org/10.1063/1.1583677  

   Simulation of the Far Infrared Spectrum of Liquid Water and Steam Along the Coexistence Curve
   from the Triple Point to the Critical Point https://doi.org/10.1007/978-94-011-0183-7_10

   Here the alpha paramters are modified to accurately reproduce the polarisabilites, computed
   using the CP2K program with the same configuration  

   The mulliken charges were used and parameterized with the help of CP2K program
*/

  if(mol.MOL[0] == 'H' && mol.MOL[1] == '2' && mol.MOL[2] == 'O')
    {
      mol.q = 0;
      mol.PPol.xx = 1.3227;mol.PPol.yy =  1.1191;mol.PPol.zz = 0.8808;
      mol.PPol.xy = 0.0018;mol.PPol.xz = -0.0003;mol.PPol.yz = 0.0006;
      mol.PPol.yx = 0.0016;mol.PPol.zx = -0.0002;mol.PPol.zy = 0.0006; 

      //mol.PPol.xy = 0.0000;mol.PPol.xz = -0.0000;mol.PPol.yz = 0.0000;
      //mol.PPol.yx = 0.0000;mol.PPol.zx = -0.0000;mol.PPol.zy = 0.0000; 

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'M' || mol.MOL[0] == 'm' ) && ( mol.MOL[1] == 'G' || mol.MOL[1] == 'g' ))
    {
      mol.q = 0;
      mol.PPol.xx = 10.090;mol.PPol.yy = 10.090;mol.PPol.zz = 10.090;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'C' || mol.MOL[0] == 'c' ) && ( mol.MOL[1] == 'L' || mol.MOL[1] == 'l' ))
    {
      mol.q = 0;
      mol.PPol.xx = 1.4300;mol.PPol.yy = 1.4300;mol.PPol.zz = 1.4300;
      mol.PPol.xy = 0.0087;mol.PPol.xz =-0.0080;mol.PPol.yz = 0.0083;
      mol.PPol.yx = 0.0087;mol.PPol.zx =-0.0080;mol.PPol.zy = 0.0083;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'N' || mol.MOL[0] == 'n' ) && ( mol.MOL[1] == 'A' || mol.MOL[1] == 'a' ))
    {
      mol.q = 0;
      mol.PPol.xx = 22.050;mol.PPol.yy = 22.050;mol.PPol.zz = 22.050;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if(mol.MOL[0] == 'S' && mol.MOL[1] == 'O' && mol.MOL[2] == '4')
    {
      mol.q = 0;
      mol.PPol.xx = 6.104;mol.PPol.yy = 6.0861;mol.PPol.zz = 6.1769;
      mol.PPol.xy =-0.013;mol.PPol.xz = 0.0011;mol.PPol.yz = 0.0808;
      mol.PPol.yx =-0.014;mol.PPol.zx = 0.0006;mol.PPol.zy = 0.0806; 

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
}


// Transforming atomic to molecular assigns dipole moment, center of mass and permanent polarisabilites to the variables
void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol)
{
  Molecular mols; 
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i;          
          // NA
          if((r[id].symbol[0] == 'N' || r[id].symbol[0] == 'n' ) && ( r[id].symbol[1] == 'A' || r[id].symbol[1] == 'a'))
            {
              mols.x = r[id].x;
              mols.y = r[id].y;
              mols.z = r[id].z;
              mols.MOL = r[id].symbol;
              mols.m = r[id].atomicmass; 
              mols.PD.x = 0 ;
              mols.PD.y = 0 ;
              mols.PD.z = 0 ;
              mols.ID.x = 0 ;
              mols.ID.y = 0 ;
              mols.ID.z = 0 ;
              parameters(mols);
              mol.push_back(mols); 
            }
          //MG
          if((r[id].symbol[0] == 'M' || r[id].symbol[0] == 'm' ) && ( r[id].symbol[1] == 'G' || r[id].symbol[1] == 'g'))
            {
              mols.x = r[id].x;
              mols.y = r[id].y;
              mols.z = r[id].z;
              mols.MOL = r[id].symbol;
              mols.m = r[id].atomicmass; 
              mols.PD.x = 0 ;
              mols.PD.y = 0 ;
              mols.PD.z = 0 ;
              mols.ID.x = 0 ;
              mols.ID.y = 0 ;
              mols.ID.z = 0 ;
              parameters(mols);
              mol.push_back(mols); 
            }
          // SO4
          if(r[id].symbol[0] == 'S' && r[id+1].symbol[0] == 'O' && r[id+2].symbol[0] == 'O' && r[id+3].symbol[0] == 'O')
            { 
              const float am_S = 32.065 * amu, am_O = 16 * amu, am_SO4 = 96.06 * amu ;
              mols.x = ( r[id-1].x * am_O + r[id].x * am_S + r[id+1].x * am_O + r[id+2].x * am_O + r[id+3].x * am_O ) / am_SO4 ;
              mols.y = ( r[id-1].y * am_O + r[id].y * am_S + r[id+1].y * am_O + r[id+2].y * am_O + r[id+3].y * am_O ) / am_SO4 ;
              mols.z = ( r[id-1].z * am_O + r[id].z * am_S + r[id+1].z * am_O + r[id+2].z * am_O + r[id+3].z * am_O ) / am_SO4 ;               
              mols.MOL = "SO4";
              mols.m = 96.06;  
              float wb_x = 0,wb_y = 0,wb_z = 0;           // bisector vector
              wb_x = min_distance(r[id-1].x - r[id].x, L[0]) + min_distance(r[id+1].x - r[id].x, L[0]) + min_distance(r[id+2].x - r[id].x, L[0]) + min_distance(r[id+3].x - r[id].x, L[0]) ;
              wb_y = min_distance(r[id-1].y - r[id].y, L[1]) + min_distance(r[id+1].y - r[id].y, L[1]) + min_distance(r[id+2].y - r[id].y, L[1]) + min_distance(r[id+3].y - r[id].y, L[1]) ;
              wb_z = min_distance(r[id-1].z - r[id].z, L[2]) + min_distance(r[id+1].z - r[id].z, L[2]) + min_distance(r[id+2].z - r[id].z, L[2]) + min_distance(r[id+3].z - r[id].z, L[2]) ;
              float sum  =  pow(wb_x * wb_x  + wb_y * wb_y + wb_z * wb_z, 0.5);
              wb_x = wb_x / sum ;
              wb_y = wb_y / sum ;
              wb_z = wb_z / sum ;
              mols.PD.x = 0  ;
              mols.PD.y = 0  ;
              mols.PD.z = 0  ;
              mols.ID.x = 0 ;
              mols.ID.y = 0 ;
              mols.ID.z = 0 ;
              parameters(mols);
              mol.push_back(mols);
            } 
          // H2O
          if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
            { 
              const float am_H = 1.00784 * amu, am_O = 15.999 * amu, am_H2O = 18.015 * amu ;
              float wb_x = 0,wb_y = 0,wb_z = 0;           // bisector vector
              float wb1_x = 0,wb1_y = 0,wb1_z = 0;           // bisector vector
              float wb2_x = 0,wb2_y = 0,wb2_z = 0;           // bisector vector
              float HH_x = 0,HH_y = 0,HH_z = 0 ;           // H-H vector
              float Pv_x = 0,Pv_y = 0,Pv_z = 0;           // vector Perpendicular to bisector and H-H vector
              mols.x = ( r[id].x * am_O + r[id+1].x * am_H + r[id+2].x * am_H ) / am_H2O ;
              mols.y = ( r[id].y * am_O + r[id+1].y * am_H + r[id+2].y * am_H ) / am_H2O ;
              mols.z = ( r[id].z * am_O + r[id+1].z * am_H + r[id+2].z * am_H ) / am_H2O ;

              // water bisector vector, Minimum image convention was applied
              wb1_x = min_distance(r[id+1].x - r[id].x, L[0]) ;
              wb1_y = min_distance(r[id+1].y - r[id].y, L[1]) ;
              wb1_z = min_distance(r[id+1].z - r[id].z, L[2]) ;
              wb2_x = min_distance(r[id+2].x - r[id].x, L[0]) ;
              wb2_y = min_distance(r[id+2].y - r[id].y, L[1]) ;
              wb2_z = min_distance(r[id+2].z - r[id].z, L[2]) ;
              float sum  =  pow(wb1_x * wb1_x  + wb1_y * wb1_y + wb1_z * wb1_z, 0.5);
              wb1_x = wb1_x / sum ;
              wb1_y = wb1_y / sum ;
              wb1_z = wb1_z / sum ;
              sum  =  pow(wb2_x * wb2_x  + wb2_y * wb2_y + wb2_z * wb2_z, 0.5);
              wb2_x = wb2_x / sum ;
              wb2_y = wb2_y / sum ;
              wb2_z = wb2_z / sum ;

              //water bisector [unit vector]
              wb_x = wb1_x + wb2_x ;
              wb_y = wb1_y + wb2_y ;
              wb_z = wb1_z + wb2_z ;
              sum  =  pow(wb_x * wb_x  + wb_y * wb_y + wb_z * wb_z, 0.5);
              wb_x = wb_x / sum ;
              wb_y = wb_y / sum ;
              wb_z = wb_z / sum ;

              // H-H vector
              HH_x = (r[id+1].x - r[id+2].x) ;
              HH_y = (r[id+1].y - r[id+2].y) ;
              HH_z = (r[id+1].z - r[id+2].z) ;
              sum  =  pow(HH_x * HH_x  + HH_y * HH_y + HH_z * HH_z, 0.5);
              HH_x = HH_x / sum ;
              HH_y = HH_y / sum ;
              HH_z = HH_z / sum ;
              // Vector perpendicular to the water bisector and H-H vector
              Pv_x = wb_y * HH_z - wb_z * HH_y; 
              Pv_y = wb_z * HH_x - wb_x * HH_z; 
              Pv_z = wb_x * HH_y - wb_y * HH_x; 
              sum  =  pow(Pv_x * Pv_x  + Pv_y * Pv_y + Pv_z * Pv_z, 0.5);
              Pv_x = Pv_x / sum ;
              Pv_y = Pv_y / sum ;
              Pv_z = Pv_z / sum ;
              // The H-H vector is slightly adjusted to form perfect orthogonality with other vectors
              HH_x = wb_y * Pv_z - wb_z * Pv_y; 
              HH_y = wb_z * Pv_x - wb_x * Pv_z; 
              HH_z = wb_x * Pv_y - wb_y * Pv_x; 
              sum  =  pow(HH_x * HH_x  + HH_y * HH_y + HH_z * HH_z, 0.5);
              HH_x = HH_x / sum ;
              HH_y = HH_y / sum ;
              HH_z = HH_z / sum ;
    
              mols.MOL = "H2O";
              mols.m = 18; 
        
              //Dipole moment [Debye]
              mols.PD.x = wb_x * Unitvectortodebye ;
              mols.PD.y = wb_y * Unitvectortodebye ;
              mols.PD.z = wb_z * Unitvectortodebye ;  
              mols.ID.x = 0 ;
              mols.ID.y = 0 ;
              mols.ID.z = 0 ;

              // Permanent polarisability [Angstrom]
              parameters(mols);
              float A_xx = mols.PPol.xx, A_yy = mols.PPol.yy, A_zz = mols.PPol.zz,
                    A_xy = mols.PPol.xy, A_xz = mols.PPol.xz, A_yz = mols.PPol.yz,
                    A_yx = mols.PPol.yx, A_zx = mols.PPol.zx, A_zy = mols.PPol.zy;     

              mols.PPol.xx =  HH_x * HH_x * A_xx + wb_x * wb_x * A_yy + Pv_x * Pv_x * A_zz
                           + HH_x * wb_x * A_xy + HH_x * Pv_x * A_xz + wb_x * Pv_x * A_yz 
                           + wb_x * HH_x * A_yx + Pv_x * HH_x * A_zx + Pv_x * wb_x * A_zy ;
              mols.PPol.yy =  HH_y * HH_y * A_xx + wb_y * wb_y * A_yy + Pv_y * Pv_y * A_zz
                           + HH_y * wb_y * A_xy + HH_y * Pv_y * A_xz + wb_y * Pv_y * A_yz 
                           + wb_y * HH_y * A_yx + Pv_y * HH_y * A_zx + Pv_y * wb_y * A_zy ;
              mols.PPol.zz =  HH_z * HH_z * A_xx + wb_z * wb_z * A_yy + Pv_z * Pv_z * A_zz
                           + HH_z * wb_z * A_xy + HH_z * Pv_z * A_xz + wb_z * Pv_z * A_yz 
                           + wb_z * HH_z * A_yx + Pv_z * HH_z * A_zx + Pv_z * wb_z * A_zy ;
              mols.PPol.xy =  HH_x * HH_y * A_xx + wb_x * wb_y * A_yy + Pv_x * Pv_y * A_zz
                           + HH_x * wb_y * A_xy + HH_x * Pv_y * A_xz + wb_x * Pv_y * A_yz 
                           + wb_x * HH_y * A_yx + Pv_x * HH_y * A_zx + Pv_x * wb_y * A_zy ;
              mols.PPol.xz =  HH_x * HH_z * A_xx + wb_x * wb_z * A_yy + Pv_x * Pv_z * A_zz
                           + HH_x * wb_z * A_xy + HH_x * Pv_z * A_xz + wb_x * Pv_z * A_yz 
                           + wb_x * HH_z * A_yx + Pv_x * HH_z * A_zx + Pv_x * wb_z * A_zy ;
              mols.PPol.yz =  HH_y * HH_z * A_xx + wb_y * wb_z * A_yy + Pv_y * Pv_z * A_zz
                           + HH_y * wb_z * A_xy + HH_y * Pv_z * A_xz + wb_y * Pv_z * A_yz 
                           + wb_y * HH_z * A_yx + Pv_y * HH_z * A_zx + Pv_y * wb_z * A_zy ;
              mols.PPol.yx =  HH_y * HH_x * A_xx + wb_y * wb_x * A_yy + Pv_y * Pv_x * A_zz
                           + HH_y * wb_x * A_xy + HH_y * Pv_x * A_xz + wb_y * Pv_x * A_yz 
                           + wb_y * HH_x * A_yx + Pv_y * HH_x * A_zx + Pv_y * wb_x * A_zy ;
              mols.PPol.zx =  HH_z * HH_x * A_xx + wb_z * wb_x * A_yy + Pv_z * Pv_x * A_zz
                           + HH_z * wb_x * A_xy + HH_z * Pv_x * A_xz + wb_z * Pv_x * A_yz 
                           + wb_z * HH_x * A_yx + Pv_z * HH_x * A_zx + Pv_z * wb_x * A_zy ;
              mols.PPol.zy =  HH_z * HH_y * A_xx + wb_z * wb_y * A_yy + Pv_z * Pv_y * A_zz
                           + HH_z * wb_y * A_xy + HH_z * Pv_y * A_xz + wb_z * Pv_y * A_yz 
                           + wb_z * HH_y * A_yx + Pv_z * HH_y * A_zx + Pv_z * wb_y * A_zy ;
              mol.push_back(mols);
            }
        }
    } 
}



void dipoletensorfield(Matrix &Tij, float rij, float x, float y, float z)
{
  float ar  = sl * rij;
  float st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
  float st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar); 
  float r3  = pow(rij, -3.0) * st1 ; 
  float r5  = 3 * pow(rij, -5.0) * st2;

  Tij.xx = r3 - r5 * x * x;
  Tij.yy = r3 - r5 * y * y;
  Tij.zz = r3 - r5 * z * z;
  Tij.xy = 0  - r5 * x * y;  
  Tij.xz = 0  - r5 * x * z;
  Tij.yz = 0  - r5 * y * z; 
  Tij.yx = 0  - r5 * y * x;  
  Tij.zx = 0  - r5 * z * x;
  Tij.zy = 0  - r5 * z * y; 
}



void copydata(const vector<Molecular> &mol, uint nsteps, uint nmol,  vector<Matrix> & TP,  vector<Vector> & TD )
{
  //copying the permanent dipole and polarisability
  for(uint t = 0; t < nsteps; ++t )
    {        
      for(uint i = 0;i < nmol;++i)
        {
          uint id = nmol*t+i; 
          TD[id].x = mol[id].PD.x ;
          TD[id].y = mol[id].PD.y ;
          TD[id].z = mol[id].PD.z ;

          TP[id].xx= mol[id].PPol.xx;TP[id].yy= mol[id].PPol.yy;TP[id].zz= mol[id].PPol.zz;
          TP[id].xy= mol[id].PPol.xy;TP[id].xz= mol[id].PPol.xz;TP[id].yz= mol[id].PPol.yz;
          TP[id].yx= mol[id].PPol.yx;TP[id].zx= mol[id].PPol.zx;TP[id].zy= mol[id].PPol.zy;
        }
    }
}

void Pol_Efield(const Matrix & A, const Vector & b, Vector & dummy)
{  
  dummy.x += polfieldtodebye * (A.xx * b.x + A.xy * b.y + A.xz * b.z );
  dummy.y += polfieldtodebye * (A.yx * b.x + A.yy * b.y + A.yz * b.z );
  dummy.z += polfieldtodebye * (A.zx * b.x + A.zy * b.y + A.zz * b.z );
}


void Pol_Ifield(const Matrix & A, const Vector & b, Vector & dummy)
{  
  dummy.x += polpointchargetodebye * (A.xx * b.x + A.xy * b.y + A.xz * b.z );
  dummy.y += polpointchargetodebye * (A.yx * b.x + A.yy * b.y + A.yz * b.z );
  dummy.z += polpointchargetodebye * (A.zx * b.x + A.zy * b.y + A.zz * b.z );
}


void Tij_dipole(const Matrix & Tij, const Vector & dipole, Vector & dummy)
{  
  dummy.x -= Tij.xx * dipole.x + Tij.xy * dipole.y + Tij.xz * dipole.z;
  dummy.y -= Tij.yx * dipole.x + Tij.yy * dipole.y + Tij.yz * dipole.z;
  dummy.z -= Tij.zx * dipole.x + Tij.zy * dipole.y + Tij.zz * dipole.z;
}


void Tij_Pol(const Matrix & Tij, const Matrix & Pol, Matrix & dummy)
{
  dummy.xx -= Tij.xx * Pol.xx + Tij.xy * Pol.yx + Tij.xz * Pol.zx;
  dummy.xy -= Tij.xx * Pol.xy + Tij.xy * Pol.yy + Tij.xz * Pol.zy;
  dummy.xz -= Tij.xx * Pol.xz + Tij.xy * Pol.yz + Tij.xz * Pol.zz;

  dummy.yx -= Tij.yx * Pol.xx + Tij.yy * Pol.yx + Tij.yz * Pol.zx;
  dummy.yy -= Tij.yx * Pol.xy + Tij.yy * Pol.yy + Tij.yz * Pol.zy;
  dummy.yz -= Tij.yx * Pol.xz + Tij.yy * Pol.yz + Tij.yz * Pol.zz;

  dummy.zx -= Tij.zx * Pol.xx + Tij.zy * Pol.yx + Tij.zz * Pol.zx;
  dummy.zy -= Tij.zx * Pol.xy + Tij.zy * Pol.yy + Tij.zz * Pol.zy;
  dummy.zz -= Tij.zx * Pol.xz + Tij.zy * Pol.yz + Tij.zz * Pol.zz;
}

