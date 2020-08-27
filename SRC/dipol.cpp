
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
  vector<Vector> Field (nmol);
  vector<Vector> rsd (nmol), zrsd (nmol), conj (nmol), vec(nmol);
  uint idi = 0, idj = 0;

  for(uint t = 0; t < nsteps;t += 1 )
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
      // (CG) approach, exactly as implemented in Tinker Molecular
      // Modeling Software
      for(uint iter = 0;iter < 10;++iter)
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

          float a = 0.0 ;
          float sum = 0.0;
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

          float b = 0.0;
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              Mat_vec(mol[idi].PPol, rsd[i], zrsd[i]);
              b = b + rsd[i].x * zrsd[i].x + rsd[i].y * zrsd[i].y + rsd[i].z * zrsd[i].z ;
            }
          b = (b/sum) ; 
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              conj[i].x = zrsd[i].x + b * conj[i].x ;
              conj[i].y = zrsd[i].y + b * conj[i].y ;
              conj[i].z = zrsd[i].z + b * conj[i].z ;
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
  uint idi = 0, idj = 0, ncell = 3;
  float x = 0, y = 0, z = 0, rij = 0, rcut = 10;
  vector<float> PB_L(3,0.0);
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
              if ((i != j && cell == 0 ) || (cell > 0))
              {
                float ar  = sl * rij;
                float st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                float st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar); 
                float r3  = pow(rij, -3.0); 
                float r5  = 3.0 * pow(rij, -5.0);

                Matrix Tij;
                Tij.xx = -r3 + r5 * x * x ;
                Tij.yy = -r3 + r5 * y * y ; 
                Tij.zz = -r3 + r5 * z * z ;
                Tij.xy =  0  + r5 * x * y ;  
                Tij.xz =  0  + r5 * x * z ;
                Tij.yz =  0  + r5 * y * z ; 
                Tij.yx =  0  + r5 * y * x ;    
                Tij.zx =  0  + r5 * z * x ;
                Tij.zy =  0  + r5 * z * y ; 

                Vector dipole;
                dipole.x = mol[idj].PD.x + mol[idj].ID.x; 
                dipole.y = mol[idj].PD.y + mol[idj].ID.y; 
                dipole.z = mol[idj].PD.z + mol[idj].ID.z; 

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
  uint idi = 0, idj = 0, ncell = 3;
  float x = 0, y = 0, z = 0, rij = 0, rcut = 10;
  vector<float> PB_L(3,0.0);
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
              if ((i != j && cell == 0 ) || (cell > 0))
              {
                float ar  = sl * rij;
                float st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                float st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar); 
                float r3  = pow(rij, -3.0); 
                float r5  = 3.0 * pow(rij, -5.0);
                float r7  = 15.0 * pow(rij, -7.0);

                // Field due to point charges and dipole (quadrupole not considered)
                float pc =  mol[idj].q;
                Field[i].x += x * ( r3 * pc * pointchargedistancetodebye ) ;
                Field[i].y += y * ( r3 * pc * pointchargedistancetodebye ) ;
                Field[i].z += z * ( r3 * pc * pointchargedistancetodebye ) ;

                // Field due to dipole moment (quadrupole moment not considered)
                Field[i].x += -r3 * mol[idj].PD.x  ;
                Field[i].y += -r3 * mol[idj].PD.y  ;
                Field[i].z += -r3 * mol[idj].PD.z  ;
              }                          
            }
        }
    }
}




void parameters(Molecular &mol)
{  
/*
   Cp2k Luber2014 Raman spectra from ab initio moleculardynamics and its application to liquid S-methyloxirane
        Putrino2002 , Anharmonic Raman Spectra in High-Pressure Ice fromAb InitioSimulations

              //rhat unit vector describing the line between i and j molecules (Imp: i->j direction)
              //http://www.physics.umd.edu/courses/Phys260/agashe/S10/notes/lecture18.pdf

//  Introductory Bioelectronics: For Engineers and Physical Scientists equaiton 3.5 rhat

   The parameter for the polarizability are generally set to axx = 1.626 angstrom,
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

   Pengyu Ren and Jay W. Ponder Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation
   https://pubs.acs.org/doi/pdf/10.1021/jp027815%2B

   IMOLECULAR POLARCIZABILITLES CALCULATED WITH A MOIXEIEID DIPOLE INTERACTiON
   AnAtomDipoleInteractionModelforMolecularOpticalProperties

    polarisabulity due to field: Flexible bent rod model with a saturatinginduced dipole moment to study theelectric linear dichroism of DNA fragments

   Extended dipole-induced dipole mechanism for generatingRaman and optical Kerr effect intensities oflow-frequency dynamics in liquids

   Polarizability response in polar solvents:Molecular-dynamics simulations ofacetonitrile and chloroform
 
   J Comput Chem. 2007 May ; 28(7): 1261–1274    Gaussian Induced Dipole Polarization Model

   Development of Polarizable Gaussian Model for MolecularMechanical Calculations I: Atomic Polarizability Parameterization ToReproduceab InitioAnisotropy

  //2018 Table of static dipole polarizabilities of theneutral elements in the periodic table

   Here the alpha paramters are modified to accurately reproduce the polarisabilites, computed
   using the CP2K program with the same configuration  

   The mulliken charges were used and parameterized with the help of CP2K program
*/

  // Atomic Polarisability [Angstrom3]
  if(mol.MOL[0] == 'H' && mol.MOL[1] == '2' && mol.MOL[2] == 'O')
    {
      mol.q = 0.0000;
      mol.PPol.xx = 1.3227;mol.PPol.yy =  1.1191;mol.PPol.zz = 0.8808;
      mol.PPol.xy = 0.0018;mol.PPol.xz = -0.0003;mol.PPol.yz = 0.0006;
      mol.PPol.yx = 0.0016;mol.PPol.zx = -0.0002;mol.PPol.zy = 0.0006;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'A' || mol.MOL[0] == 'a' ) && ( mol.MOL[1] == 'R' || mol.MOL[1] == 'r' ))
    {
      mol.q = 0;
      mol.PPol.xx = 1.6400;mol.PPol.yy = 1.6400;mol.PPol.zz = 1.6400;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'M' || mol.MOL[0] == 'm' ) && ( mol.MOL[1] == 'G' || mol.MOL[1] == 'g' ))
    {
      mol.q = 2;
      mol.PPol.xx = 0.0800;mol.PPol.yy = 0.0800;mol.PPol.zz = 0.0800;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'C' || mol.MOL[0] == 'c' ) && ( mol.MOL[1] == 'L' || mol.MOL[1] == 'l' ))
    {
      mol.q = -1;
      mol.PPol.xx = 4.0000;mol.PPol.yy = 4.0000;mol.PPol.zz = 4.0000;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if( ( mol.MOL[0] == 'N' || mol.MOL[0] == 'n' ) && ( mol.MOL[1] == 'A' || mol.MOL[1] == 'a' ))
    {
      mol.q = 0.63;
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
              parameters(mols);
              mols.PD.x = 0;
              mols.PD.y = 0;
              mols.PD.z = 0;
              mols.ID.x = 0.0 ;
              mols.ID.y = 0.0 ;
              mols.ID.z = 0.0 ;
              mol.push_back(mols); 
            }
          //Ar
          if((r[id].symbol[0] == 'A' || r[id].symbol[0] == 'a' ) && ( r[id].symbol[1] == 'R' || r[id].symbol[1] == 'r'))
            {
              mols.x = r[id].x;
              mols.y = r[id].y;
              mols.z = r[id].z;
              mols.MOL = r[id].symbol;
              mols.m = r[id].atomicmass; 
              parameters(mols);
              mols.PD.x = 0;
              mols.PD.y = 0;
              mols.PD.z = 0;
              mols.ID.x = 0.0 ;
              mols.ID.y = 0.0 ;
              mols.ID.z = 0.0 ;
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
              parameters(mols);
              mols.PD.x = 0;
              mols.PD.y = 0;
              mols.PD.z = 0;
              mols.ID.x = 0.0 ;
              mols.ID.y = 0.0 ;
              mols.ID.z = 0.0 ;
              mol.push_back(mols); 
            }
          //CL
          if((r[id].symbol[0] == 'C' || r[id].symbol[0] == 'c' ) && ( r[id].symbol[1] == 'L' || r[id].symbol[1] == 'l'))
            {
              mols.x = r[id].x;
              mols.y = r[id].y;
              mols.z = r[id].z;
              mols.MOL = r[id].symbol;
              mols.m = r[id].atomicmass; 
              parameters(mols);
              mols.PD.x = 0;
              mols.PD.y = 0;
              mols.PD.z = 0;
              mols.ID.x = 0.0 ;
              mols.ID.y = 0.0 ;
              mols.ID.z = 0.0 ;
              mol.push_back(mols); 
            }
          // SO4
          if(r[id].symbol[0] == 'S' && r[id+1].symbol[0] == 'O' && r[id+2].symbol[0] == 'O' && r[id+3].symbol[0] == 'O')
            { 
              const float am_S = 32.065 * amu, am_O = 16 * amu, am_SO4 = 96.06 * amu ;
              //COM
              mols.x = ( r[id-1].x * am_O + r[id].x * am_S + r[id+1].x * am_O + r[id+2].x * am_O + r[id+3].x * am_O ) / am_SO4 ;
              mols.y = ( r[id-1].y * am_O + r[id].y * am_S + r[id+1].y * am_O + r[id+2].y * am_O + r[id+3].y * am_O ) / am_SO4 ;
              mols.z = ( r[id-1].z * am_O + r[id].z * am_S + r[id+1].z * am_O + r[id+2].z * am_O + r[id+3].z * am_O ) / am_SO4 ;               
              mols.MOL = "SO4";
              mols.m = 96.06;  
              parameters(mols);
              mols.PD.x = 0;
              mols.PD.y = 0;
              mols.PD.z = 0;
              mols.ID.x = 0.0 ;
              mols.ID.y = 0.0 ;
              mols.ID.z = 0.0 ;
              mol.push_back(mols);
            } 
          // H2O
          if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
            { 
              const float am_H = 1.00784 * amu, am_O = 15.999 * amu, am_H2O = 18.015 * amu ;
              double wb_x = 0,wb_y = 0,wb_z = 0;           // bisector vector H2O
              double wb1_x = 0,wb1_y = 0,wb1_z = 0;        // bisector vector OH1
              double wb2_x = 0,wb2_y = 0,wb2_z = 0;        // bisector vector OH2
              double HH_x = 0,HH_y = 0,HH_z = 0 ;          // H-H vector
              double Pv_x = 0,Pv_y = 0,Pv_z = 0;           // vector Perpendicular to bisector and H-H vector
              //COM
              mols.x = ( r[id].x * am_O + r[id+1].x * am_H + r[id+2].x * am_H ) / am_H2O ;
              mols.y = ( r[id].y * am_O + r[id+1].y * am_H + r[id+2].y * am_H ) / am_H2O ;
              mols.z = ( r[id].z * am_O + r[id+1].z * am_H + r[id+2].z * am_H ) / am_H2O ;

              // water bisector vector, Minimum image convention was applied
              wb1_x = min_distance(r[id+1].x - r[id].x, L[0]) ;
              wb1_y = min_distance(r[id+1].y - r[id].y, L[1]) ;
              wb1_z = min_distance(r[id+1].z - r[id].z, L[2]) ;
              convertounivector(wb1_x, wb1_y, wb1_z);
              wb2_x = min_distance(r[id+2].x - r[id].x, L[0]) ;
              wb2_y = min_distance(r[id+2].y - r[id].y, L[1]) ;
              wb2_z = min_distance(r[id+2].z - r[id].z, L[2]) ;
              convertounivector(wb2_x, wb2_y, wb2_z);

              //water bisector [unit vector]
              wb_x = wb1_x + wb2_x ;
              wb_y = wb1_y + wb2_y ;
              wb_z = wb1_z + wb2_z ;
              convertounivector(wb_x, wb_y, wb_z);

              // H-H vector
              HH_x = (r[id+1].x - r[id+2].x) ;
              HH_y = (r[id+1].y - r[id+2].y) ;
              HH_z = (r[id+1].z - r[id+2].z) ;
              convertounivector(HH_x, HH_y, HH_z);

              // Vector perpendicular to the water bisector and H-H vector
              Pv_x = wb_y * HH_z - wb_z * HH_y; 
              Pv_y = wb_z * HH_x - wb_x * HH_z; 
              Pv_z = wb_x * HH_y - wb_y * HH_x; 
              convertounivector(Pv_x, Pv_y, Pv_z);

              // The H-H vector is slightly adjusted to form perfect orthogonality with other vectors
              HH_x = wb_y * Pv_z - wb_z * Pv_y; 
              HH_y = wb_z * Pv_x - wb_x * Pv_z; 
              HH_z = wb_x * Pv_y - wb_y * Pv_x; 
              convertounivector(HH_x, HH_y, HH_z);
    
              mols.MOL = "H2O";
              mols.m = 18; 
        
              //Dipole moment [Debye]
              mols.PD.x = wb_x * Unitvectortodebye ;
              mols.PD.y = wb_y * Unitvectortodebye ;
              mols.PD.z = wb_z * Unitvectortodebye ;  
              mols.ID.x = 0.0 ;
              mols.ID.y = 0.0 ;
              mols.ID.z = 0.0 ;

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


void replica(const vector<float> & L, uint ncell, vector<float> & PB_L, vector<Vector_int> & imageno)
{
  uint index = 0;
  PB_L[0] = L[0] * ncell ;
  PB_L[1] = L[1] * ncell ;
  PB_L[2] = L[2] * ncell ;
  for(uint cellz = 0; cellz < ncell ; ++cellz)
    {
      for(uint celly = 0; celly < ncell ; ++celly)
        {
          for(uint cellx = 0; cellx < ncell ; ++cellx)
            {
              imageno[index].x = cellx;
              imageno[index].y = celly;
              imageno[index].z = cellz;  
              index = index + 1 ;
            }
        }
    }
}



void dist(vector<Molecular> &mol, uint idi, uint idj, const vector<float> & L,  vector<float> & PB_L, vector<Vector_int> & imageno, uint index,  float &x, float & y, float &z )
{
  x = min_distance(mol[idj].x - mol[idi].x, L[0]);
  y = min_distance(mol[idj].y - mol[idi].y, L[1]);
  z = min_distance(mol[idj].z - mol[idi].z, L[2]);

  x = x + imageno[index].x * L[0];
  y = y + imageno[index].y * L[1];
  z = z + imageno[index].z * L[2];

  if(x > 0 && abs(x) > (PB_L[0] * 0.5) ) { x  -= PB_L[0] ; }
  if(y > 0 && abs(y) > (PB_L[1] * 0.5) ) { y  -= PB_L[1] ; }
  if(z > 0 && abs(z) > (PB_L[2] * 0.5) ) { z  -= PB_L[2] ; } 

  if(x < 0 && abs(x) > (PB_L[0] * 0.5) ) { x  += PB_L[0] ; }
  if(y < 0 && abs(y) > (PB_L[1] * 0.5) ) { y  += PB_L[1] ; }
  if(z < 0 && abs(z) > (PB_L[2] * 0.5) ) { z  += PB_L[2] ; }
}




