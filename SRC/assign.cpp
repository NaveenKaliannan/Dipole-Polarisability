
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
#include "../include/assign.h"
#include "../include/pbc.h"
#include "../include/velocity.h"


using namespace std;

void Assigncoordinationforpurewater(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  float x = 0, y = 0, z = 0, rij = 0;
  uint idi, idi1, idi2, idj, idj1, idj2, donor_hbond = 0, acceptor_hbond = 0;

  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          idi = natoms*t+i; 
          idi1 = natoms*t+i+1;          
          idi2 = natoms*t+i+2;    
          donor_hbond = 0; 
          acceptor_hbond = 0;
                
          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H') 
            {
              for(uint j = 0;j < natoms;++j)
                { 
                  idj = natoms*t+j; 
                  idj1 = natoms*t+j+1;          
                  idj2 = natoms*t+j+2;                  
                  if(idi != idj  && r[idj].symbol[0] == 'O' && r[idj+1].symbol[0] == 'H' && r[idj+2].symbol[0] == 'H')
                    {
                      x = min_distance(r[idj].x - r[idi].x, L[0]);
                      y = min_distance(r[idj].y - r[idi].y, L[1]);
                      z = min_distance(r[idj].z - r[idi].z, L[2]); 
                      rij = mindis(x,y,z,L); 

                      if( rij < 3.5  && rij > 0 && ( angle_btwn_3points(r,idi,idi1,idj, L) < 30 || angle_btwn_3points(r,idi,idi2,idj, L) < 30) )
                        {
                          donor_hbond += 1;

                        }      
                      if( rij < 3.5  && rij > 0 && ( angle_btwn_3points(r,idj,idj1,idi, L) < 30 || angle_btwn_3points(r,idj,idj2,idi, L) < 30) )
                        {
                          acceptor_hbond += 1;
                        } 
                    }
                }           
              r[idi].totalhbonds = donor_hbond + acceptor_hbond; 
              r[idi].totaldonorhbonds = donor_hbond; 
              r[idi].totalacceptorhbonds = acceptor_hbond; 
            }
        }
    }
}


void Assigngammaforpurewater(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  float x = 0, y = 0, z = 0, rij = 0;
  uint idi, idi1, idi2, idj, idj1, idj2;
  float drij1 = 0, drij2 = 0, arij1 = 0, arij2 = 0;

  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          idi = natoms*t+i; 
          idi1 = natoms*t+i+1;          
          idi2 = natoms*t+i+2; 
          drij1 = 0, drij2 = 0, arij1 = 0, arij2 = 0;
                
          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H' && r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].totalacceptorhbonds == 2 ) 
            {
              for(uint j = 0;j < natoms;++j)
                { 
                  idj = natoms*t+j; 
                  idj1 = natoms*t+j+1;          
                  idj2 = natoms*t+j+2;                  
                  if(idi != idj  && r[idj].symbol[0] == 'O' && r[idj+1].symbol[0] == 'H' && r[idj+2].symbol[0] == 'H')
                    {
                      x = min_distance(r[idj].x - r[idi].x, L[0]);
                      y = min_distance(r[idj].y - r[idi].y, L[1]);
                      z = min_distance(r[idj].z - r[idi].z, L[2]); 
                      rij = mindis(x,y,z,L); 

                      if( rij < 3.5  && rij > 0 &&  angle_btwn_3points(r,idi,idi1,idj, L) < 30 )
                        {
                          drij1 = rij;
                        }    
                      if( rij < 3.5  && rij > 0 &&  angle_btwn_3points(r,idi,idi2,idj, L) < 30 )
                        {
                          drij2 = rij;
                        }  
  
                      if( rij < 3.5  && rij > 0 &&  (angle_btwn_3points(r,idj,idj1,idi, L) < 30 || angle_btwn_3points(r,idj,idj2,idi, L) < 30 ) )
                        {
                          if(arij1 == 0)
                            {
                              arij1 = rij;
                            }
                          else if(arij1 != 0)
                            {
                              arij2 = rij;
                            }
                        } 
                    }
                }           
              r[idi].gamma_d = max(drij1, drij2) / min(drij1, drij2); 
              r[idi].gamma_a = max(arij1, arij2) / min(arij1, arij2); 
              /*AIMD gamma ratio value for asymmetry and symmetry environment is 1.16 and 1.02*/
            }
        }
    }
}



void AssignAtomicMass(vector<Atom> &r, uint nsteps, uint natoms)
{
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          if((r[id].symbol[0] == 'N' || r[id].symbol[0] == 'n' ) && ( r[id].symbol[1] == 'A' || r[id].symbol[1] == 'a'))
            {
              r[id].atomicmass = 22.9897;
            }
          else if((r[id].symbol[0] == 'M' || r[id].symbol[0] == 'm' ) && ( r[id].symbol[1] == 'G' || r[id].symbol[1] == 'g'))
            {
              r[id].atomicmass = 24.305;
            }
          else if((r[id].symbol[0] == 'C' || r[id].symbol[0] == 'c' ) && ( r[id].symbol[1] == 'L' || r[id].symbol[1] == 'l'))
            {
              r[id].atomicmass = 35.453;
            }
          else if(r[id].symbol[0] == 'F' || r[id].symbol[0] == 'f' )
            {
              r[id].atomicmass = 18.99;
            }
          else if(r[id].symbol[0] == 'O' || r[id].symbol[0] == 'o' )
            {
              r[id].atomicmass = 15.999;
            }
          else if(r[id].symbol[0] == 'S' || r[id].symbol[0] == 's' )
            {
              r[id].atomicmass = 32.065;
            }
          else if(r[id].symbol[0] == 'H' || r[id].symbol[0] == 'h' )
            {
              r[id].atomicmass = 1.00784;
            }
        }
    }
}


void parameters(Molecular &mol)
{  
/*
   Cp2k Luber2014 Raman spectra from ab initio moleculardynamics and its application to liquid S-methyloxirane
        Putrino2002 , Anharmonic Raman Spectra in High-Pressure Ice fromAb InitioSimulations

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

   Munir S. Skaf and Sérgio M. Vechi J. Chem. Phys 119, 2181 (2003);


   Simulation of the Far Infrared Spectrum of Liquid Water and Steam Along the Coexistence Curve
   from the Triple Point to the Critical Point https://doi.org/10.1007/978-94-011-0183-7_10

   Pengyu Ren and Jay W. Ponder Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation
   https://pubs.acs.org/doi/pdf/10.1021/jp027815%2B


   // Thole model 
   IMOLECULAR POLARCIZABILITLES CALCULATED WITH A MOIXEIEID DIPOLE INTERACTiON
   AnAtomDipoleInteractionModelforMolecularOpticalProperties


   // Exteded DID approach torri
   Extended dipole-induced dipole mechanism for generatingRaman and optical Kerr effect intensities oflow-frequency dynamics in liquids

   // DID approach 
   Polarizability response in polar solvents:Molecular-dynamics simulations ofacetonitrile and chloroform
 
   grid radii : J Comput Chem. 2007 May ; 28(7): 1261–1274    Gaussian Induced Dipole Polarization Model

   Development of Polarizable Gaussian Model for MolecularMechanical Calculations I: Atomic Polarizability Parameterization ToReproduceab InitioAnisotropy

  //2018 Table of static dipole polarizabilities of theneutral elements in the periodic table

   The mulliken charges were used and parameterized with the help of CP2K program
*/

  // Atomic Polarisability [Angstrom3]
  if(mol.MOL[0] == 'H' && mol.MOL[1] == '2' && mol.MOL[2] == 'O')
    {
      mol.q = 0.0000;
      mol.PPol.xx = 1.3227;mol.PPol.yy =  1.1191;mol.PPol.zz = 0.8808;
      //mol.PPol.xy = 0.0018;mol.PPol.xz = -0.0003;mol.PPol.yz = 0.0006;
      //mol.PPol.yx = 0.0016;mol.PPol.zx = -0.0002;mol.PPol.zy = 0.0006;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
    }
  else if(mol.MOL[0] == 'H')
    {
      mol.q = 0.25983;
      mol.PPol.xx = 0.4960;mol.PPol.yy = 0.4960;mol.PPol.zz = 0.4960;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    =-0.03859;mol.PD.y    = 0.0000;mol.PD.z    =-0.05818;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
  else if(mol.MOL[0] == 'O')
    {
      mol.q = -0.51966;
      mol.PPol.xx = 0.8370;mol.PPol.yy = 0.8370;mol.PPol.zz = 0.8370;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.14279;

      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
  else if( ( mol.MOL[0] == 'M' || mol.MOL[0] == 'm' ) && ( mol.MOL[1] == 'G' || mol.MOL[1] == 'g' ))
    {
      mol.q = +2;
      mol.PPol.xx = 0.0800;mol.PPol.yy = 0.0800;mol.PPol.zz = 0.0800;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
  else if( ( mol.MOL[0] == 'C' || mol.MOL[0] == 'c' ) && ( mol.MOL[1] == 'L' || mol.MOL[1] == 'l' ))
    {
      mol.q = -1;
      mol.PPol.xx = 4.0000;mol.PPol.yy = 4.0000;mol.PPol.zz = 4.0000;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
  else if( ( mol.MOL[0] == 'N' || mol.MOL[0] == 'n' ) && ( mol.MOL[1] == 'A' || mol.MOL[1] == 'a' ))
    {
      mol.q = 1;
      mol.PPol.xx = 0.1200;mol.PPol.yy = 0.1200;mol.PPol.zz = 0.1200;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
  else if( ( mol.MOL[0] == 'F' || mol.MOL[0] == 'f' ))
    {
      mol.q = -1;
      mol.PPol.xx = 1.3500;mol.PPol.yy = 1.3500;mol.PPol.zz = 1.3500;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
}


// Transforming atomic to atomic assigns dipole moment, center of mass and permanent polarisabilites to the variables
void TransformAtomictoAtomic(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol)
{
  Molecular mols; 
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i;  
          mols.x = r[id].x;
          mols.y = r[id].y;
          mols.z = r[id].z;
          mols.MOL = r[id].symbol;
          mols.m = r[id].atomicmass; 
          parameters(mols);
          mols.PD.x *= amutodebye;
          mols.PD.y *= amutodebye;
          mols.PD.z *= amutodebye;
          mol.push_back(mols);         
        }
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
              mols.MOL = "Na";
              mols.m = r[id].atomicmass; 
              mols.sl = 0.3900 ; 
              mols.vdwr = 2.27 ;    
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
              mols.MOL = "Mg";
              mols.m = r[id].atomicmass; 
              mols.sl = 0.1150 ;
              mols.vdwr = 1.364 ;     
              parameters(mols);
              mol.push_back(mols); 
            }
          //CL
          if((r[id].symbol[0] == 'C' || r[id].symbol[0] == 'c' ) && ( r[id].symbol[1] == 'L' || r[id].symbol[1] == 'l'))
            {
              mols.x = r[id].x;
              mols.y = r[id].y;
              mols.z = r[id].z;
              mols.MOL = r[id].symbol;
              mols.MOL = "Cl";
              mols.m = r[id].atomicmass;
              mols.sl = 0.3900 ; 
              mols.vdwr = 1.639 ;    
              parameters(mols);
              mol.push_back(mols); 
            }
          //F
          if((r[id].symbol[0] == 'F' || r[id].symbol[0] == 'f' ))
            {
              mols.x = r[id].x;
              mols.y = r[id].y;
              mols.z = r[id].z;
              mols.MOL = r[id].symbol;
              mols.MOL = "F";
              mols.m = r[id].atomicmass;
              mols.sl = 0.3900 ; 
              mols.vdwr = 1.47 ;    
              parameters(mols);
              mol.push_back(mols); 
            }
          // H2O
          if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
            { 
              const double am_H = 1.00784 * amu, am_O = 15.999 * amu, am_H2O = (2 * am_H +  am_O);
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

              mols.sl = 0.39 ; 
              mols.vdwr = 1.7 ; 

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

              /*only for water*/
              mols.totalhbonds  = r[id].totalhbonds ; 
              mols.totaldonorhbonds  = r[id].totaldonorhbonds ; 
              mols.totalacceptorhbonds  = r[id].totalacceptorhbonds ; 
              mols.gamma_d  = r[id].gamma_d ; 
              mols.gamma_a  = r[id].gamma_a ; 

              mol.push_back(mols);
            }
        }
    } 
}







