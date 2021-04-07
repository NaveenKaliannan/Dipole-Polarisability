
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
      /*iso*/
      mol.hyperpol.xxx = -0.0000;
      mol.hyperpol.yyy = -0.0000;
      mol.hyperpol.zzz = -0.0000;    
      /*xxy*/
      mol.hyperpol.xxy = -0.0000;
      mol.hyperpol.yxx = -0.0000;
      mol.hyperpol.xyx = -0.0000;
      /*xxz*/
      mol.hyperpol.xxz = -0.0000;
      mol.hyperpol.zxx = -0.0000;
      mol.hyperpol.xzx = -0.0000;
      /*yyx*/
      mol.hyperpol.yyx = -0.0000;
      mol.hyperpol.xyy = -0.0000;
      mol.hyperpol.yxy = -0.0000;
      /*yyz*/
      mol.hyperpol.yyz = -0.0000;
      mol.hyperpol.zyy = -0.0000;
      mol.hyperpol.yzy = -0.0000;
      /*zzx*/
      mol.hyperpol.zzx = -0.0000;
      mol.hyperpol.xzz = -0.0000;
      mol.hyperpol.zxz = -0.0000;
      /*zzy*/
      mol.hyperpol.zzy = -0.0000;
      mol.hyperpol.yzz = -0.0000;
      mol.hyperpol.zyz = -0.0000;
      /*xyz*/
      mol.hyperpol.xyz = -0.0000;
      mol.hyperpol.xzy = -0.0000;
      mol.hyperpol.yxz = -0.0000;
      mol.hyperpol.yzx = -0.0000;
      mol.hyperpol.zxy = -0.0000;
      mol.hyperpol.zyx = -0.0000;

      mol.shyperpol.xxxx = 0.0000;
      mol.shyperpol.yyyy = 0.0000;
      mol.shyperpol.zzzz = 0.0000;
      mol.shyperpol.xxyy = 0.0000;
      mol.shyperpol.xxzz = 0.0000;
      mol.shyperpol.yyxx = 0.0000;
      mol.shyperpol.yyzz = 0.0000;
      mol.shyperpol.zzxx = 0.0000;
      mol.shyperpol.zzyy = 0.0000;

  // Atomic Polarisability [Angstrom3]
  if(mol.MOL[0] == 'H' && mol.MOL[1] == '2' && mol.MOL[2] == 'O')
    {
      mol.q = 0.0000;
      mol.PPol.xx = 1.3227;mol.PPol.yy =  1.1191;mol.PPol.zz = 0.8808;
      mol.PPol.xy = 0.0018;mol.PPol.xz = -0.0003;mol.PPol.yz = 0.0006;
      mol.PPol.yx = 0.0016;mol.PPol.zx = -0.0002;mol.PPol.zy = 0.0006;

      //mol.PPol.xx = 1.244366;mol.PPol.yy =  1.067983;mol.PPol.zz = 0.817208;
      //mol.PPol.xy = 0.0000;mol.PPol.xz = -0.0000;mol.PPol.yz = 0.0000;
      //mol.PPol.yx = 0.0000;mol.PPol.zx = -0.0000;mol.PPol.zy = 0.0000;


      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 

      /*iso atomic unit -12.730*/
      mol.hyperpol.xxx = -0.0000;
      mol.hyperpol.yyy = -0.5282; 
      mol.hyperpol.zzz = -0.0000;    
      /*xxy  atomic unit -16.445 */
      mol.hyperpol.xxy = -0.6824; 
      mol.hyperpol.yxx = -0.6824;
      mol.hyperpol.xyx = -0.6824;
      /*xxz*/
      mol.hyperpol.xxz = -0.0000;
      mol.hyperpol.zxx = -0.0000;
      mol.hyperpol.xzx = -0.0000;
      /*yyx*/
      mol.hyperpol.yyx = -0.0000;
      mol.hyperpol.xyy = -0.0000;
      mol.hyperpol.yxy = -0.0000;
      /*yyz*/
      mol.hyperpol.yyz = -0.0000;
      mol.hyperpol.zyy = -0.0000;
      mol.hyperpol.yzy = -0.0000;
      /*zzx*/
      mol.hyperpol.zzx = -0.0000;
      mol.hyperpol.xzz = -0.0000;
      mol.hyperpol.zxz = -0.0000;
      /*zzy  atomic unit -5.5365 */
      mol.hyperpol.zzy = -0.2297;
      mol.hyperpol.yzz = -0.2297;
      mol.hyperpol.zyz = -0.2297;
      /*xyz*/
      mol.hyperpol.xyz = -0.0000;
      mol.hyperpol.xzy = -0.0000;
      mol.hyperpol.yxz = -0.0000;
      mol.hyperpol.yzx = -0.0000;
      mol.hyperpol.zxy = -0.0000;
      mol.hyperpol.zyx = -0.0000; 

      mol.shyperpol.xxxx = 269.8054*0.01162004;//3.5 times experimetal values
      mol.shyperpol.yyyy = 126.5123*0.01162004;
      mol.shyperpol.zzzz = 8.096200*0.01162004;
      mol.shyperpol.xxyy = 180.8921*0.01162004;
      mol.shyperpol.xxzz = 93.09720*0.01162004;
      mol.shyperpol.yyxx = 180.8921*0.01162004;
      mol.shyperpol.yyzz = 35.49410*0.01162004;
      mol.shyperpol.zzxx = 93.09720*0.01162004;
      mol.shyperpol.zzyy = 35.49410*0.01162004;
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
      //mol.PPol.xx = 0.041140;mol.PPol.yy = 0.041140;mol.PPol.zz = 0.041140;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;

      mol.shyperpol.xxxx = 0.9305*0.01162004;
      mol.shyperpol.yyyy = 0.9305*0.01162004;
      mol.shyperpol.zzzz = 0.9305*0.01162004;
      mol.shyperpol.xxyy = 0.3102*0.01162004;
      mol.shyperpol.xxzz = 0.3102*0.01162004;
      mol.shyperpol.yyxx = 0.3102*0.01162004;
      mol.shyperpol.yyzz = 0.3102*0.01162004;
      mol.shyperpol.zzxx = 0.3102*0.01162004;
      mol.shyperpol.zzyy = 0.3102*0.01162004;

    }
  else if( ( mol.MOL[0] == 'C' || mol.MOL[0] == 'c' ) && ( mol.MOL[1] == 'L' || mol.MOL[1] == 'l' ))
    {
      mol.q = -1;
      mol.PPol.xx = 4.0000;mol.PPol.yy = 4.0000;mol.PPol.zz = 4.0000;
      //mol.PPol.xx = 1.590515;mol.PPol.yy = 1.590515;mol.PPol.zz = 1.590515;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;

      mol.shyperpol.xxxx = -0.0661*0.01162004;
      mol.shyperpol.yyyy = -0.0661*0.01162004;
      mol.shyperpol.zzzz = -0.0661*0.01162004;
      mol.shyperpol.xxyy = -0.0210*0.01162004;
      mol.shyperpol.xxzz = -0.0210*0.01162004;
      mol.shyperpol.yyxx = -0.0210*0.01162004;
      mol.shyperpol.yyzz = -0.0210*0.01162004;
      mol.shyperpol.zzxx = -0.0210*0.01162004;
      mol.shyperpol.zzyy = -0.0210*0.01162004;

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
          //SO4
          if(r[id].symbol[0] == 'S' && r[id+1].symbol[0] == 'O' && r[id+2].symbol[0] == 'O' && r[id+3].symbol[0] == 'O')
            {
              const float am_S = 32.065 * amu, am_O = 16 * amu, am_SO4 = 96.06 * amu ;
              mols.x = ( r[id-1].x * am_O + r[id].x * am_S + r[id+1].x * am_O + r[id+2].x * am_O + r[id+3].x * am_O ) / am_SO4 ;
              mols.y = ( r[id-1].y * am_O + r[id].y * am_S + r[id+1].y * am_O + r[id+2].y * am_O + r[id+3].y * am_O ) / am_SO4 ;
              mols.z = ( r[id-1].z * am_O + r[id].z * am_S + r[id+1].z * am_O + r[id+2].z * am_O + r[id+3].z * am_O ) / am_SO4 ;              
              mols.MOL = "SO4";
              mols.m = 96.06;
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

               /*input format for DALTON 
                cout << "BASIS \ncc-pVTZ\nH2O\n\n    2    0\n        8.    1 " << endl;
                cout << "O " << (r[id].x - mols.x)  * 1.89 << "  " <<  (r[id].y - mols.y)  * 1.89 << "  " << (r[id].z - mols.z)  * 1.89 << endl;
                cout << "        1.    2"<< endl;
                cout << "H1 " << (r[id+1].x - mols.x) * 1.89  << "  " <<  (r[id+1].y - mols.y)  * 1.89  << "  " << (r[id+1].z - mols.z) * 1.89  << endl;
                cout << "H2 " << (r[id+2].x - mols.x) * 1.89  << "  " <<  (r[id+2].y - mols.y)  * 1.89 << "  " << (r[id+2].z - mols.z) * 1.89 << endl;
             */

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

              // Permanent polarisability [Angstrom3]
              parameters(mols);
              float A_xx = mols.PPol.xx, A_yy = mols.PPol.yy, A_zz = mols.PPol.zz,
                    A_xy = mols.PPol.xy, A_xz = mols.PPol.xz, A_yz = mols.PPol.yz,
                    A_yx = mols.PPol.yx, A_zx = mols.PPol.zx, A_zy = mols.PPol.zy;


//------------------------------  variable polarizability model-------------------------//
              uint hbond_cation = 0, hbond_anion = 0;   
           float rij1 = 0, rij2 = 0, rij = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*t+i;  
              uint idi1 = natoms*t+i+1;          
              uint idi2 = natoms*t+i+2;                  
              uint idj = natoms*t+j;
 
              /*first solvation shell*/
              if(r[j].symbol[0] == 'M' || r[j].symbol[0] == 'N')
                {
                  x = min_distance(r[idj].x - r[idi].x, L[0]);
                  y = min_distance(r[idj].y - r[idi].y, L[1]);
                  z = min_distance(r[idj].z - r[idi].z, L[2]); 
                  rij = mindis(x,y,z,L); 
                  if(rij < 3.2 && rij > 0)
                    {
                      hbond_cation += 1;
                    }      
                }
              if(r[j].symbol[0] == 'C' || r[j].symbol[0] == 'F')
                {
                  x = min_distance(r[idj].x - r[idi1].x, L[0]);
                  y = min_distance(r[idj].y - r[idi1].y, L[1]);
                  z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                  rij1 = mindis(x,y,z,L);

                  x = min_distance(r[idj].x - r[idi2].x, L[0]);
                  y = min_distance(r[idj].y - r[idi2].y, L[1]);
                  z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                  rij2 = mindis(x,y,z,L); 
                  if( (rij1 < 3.0 && rij1 > 0) || (rij2 < 3.0 && rij2 > 0))
                    {
                      hbond_anion += 1;
                    }      
                }
            } 

          if(hbond_cation > 0 && hbond_anion == 0)
            {
              A_xx = 1.39182, A_yy = 1.16592, A_zz = 0.9145,
              A_xy = 0.00176, A_xz =-0.00019, A_yz = 0.0002,
              A_yx = 0.00163, A_zx =-0.00038, A_zy = 0.0002;
            }
          else if(hbond_cation == 0 && hbond_anion > 0)
            {
              A_xx = 1.3680, A_yy = 1.1594, A_zz = 0.9138,
              A_xy = -0.0009, A_xz = -0.0002, A_yz = 0.0000,
              A_yx = -0.0007, A_zx = -0.0002, A_zy = 0.0000;
            }
          else if(hbond_cation > 0 && hbond_anion > 0)
            {
              A_xx = 1.3855, A_yy = 1.1604, A_zz = 0.9139,
              A_xy = 0.0050, A_xz = 0.0000, A_yz = -0.0001,
              A_yx = 0.0055, A_zx = 0.0000, A_zy = -0.0001;
            }
          else if(hbond_cation == 0 && hbond_anion == 0)
            {
              A_xx = 1.3749, A_yy = 1.1586, A_zz = 0.9138,
              A_xy = 0.0005, A_xz = 0.0000, A_yz = 0.0000,
              A_yx = 0.0005, A_zx = 0.0000, A_zy = 0.0000;
            }

//------------------------------  variable polarizability model-------------------------//
     

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

              // Permanent hyperpolarizability [Angstrom5]
              float xxx = mols.hyperpol.xxx,
                    yyy = mols.hyperpol.yyy,
                    zzz = mols.hyperpol.zzz,   

                    xxy = mols.hyperpol.xxy,
                    yxx = mols.hyperpol.yxx,
                    xyx = mols.hyperpol.xyx,

                    xxz = mols.hyperpol.xxz,
                    zxx = mols.hyperpol.zxx,
                    xzx = mols.hyperpol.xzx,

                    yyx = mols.hyperpol.yyx,
                    xyy = mols.hyperpol.xyy,
                    yxy = mols.hyperpol.yxy,
      
                    yyz = mols.hyperpol.yyz,
                    zyy = mols.hyperpol.zyy,
                    yzy = mols.hyperpol.yzy,
      
                    zzx = mols.hyperpol.zzx,
                    xzz = mols.hyperpol.xzz,
                    zxz = mols.hyperpol.zxz,
      
                    zzy = mols.hyperpol.zzy,
                    yzz = mols.hyperpol.yzz,
                    zyz = mols.hyperpol.zyz,
      
                    xyz = mols.hyperpol.xyz,
                    xzy = mols.hyperpol.xzy,
                    yxz = mols.hyperpol.yxz,
                    yzx = mols.hyperpol.yzx,
                    zxy = mols.hyperpol.zxy,
                    zyx = mols.hyperpol.zyx;  


              mols.hyperpol.xxx =  wb_x * wb_x * wb_x * yyy +
                                   HH_x * HH_x * wb_x * xxy +
                                   wb_x * HH_x * HH_x * yxx +
                                   HH_x * wb_x * HH_x * xyx +
                                   Pv_x * Pv_x * wb_x * zzy +
                                   wb_x * Pv_x * Pv_x * yzz +
                                   Pv_x * wb_x * Pv_x * zyz ;

              mols.hyperpol.yyy =  wb_y * wb_y * wb_y * yyy +
                                   HH_y * HH_y * wb_y * xxy +
                                   wb_y * HH_y * HH_y * yxx +
                                   HH_y * wb_y * HH_y * xyx +
                                   Pv_y * Pv_y * wb_y * zzy +
                                   wb_y * Pv_y * Pv_y * yzz +
                                   Pv_y * wb_y * Pv_y * zyz ;

              mols.hyperpol.zzz =  wb_z * wb_z * wb_z * yyy +
                                   HH_z * HH_z * wb_z * xxy +
                                   wb_z * HH_z * HH_z * yxx +
                                   HH_z * wb_z * HH_z * xyx +
                                   Pv_z * Pv_z * wb_z * zzy +
                                   wb_z * Pv_z * Pv_z * yzz +
                                   Pv_z * wb_z * Pv_z * zyz ;

              mols.hyperpol.xxy =  wb_x * wb_x * wb_y * yyy +
                                   HH_x * HH_x * wb_y * xxy +
                                   wb_x * HH_x * HH_y * yxx +
                                   HH_x * wb_x * HH_y * xyx +
                                   Pv_x * Pv_x * wb_y * zzy +
                                   wb_x * Pv_x * Pv_y * yzz +
                                   Pv_x * wb_x * Pv_y * zyz ;

              mols.hyperpol.yxx =  wb_y * wb_x * wb_x * yyy +
                                   HH_y * HH_x * wb_x * xxy +
                                   wb_y * HH_x * HH_x * yxx +
                                   HH_y * wb_x * HH_x * xyx +
                                   Pv_y * Pv_x * wb_x * zzy +
                                   wb_y * Pv_x * Pv_x * yzz +
                                   Pv_y * wb_x * Pv_x * zyz ;

              mols.hyperpol.xyx =  wb_x * wb_y * wb_x * yyy +
                                   HH_x * HH_y * wb_x * xxy +
                                   wb_x * HH_y * HH_x * yxx +
                                   HH_x * wb_y * HH_x * xyx +
                                   Pv_x * Pv_y * wb_x * zzy +
                                   wb_x * Pv_y * Pv_x * yzz +
                                   Pv_x * wb_y * Pv_x * zyz ;

              mols.hyperpol.xxz =  wb_x * wb_x * wb_z * yyy +
                                   HH_x * HH_x * wb_z * xxy +
                                   wb_x * HH_x * HH_z * yxx +
                                   HH_x * wb_x * HH_z * xyx +
                                   Pv_x * Pv_x * wb_z * zzy +
                                   wb_x * Pv_x * Pv_z * yzz +
                                   Pv_x * wb_x * Pv_z * zyz ;

              mols.hyperpol.zxx =  wb_z * wb_x * wb_x * yyy +
                                   HH_z * HH_x * wb_x * xxy +
                                   wb_z * HH_x * HH_x * yxx +
                                   HH_z * wb_x * HH_x * xyx +
                                   Pv_z * Pv_x * wb_x * zzy +
                                   wb_z * Pv_x * Pv_x * yzz +
                                   Pv_z * wb_x * Pv_x * zyz ;

              mols.hyperpol.xzx =  wb_x * wb_z * wb_x * yyy +
                                   HH_x * HH_z * wb_x * xxy +
                                   wb_x * HH_z * HH_x * yxx +
                                   HH_x * wb_z * HH_x * xyx +
                                   Pv_x * Pv_z * wb_x * zzy +
                                   wb_x * Pv_z * Pv_x * yzz +
                                   Pv_x * wb_z * Pv_x * zyz ;


              mols.hyperpol.yyx =  wb_y * wb_y * wb_x * yyy +
                                   HH_y * HH_y * wb_x * xxy +
                                   wb_y * HH_y * HH_x * yxx +
                                   HH_y * wb_y * HH_x * xyx +
                                   Pv_y * Pv_y * wb_x * zzy +
                                   wb_y * Pv_y * Pv_x * yzz +
                                   Pv_y * wb_y * Pv_x * zyz ;

              mols.hyperpol.xyy =  wb_x * wb_y * wb_y * yyy +
                                   HH_x * HH_y * wb_y * xxy +
                                   wb_x * HH_y * HH_y * yxx +
                                   HH_x * wb_y * HH_y * xyx +
                                   Pv_x * Pv_y * wb_y * zzy +
                                   wb_x * Pv_y * Pv_y * yzz +
                                   Pv_x * wb_y * Pv_y * zyz ;

              mols.hyperpol.yxy =  wb_y * wb_x * wb_y * yyy +
                                   HH_y * HH_x * wb_y * xxy +
                                   wb_y * HH_x * HH_y * yxx +
                                   HH_y * wb_x * HH_y * xyx +
                                   Pv_y * Pv_x * wb_y * zzy +
                                   wb_y * Pv_x * Pv_y * yzz +
                                   Pv_y * wb_x * Pv_y * zyz ;

              mols.hyperpol.yyz =  wb_y * wb_y * wb_z * yyy +
                                   HH_y * HH_y * wb_z * xxy +
                                   wb_y * HH_y * HH_z * yxx +
                                   HH_y * wb_y * HH_z * xyx +
                                   Pv_y * Pv_y * wb_z * zzy +
                                   wb_y * Pv_y * Pv_z * yzz +
                                   Pv_y * wb_y * Pv_z * zyz ;

              mols.hyperpol.zyy =  wb_z * wb_y * wb_y * yyy +
                                   HH_z * HH_y * wb_y * xxy +
                                   wb_z * HH_y * HH_y * yxx +
                                   HH_z * wb_y * HH_y * xyx +
                                   Pv_z * Pv_y * wb_y * zzy +
                                   wb_z * Pv_y * Pv_y * yzz +
                                   Pv_z * wb_y * Pv_y * zyz ;

              mols.hyperpol.yzy =  wb_y * wb_z * wb_y * yyy +
                                   HH_y * HH_z * wb_y * xxy +
                                   wb_y * HH_z * HH_y * yxx +
                                   HH_y * wb_z * HH_y * xyx +
                                   Pv_y * Pv_z * wb_y * zzy +
                                   wb_y * Pv_z * Pv_y * yzz +
                                   Pv_y * wb_z * Pv_y * zyz ;

              mols.hyperpol.zzx =  wb_z * wb_z * wb_x * yyy +
                                   HH_z * HH_z * wb_x * xxy +
                                   wb_z * HH_z * HH_x * yxx +
                                   HH_z * wb_z * HH_x * xyx +
                                   Pv_z * Pv_z * wb_x * zzy +
                                   wb_z * Pv_z * Pv_x * yzz +
                                   Pv_z * wb_z * Pv_x * zyz ;

              mols.hyperpol.xzz =  wb_x * wb_z * wb_z * yyy +
                                   HH_x * HH_z * wb_z * xxy +
                                   wb_x * HH_z * HH_z * yxx +
                                   HH_x * wb_z * HH_z * xyx +
                                   Pv_x * Pv_z * wb_z * zzy +
                                   wb_x * Pv_z * Pv_z * yzz +
                                   Pv_x * wb_z * Pv_z * zyz ;

              mols.hyperpol.zxz =  wb_z * wb_x * wb_z * yyy +
                                   HH_z * HH_x * wb_z * xxy +
                                   wb_z * HH_x * HH_z * yxx +
                                   HH_z * wb_x * HH_z * xyx +
                                   Pv_z * Pv_x * wb_z * zzy +
                                   wb_z * Pv_x * Pv_z * yzz +
                                   Pv_z * wb_x * Pv_z * zyz ;

              mols.hyperpol.zzy =  wb_z * wb_z * wb_y * yyy +
                                   HH_z * HH_z * wb_y * xxy +
                                   wb_z * HH_z * HH_y * yxx +
                                   HH_z * wb_z * HH_y * xyx +
                                   Pv_z * Pv_z * wb_y * zzy +
                                   wb_z * Pv_z * Pv_y * yzz +
                                   Pv_z * wb_z * Pv_y * zyz ;

              mols.hyperpol.yzz =  wb_y * wb_z * wb_z * yyy +
                                   HH_y * HH_z * wb_z * xxy +
                                   wb_y * HH_z * HH_z * yxx +
                                   HH_y * wb_z * HH_z * xyx +
                                   Pv_y * Pv_z * wb_z * zzy +
                                   wb_y * Pv_z * Pv_z * yzz +
                                   Pv_y * wb_z * Pv_z * zyz ;

              mols.hyperpol.zyz =  wb_z * wb_y * wb_z * yyy +
                                   HH_z * HH_y * wb_z * xxy +
                                   wb_z * HH_y * HH_z * yxx +
                                   HH_z * wb_y * HH_z * xyx +
                                   Pv_z * Pv_y * wb_z * zzy +
                                   wb_z * Pv_y * Pv_z * yzz +
                                   Pv_z * wb_y * Pv_z * zyz ;

              mols.hyperpol.xyz =  wb_x * wb_y * wb_z * yyy +
                                   HH_x * HH_y * wb_z * xxy +
                                   wb_x * HH_y * HH_z * yxx +
                                   HH_x * wb_y * HH_z * xyx +
                                   Pv_x * Pv_y * wb_z * zzy +
                                   wb_x * Pv_y * Pv_z * yzz +
                                   Pv_x * wb_y * Pv_z * zyz ;

              mols.hyperpol.xzy =  wb_x * wb_z * wb_y * yyy +
                                   HH_x * HH_z * wb_y * xxy +
                                   wb_x * HH_z * HH_y * yxx +
                                   HH_x * wb_z * HH_y * xyx +
                                   Pv_x * Pv_z * wb_y * zzy +
                                   wb_x * Pv_z * Pv_y * yzz +
                                   Pv_x * wb_z * Pv_y * zyz ;

              mols.hyperpol.yxz =  wb_y * wb_x * wb_z * yyy +
                                   HH_y * HH_x * wb_z * xxy +
                                   wb_y * HH_x * HH_z * yxx +
                                   HH_y * wb_x * HH_z * xyx +
                                   Pv_y * Pv_x * wb_z * zzy +
                                   wb_y * Pv_x * Pv_z * yzz +
                                   Pv_y * wb_x * Pv_z * zyz ;

              mols.hyperpol.yzx =  wb_y * wb_z * wb_x * yyy +
                                   HH_y * HH_z * wb_x * xxy +
                                   wb_y * HH_z * HH_x * yxx +
                                   HH_y * wb_z * HH_x * xyx +
                                   Pv_y * Pv_z * wb_x * zzy +
                                   wb_y * Pv_z * Pv_x * yzz +
                                   Pv_y * wb_z * Pv_x * zyz ;

              mols.hyperpol.zxy =  wb_z * wb_x * wb_y * yyy +
                                   HH_z * HH_x * wb_y * xxy +
                                   wb_z * HH_x * HH_y * yxx +
                                   HH_z * wb_x * HH_y * xyx +
                                   Pv_z * Pv_x * wb_y * zzy +
                                   wb_z * Pv_x * Pv_y * yzz +
                                   Pv_z * wb_x * Pv_y * zyz ;

              mols.hyperpol.zyx =  wb_z * wb_y * wb_x * yyy +
                                   HH_z * HH_y * wb_x * xxy +
                                   wb_z * HH_y * HH_x * yxx +
                                   HH_z * wb_y * HH_x * xyx +
                                   Pv_z * Pv_y * wb_x * zzy +
                                   wb_z * Pv_y * Pv_x * yzz +
                                   Pv_z * wb_y * Pv_x * zyz ;


              // Permanent second hyperpolarizability [Angstrom7]
              float xxxx = mols.shyperpol.xxxx,
                    yyyy = mols.shyperpol.yyyy,
                    zzzz = mols.shyperpol.zzzz,
                    xxyy = mols.shyperpol.xxyy,
                    xxzz = mols.shyperpol.xxzz,
                    yyxx = mols.shyperpol.yyxx,
                    yyzz = mols.shyperpol.yyzz,
                    zzxx = mols.shyperpol.zzxx,
                    zzyy = mols.shyperpol.zzyy;
 

              mols.shyperpol.xxxx = HH_x * HH_x * HH_x * HH_x * xxxx +
                                    wb_x * wb_x * wb_x * wb_x * yyyy +  
                                    Pv_x * Pv_x * Pv_x * Pv_x * zzzz +
                                    HH_x * HH_x * wb_x * wb_x * xxyy +
                                    HH_x * HH_x * Pv_x * Pv_x * xxzz +
                                    wb_x * wb_x * HH_x * HH_x * yyxx +
                                    wb_x * wb_x * Pv_x * Pv_x * yyzz +
                                    Pv_x * Pv_x * HH_x * HH_x * zzxx +
                                    Pv_x * Pv_x * wb_x * wb_x * zzyy +
                                    HH_x * wb_x * wb_x * HH_x * xxyy +
                                    HH_x * wb_x * HH_x * wb_x * xxyy +
                                    HH_x * Pv_x * Pv_x * HH_x * xxzz +
                                    HH_x * Pv_x * HH_x * Pv_x * xxzz +
                                    wb_x * HH_x * HH_x * wb_x * yyxx +  
                                    wb_x * HH_x * wb_x * HH_x * yyxx +  
                                    wb_x * Pv_x * Pv_x * wb_x * yyzz +  
                                    wb_x * Pv_x * wb_x * Pv_x * yyzz +  
                                    Pv_x * HH_x * HH_x * Pv_x * zzxx +
                                    Pv_x * HH_x * Pv_x * HH_x * zzxx +
                                    Pv_x * wb_x * wb_x * Pv_x * zzyy +
                                    Pv_x * wb_x * Pv_x * wb_x * zzyy ;


              mols.shyperpol.yyyy = HH_y * HH_y * HH_y * HH_y * xxxx +
                                    wb_y * wb_y * wb_y * wb_y * yyyy +  
                                    Pv_y * Pv_y * Pv_y * Pv_y * zzzz +
                                    HH_y * HH_y * wb_y * wb_y * xxyy +
                                    HH_y * HH_y * Pv_y * Pv_y * xxzz +
                                    wb_y * wb_y * HH_y * HH_y * yyxx +
                                    wb_y * wb_y * Pv_y * Pv_y * yyzz +
                                    Pv_y * Pv_y * HH_y * HH_y * zzxx +
                                    Pv_y * Pv_y * wb_y * wb_y * zzyy +
                                    HH_y * wb_y * wb_y * HH_y * xxyy +
                                    HH_y * wb_y * HH_y * wb_y * xxyy +
                                    HH_y * Pv_y * Pv_y * HH_y * xxzz +
                                    HH_y * Pv_y * HH_y * Pv_y * xxzz +
                                    wb_y * HH_y * HH_y * wb_y * yyxx +  
                                    wb_y * HH_y * wb_y * HH_y * yyxx +  
                                    wb_y * Pv_y * Pv_y * wb_y * yyzz +  
                                    wb_y * Pv_y * wb_y * Pv_y * yyzz +  
                                    Pv_y * HH_y * HH_y * Pv_y * zzxx +
                                    Pv_y * HH_y * Pv_y * HH_y * zzxx +
                                    Pv_y * wb_y * wb_y * Pv_y * zzyy +
                                    Pv_y * wb_y * Pv_y * wb_y * zzyy ;

              mols.shyperpol.zzzz = HH_z * HH_z * HH_z * HH_z * xxxx +
                                    wb_z * wb_z * wb_z * wb_z * yyyy +  
                                    Pv_z * Pv_z * Pv_z * Pv_z * zzzz +
                                    HH_z * HH_z * wb_z * wb_z * xxyy +
                                    HH_z * HH_z * Pv_z * Pv_z * xxzz +
                                    wb_z * wb_z * HH_z * HH_z * yyxx +
                                    wb_z * wb_z * Pv_z * Pv_z * yyzz +
                                    Pv_z * Pv_z * HH_z * HH_z * zzxx +
                                    Pv_z * Pv_z * wb_z * wb_z * zzyy +
                                    HH_z * wb_z * wb_z * HH_z * xxyy +
                                    HH_z * wb_z * HH_z * wb_z * xxyy +
                                    HH_z * Pv_z * Pv_z * HH_z * xxzz +
                                    HH_z * Pv_z * HH_z * Pv_z * xxzz +
                                    wb_z * HH_z * HH_z * wb_z * yyxx +  
                                    wb_z * HH_z * wb_z * HH_z * yyxx +  
                                    wb_z * Pv_z * Pv_z * wb_z * yyzz +  
                                    wb_z * Pv_z * wb_z * Pv_z * yyzz +  
                                    Pv_z * HH_z * HH_z * Pv_z * zzxx +
                                    Pv_z * HH_z * Pv_z * HH_z * zzxx +
                                    Pv_z * wb_z * wb_z * Pv_z * zzyy +
                                    Pv_z * wb_z * Pv_z * wb_z * zzyy ;


              mols.shyperpol.xxyy = HH_x * HH_x * HH_y * HH_y * xxxx +
                                    wb_x * wb_x * wb_y * wb_y * yyyy +  
                                    Pv_x * Pv_x * Pv_y * Pv_y * zzzz +
                                    HH_x * HH_x * wb_y * wb_y * xxyy +
                                    HH_x * HH_x * Pv_y * Pv_y * xxzz +
                                    wb_x * wb_x * HH_y * HH_y * yyxx +
                                    wb_x * wb_x * Pv_y * Pv_y * yyzz +
                                    Pv_x * Pv_x * HH_y * HH_y * zzxx +
                                    Pv_x * Pv_x * wb_y * wb_y * zzyy +
                                    HH_x * wb_x * wb_y * HH_y * xxyy +
                                    HH_x * wb_x * HH_y * wb_y * xxyy +
                                    HH_x * Pv_x * Pv_y * HH_y * xxzz +
                                    HH_x * Pv_x * HH_y * Pv_y * xxzz +
                                    wb_x * HH_x * HH_y * wb_y * yyxx +  
                                    wb_x * HH_x * wb_y * HH_y * yyxx +  
                                    wb_x * Pv_x * Pv_y * wb_y * yyzz +  
                                    wb_x * Pv_x * wb_y * Pv_y * yyzz +  
                                    Pv_x * HH_x * HH_y * Pv_y * zzxx +
                                    Pv_x * HH_x * Pv_y * HH_y * zzxx +
                                    Pv_x * wb_x * wb_y * Pv_y * zzyy +
                                    Pv_x * wb_x * Pv_y * wb_y * zzyy ;

              mols.shyperpol.xxzz = HH_x * HH_x * HH_z * HH_z * xxxx +
                                    wb_x * wb_x * wb_z * wb_z * yyyy +  
                                    Pv_x * Pv_x * Pv_z * Pv_z * zzzz +
                                    HH_x * HH_x * wb_z * wb_z * xxyy +
                                    HH_x * HH_x * Pv_z * Pv_z * xxzz +
                                    wb_x * wb_x * HH_z * HH_z * yyxx +
                                    wb_x * wb_x * Pv_z * Pv_z * yyzz +
                                    Pv_x * Pv_x * HH_z * HH_z * zzxx +
                                    Pv_x * Pv_x * wb_z * wb_z * zzyy +
                                    HH_x * wb_x * wb_z * HH_z * xxyy +
                                    HH_x * wb_x * HH_z * wb_z * xxyy +
                                    HH_x * Pv_x * Pv_z * HH_z * xxzz +
                                    HH_x * Pv_x * HH_z * Pv_z * xxzz +
                                    wb_x * HH_x * HH_z * wb_z * yyxx +  
                                    wb_x * HH_x * wb_z * HH_z * yyxx +  
                                    wb_x * Pv_x * Pv_z * wb_z * yyzz +  
                                    wb_x * Pv_x * wb_z * Pv_z * yyzz +  
                                    Pv_x * HH_x * HH_z * Pv_z * zzxx +
                                    Pv_x * HH_x * Pv_z * HH_z * zzxx +
                                    Pv_x * wb_x * wb_z * Pv_z * zzyy +
                                    Pv_x * wb_x * Pv_z * wb_z * zzyy ;

              mols.shyperpol.yyxx = HH_y * HH_y * HH_x * HH_x * xxxx +
                                    wb_y * wb_y * wb_x * wb_x * yyyy +  
                                    Pv_y * Pv_y * Pv_x * Pv_x * zzzz +
                                    HH_y * HH_y * wb_x * wb_x * xxyy +
                                    HH_y * HH_y * Pv_x * Pv_x * xxzz +
                                    wb_y * wb_y * HH_x * HH_x * yyxx +
                                    wb_y * wb_y * Pv_x * Pv_x * yyzz +
                                    Pv_y * Pv_y * HH_x * HH_x * zzxx +
                                    Pv_y * Pv_y * wb_x * wb_x * zzyy +
                                    HH_y * wb_y * wb_x * HH_x * xxyy +
                                    HH_y * wb_y * HH_x * wb_x * xxyy +
                                    HH_y * Pv_y * Pv_x * HH_x * xxzz +
                                    HH_y * Pv_y * HH_x * Pv_x * xxzz +
                                    wb_y * HH_y * HH_x * wb_x * yyxx +  
                                    wb_y * HH_y * wb_x * HH_x * yyxx +  
                                    wb_y * Pv_y * Pv_x * wb_x * yyzz +  
                                    wb_y * Pv_y * wb_x * Pv_x * yyzz +  
                                    Pv_y * HH_y * HH_x * Pv_x * zzxx +
                                    Pv_y * HH_y * Pv_x * HH_x * zzxx +
                                    Pv_y * wb_y * wb_x * Pv_x * zzyy +
                                    Pv_y * wb_y * Pv_x * wb_x * zzyy ;

              mols.shyperpol.yyzz = HH_y * HH_y * HH_z * HH_z * xxxx +
                                    wb_y * wb_y * wb_z * wb_z * yyyy +  
                                    Pv_y * Pv_y * Pv_z * Pv_z * zzzz +
                                    HH_y * HH_y * wb_z * wb_z * xxyy +
                                    HH_y * HH_y * Pv_z * Pv_z * xxzz +
                                    wb_y * wb_y * HH_z * HH_z * yyxx +
                                    wb_y * wb_y * Pv_z * Pv_z * yyzz +
                                    Pv_y * Pv_y * HH_z * HH_z * zzxx +
                                    Pv_y * Pv_y * wb_z * wb_z * zzyy +
                                    HH_y * wb_y * wb_z * HH_z * xxyy +
                                    HH_y * wb_y * HH_z * wb_z * xxyy +
                                    HH_y * Pv_y * Pv_z * HH_z * xxzz +
                                    HH_y * Pv_y * HH_z * Pv_z * xxzz +
                                    wb_y * HH_y * HH_z * wb_z * yyxx +  
                                    wb_y * HH_y * wb_z * HH_z * yyxx +  
                                    wb_y * Pv_y * Pv_z * wb_z * yyzz +  
                                    wb_y * Pv_y * wb_z * Pv_z * yyzz +  
                                    Pv_y * HH_y * HH_z * Pv_z * zzxx +
                                    Pv_y * HH_y * Pv_z * HH_z * zzxx +
                                    Pv_y * wb_y * wb_z * Pv_z * zzyy +
                                    Pv_y * wb_y * Pv_z * wb_z * zzyy ;

              mols.shyperpol.zzxx = HH_z * HH_z * HH_x * HH_x * xxxx +
                                    wb_z * wb_z * wb_x * wb_x * yyyy +  
                                    Pv_z * Pv_z * Pv_x * Pv_x * zzzz +
                                    HH_z * HH_z * wb_x * wb_x * xxyy +
                                    HH_z * HH_z * Pv_x * Pv_x * xxzz +
                                    wb_z * wb_z * HH_x * HH_x * yyxx +
                                    wb_z * wb_z * Pv_x * Pv_x * yyzz +
                                    Pv_z * Pv_z * HH_x * HH_x * zzxx +
                                    Pv_z * Pv_z * wb_x * wb_x * zzyy +
                                    HH_z * wb_z * wb_x * HH_x * xxyy +
                                    HH_z * wb_z * HH_x * wb_x * xxyy +
                                    HH_z * Pv_z * Pv_x * HH_x * xxzz +
                                    HH_z * Pv_z * HH_x * Pv_x * xxzz +
                                    wb_z * HH_z * HH_x * wb_x * yyxx +  
                                    wb_z * HH_z * wb_x * HH_x * yyxx +  
                                    wb_z * Pv_z * Pv_x * wb_x * yyzz +  
                                    wb_z * Pv_z * wb_x * Pv_x * yyzz +  
                                    Pv_z * HH_z * HH_x * Pv_x * zzxx +
                                    Pv_z * HH_z * Pv_x * HH_x * zzxx +
                                    Pv_z * wb_z * wb_x * Pv_x * zzyy +
                                    Pv_z * wb_z * Pv_x * wb_x * zzyy ;


              mols.shyperpol.zzyy = HH_z * HH_z * HH_y * HH_y * xxxx +
                                    wb_z * wb_z * wb_y * wb_y * yyyy +  
                                    Pv_z * Pv_z * Pv_y * Pv_y * zzzz +
                                    HH_z * HH_z * wb_y * wb_y * xxyy +
                                    HH_z * HH_z * Pv_y * Pv_y * xxzz +
                                    wb_z * wb_z * HH_y * HH_y * yyxx +
                                    wb_z * wb_z * Pv_y * Pv_y * yyzz +
                                    Pv_z * Pv_z * HH_y * HH_y * zzxx +
                                    Pv_z * Pv_z * wb_y * wb_y * zzyy +
                                    HH_z * wb_z * wb_y * HH_y * xxyy +
                                    HH_z * wb_z * HH_y * wb_y * xxyy +
                                    HH_z * Pv_z * Pv_y * HH_y * xxzz +
                                    HH_z * Pv_z * HH_y * Pv_y * xxzz +
                                    wb_z * HH_z * HH_y * wb_y * yyxx +  
                                    wb_z * HH_z * wb_y * HH_y * yyxx +  
                                    wb_z * Pv_z * Pv_y * wb_y * yyzz +  
                                    wb_z * Pv_z * wb_y * Pv_y * yyzz +  
                                    Pv_z * HH_z * HH_y * Pv_y * zzxx +
                                    Pv_z * HH_z * Pv_y * HH_y * zzxx +
                                    Pv_z * wb_z * wb_y * Pv_y * zzyy +
                                    Pv_z * wb_z * Pv_y * wb_y * zzyy ;

             // cout << t+1 << " " << mols.shyperpol.xxxx << " " << mols.shyperpol.yyyy << "  " << mols.shyperpol.zzzz << "  " << mols.shyperpol.xxyy << "  " << mols.shyperpol.xxzz << "  " << mols.shyperpol.yyxx << "  " << mols.shyperpol.yyzz << "  " << mols.shyperpol.zzxx << "  " << mols.shyperpol.zzyy  << "  "  <<  endl;


              //cout << t << "\n xxxx =  " << mols.shyperpol.xxxx << "\n yyyy =  " << mols.shyperpol.yyyy << "\n zzzz =  " << mols.shyperpol.zzzz << "\n xxyy =  " << mols.shyperpol.xxyy << "\n xxzz =  " << mols.shyperpol.xxzz << "\n yyxx =  " << mols.shyperpol.yyxx << "\n yyzz =  " << mols.shyperpol.yyzz << "\n zzxx =  " << mols.shyperpol.zzxx << "\n zzyy =  " << mols.shyperpol.zzyy  << "  "  <<  endl;
              


              /*
              cout << " xxx = " <<  mols.hyperpol.xxx <<  " yyy = " <<  mols.hyperpol.yyy << " zzz = " <<  mols.hyperpol.zzz <<  "\n" 
                   << " xxy = " <<  mols.hyperpol.xxy <<  " yxx = " <<  mols.hyperpol.yxx << " xyx = " <<  mols.hyperpol.xyx <<  "\n"
                   << " xxz = " <<  mols.hyperpol.xxz <<  " zxx = " <<  mols.hyperpol.zxx << " xzx = " <<  mols.hyperpol.xzx <<  "\n"
                   << " yyx = " <<  mols.hyperpol.yyx <<  " xyy = " <<  mols.hyperpol.xyy << " yxy = " <<  mols.hyperpol.yxy <<  "\n"
                   << " yyz = " <<  mols.hyperpol.yyz <<  " zyy = " <<  mols.hyperpol.zyy << " yzy = " <<  mols.hyperpol.yzy <<  "\n"
                   << " zzx = " <<  mols.hyperpol.zzx <<  " xzz = " <<  mols.hyperpol.xzz << " zxz = " <<  mols.hyperpol.zxz <<  "\n"
                   << " zzy = " <<  mols.hyperpol.zzy <<  " yzz = " <<  mols.hyperpol.yzz << " zyz = " <<  mols.hyperpol.zyz <<  "\n"
                   << " xyz = " <<  mols.hyperpol.xyz <<  " xzy = " <<  mols.hyperpol.xzy << " yxz = " <<  mols.hyperpol.yxz <<  "\n"
                   << " yzx = " <<  mols.hyperpol.yzx <<  " zxy = " <<  mols.hyperpol.zxy << " zyx = " <<  mols.hyperpol.zyx <<  "\n" << endl;
              
              cout << t << " " << mols.hyperpol.xxx << " " << mols.hyperpol.yyy << " " << mols.hyperpol.zzz << " " << mols.hyperpol.xxy << " " << mols.hyperpol.xxz << "  " <<
                           " " << mols.hyperpol.yyx << " " << mols.hyperpol.yyz << " " << mols.hyperpol.zzx << " " << mols.hyperpol.zzy << " " << mols.hyperpol.xyz << "  " << mols.hyperpol.yzx <<  endl;

              */

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







