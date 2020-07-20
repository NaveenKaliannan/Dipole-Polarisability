#include "../include/io.h"
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



using namespace std;


//Reading the xyz trajectory
void readtrajectory(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, float L)
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


// computes velocity with central difference scheme as CP2K (multiply with vAperfmstoamu converts to atomic unit )
void computevelocity(vector<Atom> &r, uint nsteps, uint natoms, float L, float dt)
{
  for(uint t = 1; t < nsteps-1; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*(t)+i; 
          uint id1 = natoms*(t-1)+i; 
          uint id2 = natoms*(t+1)+i; 
          r[id].vx = min_distance((r[id2].x - r[id1].x), L)/(2*dt) ;
          r[id].vy = min_distance((r[id2].y - r[id1].y), L)/(2*dt) ;
          r[id].vz = min_distance((r[id2].z - r[id1].z), L)/(2*dt) ;
        }
    }
    //if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition 
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


void Print(vector<Atom> &r, uint nsteps, uint natoms, float L, vector<Molecular> &mol, uint nmol, string filename, string TYPE)
{
  ofstream outfile(filename);
  if(TYPE[0] == 'A' && TYPE[1] == 'T' && TYPE[2] == 'M')
    for(uint t = 0; t < nsteps; ++t )
      { 
        outfile << natoms << endl ;
        outfile << "BOX Length " << L << "  " << L << "  " << L << "\n";
        for(uint i = 0;i < natoms;++i)
          {
            uint id = natoms*t+i; 
            outfile << r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << endl;
          }
      }
  else if(TYPE[0] == 'M' && TYPE[1] == 'O' && TYPE[2] == 'L')
    for(uint t = 0; t < nsteps; ++t )
      { 
        outfile << nmol << endl ;
        outfile << "BOX Length " << L << "  " << L << "  " << L << "\n";
        for(uint i = 0;i < nmol;++i)
          {
            uint id = nmol*t+i; 
            outfile << mol[id].MOL << " " << mol[id].x << "  " << mol[id].y << "  " << mol[id].z << endl;
          }
      }
  else if(TYPE[0] == 'D' && TYPE[1] == 'I' && TYPE[2] == 'P')
    for(uint t = 0; t < nsteps; ++t )
      { 
        outfile << nmol << endl ;
        outfile << "## Molecular name [string], Dipole moment [Debye] in x y z, Polarisability tensor [Angstrom] xx yy zz " << endl;
        for(uint i = 0;i < nmol;++i)
          {
            uint id = nmol*t+i; 
            outfile << mol[id].MOL <<  "  " << mol[id].PD_x << "  " << mol[id].PD_y << "  " << mol[id].PD_z << "  "<< mol[id].Pol_xx << "  " << mol[id].Pol_yy << "  " << mol[id].Pol_zz << endl;
          }
      }
  outfile.close();
  outfile.clear();
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

//Broken bonds due to PBC are solved here, Minimum image convention was applied
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, float L)
{
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          if(r[id].x > L || r[id].x < 0 ) { r[id].x  = min_distance(r[id].x, L); }
          if(r[id].y > L || r[id].y < 0 ) { r[id].y  = min_distance(r[id].y, L); }
          if(r[id].z > L || r[id].z < 0 ) { r[id].z  = min_distance(r[id].z, L); }

          if(r[id].x > L ) { r[id].x  -= L ; }
          if(r[id].y > L ) { r[id].y  -= L ; }
          if(r[id].z > L ) { r[id].z  -= L ; }

          if(r[id].x < 0 ) { r[id].x += L ; }
          if(r[id].y < 0 ) { r[id].y += L ; }
          if(r[id].z < 0 ) { r[id].z += L ; }
        }
    }

  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i;
 
          // Unites SO4
          if(r[id].symbol[0] == 'S')
            {
              for(uint j = 0;j < natoms;++j)
                {
                  uint id2 = natoms*t+j; 
                  if(i == j){}
                  else if(r[id2].symbol[0] == 'O')
                    {
                      float rij = mindis(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z, L);
                      if(rij < 2.25)
                        { 
                          float rij2 = norm(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z);
                          if(abs(rij2 -  rij) > 0.1)
                            {  
                             float rij3 = norm(r[id].x - r[id2].x, 0, 0);
                             if(abs(rij3) > 2.5)
                               {   
                                 if(r[id].x > r[id2].x ) { r[id2].x += L; }    
                                 else if(r[id].x < r[id2].x ) { r[id2].x -= L; }                                                                  
                               }
                             rij3 = norm(0, r[id].y - r[id2].y, 0);
                             if(abs(rij3) > 2.5)
                               {    
                                 if(r[id].y > r[id2].y ) { r[id2].y += L; }    
                                 else if(r[id].y < r[id2].y ) { r[id2].y -= L; }                                                                                                    
                               }  
                             rij3 = norm(0, 0, r[id].z - r[id2].z);
                             if(abs(rij3) > 2.5)
                               {                           
                                 if(r[id].z > r[id2].z ) { r[id2].z += L; }    
                                 else if(r[id].z < r[id2].z ) { r[id2].z -= L; }                                                                                                        
                               }                                      
                            }
                        }
                    }
                } 
            }


          // Unites H2O
          if(r[id].symbol[0] == 'O')
            {
              for(uint j = 0;j < natoms;++j)
                {
                  uint id2 = natoms*t+j; 
                  if(i == j){}
                  else if(r[id2].symbol[0] == 'H')
                    {
                      float rij = mindis(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z, L);
                      if(rij < 1.25)
                        { 
                          float rij2 = norm(r[id].x - r[id2].x, r[id].y - r[id2].y, r[id].z - r[id2].z);
                          if(abs(rij2 -  rij) > 0.1)
                            {  
                             float rij3 = norm(r[id].x - r[id2].x, 0, 0);
                             if(abs(rij3) > 2.5)
                               {   
                                 if(r[id].x > r[id2].x ) { r[id2].x += L; }    
                                 else if(r[id].x < r[id2].x ) { r[id2].x -= L; }                                                                  
                               }
                             rij3 = norm(0, r[id].y - r[id2].y, 0);
                             if(abs(rij3) > 2.5)
                               {    
                                 if(r[id].y > r[id2].y ) { r[id2].y  =  r[id2].y += L; }    
                                 else if(r[id].y < r[id2].y ) { r[id2].y -= L; }                                                                                                    
                               }  
                             rij3 = norm(0, 0, r[id].z - r[id2].z);
                             if(abs(rij3) > 2.5)
                               {                           
                                 if(r[id].z > r[id2].z ) { r[id2].z += L; }    
                                 else if(r[id].z < r[id2].z ) { r[id2].z -= L; }                                                                                                        
                               }                                      
                            }
                        }
                    }
                } 
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
      mol.Pol_xx = 1.3227;mol.Pol_yy =  1.1191;mol.Pol_zz = 0.8808;
      mol.Pol_xy = 0.0018;mol.Pol_xz = -0.0003;mol.Pol_yz = 0.0006;
      mol.Pol_yx = 0.0016;mol.Pol_zx = -0.0002;mol.Pol_zy = 0.0006; 
    }
  else if( ( mol.MOL[0] == 'M' || mol.MOL[0] == 'm' ) && ( mol.MOL[1] == 'G' || mol.MOL[1] == 'g' ))
    {
      mol.q = 0;
      mol.Pol_xx = 10.090;mol.Pol_yy = 10.090;mol.Pol_zz = 10.090;
      mol.Pol_xy = 0.0000;mol.Pol_xz = 0.0000;mol.Pol_yz = 0.0000;
      mol.Pol_yx = 0.0000;mol.Pol_zx = 0.0000;mol.Pol_zy = 0.0000;
    }
  else if( ( mol.MOL[0] == 'C' || mol.MOL[0] == 'c' ) && ( mol.MOL[1] == 'L' || mol.MOL[1] == 'l' ))
    {
      mol.q = 0;
      mol.Pol_xx = 1.4300;mol.Pol_yy = 1.4300;mol.Pol_zz = 1.4300;
      mol.Pol_xy = 0.0087;mol.Pol_xz =-0.0080;mol.Pol_yz = 0.0083;
      mol.Pol_yx = 0.0087;mol.Pol_zx =-0.0080;mol.Pol_zy = 0.0083;
    }
  else if( ( mol.MOL[0] == 'N' || mol.MOL[0] == 'n' ) && ( mol.MOL[1] == 'A' || mol.MOL[1] == 'a' ))
    {
      mol.q = 0;
      mol.Pol_xx = 22.050;mol.Pol_yy = 22.050;mol.Pol_zz = 22.050;
      mol.Pol_xy = 0.0000;mol.Pol_xz = 0.0000;mol.Pol_yz = 0.0000;
      mol.Pol_yx = 0.0000;mol.Pol_zx = 0.0000;mol.Pol_zy = 0.0000;
    }
  else if(mol.MOL[0] == 'S' && mol.MOL[1] == 'O' && mol.MOL[2] == '4')
    {
      mol.q = 0;
      mol.Pol_xx = 6.104;mol.Pol_yy = 6.0861;mol.Pol_zz = 6.1769;
      mol.Pol_xy =-0.013;mol.Pol_xz = 0.0011;mol.Pol_yz = 0.0808;
      mol.Pol_yx =-0.014;mol.Pol_zx = 0.0006;mol.Pol_zy = 0.0806; 
    }
}

void TransformAtomictoMolecular(vector<Atom> &r, uint nsteps,  uint natoms, float L, vector<Molecular> &mol, uint nmol)
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
              mols.PD_x = 0 ;
              mols.PD_y = 0 ;
              mols.PD_z = 0 ;
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
              mols.PD_x = 0 ;
              mols.PD_y = 0 ;
              mols.PD_z = 0 ;
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
              mols.PD_x = 0 ;
              mols.PD_y = 0 ;
              mols.PD_z = 0 ;
              parameters(mols);
              mol.push_back(mols);
            }
          // H2O
          if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
            {
              const float am_H = 1.00784 * amu, am_O = 15.999 * amu, am_H2O = 18.015 * amu ;
              float wb_x = 0,wb_y = 0,wb_z = 0;           // bisector vector
              float HH_x = 0,HH_y = 0,HH_z = 0 ;           // H-H vector
              float Pv_x = 0,Pv_y = 0,Pv_z = 0;           // vector Perpendicular to bisector and H-H vector
              mols.x = ( r[id].x * am_O + r[id+1].x * am_H + r[id+2].x * am_H ) / am_H2O ;
              mols.y = ( r[id].y * am_O + r[id+1].y * am_H + r[id+2].y * am_H ) / am_H2O ;
              mols.z = ( r[id].z * am_O + r[id+1].z * am_H + r[id+2].z * am_H ) / am_H2O ;

              // water bisector vector, Minimum image convention was applied
              wb_x = min_distance(r[id+1].x - r[id].x, L) + min_distance(r[id+2].x - r[id].x, L) ;
              wb_y = min_distance(r[id+1].y - r[id].y, L) + min_distance(r[id+2].y - r[id].y, L) ;
              wb_z = min_distance(r[id+1].z - r[id].z, L) + min_distance(r[id+2].z - r[id].z, L) ;
              float sum  =  pow(wb_x * wb_x  + wb_y * wb_y + wb_z * wb_z, 0.5);
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
        
              /* converting the water bisector into permanent dipole [debye]
                  3
               i =     9941, time =     3976.400, E =       -17.1856249667
               O        13.4624293660       17.2132760672        4.3695246899
               H        12.6258911614       16.7156767906        4.4333223805
               H        13.6991295972       17.4314209339        5.2307556387

               Dipole moment [Debye]
               X=   -1.02713357 Y=   -0.46771379 Z=    1.63613604     Total=      1.98763697
              // projecting the water dipole [debye] in water bisector direction using projection matrix*/
              //float d_x = -1.0271, d_y = -0.467, d_z = 1.636;
              //mols.PD_x = wb_x * wb_x * d_x + wb_x * wb_y * d_y + wb_x * wb_z * d_z;
              //mols.PD_y = wb_y * wb_x * d_x + wb_y * wb_y * d_y + wb_y * wb_z * d_z ;
              //mols.PD_z = wb_z * wb_x * d_x + wb_z * wb_y * d_y + wb_z * wb_z * d_z;
              //Dipole moment [Debye]
              mols.PD_x = wb_x * Unitvectortodebye ;
              mols.PD_y = wb_y * Unitvectortodebye ;
              mols.PD_z = wb_z * Unitvectortodebye ;

              // Permanent polarisability [Angstrom]
              parameters(mols);
              float A_xx = mols.Pol_xx, A_yy = mols.Pol_yy, A_zz = mols.Pol_zz,
                    A_xy = mols.Pol_xy, A_xz = mols.Pol_xz, A_yz = mols.Pol_yz,
                    A_yx = mols.Pol_yx, A_zx = mols.Pol_zx, A_zy = mols.Pol_zy;     

              mols.Pol_xx =  HH_x * HH_x * A_xx + wb_x * wb_x * A_yy + Pv_x * Pv_x * A_zz
                           + HH_x * wb_x * A_xy + HH_x * Pv_x * A_xz + wb_x * Pv_x * A_yz 
                           + wb_x * HH_x * A_yx + Pv_x * HH_x * A_zx + Pv_x * wb_x * A_zy ;
              mols.Pol_yy =  HH_y * HH_y * A_xx + wb_y * wb_y * A_yy + Pv_y * Pv_y * A_zz
                           + HH_y * wb_y * A_xy + HH_y * Pv_y * A_xz + wb_y * Pv_y * A_yz 
                           + wb_y * HH_y * A_yx + Pv_y * HH_y * A_zx + Pv_y * wb_y * A_zy ;
              mols.Pol_zz =  HH_z * HH_z * A_xx + wb_z * wb_z * A_yy + Pv_z * Pv_z * A_zz
                           + HH_z * wb_z * A_xy + HH_z * Pv_z * A_xz + wb_z * Pv_z * A_yz 
                           + wb_z * HH_z * A_yx + Pv_z * HH_z * A_zx + Pv_z * wb_z * A_zy ;
              mols.Pol_xy =  HH_x * HH_y * A_xx + wb_x * wb_y * A_yy + Pv_x * Pv_y * A_zz
                           + HH_x * wb_y * A_xy + HH_x * Pv_y * A_xz + wb_x * Pv_y * A_yz 
                           + wb_x * HH_y * A_yx + Pv_x * HH_y * A_zx + Pv_x * wb_y * A_zy ;
              mols.Pol_xz =  HH_x * HH_z * A_xx + wb_x * wb_z * A_yy + Pv_x * Pv_z * A_zz
                           + HH_x * wb_z * A_xy + HH_x * Pv_z * A_xz + wb_x * Pv_z * A_yz 
                           + wb_x * HH_z * A_yx + Pv_x * HH_z * A_zx + Pv_x * wb_z * A_zy ;
              mols.Pol_yz =  HH_y * HH_z * A_xx + wb_y * wb_z * A_yy + Pv_y * Pv_z * A_zz
                           + HH_y * wb_z * A_xy + HH_y * Pv_z * A_xz + wb_y * Pv_z * A_yz 
                           + wb_y * HH_z * A_yx + Pv_y * HH_z * A_zx + Pv_y * wb_z * A_zy ;
              mols.Pol_yx =  HH_y * HH_x * A_xx + wb_y * wb_x * A_yy + Pv_y * Pv_x * A_zz
                           + HH_y * wb_x * A_xy + HH_y * Pv_x * A_xz + wb_y * Pv_x * A_yz 
                           + wb_y * HH_x * A_yx + Pv_y * HH_x * A_zx + Pv_y * wb_x * A_zy ;
              mols.Pol_zx =  HH_z * HH_x * A_xx + wb_z * wb_x * A_yy + Pv_z * Pv_x * A_zz
                           + HH_z * wb_x * A_xy + HH_z * Pv_x * A_xz + wb_z * Pv_x * A_yz 
                           + wb_z * HH_x * A_yx + Pv_z * HH_x * A_zx + Pv_z * wb_x * A_zy ;
              mols.Pol_zy =  HH_z * HH_y * A_xx + wb_z * wb_y * A_yy + Pv_z * Pv_y * A_zz
                           + HH_z * wb_y * A_xy + HH_z * Pv_y * A_xz + wb_z * Pv_y * A_yz 
                           + wb_z * HH_y * A_yx + Pv_z * HH_y * A_zx + Pv_z * wb_y * A_zy ;
              mol.push_back(mols);
            }
        }
    }
}




