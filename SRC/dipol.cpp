
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
  uint idi = 0, idj = 0;
  float a = 0.0, sum = 0.0, b = 0.0, eps = 0.0;
  vector<Vector> Field (nmol);
  vector<Vector> rsd (nmol), zrsd (nmol), conj (nmol), vec(nmol);

  for(uint t = 0; t < nsteps;t += 1 )
    {
      FieldduetoPermanentMultipoles(mol, t, nmol, L, Field);
      //FieldduetoExternalField(mol, t, nmol, E,  Field);
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
      //  Modeling software package
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
            if (iter == niter - 1)  { cout << "Frame " << t <<  " is not converged " << endl;  }                                 
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
  uint idi = 0, idj = 0, ncell = 1;
  float x = 0, y = 0, z = 0, rij = 0, rcut = 40;
  float u = 0.0, ar = 0.0, st1 = 0.0, st2 = 0.0, r3 = 0.0, r5 = 0.0, pc = 0.0, pd = 0.0;
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
              if ((i != j && cell == 0 ) || (cell > 0 && rij < rcut))
              {
                mean_alpha1 = (mol[idi].PPol.xx + mol[idi].PPol.yy + mol[idi].PPol.zz ) / 3.0;
                mean_alpha2 = (mol[idj].PPol.xx + mol[idj].PPol.yy + mol[idj].PPol.zz ) / 3.0;
                u = rij / pow(mean_alpha1 * mean_alpha2, 0.16666);
                ar  = mol[idj].sl * pow(u,3);
                st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar); 
                r3  = pow(rij, -3.0)  * (st1 + 0.0); 
                r5  = 3.0 * pow(rij, -5.0)  * (st2 + 0.0);

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
  uint idi = 0, idj = 0, ncell = 1;
  float x = 0, y = 0, z = 0, rij = 0, rcut = 40;
  float u = 0.0, ar = 0.0, st1 = 0.0, st2 = 0.0, r3 = 0.0, r5 = 0.0, pc = 0.0, pd = 0.0;
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
              if ((i != j && cell == 0 ) || (cell > 0 && rij < rcut))
              {
                // S Grimme - Effect of the Damping Function in Dispersion Corrected Density Functional Theory
                // A Universal Damping Function for Empirical Dispersion Correctionon Density Functional Theory
                // vdw radii https://www.cgl.ucsf.edu/chimerax/docs/user/radii.html and Consistent van der Waals Radii for the Whole Main Group and Van der Waals Radii of Elements
                //Set of van der Waals and Coulombic Radii of Protein  Atoms for Molecular and Solvent-Accessible Surface Calculation, Packing Evaluation, and Docking
                // https://bip.weizmann.ac.il/course/structbioinfo/databases/CCDC_Mercury/appa_glossary.4.73.html Bondi, J.Phys.Chem., 68, 441, 1964
                // https://www.pamoc.it/kw_dfd.html
                // A consistent and accurate ab initioparametrization of density functionaldispersion correction (DFT-D) for the 94elements H-Pu
                // Density Functional Theory Augmented with an Empirical Dispersion Term. Interaction Energies and Geometries of 80 Noncovalent Complexes Compared with Ab Initio Quantum Mechanics Calculations
                // Damping functions in the effective fragmentpotential method
                // On the performance of molecular polarization methods. II. Water and carbon tetrachloride close to a cation
                // Gaussian and Exponeital damping function:  The polarizable point dipoles method withelectrostatic damping: Implementation on amodel system
                // Becke-Johnson  (BJ)  damping  function: A generally applicable atomic-chargedependent London dispersion correction
                // Thole damping:  Molecular and Atomic Polarizabilities: Thole’s Model Revisited
                // Thole damping: Molecular polarizabilities calculated with a modified dipole interaction
                // Empirical D3 Dispersion as a Replacement for ab Initio Dispersion Terms in Density Functional Theory-Based Symmetry-Adapted Perturbation Theory
                // atomic radii Semiempirical GGA-Type Density Functional Constructed with a Long-Range Dispersion Correction
                // Accurate Description of van der Waals Complexes by Density Functional Theory Including  Empirical Corrections
                // The polarizable point dipoles method withelectrostatic damping: Implementation on amodel system
                //Interatomic Methods for the Dispersion Energy Derived from the AdiabaticConnection Fluctuation-Dissipation Theorem
                //Note: Fermic  and Hesselmann daming function goes to zero quickly at smaller r, not suitable here. 
                // Wu and Yang damping function J. Chem. Phys. 2002, 116, 515-524. 

                //Thole damping function
                mean_alpha1 = (mol[idi].PPol.xx + mol[idi].PPol.yy + mol[idi].PPol.zz ) / 3.0;
                mean_alpha2 = (mol[idj].PPol.xx + mol[idj].PPol.yy + mol[idj].PPol.zz ) / 3.0;
                u = rij / pow(mean_alpha1 * mean_alpha2, 0.16666);
                ar  = mol[idj].sl * pow(u,3.0)  ;   
                st1 = 1.0 - (1.0 + ar + (ar*ar/2.0)) * exp(-ar);     
                st2 = 1.0 - (1.0 + ar + (ar*ar/2.0) + (pow(ar,3.0)/6.0)) * exp(-ar);                   
                //r3  = pow(rij, -3.0)  * (st1 + 0.0); 
                //r5  = 3.0 * pow(rij, -5.0)  * (st2 + 0.0);

                // enhanced thole screening length agrees with Wu and Yang damping function
                //r3  = pow(rij, -3.0)  * (st1 + 0.15); 
                //r5  = 3.0 * pow(rij, -5.0)  * (st2 + 0.15);

                // Wu and Yang damping function
                float cdamp = 3.54 ; 
                float r_vdwr = mol[idi].vdwr + mol[idj].vdwr ;
                float WY1 =  pow( 1.0  - exp( -1.0 * cdamp * pow( rij / r_vdwr ,3)) , 2);
                //r3  = pow(rij, -3.0)  * WY1 ; 
                //r5  = 3.0 * pow(rij, -5.0)  * WY1;

                r3  = pow(rij, -3.0) * WY1; 
                r5  = 3.0 * pow(rij, -5.0) * (st2 + 0.15) ;


                // Field due to point charge and dipole (quadrupole not considered)
                //point charge - http://www.physics.umd.edu/courses/Phys260/agashe/S10/notes/lecture18.pdf 
                //point dipole - http://www.physnet.org/modules/pdf_modules/m120.pdf 
                // Here we follow the standart physics convention to compare
                // with CP2K: dipole moment from - to +. hence i->j direction
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
   Yuki nagata papers
   Kenneth J Miller, Tholes, Siberstein and van Duijnen https://doi.org/10.1021/jp980221f research works
   Munir S. Skaf and Sérgio M. Vechi J. Chem. Phys 119, 2181 (2003); https://doi.org/10.1063/1.1583677  

   Simulation of the Far Infrared Spectrum of Liquid Water and Steam Along the Coexistence Curve
   from the Triple Point to the Critical Point https://doi.org/10.1007/978-94-011-0183-7_10

Pengyu Ren and Jay W. Ponder  Polarizable Atomic Multipole Water Model for Molecular Mechanics Simulation

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
      mol.q = 2;
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
      mol.q = 0.8;
      mol.PPol.xx = 0.1200;mol.PPol.yy = 0.1200;mol.PPol.zz = 0.1200;
      mol.PPol.xy = 0.0000;mol.PPol.xz = 0.0000;mol.PPol.yz = 0.0000;
      mol.PPol.yx = 0.0000;mol.PPol.zx = 0.0000;mol.PPol.zy = 0.0000;
      mol.PD.x    = 0.0000;mol.PD.y    = 0.0000;mol.PD.z    = 0.0000;
 
      mol.IPol.xx = 0.0000;mol.IPol.yy = 0.0000;mol.IPol.zz = 0.0000;
      mol.IPol.xy = 0.0000;mol.IPol.xz = 0.0000;mol.IPol.yz = 0.0000;
      mol.IPol.yx = 0.0000;mol.IPol.zx = 0.0000;mol.IPol.zy = 0.0000; 
      mol.ID.x    = 0.0000;mol.ID.y    = 0.0000;mol.ID.z    = 0.0000;
    }
  else if(mol.MOL[0] == 'S' && mol.MOL[1] == 'O' && mol.MOL[2] == '4')
    {
      mol.q = 0;
      mol.PPol.xx = 6.104;mol.PPol.yy = 6.0861;mol.PPol.zz = 6.1769;
//      mol.PPol.xy =-0.013;mol.PPol.xz = 0.0011;mol.PPol.yz = 0.0808;
//      mol.PPol.yx =-0.014;mol.PPol.zx = 0.0006;mol.PPol.zy = 0.0806; 
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
              mols.m = r[id].atomicmass;
              mols.sl = 0.3900 ; 
              mols.vdwr = 1.639 ;    
              parameters(mols);
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
              mols.sl = 0.3900 ;
              mols.vdwr = 1.75 ;   
              parameters(mols);
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

  PB_L[3] = PB_L[0] * 0.5;
  PB_L[4] = PB_L[1] * 0.5;
  PB_L[5] = PB_L[2] * 0.5;

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

  if(x > 0 && abs(x) > PB_L[3] ) { x  -= PB_L[0] ; }
  if(y > 0 && abs(y) > PB_L[4] ) { y  -= PB_L[1] ; }
  if(z > 0 && abs(z) > PB_L[5] ) { z  -= PB_L[2] ; } 

  if(x < 0 && abs(x) > PB_L[3] ) { x  += PB_L[0] ; }
  if(y < 0 && abs(y) > PB_L[4] ) { y  += PB_L[1] ; }
  if(z < 0 && abs(z) > PB_L[5] ) { z  += PB_L[2] ; }
}




