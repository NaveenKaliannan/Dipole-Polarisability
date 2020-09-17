
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
#include "../include/molar.h"

// https://chemistry.mdma.ch/hiveboard/rhodium/pdf/chemical-data/prop_aq.pdf
// Increased interfacial thickness of the NaF, NaCl and NaBr salt aqueous
// solutions probed with non-resonant surface second harmonic generation
// (SHG)
// https://aip.scitation.org/doi/am-pdf/10.1063/1.4926840
// NaF 86 page https://refubium.fu-berlin.de/bitstream/handle/fub188/4582/Rinne_Klaus_online.pdf?sequence=1&isAllowed=y
// https://refubium.fu-berlin.de/bitstream/handle/fub188/4582/Rinne_Klaus_online.pdf?sequence=1&isAllowed=y
// https://www.snc.edu/chemicalhygiene/docs/labsafety/SolutionPrep.pdf
// http://rosdok.uni-rostock.de/file/rosdok_disshab_0000000812/rosdok_derivate_0000004838/Dissertation_Riemenschneider_2012.pdf
// https://pubs.acs.org/doi/pdf/10.1021/je060220d
// A Review of Sodium Fluoride Solubility in Water
// https://www.nist.gov/system/files/documents/srd/jpcrd506.pdf
// https://link.springer.com/content/pdf/10.1007/BF02907785.pdf


using namespace std;


void EstimateCellSizenNumberofResidue()
{
  double L = 0, mol_L = 0, mol_Kg = 0, density = 0;
  double fw_water = 18.01528;
  double amu_to_kg = 1.660538921E-27;
  double N_A       = 6.023E23;
  double fw_sol = 0;
  cout << "atomic mass of sol (Example NaF 41.98817 g  per mol, mgcl2 95.211 g per mol) : " << endl;
  cin >>   fw_sol ;

  cout << "Cube Length (in angstrom)" << endl;
  cin  >> L ; 
  double V_L  = pow(L*1E-10,3) * 1000;
  double V_cm3 = V_L * 1000;

  cout << "Density" << endl;
  cin >> density ;
  cout << "Mol per liter" << endl;
  cin  >> mol_L ; 
  cout << "Mol per Kg " << endl;
  cin >> mol_Kg;
  double box_mass_gm_sol = density * V_cm3;
  double n_sol = mol_L * V_L;
  double wt_water_kg = n_sol/mol_Kg;
  double n_wat = wt_water_kg / (amu_to_kg * fw_water) ; 
  cout << n_wat << endl;
  cout << n_sol * N_A << endl; 

  double N_sol = round(n_sol * N_A) ;
  double N_wat = round(n_wat) ; cout <<" Rounded values number of sols = " << N_sol << "  number of H2O = " << N_wat << endl;

  cout << "Density  =  " << (N_sol * fw_sol + N_wat * fw_water) * amu_to_kg / V_L << endl;
  cout << "Accurate box size  =  " << pow((N_sol * fw_sol + N_wat * fw_water) * amu_to_kg * 1000 / (density * 1E-21) , 1.0/3.0) << endl;
}














