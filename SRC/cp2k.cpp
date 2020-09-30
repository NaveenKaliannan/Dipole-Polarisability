
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
#include "../include/cp2k.h"
#include "../include/io.h"

using namespace std;




/*Reads mulliken charges from cp2k output data*/
void readmullikencharges(vector<Atom> &r, uint nsteps, uint natoms, string filename)
{
  string temp;
  for(uint t = 0; t < nsteps; ++t )
    { 
      string filename = "/home/naveenk/temp/1-mgcl2/Frame" + to_string(t+1) + "/charges.out.mulliken" ;
      ifstream charges_filename(filename);
      for(uint i = 0; i <  5; ++i)
        {
          getline(charges_filename, temp);
        }
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          charges_filename >> temp >> temp >> temp >> temp >> r[id].charge; 
          //cout << t <<  "  " << i << "  " << r[id].symbol << "  " <<  r[id].charge << endl;
        }
      charges_filename.close();
      charges_filename.clear();
    }
}








