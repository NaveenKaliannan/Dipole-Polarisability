
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
#include "../include/tinker.h"


using namespace std;



//Reading the xyz trajectory
void readtrajectory_tinker(vector<Atom> &r, uint nsteps, uint natoms, string xyzfilename, const vector<float> & L)
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
          xyzfile >> temp >> r[id].symbol ;
          if((r[id].symbol[0] == 'N' || r[id].symbol[0] == 'n' ) && ( r[id].symbol[1] == 'A' || r[id].symbol[1] == 'a'))
            {
              xyzfile >> r[id].x  >> r[id].y  >> r[id].z >> temp;          
            }
          if((r[id].symbol[0] == 'M' || r[id].symbol[0] == 'm' ) && ( r[id].symbol[1] == 'G' || r[id].symbol[1] == 'g'))
            {
              xyzfile >> r[id].x  >> r[id].y  >> r[id].z >> temp;          
            }
          if((r[id].symbol[0] == 'C' || r[id].symbol[0] == 'c' ) && ( r[id].symbol[1] == 'L' || r[id].symbol[1] == 'l'))
            {
              xyzfile >> r[id].x  >> r[id].y  >> r[id].z >> temp;          
            }
          if((r[id].symbol[0] == 'F' || r[id].symbol[0] == 'f' ))
            {
              xyzfile >> r[id].x  >> r[id].y  >> r[id].z >> temp;          
            }
          if((r[id].symbol[0] == 'O' || r[id].symbol[0] == 'o' ))
            {
              xyzfile >> r[id].x  >> r[id].y  >> r[id].z >> temp >> temp >> temp;          
            }
          if((r[id].symbol[0] == 'H' || r[id].symbol[0] == 'h' ))
            {
              xyzfile >> r[id].x  >> r[id].y  >> r[id].z >> temp >> temp ;          
            }
          //cout << r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << endl;
        }
    }
  xyzfile.close();
  xyzfile.clear();
}


void Print_tinker(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, vector<Molecular> &mol, uint nmol, float dt, string filename, string TYPE)
{
  ofstream outfile(filename);
  for(uint t = 0; t < nsteps; ++t )
    { 
      if(TYPE[0] == 'T' && TYPE[1] == 'I' && TYPE[2] == 'N' && TYPE[3] == 'K' && TYPE[4] == 'E' && TYPE[5] == 'R')
        {
          outfile << "   " << natoms  <<  "  TINKER Format" << endl ;
          for(uint i = 0;i < natoms;++i)
            {
              uint id = natoms*t+i; 
              if((r[id].symbol[0] == 'N' || r[id].symbol[0] == 'n' ) && ( r[id].symbol[1] == 'A' || r[id].symbol[1] == 'a'))
                {
                  outfile << "   " << i+1 << "  " <<  r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << "    106"  << endl;
                }
              if((r[id].symbol[0] == 'M' || r[id].symbol[0] == 'm' ) && ( r[id].symbol[1] == 'G' || r[id].symbol[1] == 'g'))
                {
                  outfile << "   " << i+1 << "  " <<  r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << "    111"  << endl;
                }
              if((r[id].symbol[0] == 'C' || r[id].symbol[0] == 'c' ) && ( r[id].symbol[1] == 'L' || r[id].symbol[1] == 'l'))
                {
                  outfile << "   " << i+1 << "  " <<  r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << "    115"  << endl;
                }
              if((r[id].symbol[0] == 'F' || r[id].symbol[0] == 'f' ))
                {
                  outfile << "   " << i+1 << "  " <<  r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << "    114"  << endl;
                }
              if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
                {
                  outfile << "   " << i+1 << "  " <<  r[id].symbol <<  "  " << r[id].x << "  " << r[id].y << "  " << r[id].z << "    103 " << i+2 << "  "<< i+3 << endl;
                  outfile << "   " << i+2 << "  " <<  r[id+1].symbol <<  "  " << r[id+1].x << "  " << r[id+1].y << "  " << r[id+1].z << "    104 " << i+1 << endl;
                  outfile << "   " << i+3 << "  " <<  r[id+2].symbol <<  "  " << r[id+2].x << "  " << r[id+2].y << "  " << r[id+2].z << "    104 " << i+1 << endl;
                }
            }
        }
    } 
  outfile.close();
  outfile.clear();  
}
