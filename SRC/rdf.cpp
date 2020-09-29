
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
#include "../include/rdf.h"
#include "../include/dipol.h"

using namespace std;


void Printrdf(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  float x = 0, y = 0, z = 0, rij = 0, V = L[0] * L[1] * L[2], len = max(L[0], L[1]), RDF_h = 0.1, count = 0 ;
  uint idi = 0, idj = 0, ncell = 0, count_atoms = 0;
  len = max(len, L[2]);
  len = len/2.0;
  uint RDF_size = len*100;
  vector<float> RDF(RDF_size,0.0), rad(RDF_size,0.0);
  for(uint i = 1;i < RDF_size ;i++)
    {
      rad[i]          = i * 0.01;
    }  
  float no_of_residue_first  = 0;
  float no_of_residue_second = 0; 
  cout << "g \u03B1 - \u03B2 (rij)" << endl; ;
  cout << "Enter number of \u03B1 species : " ;
  cin >>  no_of_residue_first ;
  cout << "Enter number of \u03B2 species : " ;
  cin >>  no_of_residue_second ;

  for(uint t = 0; t < nsteps; t += 5 )  
    {
      cout << t << endl;
      count += 1;
      for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
        {
          count_atoms = 0;
          for(uint i = 0;i < natoms-1;++i)
            {
              idi = natoms*t+i;          
              for(uint j = i+1;j < natoms;++j)
                {
                  idj = natoms*t+j;
                  if (r[idi].symbol[0] == 'F' && r[idj].symbol[0] == 'H' && i != j)
                    {                                       
                      x = min_distance(r[idj].x - r[idi].x, L[0]);
                      y = min_distance(r[idj].y - r[idi].y, L[1]);
                      z = min_distance(r[idj].z - r[idi].z, L[2]); 
                      rij = mindis(x,y,z,L);       
                      if (rij <= rad[RDF_i] + (RDF_h/2.0) && rij > rad[RDF_i] - (RDF_h/2.0))
                        {
                          count_atoms += 1;  
                        }
                    }
                }

            }
          RDF[RDF_i]  += count_atoms  / ( ( no_of_residue_second / V) * 4.0 * PI * rad[RDF_i] * rad[RDF_i] *  RDF_h);
        }
    }

  // consine thetha per water molecule
  ofstream outfile(filename);
  for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
    {
      outfile << rad[RDF_i]  << "  " << RDF[RDF_i] / (count * no_of_residue_first )  << endl;
    }
  outfile.close();
  outfile.clear();
}















