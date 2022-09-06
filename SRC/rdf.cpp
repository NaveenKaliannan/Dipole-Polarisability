
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

  for(uint t = 0; t < nsteps; t += deltat )  
    {
      cout << t << endl;
      count += 1;
      for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
        {
          count_atoms = 0;
          for(uint i = 0;i < natoms;++i)
            {
              idi = natoms*t+i;          
              for(uint j = 0;j < natoms;++j)
                {
                  idj = natoms*t+j;
                  if (r[idi].symbol[0] == 'M' && r[idj].symbol[0] == 'O' && i != j)
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



void population_hbonds(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<double> vec1(nsteps,0);
  vector<double> vec2(nsteps,0);
  vector<double> vec3(nsteps,0);
  vector<double> vec4(nsteps,0);
  vector<double> vec5(nsteps,0);

  for(uint t = 0; t < nsteps; ++t )
    {
      for(uint i = 0;i < natoms;++i)
        {
          uint idi = natoms*t+i;
          uint idi1 = natoms*t+i+1;
          uint idi2 = natoms*t+i+2;

          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H')
            {
              if( r[idi].totalhbonds  == 1) 
                {   
                  vec1[t] += 1;          
                }            
              if(r[idi].totalhbonds == 2) 
                {   
                  vec2[t] += 1;            
                }
              if(r[idi].totalhbonds == 3) 
                {    
                  vec3[t] += 1;                        
                }
              if(r[idi].totalhbonds == 4) 
                {     
                  vec4[t] += 1;                       
                }
              if(r[idi].totalhbonds == 5) 
                {     
                  vec5[t] += 1;                       
                }  
            }
        }
    }

  ofstream outfile_1(filename);
  for(uint t = 0;t < nsteps;++t)
    {
      outfile_1 << t*dt << "   " << vec1[t] * 100 / (vec1[t] + vec2[t] + vec3[t] + vec4[t] + vec5[t] )
                        << "   " << vec2[t] * 100 / (vec1[t] + vec2[t] + vec3[t] + vec4[t] + vec5[t] )
                        << "   " << vec3[t] * 100 / (vec1[t] + vec2[t] + vec3[t] + vec4[t] + vec5[t] ) 
                        << "   " << vec4[t] * 100 / (vec1[t] + vec2[t] + vec3[t] + vec4[t] + vec5[t] )
                        << "   " << vec5[t] * 100 / (vec1[t] + vec2[t] + vec3[t] + vec4[t] + vec5[t] ) << "\n";
    }
  outfile_1.close();
  outfile_1.clear();

}












