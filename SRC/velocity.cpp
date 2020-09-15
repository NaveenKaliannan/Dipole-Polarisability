
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
#include "../include/velocity.h"


using namespace std;




// computes atomic velocity with central difference scheme as CP2K (multiply with vAperfmstoamu converts to atomic unit )
void computeatomicvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  for(uint t = 1; t < nsteps-1; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*(t)+i; 
          uint id1 = natoms*(t-1)+i; 
          uint id2 = natoms*(t+1)+i; 
          r[id].vx = min_distance((r[id2].x - r[id1].x), L[0])/(2*dt) ;
          r[id].vy = min_distance((r[id2].y - r[id1].y), L[1])/(2*dt) ;
          r[id].vz = min_distance((r[id2].z - r[id1].z), L[2])/(2*dt) ;
        }
    }
    //if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition 
}


void PrintKEs(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_KE (nsteps, 0.0), trans_KE (nsteps, 0.0), rot_KE (nsteps, 0.0);
  float count = 0;
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  double com[3]; // com-center of mass velocity
  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          count += 1;
          for(uint t = 1; t < nsteps-1; ++t )
            {
              uint id = natoms*t+i;
              com[0] = r[id].vx * (am_O/am_H2O) +  r[id+1].vx * (am_H/am_H2O) +  r[id+2].vx * (am_H/am_H2O) ;
              com[1] = r[id].vy * (am_O/am_H2O) +  r[id+1].vy * (am_H/am_H2O) +  r[id+2].vy * (am_H/am_H2O) ;
              com[2] = r[id].vz * (am_O/am_H2O) +  r[id+1].vz * (am_H/am_H2O) +  r[id+2].vz * (am_H/am_H2O) ;
              trans_KE[t] +=  jtohartree * 0.5 * am_H2O *  norm2(com[0],com[1],com[2]) ;
              total_KE[t]  += jtohartree * 0.5 * am_O * norm2(  r[id].vx,  r[id].vy,  r[id].vz) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+1].vx,r[id+1].vy,r[id+1].vz) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+2].vx,r[id+2].vy,r[id+2].vz) ;

              rot_KE[t]    += jtohartree * 0.5 * am_O * norm2(  r[id].vx - com[0],  r[id].vy - com[1],  r[id].vz - com[2]) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+1].vx - com[0],r[id+1].vy - com[1],r[id+1].vz - com[2]) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+2].vx - com[0],r[id+2].vy - com[1],r[id+2].vz - com[2]) ;
            }
        }
    }

  // Total KE is in Atomic Unit in order to compare with CP2K data
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1; ++t )
    {
      outfile << t*dt << "  " << total_KE[t]  << "  " << trans_KE[t] / total_KE[t] << "  " << rot_KE[t] / total_KE[t]  << endl;
    }
  outfile.close();
  outfile.clear();
}


void Printcosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> cosine (nsteps, 0.0);
  float count = 0;
  double vec[3]; 
  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          count += 1;
          for(uint t = 1; t < nsteps-1; ++t )
            {
              uint id = natoms*t+i;
              vec[0] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].x - r[id].x)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].x - r[id].x) ;

              vec[1] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].y - r[id].y)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].y - r[id].y) ;

              vec[2] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].z - r[id].z)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].z - r[id].z) ;

              cosine[t] += vec[0] / pow(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2],0.5); 
            }
        }
    }

  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1; ++t )
    {
      outfile << t*dt << "  " << cosine[t] / count  << endl;
    }
  outfile.close();
  outfile.clear();
}















