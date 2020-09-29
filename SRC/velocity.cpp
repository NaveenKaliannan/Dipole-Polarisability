
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


/* computes atomic velocity with central difference scheme as CP2K (multiply with vAperfmstoamu converts to atomic unit )*/
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
    /*if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition */
}


void PrintKEs(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_KE (nsteps, 0.0), trans_KE (nsteps, 0.0), rot_KE (nsteps, 0.0);
  vector<float> total_KE_cation (nsteps, 0.0), trans_KE_cation (nsteps, 0.0), rot_KE_cation (nsteps, 0.0);
  vector<float> total_KE_anion (nsteps, 0.0), trans_KE_anion (nsteps, 0.0), rot_KE_anion (nsteps, 0.0);
  vector<float> total_KE_both (nsteps, 0.0), trans_KE_both (nsteps, 0.0), rot_KE_both (nsteps, 0.0);
  vector<float> total_KE_rem (nsteps, 0.0), trans_KE_rem (nsteps, 0.0), rot_KE_rem (nsteps, 0.0);

  float rij1 = 0, rij2 = 0, rij = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_both = 0, count_rem = 0;
  uint hbond_cation = 0, hbond_anion = 0;
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  double com[3]; 


  computeatomicvelocity(r, nsteps, natoms, L, dt);
  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {

          hbond_cation = 0, hbond_anion = 0;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*3500+i;  
              uint idi1 = natoms*3500+i+1;          
              uint idi2 = natoms*3500+i+2;                  
              uint idj = natoms*3500+j;
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


          count += 1;
          if(hbond_cation == 0 && hbond_anion == 0)
            {
              count_rem   += 1; 
            }
          else if(hbond_cation > 0 && hbond_anion > 0)
            {
              count_both  += 1; 
            }
          else if(hbond_cation == 0 && hbond_anion > 0)
            {
              count_anion += 1; 
            }
          else if(hbond_cation > 0 && hbond_anion == 0)
            {
              count_cation += 1; 
            }

          for(uint t = 1; t < nsteps-1; ++t )
            {
              uint id = natoms*t+i;
              com[0] = r[id].vx * (am_O/am_H2O) +  r[id+1].vx * (am_H/am_H2O) +  r[id+2].vx * (am_H/am_H2O) ;
              com[1] = r[id].vy * (am_O/am_H2O) +  r[id+1].vy * (am_H/am_H2O) +  r[id+2].vy * (am_H/am_H2O) ;
              com[2] = r[id].vz * (am_O/am_H2O) +  r[id+1].vz * (am_H/am_H2O) +  r[id+2].vz * (am_H/am_H2O) ;

              temp1         = jtohartree * 0.5 * am_H2O *  norm2(com[0],com[1],com[2]) ;
              temp2         = jtohartree * 0.5 * am_O * norm2(  r[id].vx,  r[id].vy,  r[id].vz) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+1].vx,r[id+1].vy,r[id+1].vz) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+2].vx,r[id+2].vy,r[id+2].vz) ;

              temp3         = jtohartree * 0.5 * am_O * norm2(  r[id].vx - com[0],  r[id].vy - com[1],  r[id].vz - com[2]) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+1].vx - com[0],r[id+1].vy - com[1],r[id+1].vz - com[2]) 
                            + jtohartree * 0.5 * am_H * norm2(r[id+2].vx - com[0],r[id+2].vy - com[1],r[id+2].vz - com[2]) ;

              trans_KE[t]  += temp1 ;
              total_KE[t]  += temp2 ;
              rot_KE[t]    += temp3 ;

              if(hbond_cation == 0 && hbond_anion == 0)
                {
                  trans_KE_rem[t]  += temp1 ;
                  total_KE_rem[t]  += temp2 ;
                  rot_KE_rem[t]    += temp3 ;
                }
              else if(hbond_cation > 0 && hbond_anion > 0)
                {
                  trans_KE_both[t]  += temp1 ;
                  total_KE_both[t]  += temp2 ;
                  rot_KE_both[t]    += temp3 ;
                }
              else if(hbond_cation == 0 && hbond_anion > 0)
                {
                  trans_KE_anion[t]  += temp1 ;
                  total_KE_anion[t]  += temp2 ;
                  rot_KE_anion[t]    += temp3 ;
                }
              else if(hbond_cation > 0 && hbond_anion == 0)
                {
                  trans_KE_cation[t]  += temp1 ;
                  total_KE_cation[t]  += temp2 ;
                  rot_KE_cation[t]    += temp3 ;
                }


            }
        }
    }

  /* Total KE is in Atomic Unit in order to compare with CP2K data */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1; ++t )
    {
      outfile << t*dt << "  " << total_KE[t]       << "  " << trans_KE[t] / total_KE[t]               << "  " << rot_KE[t] / total_KE[t]               << "   " << 
                                total_KE_cation[t] << "  " << trans_KE_cation[t] / total_KE_cation[t] << "  " << rot_KE_cation[t] / total_KE_cation[t] << "   " <<
                                total_KE_anion[t]  << "  " << trans_KE_anion[t] / total_KE_anion[t]   << "  " << rot_KE_anion[t] / total_KE_anion[t]   << "   " <<
                                total_KE_both[t]   << "  " << trans_KE_both[t] / total_KE_both[t]     << "  " << rot_KE_both[t] / total_KE_both[t]     << "   " <<
                                total_KE_rem[t]    << "  " << trans_KE_rem[t] / total_KE_rem[t]       << "  " << rot_KE_rem[t] / total_KE_rem[t]       << endl;
    }
  outfile.close();
  outfile.clear();
}


void Printcosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> cosine (nsteps, 0.0), cosine_square (nsteps, 0.0);
  vector<float> cosine_cation (nsteps, 0.0), cosine_anion (nsteps, 0.0);
  vector<float> cosine_both (nsteps, 0.0), cosine_rem (nsteps, 0.0);

  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_both = 0, count_rem = 0;
  uint hbond_cation = 0, hbond_anion = 0;
  double vec[3]; 
  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {

          hbond_cation = 0, hbond_anion = 0;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*3500+i;  
              uint idi1 = natoms*3500+i+1;          
              uint idi2 = natoms*3500+i+2;                  
              uint idj = natoms*3500+j;
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

          count += 1;
          if(hbond_cation == 0 && hbond_anion == 0)
            {
              count_rem   += 1; 
            }
          else if(hbond_cation > 0 && hbond_anion > 0)
            {
              count_both  += 1; 
            }
          else if(hbond_cation == 0 && hbond_anion > 0)
            {
              count_anion += 1; 
            }
          else if(hbond_cation > 0 && hbond_anion == 0)
            {
              count_cation += 1; 
            }

          for(uint t = 1; t < nsteps-1; ++t )
            {
              uint id = natoms*t+i;
              vec[0] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].x - r[id].x)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].x - r[id].x) ;

              vec[1] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].y - r[id].y)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].y - r[id].y) ;

              vec[2] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].z - r[id].z)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].z - r[id].z) ;

              temp              =  vec[0] / pow(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2],0.5); 
              cosine[t]        += temp; 
              cosine_square[t] += pow(temp,2.0);
              if(hbond_cation == 0 && hbond_anion == 0)
                {
                  cosine_rem[t]   += temp; 
                }
              else if(hbond_cation > 0 && hbond_anion > 0)
                {
                  cosine_both[t]  += temp; 
                }
              else if(hbond_cation == 0 && hbond_anion > 0)
                {
                  cosine_anion[t] += temp; 
                }
              else if(hbond_cation > 0 && hbond_anion == 0)
                {
                  cosine_cation[t] += temp; 
                }
            }
        }
    }

  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1; ++t )
    {
      outfile << t*dt << "  " << cosine[t]        / count << "  " << 
                                 cosine_square[t] / count << "  " << 
                                 cosine_cation[t] / count_cation << "  " << 
                                 cosine_anion[t]  / count_anion << "  " << 
                                 cosine_both[t]   / count_both << "  " <<  
                                 cosine_rem[t]    / count_rem  << endl;
    }
  outfile.close();
  outfile.clear();
}















