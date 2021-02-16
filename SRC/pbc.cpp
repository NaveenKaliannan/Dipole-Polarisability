
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
#include "../include/pbc.h"
#include "../include/dipol.h"

using namespace std;

/*Calculates the distance between two atoms in different unit cell*/
void dist(vector<Molecular> &mol, uint idi, uint idj, const vector<float> & L,  vector<float> & PB_L, vector<Vector_int> & imageno, uint index,  double &x, double & y, double &z )
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



/*Replicates the center simulation box in all the directions*/
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



/* Broken bonds due to PBC are solved here, Minimum image convention was applied to bring them together */
void BringintoBox(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L)
{
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          if(r[id].x > L[0] || r[id].x < 0 ) { r[id].x  = min_distance(r[id].x, L[0]); }
          if(r[id].y > L[1] || r[id].y < 0 ) { r[id].y  = min_distance(r[id].y, L[1]); }
          if(r[id].z > L[2] || r[id].z < 0 ) { r[id].z  = min_distance(r[id].z, L[2]); }

          if(r[id].x > L[0] ) { r[id].x  -= L[0] ; }
          if(r[id].y > L[1] ) { r[id].y  -= L[1] ; }
          if(r[id].z > L[2] ) { r[id].z  -= L[2] ; }

          if(r[id].x < 0 ) { r[id].x += L[0] ; }
          if(r[id].y < 0 ) { r[id].y += L[1] ; }
          if(r[id].z < 0 ) { r[id].z += L[2] ; }
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
                                 if(r[id].x > r[id2].x ) { r[id2].x += L[0]; }    
                                 else if(r[id].x < r[id2].x ) { r[id2].x -= L[0]; }                                                                  
                               }
                             rij3 = norm(0, r[id].y - r[id2].y, 0);
                             if(abs(rij3) > 2.5)
                               {    
                                 if(r[id].y > r[id2].y ) { r[id2].y += L[1]; }    
                                 else if(r[id].y < r[id2].y ) { r[id2].y -= L[1]; }                                                                                                    
                               }  
                             rij3 = norm(0, 0, r[id].z - r[id2].z);
                             if(abs(rij3) > 2.5)
                               {                           
                                 if(r[id].z > r[id2].z ) { r[id2].z += L[2]; }    
                                 else if(r[id].z < r[id2].z ) { r[id2].z -= L[2]; }                                                                                                        
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
                                 if(r[id].x > r[id2].x ) { r[id2].x += L[0]; }    
                                 else if(r[id].x < r[id2].x ) { r[id2].x -= L[0]; }                                                                  
                               }
                             rij3 = norm(0, r[id].y - r[id2].y, 0);
                             if(abs(rij3) > 2.5)
                               {    
                                 if(r[id].y > r[id2].y ) { r[id2].y  =  r[id2].y += L[1]; }    
                                 else if(r[id].y < r[id2].y ) { r[id2].y -= L[1]; }                                                                                                    
                               }  
                             rij3 = norm(0, 0, r[id].z - r[id2].z);
                             if(abs(rij3) > 2.5)
                               {                           
                                 if(r[id].z > r[id2].z ) { r[id2].z += L[2]; }    
                                 else if(r[id].z < r[id2].z ) { r[id2].z -= L[2]; }                                                                                                        
                               }                                      
                            }
                        }
                    }
                } 
            }
        }
    }
}










