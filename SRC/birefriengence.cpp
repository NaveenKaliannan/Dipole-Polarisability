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
#include "../include/birefriengence.h"



void Print_birefriengenceT_purewater(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0;

  for(uint i = 0;i < nmol;++i)
    {
      if(mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O') 
        {
          uint idi = nmol*3500+i;  
          /*counting*/
          count += 1;
          if(mol[idi].totalhbonds == 2)
            {
              count_hbond2   += 1; 
            }
          if(mol[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(mol[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.16 && mol[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            }


          for(uint t = 0; t < nsteps; t += deltat)
            {
              uint id = nmol*t+i;
              temp          =  (mol[id].PPol.xx + mol[id].IPol.xx) - 0.5 * (mol[id].PPol.yy + mol[id].IPol.yy + mol[id].PPol.zz + mol[id].IPol.zz )  ; 
              total_OBF[t] += temp ;

              if(mol[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02)
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.16 && mol[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }
            }
        } 
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 0; t < nsteps; t += deltat)
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}


void Print_birefriengenceP_purewater(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0;

  for(uint i = 0;i < nmol;++i)
    {
      if(mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O') 
        {
          uint idi = nmol*3500+i;  
          /*counting*/
          count += 1;
          if(mol[idi].totalhbonds == 2)
            {
              count_hbond2   += 1; 
            }
          if(mol[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(mol[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.16 && mol[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            }


          for(uint t = 0; t < nsteps; t += deltat)
            {
              uint id = nmol*t+i;
              temp          =  (mol[id].PPol.xx) - 0.5 * (mol[id].PPol.yy + mol[id].PPol.zz )  ; 
              total_OBF[t] += temp ;

              if(mol[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02)
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.16 && mol[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }
            }
        } 
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 0; t < nsteps; t += deltat)
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}


void Print_birefriengenceI_purewater(vector<Molecular> &mol, uint nsteps, uint nmol, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0;

  for(uint i = 0;i < nmol;++i)
    {
      if(mol[i].MOL[0] == 'H' && mol[i].MOL[1] == '2' && mol[i].MOL[2] == 'O') 
        {
          uint idi = nmol*3500+i;  
          /*counting*/
          count += 1;
          if(mol[idi].totalhbonds == 2)
            {
              count_hbond2   += 1; 
            }
          if(mol[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(mol[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.16 && mol[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            }


          for(uint t = 0; t < nsteps; t += deltat)
            {
              uint id = nmol*t+i;
              temp          =  (mol[id].IPol.xx) - 0.5 * (mol[id].IPol.yy  + mol[id].IPol.zz )  ; 
              total_OBF[t] += temp ;

              if(mol[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02)
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.16 && mol[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d < 1.02 && mol[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(mol[idi].totalhbonds == 4 && mol[idi].totaldonorhbonds == 2 && mol[idi].gamma_d > 1.16 && mol[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }
            }
        } 
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 0; t < nsteps; t += deltat)
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}



void Print_Cosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0;
  double vec[3]; 

  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          uint idi = natoms*3500+i;  
          /*counting*/
          count += 1;
          if(r[idi].totalhbonds == 2)
            {
              count_hbond2   += 1;
            }
          if(r[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(r[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            } 

          for(uint t = 1; t < nsteps-1;  t += deltat )
            {
              uint id = natoms*t+i; 

              /* cosine and cosine square theta */
              vec[0] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].x - r[id].x)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].x - r[id].x) ;

              vec[1] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].y - r[id].y)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].y - r[id].y) ;

              vec[2] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].z - r[id].z)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].z - r[id].z) ;

              temp              =  vec[0] / pow(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2],0.5); 
              total_OBF[t] += temp ; 
              if(r[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(r[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)    
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }    
            }
        }
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}

void Print_Cosine2(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0;
  double vec[3]; 

  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          uint idi = natoms*3500+i;  
          /*counting*/
          count += 1;
          if(r[idi].totalhbonds == 2)
            {
              count_hbond2   += 1; 
            }
          if(r[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(r[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            }

          for(uint t = 1; t < nsteps-1;  t += deltat )
            {
              uint id = natoms*t+i;

              /* cosine and cosine square theta */
              vec[0] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].x - r[id].x)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].x - r[id].x) ;

              vec[1] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].y - r[id].y)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].y - r[id].y) ;

              vec[2] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].z - r[id].z)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].z - r[id].z) ;

              temp              =  vec[0] / pow(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2],0.5); 
              total_OBF[t]      += pow(temp,2.0) ;

              if(r[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(r[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)    
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }    
            }
        }
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}


void Print_KEtrans(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0, temp1 = 0, temp2 = 0, temp3 = 0;
  double com[3]; 
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;

  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          uint idi = natoms*3500+i;  
          /*counting*/
          count += 1;
          if(r[idi].totalhbonds == 2)
            {
              count_hbond2   += 1; 
            }
          if(r[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(r[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            }

          for(uint t = 1; t < nsteps-1;  t += deltat )
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

              temp         = temp1/temp2 ;
              total_OBF[t]      += temp ;

              if(r[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(r[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)    
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }    
            }
        }
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}


void Print_KErot(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_OBF        (nsteps, 0.0),
                total_OBF_hbond2 (nsteps, 0.0),
                total_OBF_hbond3  (nsteps, 0.0),
                total_OBF_hbond4   (nsteps, 0.0),
                total_OBF_hbond4_lgd   (nsteps, 0.0), 
                total_OBF_hbond4_mgd   (nsteps, 0.0), 
                total_OBF_hbond4_hgd   (nsteps, 0.0), 
                total_OBF_hbond4_lgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdlga   (nsteps, 0.0), 
                total_OBF_hbond4_lgdhga   (nsteps, 0.0), 
                total_OBF_hbond4_hgdhga    (nsteps, 0.0);

  float count = 0, count_hbond2 = 0, count_hbond3 = 0, count_hbond4 = 0, 
        count_hbond4_lgd = 0, count_hbond4_hgd = 0, count_hbond4_mgd = 0, 
        count_hbond4_lgdlga = 0, count_hbond4_hgdlga = 0, count_hbond4_lgdhga = 0, count_hbond4_hgdhga = 0;

  float temp = 0, temp1 = 0, temp2 = 0, temp3 = 0;
  double com[3]; 
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;

  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          uint idi = natoms*3500+i;  
          /*counting*/
          count += 1;
          if(r[idi].totalhbonds == 2)
            {
              count_hbond2   += 1; 
            }
          if(r[idi].totalhbonds == 3)
            {
              count_hbond3   += 1; 
            }
          if(r[idi].totalhbonds == 4)
            {
              count_hbond4   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)
            {
              count_hbond4_lgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
            {
              count_hbond4_hgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
            {
              count_hbond4_mgd   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_lgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
            {
              count_hbond4_hgdlga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_lgdhga   += 1; 
            }
          if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
            {
              count_hbond4_hgdhga   += 1; 
            }

          for(uint t = 1; t < nsteps-1;  t += deltat )
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

              temp         = temp3/temp2 ;
              total_OBF[t]      += temp ;

              if(r[idi].totalhbonds == 2)
                {
                  total_OBF_hbond2[t]   += temp; 
                }
              if(r[idi].totalhbonds == 3)
                {
                  total_OBF_hbond3[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4)
                {
                  total_OBF_hbond4[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02)    
                {
                  total_OBF_hbond4_lgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16)
                {
                  total_OBF_hbond4_hgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.16 && r[idi].gamma_d > 1.02)
                {
                  total_OBF_hbond4_mgd[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_lgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a < 1.02)
                {
                  total_OBF_hbond4_hgdlga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d < 1.02 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_lgdhga[t]   += temp; 
                }
              if(r[idi].totalhbonds == 4 && r[idi].totaldonorhbonds == 2 && r[idi].gamma_d > 1.16 && r[idi].gamma_a > 1.16)
                {
                  total_OBF_hbond4_hgdhga[t]   += temp; 
                }    
            }
        }
    }

  if(count == 0)               {count = 1;}
  if(count_hbond2 == 0)        {count_hbond2 = 1;}
  if(count_hbond3 == 0)        {count_hbond3 = 1;}
  if(count_hbond4 == 0)        {count_hbond4 = 1;}
  if(count_hbond4_lgd == 0)    {count_hbond4_lgd = 1;}
  if(count_hbond4_hgd == 0)    {count_hbond4_hgd = 1;}
  if(count_hbond4_mgd == 0)    {count_hbond4_mgd = 1;}
  if(count_hbond4_lgdlga == 0) {count_hbond4_lgdlga = 1;}
  if(count_hbond4_hgdlga == 0) {count_hbond4_hgdlga = 1;}
  if(count_hbond4_lgdhga == 0) {count_hbond4_lgdhga = 1;}
  if(count_hbond4_hgdhga == 0) {count_hbond4_hgdhga = 1;}


  /* Optical birefringence is in angstrom3 and normalized to per molecule */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      outfile << t*dt << "  " << total_OBF[t]               / count        << "  " << 
                                 total_OBF_hbond2[t]        / count_hbond2 << "  " << 
                                 total_OBF_hbond3[t]        / count_hbond3  << "  " << 
                                 total_OBF_hbond4[t]        / count_hbond4   << "  " << 
                                 total_OBF_hbond4_lgd[t]    / count_hbond4_lgd   << "  " << 
                                 total_OBF_hbond4_mgd[t]    / count_hbond4_mgd  << "  " << 
                                 total_OBF_hbond4_hgd[t]    / count_hbond4_hgd   << "  " << 
                                 total_OBF_hbond4_lgdlga[t] / count_hbond4_lgdlga   << "  " << 
                                 total_OBF_hbond4_hgdlga[t] / count_hbond4_hgdlga  << "  " << 
                                 total_OBF_hbond4_lgdhga[t] / count_hbond4_lgdhga   << "  " << 
                                 total_OBF_hbond4_hgdhga[t] / count_hbond4_hgdhga    << "  " <<  endl;
    }
  outfile.close();
  outfile.clear();
}


