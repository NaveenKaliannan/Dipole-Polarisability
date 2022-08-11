#include <stdlib.h>
#include<cstdlib>
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

bool isNaN(double x) { 
  return x != x;
}


/* computes atomic velocity via central difference scheme as CP2K (multiply with vAperfmstoamu converts to atomic unit )*/
void computeatomicvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  for(uint t = 1; t < nsteps-1;  t += deltat )
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


void computecomposition(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      for(uint i = 0;i < natoms;++i)
        {      
          uint idi = natoms*t+i;
          uint idi1 = natoms*t+i+1;
          uint idi2 = natoms*t+i+2;

          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H')
            {
             r[idi].x -=  L[0] * round(r[idi].x/L[0]);
             r[idi].y -=  L[1] * round(r[idi].y/L[1]);
             r[idi].z -=  L[2] * round(r[idi].z/L[2]);

             r[idi+1].x -=  L[0] * round(r[idi+1].x/L[0]);
             r[idi+1].y -=  L[1] * round(r[idi+1].y/L[1]);
             r[idi+1].z -=  L[2] * round(r[idi+1].z/L[2]);

             r[idi+2].x -=  L[0] * round(r[idi+2].x/L[0]);
             r[idi+2].y -=  L[1] * round(r[idi+2].y/L[1]);
             r[idi+2].z -=  L[2] * round(r[idi+2].z/L[2]);

             r[idi].comx = r[idi].x * (am_O/am_H2O) +  r[idi+1].x * (am_H/am_H2O) +  r[idi+2].x * (am_H/am_H2O) ;
             r[idi].comy = r[idi].y * (am_O/am_H2O) +  r[idi+1].y * (am_H/am_H2O) +  r[idi+2].y * (am_H/am_H2O) ;
             r[idi].comz = r[idi].z * (am_O/am_H2O) +  r[idi+1].z * (am_H/am_H2O) +  r[idi+2].z * (am_H/am_H2O) ;
            }
        }
    }
    /*if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition */
}

void computepositionmolframe(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      for(uint i = 0;i < natoms;i += 3)
        {      
          uint idi = natoms*t+i;
          uint idi1 = natoms*t+i+1;
          uint idi2 = natoms*t+i+2;

          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H')
            {
             r[idi+1].x -= r[idi].comx ;
             r[idi+1].y -= r[idi].comy ;
             r[idi+1].z -= r[idi].comz ;

             r[idi+2].x -= r[idi].comx ;
             r[idi+2].y -= r[idi].comy ;
             r[idi+2].z -= r[idi].comz ;

             r[idi].x -= r[idi].comx ;
             r[idi].y -= r[idi].comy ;
             r[idi].z -= r[idi].comz ;
            }
        }
    }
    /*if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition */
}

void computecomvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      for(uint i = 0;i < natoms;++i)
        {     
          uint idi = natoms*t+i;
          uint idi1 = natoms*t+i+1;
          uint idi2 = natoms*t+i+2;

          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H') 
            {
              r[idi].comvx = r[idi].vx * (am_O/am_H2O) +  r[idi+1].vx * (am_H/am_H2O) +  r[idi+2].vx * (am_H/am_H2O) ;
              r[idi].comvy = r[idi].vy * (am_O/am_H2O) +  r[idi+1].vy * (am_H/am_H2O) +  r[idi+2].vy * (am_H/am_H2O) ;
              r[idi].comvz = r[idi].vz * (am_O/am_H2O) +  r[idi+1].vz * (am_H/am_H2O) +  r[idi+2].vz * (am_H/am_H2O) ;
            }
        }
    }
    /*if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition */
}


void computeangularvelocity(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{
  computecomposition(r, nsteps, natoms, L, dt);
  computepositionmolframe(r, nsteps, natoms, L, dt);
  computecomposition(r, nsteps, natoms, L, dt);
  computecomvelocity(r, nsteps, natoms, L, dt);
  double vect_A[3], vect_B[3], normr ;
  double wb_x = 0,wb_y = 0,wb_z = 0;           // bisector vector H2O
  double wb1_x = 0,wb1_y = 0,wb1_z = 0;        // bisector vector OH1
  double wb2_x = 0,wb2_y = 0,wb2_z = 0;        // bisector vector OH2
  double HH_x = 0,HH_y = 0,HH_z = 0 ;          // H-H vector
  double Pv_x = 0,Pv_y = 0,Pv_z = 0;           // vector Perpendicular to bisector and H-H vector
  double v1,v2,v3;

  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      for(uint i = 0;i < natoms;++i)
        {
          uint idi = natoms*t+i;
          uint idi1 = natoms*t+i+1;
          uint idi2 = natoms*t+i+2;

          if(r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' && r[idi+2].symbol[0] == 'H')
            {
              vect_A[0] = r[idi].x - r[idi].comx ;
              vect_A[1] = r[idi].y - r[idi].comy ;
              vect_A[2] = r[idi].z - r[idi].comz ;
              vect_B[0] = r[idi].vx - r[idi].comvx ;
              vect_B[1] = r[idi].vy - r[idi].comvy ;
              vect_B[2] = r[idi].vz - r[idi].comvz ;
              normr =  pow(vect_A[0],2) + pow(vect_A[1],2) + pow(vect_A[2],2) ;
              r[idi].angvx = (vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]) / normr  ;
              r[idi].angvy = (vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]) / normr  ;
              r[idi].angvz = (vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]) / normr  ;

              vect_A[0] = r[idi+1].x - r[idi].comx  ;
              vect_A[1] = r[idi+1].y - r[idi].comy ;
              vect_A[2] = r[idi+1].z - r[idi].comz ;
              vect_B[0] = r[idi+1].vx - r[idi].comvx ;
              vect_B[1] = r[idi+1].vy - r[idi].comvy ;
              vect_B[2] = r[idi+1].vz - r[idi].comvz;
              normr =  pow(vect_A[0],2) + pow(vect_A[1],2) + pow(vect_A[2],2) ;
              r[idi].angvx += (vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]) / normr  ;
              r[idi].angvy += (vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]) / normr  ;
              r[idi].angvz += (vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]) / normr  ;

              vect_A[0] = r[idi+2].x - r[idi].comx ;
              vect_A[1] = r[idi+2].y - r[idi].comy ;
              vect_A[2] = r[idi+2].z - r[idi].comz ;
              vect_B[0] = r[idi+2].vx - r[idi].comvx ;
              vect_B[1] = r[idi+2].vy - r[idi].comvy ;
              vect_B[2] = r[idi+2].vz - r[idi].comvz ;
              normr =  pow(vect_A[0],2) + pow(vect_A[1],2) + pow(vect_A[2],2) ;
              r[idi].angvx += (vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]) / normr  ;
              r[idi].angvy += (vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]) / normr  ;
              r[idi].angvz += (vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]) / normr  ;

              wb1_x = min_distance(r[idi+1].x - r[idi].x, L[0]) ;
              wb1_y = min_distance(r[idi+1].y - r[idi].y, L[1]) ;
              wb1_z = min_distance(r[idi+1].z - r[idi].z, L[2]) ;
              convertounivector(wb1_x, wb1_y, wb1_z);
              wb2_x = min_distance(r[idi+2].x - r[idi].x, L[0]) ;
              wb2_y = min_distance(r[idi+2].y - r[idi].y, L[1]) ;
              wb2_z = min_distance(r[idi+2].z - r[idi].z, L[2]) ;
              convertounivector(wb2_x, wb2_y, wb2_z);

              wb_x = wb1_x + wb2_x ;
              wb_y = wb1_y + wb2_y ;
              wb_z = wb1_z + wb2_z ;
              convertounivector(wb_x, wb_y, wb_z);

              HH_x = (r[idi+1].x - r[idi+2].x) ;
              HH_y = (r[idi+1].y - r[idi+2].y) ;
              HH_z = (r[idi+1].z - r[idi+2].z) ;
              convertounivector(HH_x, HH_y, HH_z);

              Pv_x = wb_y * HH_z - wb_z * HH_y;
              Pv_y = wb_z * HH_x - wb_x * HH_z;
              Pv_z = wb_x * HH_y - wb_y * HH_x;
              convertounivector(Pv_x, Pv_y, Pv_z);

              HH_x = wb_y * Pv_z - wb_z * Pv_y;
              HH_y = wb_z * Pv_x - wb_x * Pv_z;
              HH_z = wb_x * Pv_y - wb_y * Pv_x;
              convertounivector(HH_x, HH_y, HH_z);

//              for(unsigned int i_ = 0; i_ < 3 ; i_ = i_ + 1)
//              {
//                uint idi_ = natoms*t+(i+i_);
//                v1 = r[idi_].x ;
//                v2 = r[idi_].y ;
//                v3 = r[idi_].z ;
//                r[idi_].x = HH_x * v1 + HH_y * v2  + HH_z * v3 ;
//                r[idi_].y = wb_x * v1 + wb_y * v2  + wb_z * v3 ;
//                r[idi_].z = Pv_x * v1 + Pv_y * v2  + Pv_z * v3 ;
//              
//                v1 = r[idi_].vx ;
//                v2 = r[idi_].vy ;
//                v3 = r[idi_].vz ;
//                r[idi_].vx = HH_x * v1 + HH_y * v2  + HH_z * v3 ;
//                r[idi_].vy = wb_x * v1 + wb_y * v2  + wb_z * v3 ;
//                r[idi_].vz = Pv_x * v1 + Pv_y * v2  + Pv_z * v3 ;
//              }
//

              v1 = r[idi].angvx ;
              v2 = r[idi].angvy ;
              v3 = r[idi].angvz ;
              r[idi].angvx = HH_x * v1 + HH_y * v2  + HH_z * v3 ;
              r[idi].angvy = wb_x * v1 + wb_y * v2  + wb_z * v3 ;
              r[idi].angvz = Pv_x * v1 + Pv_y * v2  + Pv_z * v3 ;

              v1 = r[idi].comvx ;
              v2 = r[idi].comvy ;
              v3 = r[idi].comvz ;
              r[idi].comvx = HH_x * v1 + HH_y * v2  + HH_z * v3 ;
              r[idi].comvy = wb_x * v1 + wb_y * v2  + wb_z * v3 ;
              r[idi].comvz = Pv_x * v1 + Pv_y * v2  + Pv_z * v3 ;
            }
        }
    }
    /*if(abs(r[id].vx * vAperfmstoamu  - v[id].x) > 1.E-4 ) checking condition */
}


void Print_trans_rot_transrot_cf(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  uint tcfl=2000;
  vector<double> comacf(tcfl,0.0);
  vector<double> angacf(tcfl,0.0);
  vector<double> ccfxx(tcfl,0.0);
  vector<double> ccfyy(tcfl,0.0);
  vector<double> ccfzz(tcfl,0.0);
  vector<double> ccfxy(tcfl,0.0);
  vector<double> ccfxz(tcfl,0.0);
  vector<double> ccfyx(tcfl,0.0);
  vector<double> ccfyz(tcfl,0.0);
  vector<double> ccfzx(tcfl,0.0);
  vector<double> ccfzy(tcfl,0.0);
  unsigned int from = 1, to = nsteps-10-tcfl ;

  // compute atomic velocity using CD scheme
  computeatomicvelocity(r, nsteps, natoms, L, dt); 
  computecomvelocity(r, nsteps, natoms, L, dt);

  // compute angular velocity and translational velocity in molecular frame
  computeangularvelocity(r, nsteps, natoms, L, dt);

  double mean_ = 0;
  for(unsigned int t = from;t < to; t += 100)
    { 
      for(unsigned int i = 0;i < natoms; i+=3)
        { 
              uint id = natoms*t+i;
              mean_ += 1; 
              for(unsigned int t_ = 0;t_ < tcfl;++t_)
                {
                  uint id_t = natoms*(t+t_)+i;

                  comacf[t_] += ( r[id].comvx * r[id_t].comvx   + r[id].comvy * r[id_t].comvy   + r[id].comvz * r[id_t].comvz  ) / 
                                ( r[id].comvx * r[id].comvx + r[id].comvy * r[id].comvy + r[id].comvz * r[id].comvz )  ; 

                  angacf[t_] += ( r[id].angvx * r[id_t].angvx   + r[id].angvy * r[id_t].angvy   + r[id].angvz * r[id_t].angvz  ) /
                                ( r[id].angvx * r[id].angvx + r[id].angvy * r[id].angvy + r[id].angvz * r[id].angvz )  ;

                  ccfxx[t_] +=  (r[id].angvx * r[id_t].comvx )/(abs(r[id].angvx) * abs(r[id].comvx) ) ; 
                  ccfyy[t_] +=  (r[id].angvy * r[id_t].comvy )/(abs(r[id].angvy) * abs(r[id].comvy) ) ; 
                  ccfzz[t_] +=  (r[id].angvz * r[id_t].comvz )/(abs(r[id].angvz) * abs(r[id].comvz) ) ; 

                  ccfxy[t_] +=  (r[id].angvx * r[id_t].comvy )/(abs(r[id].angvx) * abs(r[id].comvy) ) ; 
                  ccfxz[t_] +=  (r[id].angvx * r[id_t].comvz )/(abs(r[id].angvx) * abs(r[id].comvz) ) ; 

                  ccfyx[t_] +=  (r[id].angvy * r[id_t].comvx )/(abs(r[id].angvy) * abs(r[id].comvx) ) ; 
                  ccfyz[t_] +=  (r[id].angvy * r[id_t].comvz )/(abs(r[id].angvy) * abs(r[id].comvz) ) ; 

                  ccfzx[t_] +=  (r[id].angvz * r[id_t].comvx )/(abs(r[id].angvz) * abs(r[id].comvx) ) ; 
                  ccfzy[t_] +=  (r[id].angvz * r[id_t].comvy )/(abs(r[id].angvz) * abs(r[id].comvy) ) ; 
                }
         } 
    }


  ofstream outfile(filename);
  for(uint t = 1; t < tcfl-1;    t += 1 )
    {
      outfile << t*dt << "  "<< comacf[t]  / mean_ << "  " << 
                                angacf[t] / mean_ << "  " << 
                                ccfxx[t] / mean_ << "  " << 
                                ccfyy[t] / mean_ << "  " << 
                                ccfzz[t] / mean_ << "  " << 
                                ccfxy[t] / mean_ << "  " << 
                                ccfxz[t] / mean_ << "  " << 
                                ccfyx[t] / mean_ << "  " << 
                                ccfyz[t] / mean_ << "  " << 
                                ccfzx[t] / mean_ << "  " << 
                                ccfzy[t]    / mean_ << "  " << endl;
    }
  outfile.close();
  outfile.clear();
}



/* Finds the angle */
float angle_btwn_3points(const vector<Atom> &r, uint i, uint j1, uint j2, const vector<float> & L )
{
  float a = L[0], b = L[1], c = L[2];
  float x1,x2;
  float y1,y2;
  float z1,z2;

  x1 = r[j1].x - r[i].x;
  x2 = r[i].x - r[j2].x;
  y1 = r[j1].y - r[i].y;
  y2 = r[i].y - r[j2].y;
  z1 = r[j1].z - r[i].z;
  z2 = r[i].z - r[j2].z; 


  x1 = x1 - a * round(x1/a);
  x2 = x2 - a * round(x2/a);
  y1 = y1 - b * round(y1/b);
  y2 = y2 - b * round(y2/b);
  z1 = z1 - c * round(z1/c);
  z2 = z2 - c * round(z2/c);

  float top =  (x1*x2+y1*y2+z1*z2) - 0.0001 ;
  float bot =  norm(x1,y1,z1) * norm(x2,y2,z2);
  if( top == bot )
    {
      return 180; 
    }
  else
    {
      return  180 - (acos(top/bot) * 57.296);
    }
}

void PrintKEion(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_KE (nsteps, 0.0), trans_KE (nsteps, 0.0), rot_KE (nsteps, 0.0),
                total_KE_cation (nsteps, 0.0), trans_KE_cation (nsteps, 0.0), rot_KE_cation (nsteps, 0.0),
                total_KE_anion (nsteps, 0.0), trans_KE_anion (nsteps, 0.0), rot_KE_anion (nsteps, 0.0);

  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_cation2 = 0, count_anion2 = 0, count_both = 0, count_rem = 0, count_ions = 0, count_oxygen = 0;
  uint hbond_cation = 0, hbond_anion = 0, hbond_cation2 = 0, hbond_anion2 = 0, hbond_oxygen = 0;
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  double com[3]; 
  double vec[3]; 

  computeatomicvelocity(r, nsteps, natoms, L, dt); 
  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'M' || r[i].symbol[0] == 'N' || r[i].symbol[0] == 'S' || r[i].symbol[0] == 'C' || r[i].symbol[0] == 'F' ) 
        {

          count += 1;
          if(r[i].symbol[0] == 'M' || r[i].symbol[0] == 'N')
            {
              count_cation   += 1; 
            }
          else if(r[i].symbol[0] == 'S' || r[i].symbol[0] == 'C' || r[i].symbol[0] == 'F')
            {
              count_anion  += 1; 
            }

          for(uint t = 1; t < nsteps-1;  t += deltat )
            {
              uint id = natoms*t+i;

              /* translational and rotation KEs */
              com[0] = r[id].vx ;
              com[1] = r[id].vy ;
              com[2] = r[id].vz ;

          if(r[i].symbol[0] == 'M')
            {
              temp1         = jtohartree * 0.5 * 24.305 * amu *  norm2(com[0],com[1],com[2]) ;
              temp2         = jtohartree * 0.5 * 24.305 * amu * norm2(  r[id].vx,  r[id].vy,  r[id].vz) ;
            }
          if( r[i].symbol[0] == 'N')
            {
              temp1         = jtohartree * 0.5 *11 * amu *  norm2(com[0],com[1],com[2]) ;
              temp2         = jtohartree * 0.5 * 11 * amu * norm2(  r[id].vx,  r[id].vy,  r[id].vz) ;
            }
          if(r[i].symbol[0] == 'C')
            {
              temp1         = jtohartree * 0.5 * 35.453 * amu *  norm2(com[0],com[1],com[2]) ;
              temp2         = jtohartree * 0.5 * 35.453 * amu * norm2(  r[id].vx,  r[id].vy,  r[id].vz) ;
            }
          if(r[i].symbol[0] == 'S')
            {
              temp1         = jtohartree * 0.5 * 96.06 * amu *  norm2(com[0],com[1],com[2]) ;
              temp2         = jtohartree * 0.5 * 96.06 * amu * norm2(  r[id].vx,  r[id].vy,  r[id].vz) ;
            }



              trans_KE[t]  += temp1 ;
              total_KE[t]  += temp2 ;

               if(r[i].symbol[0] == 'M' || r[i].symbol[0] == 'N')
            {
                  trans_KE_cation[t]  += temp1 ;
                  total_KE_cation[t]  += temp2 ;
            }
          else if(r[i].symbol[0] == 'S' || r[i].symbol[0] == 'C' || r[i].symbol[0] == 'F')
            {
                  trans_KE_anion[t]  += temp1 ;
                  total_KE_anion[t]  += temp2 ;
            }

            }
        }
    }

  if(count == 0)        {count = 1;}
  if(count_rem == 0)    {count_rem = 1;}
  if(count_both == 0)   {count_both = 1;}
  if(count_anion == 0)  {count_anion = 1;}
  if(count_cation == 0) {count_cation = 1;}
  if(count_anion2 == 0)  {count_anion2 = 1;}
  if(count_cation2 == 0) {count_cation2 = 1;}
  if(count_ions == 0)   {count_ions = 1;}
  if(count_oxygen == 0) {count_oxygen = 1;}

  /* Total KE is in Atomic Unit in order to compare with CP2K data */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      if(total_KE[t] == 0)        {total_KE[t]        = 1;}
      if(total_KE_cation[t] == 0) {total_KE_cation[t] = 1;}
      if(total_KE_anion[t] == 0)  {total_KE_anion[t]  = 1;}


      outfile << t*dt << "  " << total_KE[t]       << "  " << trans_KE[t] / count              << "  " << rot_KE[t] / count            << "   " << 
                                total_KE_cation[t] << "  " << trans_KE_cation[t] / count_cation << "  " << rot_KE_cation[t] / count_cation << "   " <<
                                total_KE_anion[t]  << "  " << trans_KE_anion[t] / count_anion   << "  " << rot_KE_anion[t] / count_anion   << "   "  << endl;
    }
  outfile.close();
  outfile.clear();
}





void PrintKEnCosine(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  vector<float> total_KE (nsteps, 0.0), trans_KE (nsteps, 0.0), rot_KE (nsteps, 0.0),
                total_KE_cation (nsteps, 0.0), trans_KE_cation (nsteps, 0.0), rot_KE_cation (nsteps, 0.0),
                total_KE_anion (nsteps, 0.0), trans_KE_anion (nsteps, 0.0), rot_KE_anion (nsteps, 0.0),
                total_KE_cation2 (nsteps, 0.0), trans_KE_cation2 (nsteps, 0.0), rot_KE_cation2 (nsteps, 0.0),
                total_KE_anion2 (nsteps, 0.0), trans_KE_anion2 (nsteps, 0.0), rot_KE_anion2 (nsteps, 0.0),
                total_KE_both (nsteps, 0.0), trans_KE_both (nsteps, 0.0), rot_KE_both (nsteps, 0.0),
                total_KE_rem (nsteps, 0.0), trans_KE_rem (nsteps, 0.0), rot_KE_rem (nsteps, 0.0),
                total_KE_ions (nsteps, 0.0), trans_KE_ions (nsteps, 0.0), rot_KE_ions (nsteps, 0.0),
                total_KE_oxygen (nsteps, 0.0), trans_KE_oxygen (nsteps, 0.0), rot_KE_oxygen (nsteps, 0.0),
                cosine (nsteps, 0.0),        cosine2 (nsteps, 0.0),
                cosine_cation (nsteps, 0.0), cosine2_cation (nsteps, 0.0),
                cosine_anion (nsteps, 0.0),  cosine2_anion (nsteps, 0.0),
                cosine_both (nsteps, 0.0),   cosine2_both (nsteps, 0.0),
                cosine_rem (nsteps, 0.0),    cosine2_rem (nsteps, 0.0),
                cosine_ions (nsteps, 0.0),   cosine2_ions (nsteps, 0.0),
                cosine_oxygen (nsteps, 0.0), cosine2_oxygen (nsteps, 0.0);

  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_cation2 = 0, count_anion2 = 0, count_both = 0, count_rem = 0, count_ions = 0, count_oxygen = 0;
  uint hbond_cation = 0, hbond_anion = 0, hbond_cation2 = 0, hbond_anion2 = 0, hbond_oxygen = 0;
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  double com[3]; 
  double vec[3]; 

  computeatomicvelocity(r, nsteps, natoms, L, dt); 
  for(uint i = 0;i < natoms;++i)
    { 
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {

          hbond_cation = 0, hbond_anion = 0, hbond_oxygen = 0;
          hbond_cation2 = 0, hbond_anion2 = 0;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*0+i;  
              uint idi1 = natoms*0+i+1;          
              uint idi2 = natoms*0+i+2;                  
              uint idj = natoms*0+j;
 
              /*first solvation shell*/
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
              if(r[j].symbol[0] == 'S')
                {
                  x = min_distance(r[idj].x - r[idi1].x, L[0]);
                  y = min_distance(r[idj].y - r[idi1].y, L[1]);
                  z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                  rij1 = mindis(x,y,z,L);

                  x = min_distance(r[idj].x - r[idi2].x, L[0]);
                  y = min_distance(r[idj].y - r[idi2].y, L[1]);
                  z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                  rij2 = mindis(x,y,z,L); 
                  if( (rij1 < 3.75 && rij1 > 0) || (rij2 < 3.75 && rij2 > 0))
                    {
                      hbond_anion += 1; 
                    }      
                }

              /*second solvation shell*/
              if(r[j].symbol[0] == 'M' || r[j].symbol[0] == 'N')
                {
                  x = min_distance(r[idj].x - r[idi].x, L[0]);
                  y = min_distance(r[idj].y - r[idi].y, L[1]);
                  z = min_distance(r[idj].z - r[idi].z, L[2]); 
                  rij = mindis(x,y,z,L); 
                  if(rij < 6 && rij > 3.2)
                    {
                      hbond_cation2 += 1;
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
                  if( (rij1 < 4.2 && rij1 > 3.0) || (rij2 < 4.2 && rij2 > 3.0))
                    {
                      hbond_anion2 += 1;
                    }      
                }
              if(r[j].symbol[0] == 'S')
                {
                  x = min_distance(r[idj].x - r[idi1].x, L[0]);
                  y = min_distance(r[idj].y - r[idi1].y, L[1]);
                  z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                  rij1 = mindis(x,y,z,L);

                  x = min_distance(r[idj].x - r[idi2].x, L[0]);
                  y = min_distance(r[idj].y - r[idi2].y, L[1]);
                  z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                  rij2 = mindis(x,y,z,L);  
                  if( (rij1 < 6 && rij1 > 3.75) || (rij2 < 6 && rij2 > 3.75))
                    {
                      hbond_anion2 += 1;
                    }      
                }

              if(i != j  && r[j].symbol[0] == 'O' && r[j+1].symbol[0] == 'H' && r[j+2].symbol[0] == 'H')
                {
                  x = min_distance(r[idj].x - r[idi].x, L[0]);
                  y = min_distance(r[idj].y - r[idi].y, L[1]);
                  z = min_distance(r[idj].z - r[idi].z, L[2]); 
                  rij = mindis(x,y,z,L); 
                  if( rij < 3.5  && rij > 0 && ( angle_btwn_3points(r,idi,idi1,idj, L) < 30 || angle_btwn_3points(r,idi,idi2,idj, L) < 30) )
                    {
                      hbond_oxygen += 1;
                    }      
                }
            }

          /*counting first solvation shell*/
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

          if(hbond_cation > 0 || hbond_anion > 0  || (hbond_cation + hbond_anion) > 0 )
            {
              count_ions += 1; 
            }


          /*counting second solvation shell*/
          if(hbond_cation2 == 0 && hbond_anion2 > 0)
            {
              count_anion2 += 1; 
            }
          if(hbond_cation2 > 0 && hbond_anion2 == 0)
            {
              count_cation2 += 1; 
            }

          if(hbond_oxygen > 0 && hbond_cation == 0 && hbond_anion == 0)
            {
              count_oxygen += 1; 
            }

          for(uint t = 1; t < nsteps-1;  t += deltat )
            {
              uint id = natoms*t+i;

              /* translational and rotation KEs */
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

              /* cosine and cosine square theta */
              vec[0] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].x - r[id].x)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].x - r[id].x) ;

              vec[1] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].y - r[id].y)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].y - r[id].y) ;

              vec[2] =  mindis(r[id].x - r[id+1].x, r[id].y - r[id+1].y, r[id].z - r[id+1].z, L) * (r[id+1].z - r[id].z)
                      + mindis(r[id].x - r[id+2].x, r[id].y - r[id+2].y, r[id].z - r[id+2].z, L) * (r[id+2].z - r[id].z) ;

              temp              =  vec[0] / pow(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2],0.5); 
              cosine[t]        += temp; 
              cosine2[t]       += pow(temp,2.0);

              /*first solvation shell*/
              if(hbond_cation == 0 && hbond_anion == 0)
                {
                  trans_KE_rem[t]  += temp1 ;
                  total_KE_rem[t]  += temp2 ;
                  rot_KE_rem[t]    += temp3 ;
                  cosine_rem[t]    += temp ; 
                  cosine2_rem[t]   += pow(temp,2.0) ; 
                }
              else if(hbond_cation > 0 && hbond_anion > 0)
                {
                  trans_KE_both[t]  += temp1 ;
                  total_KE_both[t]  += temp2 ;
                  rot_KE_both[t]    += temp3 ;
                  cosine_both[t]    += temp ; 
                  cosine2_both[t]   += pow(temp,2.0) ; 
                }
              else if(hbond_cation == 0 && hbond_anion > 0)
                {
                  trans_KE_anion[t]  += temp1 ;
                  total_KE_anion[t]  += temp2 ;
                  rot_KE_anion[t]    += temp3 ;
                  cosine_anion[t]    += temp ; 
                  cosine2_anion[t]   += pow(temp,2.0) ; 
                }
              else if(hbond_cation > 0 && hbond_anion == 0)
                {
                  trans_KE_cation[t]  += temp1 ;
                  total_KE_cation[t]  += temp2 ;
                  rot_KE_cation[t]    += temp3 ;
                  cosine_cation[t]    += temp ; 
                  cosine2_cation[t]   += pow(temp,2.0) ; 
                }

              /*second solvation shell*/
              if(hbond_cation2 == 0 && hbond_anion2 > 0)
                {
                  trans_KE_anion2[t]  += temp1 ;
                  total_KE_anion2[t]  += temp2 ;
                  rot_KE_anion2[t]    += temp3 ;
                }
              else if(hbond_cation2 > 0 && hbond_anion2 == 0)
                {
                  trans_KE_cation2[t]  += temp1 ;
                  total_KE_cation2[t]  += temp2 ;
                  rot_KE_cation2[t]    += temp3 ;
                }


              if(hbond_cation > 0 || hbond_anion > 0  || (hbond_cation + hbond_anion) > 0 )
                {
                  trans_KE_ions[t]  += temp1 ;
                  total_KE_ions[t]  += temp2 ;
                  rot_KE_ions[t]    += temp3 ;
                  cosine_ions[t]    += temp ; 
                  cosine2_ions[t]   += pow(temp,2.0)  ; 
                }

              if(hbond_oxygen > 0 && hbond_cation == 0 && hbond_anion == 0)
                {
                  trans_KE_oxygen[t]  += temp1 ;
                  total_KE_oxygen[t]  += temp2 ;
                  rot_KE_oxygen[t]    += temp3 ;
                  cosine_oxygen[t]    += temp ;
                  cosine2_oxygen[t]   += pow(temp,2.0)  ;  
                }

            }
        }
    }

  if(count == 0)        {count = 1;}
  if(count_rem == 0)    {count_rem = 1;}
  if(count_both == 0)   {count_both = 1;}
  if(count_anion == 0)  {count_anion = 1;}
  if(count_cation == 0) {count_cation = 1;}
  if(count_anion2 == 0)  {count_anion2 = 1;}
  if(count_cation2 == 0) {count_cation2 = 1;}
  if(count_ions == 0)   {count_ions = 1;}
  if(count_oxygen == 0) {count_oxygen = 1;}

  /* Total KE is in Atomic Unit in order to compare with CP2K data */
  ofstream outfile(filename);
  for(uint t = 1; t < nsteps-1;  t += deltat )
    {
      if(total_KE[t] == 0)        {total_KE[t]        = 1;}
      if(total_KE_cation[t] == 0) {total_KE_cation[t] = 1;}
      if(total_KE_anion[t] == 0)  {total_KE_anion[t]  = 1;}
      if(total_KE_cation2[t] == 0) {total_KE_cation2[t] = 1;}
      if(total_KE_anion2[t] == 0)  {total_KE_anion2[t]  = 1;}
      if(total_KE_both[t] == 0)   {total_KE_both[t]   = 1;}
      if(total_KE_ions[t] == 0)   {total_KE_ions[t]   = 1;}
      if(total_KE_rem[t] == 0)    {total_KE_rem[t]    = 1;}
      if(total_KE_oxygen[t] == 0) {total_KE_oxygen[t]    = 1;}

      outfile << t*dt << "  " << total_KE[t]       << "  " << trans_KE[t] / total_KE[t]               << "  " << rot_KE[t] / total_KE[t]               << "   " << 
                                total_KE_cation[t] << "  " << trans_KE_cation[t] / total_KE_cation[t] << "  " << rot_KE_cation[t] / total_KE_cation[t] << "   " <<
                                total_KE_anion[t]  << "  " << trans_KE_anion[t] / total_KE_anion[t]   << "  " << rot_KE_anion[t] / total_KE_anion[t]   << "   " <<

                                total_KE_cation2[t] << "  " << trans_KE_cation2[t] / total_KE_cation2[t] << "  " << rot_KE_cation2[t] / total_KE_cation2[t] << "   " <<
                                total_KE_anion2[t]  << "  " << trans_KE_anion2[t] / total_KE_anion2[t]   << "  " << rot_KE_anion2[t] / total_KE_anion2[t]   << "   " <<

                                total_KE_both[t]   << "  " << trans_KE_both[t] / total_KE_both[t]     << "  " << rot_KE_both[t] / total_KE_both[t]     << "   " <<
                                total_KE_ions[t]   << "  " << trans_KE_ions[t] / total_KE_ions[t]     << "  " << rot_KE_ions[t] / total_KE_ions[t]     << "   " <<
                                total_KE_oxygen[t] << "  " << trans_KE_oxygen[t] / total_KE_oxygen[t] << "  " << rot_KE_oxygen[t] / total_KE_oxygen[t] << "   " <<
                                total_KE_rem[t]    << "  " << trans_KE_rem[t] / total_KE_rem[t]       << "  " << rot_KE_rem[t] / total_KE_rem[t]       << "  " << 
                                cosine[t]        / count << "  " << 
                                cosine2[t]       / count << "  " << 
                                cosine_cation[t] / count_cation << "  " << 
                                cosine2_cation[t]/ count_cation << "  " << 
                                cosine_anion[t]  / count_anion << "  " << 
                                cosine2_anion[t] / count_anion << "  " << 
                                cosine_both[t]   / count_both << "  " <<  
                                cosine2_both[t]  / count_both << "  " <<  
                                cosine_ions[t]   / count_ions << "  " <<  
                                cosine2_ions[t]  / count_ions << "  " <<  
                                cosine_oxygen[t] / count_oxygen << "  " << 
                                cosine2_oxygen[t]/ count_oxygen << "  " <<   
                                cosine_rem[t]    / count_rem     << "  " <<
                                cosine2_rem[t]   / count_rem   << endl;
    }
  outfile.close();
  outfile.clear();
}



void Printwaterioncoordination(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_both = 0, count_rem = 0, count_ions = 0;
  float hbond_cation = 0, hbond_anion = 0;
  double vec[3];  

  for(uint t = 1; t < nsteps-1; ++t )
    {

      for(uint j = 0;j < natoms;++j)
        { 
          if(r[j].symbol[0] == 'M' || r[j].symbol[0] == 'N')
            {
              count_cation += 1; 
              for(uint i = 0;i < natoms;++i)
                { 
                  if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
                    {

                      uint idj = natoms*t+j;
                      uint idi = natoms*t+i;  
                      uint idi1 = natoms*t+i+1;          
                      uint idi2 = natoms*t+i+2;                  
                      x = min_distance(r[idj].x - r[idi].x, L[0]);
                      y = min_distance(r[idj].y - r[idi].y, L[1]);
                      z = min_distance(r[idj].z - r[idi].z, L[2]); 
                      rij = mindis(x,y,z,L); 
                      if(rij < 3.0 && rij > 0)
                        {
                          hbond_cation += 1; 
                        }  
                    } 
                }
            }

          if(r[j].symbol[0] == 'C' || r[j].symbol[0] == 'F')
            {
              count_anion += 1;
              for(uint i = 0;i < natoms;++i)
                { 
                  if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
                    {
                      uint idj = natoms*t+j;
                      uint idi = natoms*t+i;  
                      uint idi1 = natoms*t+i+1;          
                      uint idi2 = natoms*t+i+2;                  
 
                      x = min_distance(r[idj].x - r[idi1].x, L[0]);
                      y = min_distance(r[idj].y - r[idi1].y, L[1]);
                      z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                      rij1 = mindis(x,y,z,L);

                      x = min_distance(r[idj].x - r[idi2].x, L[0]);
                      y = min_distance(r[idj].y - r[idi2].y, L[1]);
                      z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                      rij2 = mindis(x,y,z,L); 
                      if( (rij1 < 3.1 && rij1 > 0) || (rij2 < 3.1 && rij2 > 0))
                        {
                          hbond_anion += 1;
                        }  
                    } 
                } 
            }
        } 
    }

      cout << "number Of H2O aroud cation  = " << hbond_cation / count_cation << "\n"  
           << "number Of H2O aroud anion   = " << hbond_anion / count_anion << "\n"  ;                            
}

void classifywater(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt)
{

  uint idi = 0, idj = 0;
  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;;
  uint hbond_cation = 0, hbond_anion = 0, hbond_cation2 = 0, hbond_anion2 = 0, hbond_oxygen = 0;
 
  for(uint t = 0; t < nsteps; t += deltat )
    {
  for(uint i = 0;i < natoms;++i)
    { 
      uint idi = natoms*t+i;  
      uint idi1 = natoms*t+i+1;          
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          hbond_cation = 0, hbond_anion = 0, hbond_oxygen = 0;
          hbond_cation2 = 0, hbond_anion2 = 0;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi2 = natoms*t+i+2;                  
              uint idj = natoms*t+j;
 
              /*first solvation shell*/
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
              if(r[j].symbol[0] == 'S')
                {
                  x = min_distance(r[idj].x - r[idi1].x, L[0]);
                  y = min_distance(r[idj].y - r[idi1].y, L[1]);
                  z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                  rij1 = mindis(x,y,z,L);

                  x = min_distance(r[idj].x - r[idi2].x, L[0]);
                  y = min_distance(r[idj].y - r[idi2].y, L[1]);
                  z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                  rij2 = mindis(x,y,z,L); 
                  if( (rij1 < 3.75 && rij1 > 0) || (rij2 < 3.75 && rij2 > 0))
                    {
                      hbond_anion += 1; 
                    }      
                }
            }

          if(hbond_cation == 0 && hbond_anion == 0)
            {
              r[idi].watertype = "remwater"; 
            }
          else if(hbond_cation > 0 && hbond_anion > 0)
            {
              r[idi].watertype = "waterbetweencationandanion"; 
            }
          else if(hbond_cation == 0 && hbond_anion > 0)
            {
              r[idi].watertype = "wateraroundanion"; 
            }
          else if(hbond_cation > 0 && hbond_anion == 0)
            {
              r[idi].watertype = "wateraroundcation"; 
            }
        }
    }
    }
}


void PrintOOdistance1(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
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
  //cout << "g \u03B1 - \u03B2 (rij)" << endl; ;
  //cout << "Enter number of \u03B1 species : " ;
  //cin >>  no_of_residue_first ;
  //cout << "Enter number of \u03B2 species : " ;
  //cin >>  no_of_residue_second ;

  string whichtype = "waterbetweencationandanion";
  for(uint t = 0; t < nsteps; t += deltat )  
    {
      cout << t << endl;
      count += 1;      
      for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
        {
          count_atoms = 0;
          no_of_residue_first = 0;        
          for(uint i = 0;i < natoms;++i)
            {
              idi = natoms*t+i;  
              if(r[idi].watertype == whichtype && r[idi].symbol[0] == 'O' && r[idi+1].symbol[0] == 'H' ) {no_of_residue_first+=1;}
              no_of_residue_second = 0;        
              for(uint j = 0;j < natoms;++j)
                {
                  idj = natoms*t+j;
                  if(r[idj].watertype != whichtype  && r[idj].symbol[0] == 'O' && r[idj+1].symbol[0] == 'H'  ) {no_of_residue_second+=1;}
                  if (r[idi].symbol[0] == 'O' && r[idj].symbol[0] == 'O' &&
                      r[idi+1].symbol[0] == 'H' && r[idj+1].symbol[0] == 'H' && 
                      r[idi+2].symbol[0] == 'H' && r[idj+2].symbol[0] == 'H'  &&
                       i != j  &&
                       r[idi].watertype == whichtype &&
                       r[idj].watertype != whichtype  )
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
          RDF[RDF_i]  += count_atoms;
        }
    }

  // consine thetha per water molecule
  ofstream outfile(whichtype);
  for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
    {
      float bulkdensity = no_of_residue_second/V ; 
      outfile << rad[RDF_i]  << "  " << RDF[RDF_i] / (bulkdensity * 4.0 * PI * rad[RDF_i] * rad[RDF_i] *  RDF_h * count * no_of_residue_first)  << endl;
    }
  outfile.close();
  outfile.clear();
}



void Printtrans_rot_ccfn(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{

  uint tcfl = 2000; 
  uint time_cond = 1; 
  vector<float> total_ccfn (tcfl, 0.0), 
                ccfn_cation (tcfl, 0.0),
                ccfn_anion (tcfl, 0.0), 
                ccfn_both (tcfl, 0.0), 
                ccfn_rem (tcfl, 0.0);

  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_both = 0, count_rem = 0;
  uint hbond_cation = 0, hbond_anion = 0;
  float am_H = 1 * amu, am_O = 16 * amu, am_H2O = 18 * amu ;
  double comv_t0[3]; 
  double comr_t0[3]; 
  double comv_t[3]; 
  double comr_t[3]; 
  double angv_t0[3]; 
  int count_tcfl=0;


  computeatomicvelocity(r, nsteps, natoms, L, dt); 
  for(uint i = 0;i < natoms;++i)
    { cout << i << endl;
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          hbond_cation = 0, hbond_anion = 0;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*time_cond+i;  
              uint idi1 = natoms*time_cond+i+1;          
              uint idi2 = natoms*time_cond+i+2;                  
              uint idj = natoms*time_cond+j;
 
              /*first solvation shell*/
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

          /*counting first solvation shell*/
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

          count_tcfl=0; 
          for(uint time_cond = 1; time_cond < nsteps - tcfl - 2000 ;  time_cond += 100 ){
          count_tcfl+=1;
          uint indexstart = 0;
      
              uint id_t0 = natoms*time_cond+i;
              comv_t0[0] = r[id_t0].vx * (am_O/am_H2O) +  r[id_t0+1].vx * (am_H/am_H2O) +  r[id_t0+2].vx * (am_H/am_H2O) ;
              comv_t0[1] = r[id_t0].vy * (am_O/am_H2O) +  r[id_t0+1].vy * (am_H/am_H2O) +  r[id_t0+2].vy * (am_H/am_H2O) ;
              comv_t0[2] = r[id_t0].vz * (am_O/am_H2O) +  r[id_t0+1].vz * (am_H/am_H2O) +  r[id_t0+2].vz * (am_H/am_H2O) ;


              comr_t0[0] = r[id_t0].x * (am_O/am_H2O) +  r[id_t0+1].x * (am_H/am_H2O) +  r[id_t0+2].x * (am_H/am_H2O) ;
              comr_t0[1] = r[id_t0].y * (am_O/am_H2O) +  r[id_t0+1].y * (am_H/am_H2O) +  r[id_t0+2].y * (am_H/am_H2O) ;
              comr_t0[2] = r[id_t0].z * (am_O/am_H2O) +  r[id_t0+1].z * (am_H/am_H2O) +  r[id_t0+2].z * (am_H/am_H2O) ;

          for(uint t = time_cond; t < (tcfl + time_cond) - 1 ;  t += 1 )
            {

              uint id = natoms*t+i;
              comv_t[0] = r[id].vx * (am_O/am_H2O) +  r[id+1].vx * (am_H/am_H2O) +  r[id+2].vx * (am_H/am_H2O) ;
              comv_t[1] = r[id].vy * (am_O/am_H2O) +  r[id+1].vy * (am_H/am_H2O) +  r[id+2].vy * (am_H/am_H2O) ;
              comv_t[2] = r[id].vz * (am_O/am_H2O) +  r[id+1].vz * (am_H/am_H2O) +  r[id+2].vz * (am_H/am_H2O) ;


              comr_t[0] = r[id].x * (am_O/am_H2O) +  r[id+1].x * (am_H/am_H2O) +  r[id+2].x * (am_H/am_H2O) ;
              comr_t[1] = r[id].y * (am_O/am_H2O) +  r[id+1].y * (am_H/am_H2O) +  r[id+2].y * (am_H/am_H2O) ;
              comr_t[2] = r[id].z * (am_O/am_H2O) +  r[id+1].z * (am_H/am_H2O) +  r[id+2].z * (am_H/am_H2O) ;
            

              /*angular velocity*/
              double vect_A[3], vect_B[3], normr ;
              vect_A[0] = r[id_t0].x - comr_t0[0] ; 
              vect_A[1] = r[id_t0].y - comr_t0[1] ; 
              vect_A[2] = r[id_t0].z - comr_t0[2] ;
              vect_B[0] = r[id_t0].vx - comv_t0[0]; 
              vect_B[1] = r[id_t0].vy - comv_t0[1]; 
              vect_B[2] = r[id_t0].vz - comv_t0[2] ;
              normr =  pow(vect_A[0],2) + pow(vect_A[1],2) + pow(vect_A[2],2) ; 
              angv_t0[0] = (vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]) / normr ;
              angv_t0[1] = (vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]) / normr ;
              angv_t0[2] = (vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]) / normr ;

              vect_A[0] = r[id_t0+1].x - comr_t0[0] ; 
              vect_A[1] = r[id_t0+1].y - comr_t0[1] ; 
              vect_A[2] = r[id_t0+1].z - comr_t0[2] ;
              vect_B[0] = r[id_t0+1].vx - comv_t0[0]; 
              vect_B[1] = r[id_t0+1].vy - comv_t0[1]; 
              vect_B[2] = r[id_t0+1].vz - comv_t0[2] ;
              normr =  pow(vect_A[0],2) + pow(vect_A[1],2) + pow(vect_A[2],2) ;
              angv_t0[0] += (vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]) / normr ;
              angv_t0[1] += (vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]) / normr ;
              angv_t0[2] += (vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]) / normr ;

              vect_A[0] = r[id_t0+2].x - comr_t0[0] ; 
              vect_A[1] = r[id_t0+2].y - comr_t0[1] ; 
              vect_A[2] = r[id_t0+2].z - comr_t0[2] ;
              vect_B[0] = r[id_t0+2].vx - comv_t0[0]; 
              vect_B[1] = r[id_t0+2].vy - comv_t0[1]; 
              vect_B[2] = r[id_t0+2].vz - comv_t0[2] ;
              normr =  pow(vect_A[0],2) + pow(vect_A[1],2) + pow(vect_A[2],2) ;
              angv_t0[0] += (vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1]) / normr ;
              angv_t0[1] += (vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2]) / normr ;
              angv_t0[2] += (vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0]) / normr ;


              double a1, a2, a3 ;
              a1 = angv_t0[2] ;
              a2 = comv_t0[1] ;
              a3 = comv_t[1] ; 

              temp         = (comv_t0[0]*comv_t[0]+comv_t0[1]*comv_t[1]+comv_t0[2]*comv_t[2]) / (comv_t0[0]*comv_t0[0]+comv_t0[1]*comv_t0[1]+comv_t0[2]*comv_t0[2]);

              //OH vector decay
              if(1){
                 float a = L[0], b = L[1], c = L[2] ; 
                 float x1,x2,y1,y2,z1,z2,top,bot;

                x1 = r[id_t0+1].x - r[id_t0+0].x ;
                x2 = r[id_t0+0].x - r[id+1].x;

                y1 = r[id_t0+1].y - r[id_t0+0].y;
                y2 = r[id_t0+0].y - r[id+1].y ;

                z1 = r[id_t0+1].z - r[id_t0+0].z;
                z2 = r[id_t0+0].z - r[id+1].z; 

                x1 = x1 - a * round(x1/a);
                x2 = x2 - a * round(x2/a);
                y1 = y1 - b * round(y1/b);
                y2 = y2 - b * round(y2/b);
                z1 = z1 - c * round(z1/c);
                z2 = z2 - c * round(z2/c);

                top =  (x1*x2+y1*y2+z1*z2) ;// - 0.0001 ;
                bot =  norm(x1,y1,z1) * norm(x2,y2,z2);
                float ang = 0;
                if( top == bot )
                {
                    ang =  180; 
                 }
                else
                {
                 ang =   180 - (acos(top/bot) * 57.296);
                }

               if(isnan(ang)) {ang = 180;}
               temp   = cos(ang * PI / 180.0); //cout << vec[t] << endl; 
              }
  
              //temp         = (a1 * a3 ) / (pow(pow(a1,2),0.5) * pow(pow(a2,2),0.5));
              //cout << r[id+0].vx << "  " << r[id+0].vy << "  " << r[id+0].vz << endl;
              //cout << a1 << " " << a2 << "  " << a3 <<   "  " <<  a1*a3 << "  " << (pow(pow(a1,2),0.5) * pow(pow(a2,2),0.5))  << "  " << temp << endl;

              total_ccfn[indexstart]  += temp ;

              /*first solvation shell*/
              if(hbond_cation == 0 && hbond_anion == 0)
                {
                  ccfn_rem[indexstart]  += temp ;
                }
              else if(hbond_cation > 0 && hbond_anion > 0)
                {
                  ccfn_both[indexstart]  += temp ;
                }
              else if(hbond_cation == 0 && hbond_anion > 0)
                {
                  ccfn_anion[indexstart]  += temp ;
                }
              else if(hbond_cation > 0 && hbond_anion == 0)
                {
                  ccfn_cation[indexstart]  += temp ;
                }
              indexstart += 1;
            }}
        }
    }

  if(count == 0)        {count = 1;}
  if(count_rem == 0)    {count_rem = 1;}
  if(count_both == 0)   {count_both = 1;}
  if(count_anion == 0)  {count_anion = 1;}
  if(count_cation == 0) {count_cation = 1;}

  ofstream outfile(filename);
  for(uint t = 1; t < tcfl-1;    t += 1 )
    {
      if(total_ccfn[t] == 0)    {total_ccfn[t]        = 1;}
      if(ccfn_rem[t] == 0)      {ccfn_rem[t] = 1;}
      if(ccfn_both[t] == 0)     {ccfn_both[t]  = 1;}
      if(ccfn_anion[t] == 0)    {ccfn_anion[t] = 1;}
      if(ccfn_cation[t] == 0)   {ccfn_cation[t]  = 1;}

      outfile << t*dt << "  "<< total_ccfn[t]  / (count * count_tcfl) << "  " << 
                                ccfn_cation[t] / count_cation << "  " << 
                                ccfn_anion[t]  / count_anion << "  " << 
                                ccfn_both[t]   / count_both << "  " << 
                                ccfn_rem[t]    / count_rem << "  " << endl;
    }
  outfile.close();
  outfile.clear();
}



void PrintOOdistance(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string filename)
{
  float V = L[0] * L[1] * L[2], len = max(L[0], L[1]), RDF_h = 0.1 ;
  uint idi = 0, idj = 0, ncell = 0, count_atoms = 0;
  len = max(len, L[2]);
  len = len/2.0;
  uint RDF_size = len*100;
  vector<float> cation(RDF_size,0.0), anion(RDF_size,0.0), both(RDF_size,0.0), rem(RDF_size,0.0), total(RDF_size,0.0), rad(RDF_size,0.0);
  for(uint i = 1;i < RDF_size ;i++)
    {
      rad[i]          = i * 0.01;
    }   

  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;
  float count = 0, count_cation = 0, count_anion = 0, count_cation2 = 0, count_anion2 = 0, count_both = 0, count_rem = 0, count_ions = 0, count_oxygen = 0;
  uint hbond_cation = 0, hbond_anion = 0, hbond_cation2 = 0, hbond_anion2 = 0, hbond_oxygen = 0;
 
  for(uint i = 0;i < natoms;++i)
    { 
      //cout << i << endl;
      if(r[i].symbol[0] == 'O' && r[i+1].symbol[0] == 'H' && r[i+2].symbol[0] == 'H') 
        {
          hbond_cation = 0, hbond_anion = 0, hbond_oxygen = 0;
          hbond_cation2 = 0, hbond_anion2 = 0;
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*0+i;  
              uint idi1 = natoms*0+i+1;          
              uint idi2 = natoms*0+i+2;                  
              uint idj = natoms*0+j;
 
              /*first solvation shell*/
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
              if(r[j].symbol[0] == 'S')
                {
                  x = min_distance(r[idj].x - r[idi1].x, L[0]);
                  y = min_distance(r[idj].y - r[idi1].y, L[1]);
                  z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                  rij1 = mindis(x,y,z,L);

                  x = min_distance(r[idj].x - r[idi2].x, L[0]);
                  y = min_distance(r[idj].y - r[idi2].y, L[1]);
                  z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                  rij2 = mindis(x,y,z,L); 
                  if( (rij1 < 3.75 && rij1 > 0) || (rij2 < 3.75 && rij2 > 0))
                    {
                      hbond_anion += 1; 
                    }      
                }
            }

          /*counting first solvation shell*/
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


      for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
        {
          if(rad[RDF_i] > 2.2 && rad[RDF_i] < 3.5){
          temp = 0;
          for(uint t = 1; t < nsteps-1;  t += deltat)
            {
              idi = natoms*t+i;          
              for(uint j = 0;j < natoms;++j)
                {
                  idj = natoms*t+j; 
                  // donating + accepting
                  if (r[idi].symbol[0] == 'O' && r[idj].symbol[0] == 'O' && r[idj+1].symbol[0] == 'H' && i != j)
                    {                  
                      x = min_distance(r[idj].x - r[idi].x, L[0]);
                      y = min_distance(r[idj].y - r[idi].y, L[1]);
                      z = min_distance(r[idj].z - r[idi].z, L[2]); 
                      rij = mindis(x,y,z,L);     
  
                      if (rij <= rad[RDF_i] + (RDF_h/2.0) && rij > rad[RDF_i] - (RDF_h/2.0) && rij < 3.5  && rij > 0 && ( angle_btwn_3points(r,idi,idi+1,idj, L) < 30 || angle_btwn_3points(r,idi,idi+2,idj, L) < 30 || angle_btwn_3points(r,idj,idj+1,idi, L) < 30 || angle_btwn_3points(r,idj,idj+2,idi, L) < 30)  )
                        {
                          temp += 1; // cout << rij << endl;
                        }
                    }
                }
            }}

                  total[RDF_i]  += temp    ; 

              /*first solvation shell*/
              if(hbond_cation == 0 && hbond_anion == 0)
                {
                  rem[RDF_i]  += temp ;
                }
              else if(hbond_cation > 0 && hbond_anion > 0)
                {
                  both[RDF_i]  += temp ;
                }
              else if(hbond_cation == 0 && hbond_anion > 0)
                {
                  anion[RDF_i]  += temp ;
                }
              else if(hbond_cation > 0 && hbond_anion == 0)
                {
                  cation[RDF_i]  += temp ;
                }
        }
        } //if condition ends
    }

  if(count == 0)        {count = 1;}
  if(count_rem == 0)    {count_rem = 1;}
  if(count_both == 0)   {count_both = 1;}
  if(count_anion == 0)  {count_anion = 1;}
  if(count_cation == 0) {count_cation = 1;}
  if(count_anion2 == 0)  {count_anion2 = 1;}
  if(count_cation2 == 0) {count_cation2 = 1;}
  if(count_ions == 0)   {count_ions = 1;}
  if(count_oxygen == 0) {count_oxygen = 1;}

  double sum1=0;
  double sum2=0;
  double sum3=0;
  double sum4=0;
  double sum5=0;


  for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
    {
      sum1 +=total[RDF_i] / count ;
      sum2 +=cation[RDF_i] / count_cation;
      sum3 +=anion[RDF_i] / count_anion;
      sum4 +=both[RDF_i] / count_both;
      sum5 +=rem[RDF_i] / count_rem;
    }

  if(sum1 == 0)        {sum1 = 1;}
  if(sum2 == 0)        {sum2 = 1;}
  if(sum3 == 0)        {sum3 = 1;}
  if(sum4 == 0)        {sum4 = 1;}
  if(sum5 == 0)        {sum5 = 1;}


  ofstream outfile(filename);
  for(uint RDF_i = 1;RDF_i < RDF_size - 1; ++RDF_i)
    {
      double bulkden = 1 ;
      outfile << rad[RDF_i]  << "  " << total[RDF_i] / (count * sum1 * bulkden )  << "  " 
                                     << cation[RDF_i] / (count_cation * sum2* bulkden)  << "  " 
                                     << anion[RDF_i] / (count_anion * sum3 * bulkden )  << "  " 
                                     << both[RDF_i]  /(count_both * sum4 * bulkden )  << "  " 
                                     << rem[RDF_i] / (count_rem * sum5 * bulkden )  << endl;
    }
  outfile.close();
  outfile.clear();
}

