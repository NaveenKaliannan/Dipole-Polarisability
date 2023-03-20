
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

// Constants


double min_distance(double r, float l){ return r - l * round(r/l);}
double norm(double x,double y,double z){  return sqrt(x*x+y*y+z*z);}
double mindis(double dx,double dy,double dz, const vector<float> & L)
{ 
  return norm(dx - L[0] * round(dx/L[0]),dy - L[1] * round(dy/L[1]), dz - L[2] * round(dz/L[2]));
}


double norm2(double x,double y,double z)
{
  return x*x+y*y+z*z;
}


void convertounivector(double &b_x,double &b_y,double &b_z)
{
  double sum  =  pow(b_x * b_x  + b_y * b_y + b_z * b_z, 0.5);
  b_x = b_x / sum ;
  b_y = b_y / sum ;
  b_z = b_z / sum ;
}


void FFT(vector<double> &vvacf,vector<double> & vvacf_r_fft,vector<double> & vvacf_i_fft, float tcfl, float dt)
{
  for(unsigned int freq = 0; freq < Nfreq;++freq)
    {
      for(unsigned int t_ = 0;t_ < tcfl;++t_)
        {
          vvacf_r_fft[freq] = vvacf_r_fft[freq] + vvacf[t_] * cos(2.0*PI*light_vel*t_*dt*freq) * dt;
          vvacf_i_fft[freq] = vvacf_i_fft[freq] + vvacf[t_] * sin(2.0*PI*light_vel*t_*dt*freq) * dt;
        }
    }

}

