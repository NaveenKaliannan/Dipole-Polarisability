
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


double min_distance(double r, float l){ return r - l * round(r/l);}
double mindis(double dx,double dy,double dz, const vector<float> & L){  return norm(dx - L[0] * round(dx/L[0]),dy - L[1] * round(dy/L[1]), dz - L[2] * round(dz/L[2]));}
double norm(double x,double y,double z){  return sqrt(x*x+y*y+z*z);}

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

