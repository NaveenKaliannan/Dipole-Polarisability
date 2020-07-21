
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

float min_distance(float r, float l){ return r - l * round(r/l);}
float mindis(float dx,float dy,float dz, const vector<float> & L){  return norm(dx - L[0] * round(dx/L[0]),dy - L[1] * round(dy/L[1]), dz - L[2] * round(dz/L[2]));}
float norm(float x,float y,float z){  return sqrt(x*x+y*y+z*z);}
