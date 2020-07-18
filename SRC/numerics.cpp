#include "../include/numerics.h"
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

float min_distance(float r, float l){ return r - l * round(r/l);}
float mindis(float dx,float dy,float dz,float L){  return norm(dx - L * round(dx/L),dy - L * round(dy/L), dz - L * round(dz/L));}
float norm(float x,float y,float z){  return sqrt(x*x+y*y+z*z);}
