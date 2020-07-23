#ifndef NUMERICS_H
#define NUMERICS_H

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

#include "../include/io.h"

using namespace std;



float min_distance(float r, float l);
float mindis(float dx,float dy,float dz, const vector<float> & L);
float norm(float x,float y,float z);



#endif
