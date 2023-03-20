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



double min_distance(double r, float l);
double mindis(double dx,double dy,double dz, const vector<float> & L);
double norm(double x,double y,double z);
double norm2(double x,double y,double z);
void convertounivector(double &b_x,double &b_y,double &b_z);
void FFT(vector<double> &vvacf,vector<double> & vvacf_r_fft,vector<double> & vvacf_i_fft, float tcfl, float dt);


#endif
