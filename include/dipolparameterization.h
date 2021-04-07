#ifndef DIPOLPARAMETERIZATION_H
#define DIPOLPARAMETERIZATION_H

#include<vector>
#include<fstream>
#include "../include/numerics.h"
#include "../include/io.h"

using namespace std;


void printsinglewatercoordinates(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L) ;
void parameterizationpolarizability(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L) ;

#endif
