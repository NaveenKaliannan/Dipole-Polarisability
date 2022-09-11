
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
#include "../include/cp2k.h"
#include "../include/io.h"

using namespace std;




/*Reads mulliken charges from cp2k output data*/
void readmullikencharges(vector<Atom> &r, uint nsteps, uint natoms, string filename)
{
  string temp;
  for(uint t = 0; t < nsteps; ++t )
    { 
      string filename = "/home/naveenk/temp/1-mgcl2/Frame" + to_string(t+1) + "/charges.out.mulliken" ;
      ifstream charges_filename(filename);
      for(uint i = 0; i <  5; ++i)
        {
          getline(charges_filename, temp);
        }
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i; 
          charges_filename >> temp >> temp >> temp >> temp >> r[id].charge; 
          //cout << t <<  "  " << i << "  " << r[id].symbol << "  " <<  r[id].charge << endl;
        }
      charges_filename.close();
      charges_filename.clear();
    }
}


double mindis(double dx,double dy,double dz,double a,double b,double c)
{
  return norm(dx - a * round(dx/a),dy - b * round(dy/b),dz - c * round(dz/c));
}

double dist(const vector<Atom> &r, uint i, uint j, float L)
{
  return mindis(r[i].x - r[j].x,
                r[i].y - r[j].y,
                r[i].z - r[j].z, L,L,L);
}


float angle_btwn_3points(const vector<Atom> &r, uint i, uint j1, uint j2, uint t, float L, uint N)
{
  float x1,x2;
  float y1,y2;
  float z1,z2;


  x1 = r[N*t+(j1+0)].x - r[N*t+(i+0) ].x;
  x2 = r[N*t+(i+0) ].x - r[N*t+(j2+0)].x;
  y1 = r[N*t+(j1+0)].y - r[N*t+(i+0) ].y;
  y2 = r[N*t+(i+0) ].y - r[N*t+(j2+0)].y;
  z1 = r[N*t+(j1+0)].z - r[N*t+(i+0) ].z;
  z2 = r[N*t+(i+0) ].z - r[N*t+(j2+0)].z;
  float a = L, b = L, c = L;

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

uint check_Hbondd(const vector<Atom> &r, uint i, uint ii, uint t, uint Natoms, float L)
{
  float OOdistance = mindis(r[Natoms*t+(i+0)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+0)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+0)].z - r[Natoms*t+(ii+0)].z, L,L,L);

  if(OOdistance < 3.5 && angle_btwn_3points(r,i,i+1,ii, t, L, Natoms) < 30) //Hbond creteria
    {
      return 1;
    }
  else if(OOdistance < 3.5 && angle_btwn_3points(r,i,i+2,ii, t, L,Natoms) < 30) //Hbond creteria
    {
      return 1;
    }
  else
    {
      return 0;
    }
}


uint check_Hbonda(const vector<Atom> &r, uint i, uint ii, uint t, uint Natoms, float L)
{
  float OOdistance = mindis(r[Natoms*t+(ii+0)].x - r[Natoms*t+(i+0)].x,
                             r[Natoms*t+(ii+0)].y - r[Natoms*t+(i+0)].y,
                             r[Natoms*t+(ii+0)].z - r[Natoms*t+(i+0)].z, L,L,L);

  if(OOdistance < 3.5 && angle_btwn_3points(r,ii,ii+1,i, t, L, Natoms) < 30) //Hbond creteria
    {
      return 1;
    }
  else if(OOdistance < 3.5 && angle_btwn_3points(r,ii,ii+2,i, t, L,Natoms) < 30) //Hbond creteria
    {
      return 1;
    }
  else
    {
      return 0;
    }
}


float Hbondangled(const vector<Atom> &r, uint i, uint ii, uint t, uint Natoms, float L)
{
  float OOdistance = mindis(r[Natoms*t+(i+0)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+0)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+0)].z - r[Natoms*t+(ii+0)].z, L,L,L);

  float OH1distance = mindis(r[Natoms*t+(i+1)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+1)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+1)].z - r[Natoms*t+(ii+0)].z, L,L,L);

  float OH2distance = mindis(r[Natoms*t+(i+2)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+2)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+2)].z - r[Natoms*t+(ii+0)].z, L,L,L);

  if(OH1distance < OH2distance) //Hbond creteria
    {
      return angle_btwn_3points(r,i,i+1,ii, t, L, Natoms);
    }
  else
    {
      return angle_btwn_3points(r,i,i+2,ii, t, L,Natoms);
    }
}


float Hbondanglea(const vector<Atom> &r, uint i, uint ii, uint t, uint Natoms, float L)
{
  float OOdistance = mindis(r[Natoms*t+(ii+0)].x - r[Natoms*t+(i+0)].x,
                             r[Natoms*t+(ii+0)].y - r[Natoms*t+(i+0)].y,
                             r[Natoms*t+(ii+0)].z - r[Natoms*t+(i+0)].z, L,L,L);

    if(OOdistance < 3.5 && angle_btwn_3points(r,ii,ii+1,i, t, L, Natoms) < 30) //Hbond creteria
    {
      return angle_btwn_3points(r,ii,ii+1,i, t, L, Natoms);
    }
    else if(OOdistance < 3.5 && angle_btwn_3points(r,ii,ii+2,i, t, L,Natoms) < 30) //Hbond creteria
    {
      return angle_btwn_3points(r,ii,ii+2,i, t, L,Natoms);
    }
  else
    {
      return 0;
    }
}

float OOdist_donor(const vector<Atom> &r, uint i, uint ii, uint t, uint Natoms, float L)
{
  float OOdistance = mindis(r[Natoms*t+(i+0)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+0)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+0)].z - r[Natoms*t+(ii+0)].z, L,L,L);

      return OOdistance;
}



float OHdist_donor(const vector<Atom> &r, uint i, uint ii, uint t, uint Natoms, float L)
{
  float OH1distance = mindis(r[Natoms*t+(i+1)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+1)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+1)].z - r[Natoms*t+(ii+0)].z, L,L,L);

  float OH2distance = mindis(r[Natoms*t+(i+2)].x - r[Natoms*t+(ii+0)].x,
                             r[Natoms*t+(i+2)].y - r[Natoms*t+(ii+0)].y,
                             r[Natoms*t+(i+2)].z - r[Natoms*t+(ii+0)].z, L,L,L);

      return min(OH1distance, OH2distance);
}


void Print_ALMO_data(vector<Atom> &r, uint nsteps, uint natoms, const vector<float> & L, float dt, string donorfile, string acceptorfile, string filename)
{
  uint deltat_ = 5;
  uint id = 0;
  string temp; stringstream oneline;
  for(uint t = 5;t < nsteps; t = t + deltat_)
    {
      string gammaDfilename = donorfile;
      gammaDfilename.append(to_string(t));
      gammaDfilename.append(".dat");
      ifstream gammaDfile(gammaDfilename);
      getline(gammaDfile, temp); 

      string gammaAfilename = acceptorfile;
      gammaAfilename.append(to_string(t));
      gammaAfilename.append(".dat");
      ifstream gammaAfile(gammaAfilename);
      getline(gammaAfile, temp);     
      for(uint i = 0;i < natoms;i = i + 3)
        {
          id = natoms*t+i; 
          gammaDfile >> temp  >>  r[id].D1 >>  r[id].ED1  >>  r[id].CTD1  >>  r[id].D2  >>  r[id].ED2 >>  r[id].CTD2  >> r[id].gamma_d ; 
          gammaAfile >> temp  >>  r[id].A1 >>  r[id].EA1  >>  r[id].CTA1  >>  r[id].A2  >>  r[id].EA2 >>  r[id].CTA2  >> r[id].gamma_a ; 
        }
    }

  ofstream outfile(filename);
  for(uint t = 5;t < nsteps; t = t + deltat_)
    {
      float tot = 0;
      float strong_d = 0;
      float secondstrong_d = 0;     
      uint mean = 0, mean1 = 0;
      for(uint i = 0;i < natoms;i = i + 3)
        {
          id = natoms*t+i; 
          uint idt = natoms*t+i;  

          uint NHbond = 0;
          NHbond += check_Hbondd(r,i, r[idt].D1, t, natoms, L[0]); // D to A
          NHbond += check_Hbondd(r,i, r[idt].D2, t, natoms, L[0]); // D to A
          NHbond += check_Hbonda(r,i, r[idt].A1, t, natoms, L[0]); // A to D
          NHbond += check_Hbonda(r,i, r[idt].A2, t, natoms, L[0]); // A to D
          if(NHbond == 4)
	  {
            mean += 4;
	    tot +=  r[id].ED1 * 2625 + r[id].ED2 * 2625 ;
	    tot +=  r[id].EA1 * 2625 + r[id].EA2 * 2625 ;
	  }
         if(check_Hbondd(r,i, r[idt].D1, t, natoms, L[0]) and check_Hbondd(r,i, r[idt].D2, t, natoms, L[0]))
	  {
            mean1 += 1;
	    strong_d +=  r[id].ED1 * 2625 ;
	    secondstrong_d += r[id].ED2 * 2625 ;
	  }
        }
      outfile << t*dt << "  "<< tot/mean << "  " << strong_d/mean1 << "  " << secondstrong_d/mean1  << endl;
    }
  outfile.close();
  outfile.clear();
}


