
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>
#include <vector> 
#include <numeric>
#include <sstream>     

using namespace std;


int main ( int argc, char** argv )
{
  float t = atof(argv[1]);

  ifstream dipole_filename(argv[2]);  
  string temp;
  double mu_x = 0, mu_y = 0, mu_z = 0;
  for(unsigned int i = 0; i <  57; ++i) ///57 24
    {
    dipole_filename >> temp; 
    }
  dipole_filename >> mu_x >> temp >> mu_y >> temp >> mu_z ;
  cout << ((t) - 1.0) / 1000.0 << "  " <<  mu_x << "  " << mu_y << "   " << mu_z << "  " ;

  ifstream alpha_filename(argv[3]);  
  double alpha_xx = 0, alpha_yy = 0, alpha_zz = 0;
  double alpha_xy = 0, alpha_xz = 0, alpha_yz = 0;
  double alpha_yx = 0, alpha_zx = 0, alpha_zy = 0;
  for(unsigned int i = 0; i <  20; ++i)
    {
    alpha_filename >> temp;
    }
  alpha_filename >> alpha_xx >> alpha_yy >> alpha_zz >> temp ;
  alpha_filename >> alpha_xy >> alpha_xz >> alpha_yz >> temp ;
  alpha_filename >> alpha_yx >> alpha_zx >> alpha_zy >> temp ;
  cout << alpha_xx << "  " << alpha_yy << "   " << alpha_zz << "  " ;
  cout << alpha_xy << "  " << alpha_xz << "   " << alpha_yz << "  " ;
  cout << alpha_yx << "  " << alpha_zx << "   " << alpha_zy << "  " ;

  ifstream charges_filename(argv[4]);  
  int N = 3;
  double q[N]; 
  for(unsigned int i = 0; i <  12; ++i)
    {
     charges_filename >> temp; 
    }
  for(unsigned int i = 0; i <  N; ++i)
    {
     charges_filename >> temp >> temp >> temp >> temp >> q[i];
    }

  for(unsigned int i = 0; i <  N; ++i)
    {
     //cout << q[i] << "  " ;
    }

  cout << t << endl;
  
  return 0;
}







