
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
  cout << (t * 1.0 ) - 1.0 << "  " <<  mu_x << "  " << mu_y << "   " << mu_z << "  " ;

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


/*
  for(uint t = 0; t < nsteps;t += 1 )
    {  
      for(uint i = 0;i < nmol;++i)
        {
          idi = nmol*t+i; 
          init_vector_zero(dummyv);  
          init_matrix_zero(dummyM);  
          Pol_Efield(mol[idi].PPol, E[t], dummyv);
          for(uint j = 0;j < nmol;++j)
            {
              idj = nmol*t+j; 
              dist(mol, idi, idj, L, x, y, z );
              rij = mindis(x,y,z,L); 
              if(i == j){}
              else if (i != j && rij < rcut )
                {
                  float rhat_x = x, rhat_y = y, rhat_z = z;          
                  convertounivector(rhat_x, rhat_y, rhat_z);
                  Eion.x = mol[idj].q * rhat_x * pow(rij, -2.0) ;
                  Eion.y = mol[idj].q * rhat_y * pow(rij, -2.0) ;
                  Eion.z = mol[idj].q * rhat_z * pow(rij, -2.0) ; 
                  Pol_Ifield(mol[idi].PPol, Eion, dummyv); 

                  iontensorfield(Tij_ion, mol[idj].q, rij, x, y, z); 
                     Tij_dipole(Tij_ion, TD[idj], dummyv);
                       Tij_Pol(Tij_ion, TP[idj], dummyM);    
                }
            }  
          mol[idi].ID.x += dummyv.x ;
          mol[idi].ID.y += dummyv.y ; 
          mol[idi].ID.z += dummyv.z ;

          mol[idi].IPol.xx += dummyM.xx;mol[idi].IPol.yy += dummyM.yy;mol[idi].IPol.zz += dummyM.zz;
          mol[idi].IPol.xy += dummyM.xy;mol[idi].IPol.xz += dummyM.xz;mol[idi].IPol.yz += dummyM.yz;
          mol[idi].IPol.yx += dummyM.yx;mol[idi].IPol.zx += dummyM.zx;mol[idi].IPol.zy += dummyM.zy;
          //TD[idi].x += dummyv.x ; TD[idi].y += dummyv.y ; TD[idi].z += dummyv.z ;
        }

      for(uint iter = 0; iter < niter ; ++iter)
        {       
          init_Vector_zero(dummyTD, 1, nmol); 
          init_Matrix_zero(dummyTP, 1, nmol);
          init_vector_zero(dummyv);  
          init_matrix_zero(dummyM);  
          for(uint i = 0;i < nmol;++i)
            {
              idi = nmol*t+i; 
              for(uint j = 0;j < nmol;++j)
                {
                  idj = nmol*t+j; 
                  dist(mol, idi, idj, L, x, y, z );
                  rij = mindis(x,y,z,L); 
                  if(i == j){}
                  else if (i != j && rij < rcut  )
                  {
                    //Dipole induced interactions
                    dipoletensorfield(Tij, rij, x, y, z);
                    Tij_dipole(Tij, TD[idj], dummyTD[i]);
                       Tij_Pol(Tij, TP[idj], dummyTP[i]); 
                  }

                }
            }
          for(uint i = 0;i < nmol;++i)
            {
              uint idi = nmol*t+i; 

              Mat_vec(mol[idi].PPol, dummyTD[i], TD[idi]);
              Mat_Mat(mol[idi].PPol, dummyTP[i], TP[idi]);

              mol[idi].ID.x += TD[idi].x ;
              mol[idi].ID.y += TD[idi].y ;
              mol[idi].ID.z += TD[idi].z ;

              mol[idi].IPol.xx += TP[idi].xx;mol[idi].IPol.yy += TP[idi].yy;mol[idi].IPol.zz += TP[idi].zz;
              mol[idi].IPol.xy += TP[idi].xy;mol[idi].IPol.xz += TP[idi].xz;mol[idi].IPol.yz += TP[idi].yz;
              mol[idi].IPol.yx += TP[idi].yx;mol[idi].IPol.zx += TP[idi].zx;mol[idi].IPol.zy += TP[idi].zy;
	    }
        } 
    }*/







