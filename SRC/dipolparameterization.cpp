
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
#include "../include/dipolparameterization.h"
#include "../include/dipol.h"

using namespace std;

#define M 9
  
// Function to get cofactor of A[p][q] in temp[][]. n is current 
// dimension of A[][] 
void getCofactor(double A[M][M], double temp[M][M], int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 
  

double determinant(double A[M][M], int n) 
{ 
    double D = 0; // Initialize result 
  
    //  Base case : if matrix contains single element 
    if (n == 1) 
        return A[0][0]; 
  
    double temp[M][M]; // To store cofactors 
  
    int sign = 1;  // To store sign multiplier 
  
     // Iterate for each element of first row 
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n); 
        D += sign * A[0][f] * determinant(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
} 
  
// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(double A[M][M],double adj[M][M]) 
{ 
    if (M == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
  
    // temp is used to store cofactors of A[][] 
    int sign = 1;

    double temp[M][M]; // To store cofactors  
  
    for (int i=0; i<M; i++) 
    { 
        for (int j=0; j<M; j++) 
        { 
            // Get cofactor of A[i][j] 
            getCofactor(A, temp, i, j, M); 
  
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j][i] = (sign)*(determinant(temp, M-1)); 
        } 
    } 
} 

template<class T> 
void display(T A[M][M]) 
{ 
    for (int i=0; i<M; i++) 
    { 
        for (int j=0; j<M; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
} 

 
// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(double A[M][M], double inverse[M][M]) 
{ 
    // Find determinant of A[][] 
    double det = determinant(A, M); 
    if (det == 0) 
    { 
        cout << "Singular matrix, can't find its inverse  " << endl; 
        return false; 
    } 

    // Find adjoint 
    double adj[M][M]; 
    adjoint(A, adj); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<M; i++) 
        for (int j=0; j<M; j++) 
            inverse[i][j] = adj[i][j]/double(det); 
  
    return true; 
} 


void parameterizationpolarizability(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L)
{
  for(uint t = 0; t < nsteps; ++t )
    { 
      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i;
          if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
            {      
              const double am_H = 1.00784 * amu, am_O = 15.999 * amu, am_H2O = (2 * am_H +  am_O);
              double wb_x = 0,wb_y = 0,wb_z = 0;           // bisector vector H2O
              double wb1_x = 0,wb1_y = 0,wb1_z = 0;        // bisector vector OH1
              double wb2_x = 0,wb2_y = 0,wb2_z = 0;        // bisector vector OH2
              double HH_x = 0,HH_y = 0,HH_z = 0 ;          // H-H vector
              double Pv_x = 0,Pv_y = 0,Pv_z = 0;           // vector Perpendicular to bisector and H-H vector

              // water bisector vector, Minimum image convention was applied
              wb1_x = min_distance(r[id+1].x - r[id].x, L[0]) ;
              wb1_y = min_distance(r[id+1].y - r[id].y, L[1]) ;
              wb1_z = min_distance(r[id+1].z - r[id].z, L[2]) ;
              convertounivector(wb1_x, wb1_y, wb1_z);
              wb2_x = min_distance(r[id+2].x - r[id].x, L[0]) ;
              wb2_y = min_distance(r[id+2].y - r[id].y, L[1]) ;
              wb2_z = min_distance(r[id+2].z - r[id].z, L[2]) ;
              convertounivector(wb2_x, wb2_y, wb2_z);

              //water bisector [unit vector]
              wb_x = wb1_x + wb2_x ;
              wb_y = wb1_y + wb2_y ;
              wb_z = wb1_z + wb2_z ;
              convertounivector(wb_x, wb_y, wb_z);

              // H-H vector
              HH_x = (r[id+1].x - r[id+2].x) ;
              HH_y = (r[id+1].y - r[id+2].y) ;
              HH_z = (r[id+1].z - r[id+2].z) ;
              convertounivector(HH_x, HH_y, HH_z);

              // Vector perpendicular to the water bisector and H-H vector
              Pv_x = wb_y * HH_z - wb_z * HH_y; 
              Pv_y = wb_z * HH_x - wb_x * HH_z; 
              Pv_z = wb_x * HH_y - wb_y * HH_x; 
              convertounivector(Pv_x, Pv_y, Pv_z);

              // The H-H vector is slightly adjusted to form perfect orthogonality with other vectors
              HH_x = wb_y * Pv_z - wb_z * Pv_y; 
              HH_y = wb_z * Pv_x - wb_x * Pv_z; 
              HH_z = wb_x * Pv_y - wb_y * Pv_x; 
              convertounivector(HH_x, HH_y, HH_z);

          double A[M][M] = {{HH_x * HH_x, wb_x * wb_x, Pv_x * Pv_x,
                             HH_x * wb_x, HH_x * Pv_x, wb_x * Pv_x,
                             wb_x * HH_x, Pv_x * HH_x, Pv_x * wb_x },

                            {HH_y * HH_y, wb_y * wb_y, Pv_y * Pv_y,
                             HH_y * wb_y, HH_y * Pv_y, wb_y * Pv_y,
                             wb_y * HH_y, Pv_y * HH_y, Pv_y * wb_y },

                            {HH_z * HH_z, wb_z * wb_z, Pv_z * Pv_z,
                             HH_z * wb_z, HH_z * Pv_z, wb_z * Pv_z,
                             wb_z * HH_z, Pv_z * HH_z, Pv_z * wb_z },

                            {HH_x * HH_y, wb_x * wb_y, Pv_x * Pv_y,
                             HH_x * wb_y, HH_x * Pv_y, wb_x * Pv_y,
                             wb_x * HH_y, Pv_x * HH_y, Pv_x * wb_y },

                            {HH_x * HH_z, wb_x * wb_z, Pv_x * Pv_z,
                             HH_x * wb_z, HH_x * Pv_z, wb_x * Pv_z,
                             wb_x * HH_z, Pv_x * HH_z, Pv_x * wb_z },

                            {HH_y * HH_z, wb_y * wb_z, Pv_y * Pv_z,
                             HH_y * wb_z, HH_y * Pv_z, wb_y * Pv_z,
                             wb_y * HH_z, Pv_y * HH_z, Pv_y * wb_z },

                            {HH_y * HH_x, wb_y * wb_x, Pv_y * Pv_x,
                             HH_y * wb_x, HH_y * Pv_x, wb_y * Pv_x,
                             wb_y * HH_x, Pv_y * HH_x, Pv_y * wb_x },

                            {HH_z * HH_x, wb_z * wb_x, Pv_z * Pv_x,
                             HH_z * wb_x, HH_z * Pv_x, wb_z * Pv_x,
                             wb_z * HH_x, Pv_z * HH_x, Pv_z * wb_x },

                            {HH_z * HH_y, wb_z * wb_y, Pv_z * Pv_y,
                             HH_z * wb_y, HH_z * Pv_y, wb_z * Pv_y,
                             wb_z * HH_y, Pv_z * HH_y, Pv_z * wb_y }
                             };

          double inv_A[M][M] = {{0, 0, 0, 0, 0, 0, 0, 0, 0}, 
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                {0, 0, 0, 0, 0, 0, 0, 0, 0},};

          inverse(A, inv_A); 
          float L_w[9]; // cp2k polarisability to parameterize the constants 
          string temp ; 

          ifstream cp2kfile("../output/polar" + to_string(t+1) );
          for(unsigned int i = 0; i <  20; ++i)
            {
              cp2kfile >> temp;
            }
          cp2kfile >> L_w[0] >> L_w[1]  >> L_w[2]  >> temp; 
          cp2kfile >> L_w[3] >> L_w[4]  >> L_w[5]  >> temp; 
          cp2kfile >> L_w[6] >> L_w[7]  >> L_w[8] ; 
          cp2kfile.close();

              float A_xx = 0, A_yy = 0, A_zz = 0,
                    A_xy = 0, A_xz = 0, A_yz = 0,
                    A_yx = 0, A_zx = 0, A_zy = 0;     

          //cout <<  L_w[0] << " " << L_w[1] << " " << L_w[2] << "\n" << L_w[3] << " " << L_w[4] << "  " << L_w[5] << "\n" << L_w[6] << " " << L_w[7] <<  " " << L_w[8] <<  endl;      
          for(unsigned int i = 0; i < M ; ++i)
            {
              A_xx += inv_A[0][i] * L_w[i] ;
              A_yy += inv_A[1][i] * L_w[i] ;
              A_zz += inv_A[2][i] * L_w[i] ; 
              A_xy += inv_A[3][i] * L_w[i] ; 
              A_xz += inv_A[4][i] * L_w[i] ; 
              A_yz += inv_A[5][i] * L_w[i] ;  
              A_yx += inv_A[6][i] * L_w[i] ;   
              A_zx += inv_A[7][i] * L_w[i] ; 
              A_zy += inv_A[8][i] * L_w[i] ; 
            } 

          ofstream outfile("polarizability" + to_string(t+1));       
          outfile <<  A_xx << "\n" << A_yy << "\n" << A_zz << "  \n" << A_xy << "  \n" << A_xz << "  \n" << A_yz << "  \n" << A_yx << "  \n" << A_zx <<  "  \n" << A_zy <<  endl; 
          outfile.close();  
            }
        }
    } 
   
}



void printsinglewatercoordinates(vector<Atom> &r, uint nsteps,  uint natoms, const vector<float> & L)
{

  float rij1 = 0, rij2 = 0, rij = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, x = 0, y = 0, z = 0 ;
  for(uint t = 0; t < nsteps; t = t + 1000 )
    { 
      //pure water
      /*
      for(uint i = 0;i < 128;++i)
        {
          uint id = natoms*t+i;
          uint id1 = natoms*t+i*2+128 ;
          uint id2 = natoms*t+i*2+129 ;
          if(r[id].symbol[0] == 'O' && r[id1].symbol[0] == 'H' && r[id2].symbol[0] == 'H')
            { 
              cout << "3" << endl;
              cout << endl;
              cout << "O " << r[id].x     << "  " <<  r[id].y   << "  " << r[id].z   << endl;
              cout << "H " << r[id1].x   << "  " <<  r[id1].y << "  " << r[id1].z << endl;
              cout << "H " << r[id2].x   << "  " <<  r[id2].y << "  " << r[id2].z << endl;
            }
        } */

      for(uint i = 0;i < natoms;++i)
        {
          uint id = natoms*t+i;  
          uint hbond_cation = 0, hbond_anion = 0;    
          if(r[id].symbol[0] == 'O' && r[id+1].symbol[0] == 'H' && r[id+2].symbol[0] == 'H')
            { 
          for(uint j = 0;j < natoms;++j)
            { 
              uint idi = natoms*t+i;  
              uint idi1 = natoms*t+i+1;          
              uint idi2 = natoms*t+i+2;                  
              uint idj = natoms*t+j;
 
              /*first solvation shell*/
              if(r[j].symbol[0] == 'M' || r[j].symbol[0] == 'N')
                {
                  x = min_distance(r[idj].x - r[idi].x, L[0]);
                  y = min_distance(r[idj].y - r[idi].y, L[1]);
                  z = min_distance(r[idj].z - r[idi].z, L[2]); 
                  rij = mindis(x,y,z,L); 
                  if(rij < 3.2 && rij > 0)
                    {
                      hbond_cation += 1;
                    }      
                }
              if(r[j].symbol[0] == 'C' || r[j].symbol[0] == 'F')
                {
                  x = min_distance(r[idj].x - r[idi1].x, L[0]);
                  y = min_distance(r[idj].y - r[idi1].y, L[1]);
                  z = min_distance(r[idj].z - r[idi1].z, L[2]); 
                  rij1 = mindis(x,y,z,L);

                  x = min_distance(r[idj].x - r[idi2].x, L[0]);
                  y = min_distance(r[idj].y - r[idi2].y, L[1]);
                  z = min_distance(r[idj].z - r[idi2].z, L[2]); 
                  rij2 = mindis(x,y,z,L); 
                  if( (rij1 < 3.0 && rij1 > 0) || (rij2 < 3.0 && rij2 > 0))
                    {
                      hbond_anion += 1;
                    }      
                }
            }

          if(hbond_cation == 0 && hbond_anion == 0)
            {
              cout << "3" << endl;
              cout << endl;
              cout << "O " << r[id].x     << "  " <<  r[id].y   << "  " << r[id].z   << endl;
              cout << "H " << r[id+1].x   << "  " <<  r[id+1].y << "  " << r[id+1].z << endl;
              cout << "H " << r[id+2].x   << "  " <<  r[id+2].y << "  " << r[id+2].z << endl;
            }                  
            }
        }

    }
}


















