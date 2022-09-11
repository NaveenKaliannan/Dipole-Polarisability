#include<iostream>
#include<cmath>
#include<fstream>
#include<typeinfo>
#include<cstdlib>
#include<string>

using namespace std;

int main(int argc, char* argv[])
{

  unsigned int no_of_files = atoi(argv[1]);
  unsigned int no_of_coloumn = atoi(argv[2]);
  unsigned int no_of_lines = atoi(argv[3]);
  string inputfilename = argv[4];
  string outputfilename = argv[5];

  double vec[no_of_coloumn][no_of_lines];
  for( unsigned int l = 0;l < no_of_lines; ++l )
    {
      for(unsigned int c = 0;c < no_of_coloumn; ++c)
        {
          vec[c][l] = 0;
        }
    } 

  for(unsigned int i = 0;i < no_of_files;++i)
    {
      string filename = inputfilename + to_string(i+1);// + ".dat" ;
      ifstream infile(filename);
      for( unsigned int l = 0;l < no_of_lines; ++l )
        {
          for(unsigned int c = 0;c < no_of_coloumn; ++c)
            {
              double temp;
              infile >> temp; 
              vec[c][l] += temp;// if(l==0 && c==7){cout << temp << endl; }
            } 
        }      
      infile.close();
      infile.clear();
    } 

  ofstream outfile(outputfilename);
  for( unsigned int l = 0;l < no_of_lines; ++l )
    {
      for(unsigned int c = 0;c < no_of_coloumn; ++c)
        {
          outfile << vec[c][l] / no_of_files  << " "; 
        }
          outfile << "\n";
    }  
  outfile.close();
  outfile.clear();


 
  return 0;
}
