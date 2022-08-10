export OMP_NUM_THREADS=1
export PATH=/home/naveenk/my_programs/cmake-3.20.1/bin:$PATH

mkdir build
cp CMakeLists.txt build/CMakeLists.txt
cd build
cmake . 
make
cd ..
