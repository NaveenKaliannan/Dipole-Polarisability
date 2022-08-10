# Minimum version for CMake
cmake_minimum_required(VERSION 2.8.7)

# Project name
project (exe)

##find_package(MPI REQUIRED)
add_compile_options(-std=c++11)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

set(CMAKE_CXX_FLAGS "-lstdc++")
set(CMAKE_CXX_FLAGS "-static-libstdc++")
set(CMAKE_CXX_FLAGS "-static-libgcc")



#Brings the headers
include_directories(../include)
##include_directories(${MPI_INCLUDE_PATH})

# Source files are added
file(GLOB SOURCES "../SRC/*.cpp") 
add_executable(exe ${SOURCES})


##target_link_libraries(exe ${MPI_LIBRARIES})
target_link_libraries(exe )

