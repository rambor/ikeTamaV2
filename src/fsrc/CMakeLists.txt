cmake_minimum_required(VERSION 3.8)
project(fSAXS)
set(VERSION 1.0)

ENABLE_LANGUAGE(Fortran)
set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS}")
MESSAGE("${CMAKE_Fortran_FLAGS}")
#
# use Intel mkl_vml.fi library
# need to specify its location
#
#set(MKL_DIR "/opt/intel/compilers_and_libraries_2019.4.233/mac/mkl/include")

find_path(MKL_LIB mkl_vml.fi /opt/intel/compilers_and_libraries_2019/mac/mkl/include)
link_directories(/Users/xos81802/Downloads/lapack-3.8.0)
#find_path(MKL_MOD lapack95.mod /opt/intel/parallel_studio_xe_2019/compilers_and_libraries_2019/mac/mkl/include/intel64/lp64)
MESSAGE("** SETTING  MKL_LIB=${MKL_LIB}")
include_directories(fSAXS PRIVATE ${MKL_LIB} ${MKL_MOD})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(CMAKE_Fortran_COMPILER MATCHES "gfortran*")
    #    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -std=f2003 -fimplicit-none")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -O0 -g3 -fbounds-check -I/opt/intel/mkl/include")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops")
endif()

set(SOURCE_FILES PlmBar.f90 sphbessel_alt.f90)

add_library(fSAXS STATIC ${SOURCE_FILES})
target_link_libraries (fSAXS lapack)