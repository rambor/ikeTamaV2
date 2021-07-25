cmake_minimum_required(VERSION 3.8)
project(ikeTama)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG -O3 ")
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -msse4.1 -mavx -I${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty")
MESSAGE("** SYSTEM : ${CMAKE_SYSTEM_NAME}")


set(CMAKE_VERBOSE_MAKEFILE ON )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

SET (QHULL_INCLUDE_DIRS ~/src_libs/qhull-2015.2/src ~/src_libs/qhull-2015.2/src/libqhull ~/src_libs/qhull-2015.2/src/libqhullcpp)
include_directories(${QHULL_INCLUDE_DIRS})
#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)

# if you have boost in a non-standard directory, like you home directory, set the following:
#/Users/robertrambo/boost/boost_1_57_0
#/Users/robertrambo/boost/boost_1_57_0/stage/lib

# Boost root variable set as environment variable in CLION Preferences
MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
SET(BOOST_ROOT ${BOOSTROOT})
SET(BOOST_LIBRARYDIR ${BOOSTROOT}/stage/lib)
SET(BOOST_INCLUDEDIR ${BOOSTROOT})
include_directories(${BOOST_INCLUDEDIR})

set(BOOST_FILESYSTEM ${BOOSTROOT}/stage/lib/libboost_filesystem.a )
set(BOOST_PO ${BOOSTROOT}/stage/lib/libboost_program_options.a )
set(BOOST_RANDOM ${BOOSTROOT}/stage/lib/libboost_random.a )
set(BOOST_REGEX ${BOOSTROOT}/stage/lib/libboost_regex.a )
set(BOOST_SYSTEM ${BOOSTROOT}/stage/lib/libboost_system.a )


MESSAGE("** BOOST INCLUDE => ${BOOST_INCLUDEDIR}")
MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")
MESSAGE("** gFORT => ${CMAKE_Fortran_COMPILER}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${BOOSTROOT}/src/fsrc)
add_subdirectory(src/fsrc)

MESSAGE("** CURRENT SOURCE DIRECTORY MAIN : ${CMAKE_CURRENT_SOURCE_DIR}")

# download and compile QHULL as static lib
# use the same compiler for all programs
# change paths below to point to static libs
find_library(QHULL_LIB libqhullstatic.a PATHS ~/src_libs/qhull-2015.2/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/src_libs/qhull-2015.2/lib)

find_library(SASTools REQUIRED)

find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
link_directories(${OpenCL_LIBRARY})
if (OpenCL_FOUND)
    MESSAGE("** OPENCL FOUND AT => ${OpenCL_INCLUDE_DIR}")
    MESSAGE("** OPENCL FOUND AT => ${OpenCL_INCLUDE_DIRS}")
endif()

set(SOURCE_FILES main.cpp)
#add_subdirectory(src/utils)
include(src/CMakeLists.txt)
include(src/utils/CMakeLists.txt)
include(src/EulerTour/CMakeLists.txt)
include(src/Annealer/CMakeLists.txt)
include(src/SHEMapper/CMakeLists.txt)
MESSAGE("** SOURCE FILES => ${SOURCE_FILES}")

add_executable(ikeTama ${SOURCE_FILES})


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE("** with ${CMAKE_CXX_COMPILER_ID}")
    SET(CMAKE_EXE_LINKER_FLAGS "-lgfortran" )
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("** LINKING WILL BE PERFORMED WITH -lgfortran")
    SET(CMAKE_EXE_LINKER_FLAGS "-L${GFORTRAN_LIB} -lgfortran" )
    # if using icc and gfortran
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

    # if using MS Visual Studio and gfortran
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

endif()

#find_package(OpenMP)
set (OPENMP_FOUND FALSE)
if(OPENMP_FOUND)
    MESSAGE("** openMP ")
    target_link_libraries(ikeTama
            PUBLIC
            ${OpenMP_CXX_FLAGS}
            ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM}
            ${QHULL_LIB} ${QHULLCPP_LIB} SASTools fSAXS
            )

else()
    target_link_libraries(ikeTama
            PUBLIC
            ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM}
            OpenCL::OpenCL ${QHULL_LIB} ${QHULLCPP_LIB} SASTools fSAXS
            )
endif()




# on install, create distribution directory and bin and lib subdirectores
# copy dylib to lib and link binaries to it using rpath
enable_testing()
add_subdirectory("${PROJECT_SOURCE_DIR}/lib")