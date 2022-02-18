cmake_minimum_required(VERSION 3.8)
project(ikeTama)

set(CMAKE_CXX_FLAGS_Release "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG -O3 ")
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -pthread -ffast-math -fopenmp -I${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty")
MESSAGE("** SYSTEM : ${CMAKE_SYSTEM_NAME}")


set(CMAKE_VERBOSE_MAKEFILE ON )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

SET (QHULL_INCLUDE_DIRS ~/src_libs/qhull-master/src ~/src_libs/qhull-master/src/libqhull ~/src_libs/qhull-master/src/libqhullcpp)
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
find_package(Boost)
MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")
MESSAGE("** BOOST LIB DIRECTORY: ${Boost_LIBRARY_DIRS}")
include_directories(${Boost_INCLUDE_DIR})

set(BOOST_FILESYSTEM ${Boost_LIBRARY_DIRS}/libboost_filesystem.a )
set(BOOST_PO ${Boost_LIBRARY_DIRS}/libboost_program_options.a )
set(BOOST_RANDOM ${Boost_LIBRARY_DIRS}/libboost_random.a )
set(BOOST_REGEX ${Boost_LIBRARY_DIRS}/libboost_regex.a )
set(BOOST_SYSTEM ${Boost_LIBRARY_DIRS}/libboost_system.a )

MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")
MESSAGE("** gFORT => ${CMAKE_Fortran_COMPILER}")

MESSAGE("** CURRENT SOURCE DIRECTORY MAIN : ${CMAKE_CURRENT_SOURCE_DIR}")

SET (SASTOOLS_INCLUDE_DIRS $ENV{HOME}/usr/local/include)
SET(USER_LOCAL_DIR $ENV{HOME}/usr/local)
include_directories(${SASTOOLS_INCLUDE_DIRS})

# download and compile QHULL as static lib
# use the same compiler for all programs
# change paths below to point to static libs
find_library(QHULL_LIB libqhullstatic.a PATHS ~/src_libs/qhull-master/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/src_libs/qhull-master/lib)
# download and compile SASTools as a static library
find_library(SASTools_LIB libSASTools.a PATHS ~/usr/local/lib ~/lib)

#set_target_properties(SASTools PROPERTIES
#        IMPORTED_LOCATION "${USER_LOCAL_DIR}/lib/libSASTools.a"
#        INTERFACE_INCLUDE_DIRECTORIES "${USER_LOCAL_DIR}/include"
#        )

set(SOURCE_FILES main.cpp)

include(src/CMakeLists.txt)
include(src/utils/CMakeLists.txt)
include(src/EulerTour/CMakeLists.txt)
include(src/Annealer/CMakeLists.txt)
include(src/SHEMapper/CMakeLists.txt)

MESSAGE("** SOURCE FILES => ${SOURCE_FILES}")

add_executable(ikeTama ${SOURCE_FILES})

target_link_libraries(ikeTama
        PUBLIC
        ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM}
        ${QHULL_LIB} ${QHULLCPP_LIB} ${SASTools_LIB}
        )


# on install, create distribution directory and bin and lib subdirectores
# copy dylib to lib and link binaries to it using rpath
enable_testing()
add_subdirectory("${PROJECT_SOURCE_DIR}/lib")