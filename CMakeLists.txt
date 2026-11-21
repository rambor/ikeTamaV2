cmake_minimum_required(VERSION 3.8)
project(ikeTama)

set(CMAKE_CXX_FLAGS_Release "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG -O3 ")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -mtune=native -pthread -ffast-math -fopenmp -I${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty")
MESSAGE("** SYSTEM : ${CMAKE_SYSTEM_NAME}")

set(CMAKE_VERBOSE_MAKEFILE ON )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

# SASTOOLS is installed in home directory
SET(SASTOOLS_INCLUDE_DIRS $ENV{HOME}/usr/local/include)
SET(SASTOOLS_LIB_DIR $ENV{HOME}/usr/local/lib)
# set location of QHULL libs (same, installed in home directory but build directory)
set(QHULL_LOCALE "$ENV{HOME}/libs/qhull-master")
SET(QHULL_INCLUDE_DIRS ${QHULL_LOCALE}/src ${QHULL_LOCALE}/src/libqhull ${QHULL_LOCALE}/src/libqhullcpp)

#
# BOOST SECTION
#
# set Boost directory
#set(BOOST_REQUIRED_COMPONENTS
#        filesystem
#        random
#        program_options
#        system
#        regex
#)

#
#if(POLICY CMP0167)
#        cmake_policy(SET CMP0167 OLD)
#endif()

#set(Boost_DIR ${BOOSTROOT})
#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_STATIC_RUNTIME ON)
#set(Boost_NO_SYSTEM_PATHS ON)
#set(Boost_NO_BOOST_CMAKE ON)
# Boost root variable set as environment variable in CLION Preferences

#
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${BOOSTROOT})
#MESSAGE("** BOOST INCLUDE => ${Boost_INCLUDE_DIR}")
#
#find_package(Boost CONFIG REQUIRED COMPONENTS ${BOOST_REQUIRED_COMPONENTS})
#
#IF (Boost_FOUND)
#        MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
#        INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
#        ADD_DEFINITIONS( "-DHAS_BOOST" )
#ENDIF()

# BOOST AND CMAKE ARE NUTS - lots of problems, best to directly hardcord, so specify in environment variable as follows
# I keep a built boost library in my home directory
# -DBOOSTROOT=/Users/xos81802/libs/boost
# -DBOOST_ROOT=/Users/xos81802/libs/boost
# -DBOOST_LIBRARYDIR=/Users/xos81802/libs/boost/lib

#MESSAGE("** BOOST INCLUDE DIRECTORY: ${BOOST_LIBRARIES}")
#MESSAGE("** BOOST LIB DIRECTORY: ${Boost_LIBRARY_DIRS}")
#include_directories(${Boost_INCLUDE_DIR})

set(BOOST_FILESYSTEM_LIB ${BOOSTROOT}/lib/libboost_filesystem.a )
set(BOOST_PO_LIB ${BOOSTROOT}/lib/libboost_program_options.a )
set(BOOST_RANDOM_LIB ${BOOSTROOT}/lib/libboost_random.a )
set(BOOST_REGEX_LIB ${BOOSTROOT}/lib/libboost_regex.a )
set(BOOST_SYSTEM_LIB ${BOOSTROOT}/lib/libboost_system.a )
#
#MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
#MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
#MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
#MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
#MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
#MESSAGE("** BOOST INCLUDE DIRECTORY: ${BOOST_INCLUDEDIR}")
#MESSAGE("** CURRENT SOURCE DIRECTORY MAIN : ${CMAKE_CURRENT_SOURCE_DIR}")


include_directories(${SASTOOLS_INCLUDE_DIRS})
include_directories(${BOOST_INCLUDEDIR})
include_directories(${QHULL_INCLUDE_DIRS})

# download and compile QHULL as static lib http://www.qhull.org
# use the same compiler for all programs
# change paths below to point to static libs

find_library(QHULL_LIB libqhullstatic.a PATHS ${QHULL_LOCALE}/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ${QHULL_LOCALE}/lib)

# download and compile SASTools as a static library
find_library(SASTools_LIB libSASTools.a PATHS ${SASTOOLS_LIB_DIR})

set(SOURCE_FILES main.cpp)

include(src/CMakeLists.txt)
include(src/utils/CMakeLists.txt)
include(src/EulerTour/CMakeLists.txt)
include(src/Annealer/CMakeLists.txt)

MESSAGE("** SOURCE FILES => ${SOURCE_FILES}")

add_executable(ikeTama ${SOURCE_FILES})

target_link_libraries(ikeTama
        PUBLIC
        ${BOOST_FILESYSTEM_LIB} ${BOOST_PO_LIB} ${BOOST_RANDOM_LIB} ${BOOST_REGEX_LIB} ${BOOST_SYSTEM_LIB}
        ${QHULL_LIB} ${QHULLCPP_LIB} ${SASTools_LIB}
        )

#target_link_libraries(ikeTama
#        PUBLIC
#        Boost::filesystem Boost::random Boost::program_options Boost::regex Boost::system
#        ${QHULL_LIB} ${QHULLCPP_LIB} ${SASTools_LIB}
#)


# on install, create distribution directory and bin and lib subdirectores
# copy dylib to lib and link binaries to it using rpath
enable_testing()
add_subdirectory("${PROJECT_SOURCE_DIR}/lib")