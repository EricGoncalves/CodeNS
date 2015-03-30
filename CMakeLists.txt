cmake_minimum_required(VERSION 2.8)
project(CodeNS Fortran)

#set(BUILD_SHARED_LIBS OFF)
#set(LINK_SEARCH_START_STATIC ON)
#set(LINK_SEARCH_END_STATIC ON)

# ifort compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} )
set(CMAKE_Fortran_FLAGS "-fPIC -static -g  -traceback -u -override-limits") # COMMON FLAGS
set(CMAKE_Fortran_FLAGS_RELEASE " -xHOST -O3 -no-prec-div -fp-model fast=2 -fno-alias -init=zero -opt-subscript-in-range ")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -C -debug all -warn all,nounused -fpe0 -ftrapuv -fstack-protector-all  -init=snan -WB -fp-stack-check -u -gen-interfaces ") # -stand f08 -warn nodeclarations
set(DOUBLE "-r8")
set(OPENMP "-openmp -parallel -threads")
set(OPTRPT "-qopt-report=5 -guide -diag-enable=openmp,vec,par")
endif ()

# gfortran compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} )
set(CMAKE_Fortran_FLAGS "-fPIC -ffree-line-length-none -g -fbacktrace  -fimplicit-none ")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funsafe-math-optimizations -finit-local-zero") #-ftree-loop-distribution -ftree-loop-im
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g3 -Wall -ffpe-trap=zero,overflow,invalid -Wextra  -Wimplicit-interface -Wno-unused-parameter -Wno-unused ") #  -Wno-unused-dummy-argument -fcheck=all -pedantic -std=f2008 -finit-real=snan
set(DOUBLE "-fdefault-real-8 -fdefault-double-8")
set(OPENMP "-fopenmp")
endif ()

set(PROFILE "-p -pg")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release)
endif()
option(openmp "Enable OpenMP" OFF)
option(profile "Enable instrumentation of the code for profiling purpose" OFF)
option(report "Enable report of optimizations" OFF)

add_definitions(${DOUBLE})

if( openmp )
add_definitions(${OPENMP})
endif()
if( profile )
add_definitions(${PROFILE})
endif()
if( report )
add_definitions(${OPTRPT})
endif()
if (${CMAKE_BUILD_TYPE} MATCHES "Release")
set(CMAKE_EXE_LINK_FLAGS ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
set(CMAKE_EXE_LINK_FLAGS ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

FILE(GLOB SRC_FILES "*.f90")
find_library(M_LIB m)

add_executable(solver_air ${SRC_FILES})
#target_link_libraries(solver_air ${M_LIB})

#set_target_properties(solver_air PROPERTIES LINK_SEARCH_START_STATIC ON)
#set_target_properties(solver_air PROPERTIES LINK_SEARCH_END_STATIC ON)
#set_property(TARGET solver_air PROPERTY LINK_SEARCH_START_STATIC ON)
#set_property(TARGET solver_air PROPERTY LINK_SEARCH_END_STATIC ON)

#if (${CMAKE_BUILD_TYPE} MATCHES "Release")
#set_target_properties(solver_air PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
#set_property(TARGET solver_air PROPERTY LINK_INTERPROCEDURAL_OPTIMIZATIONC ON)
#endif()

