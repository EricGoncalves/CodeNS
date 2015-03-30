cmake_minimum_required(VERSION 2.8)
project(CodeNS Fortran)

site_name(hostname)
exec_program( ${CMAKE_Fortran_COMPILER} ARGS --version OUTPUT_VARIABLE _COMPILER_output)
string(REGEX REPLACE ".* ([0-9]*\\.[0-9]*\\.[0-9]*).*" "\\1" COMPILER_VERSION ${_COMPILER_output})



### Definition des options de compilation ###

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release)
endif()
option(openmp  "Enable OpenMP" OFF)
option(profile "Enable instrumentation of the code for profiling purpose" OFF)
option(report  "Report optimizations instead of compilation" OFF)
option(static  "Enable static linking" ON)
option(lto     "Enable link-time optimizations" ON)



### Definition des flags de compilation suivant le compilateur ###

# ifort compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "-fPIC -g  -traceback -u -override-limits -static")
  set(CMAKE_Fortran_FLAGS_RELEASE " -xHOST -O3 -no-prec-div -fp-model fast=2 -fno-alias -opt-subscript-in-range ") #  -init=zero
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -C -debug all -warn all,nounused -fpe0 -ftrapuv -fstack-protector-all  -init=snan -WB -fp-stack-check -u -gen-interfaces ") # -stand f08 -warn nodeclarations
  set(DOUBLE_FLAGS "-r8")
  set(OPTRPT_FLAGS "-qopt-report=5 -guide -diag-enable=openmp,vec,par")
  set(LTO_FLAGS "ipo")
  set(OPENMP_FLAGS "-openmp -parallel -threads")
endif ()

# gfortran compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-fPIC -ffree-line-length-none -g -fbacktrace  -fimplicit-none ")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funsafe-math-optimizations -finit-local-zero -ftree-loop-distribution -ftree-loop-im") #
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g3 -Wall -ffpe-trap=zero,overflow,invalid -Wextra  -Wimplicit-interface -Wno-unused-parameter -Wno-unused ") #  -Wno-unused-dummy-argument -fcheck=all -pedantic -std=f2008 -finit-real=snan
  set(DOUBLE_FLAGS "-fdefault-real-8 -fdefault-double-8")
  set(OPENMP_FLAGS "-fopenmp")
  set(OPTRPT_FLAGS "-ftree-vectorizer-verbose")
  set(LTO_FLAGS "-fwhole-program -flto") # -combine
endif ()

set(PROFILE_FLAGS "-p -pg")



### Definitions de quelques contraintes ###

if(openmp AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND COMPILER_VERSION VERSION_LESS "13.0.0")
    message(WARNING "Insufficient IFort version for OpenMP\n Deactivating it for you")
    set(openmp OFF)
    set(OPENMP_FLAGS "")
endif()
if(openmp AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND COMPILER_VERSION VERSION_LESS "4.9.0")
    message(WARNING "Insufficient GFortran version for OpenMP\n Deactivating it for you")
    set(openmp OFF)
    set(OPENMP_FLAGS "")
endif()
if(lto AND ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(WARNING "Deactivating lto for debug build")
    set(lto OFF)
    set(LTO_FLAGS "")
endif()
if(report AND ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(WARNING "Deactivating optimization report for debug build")
    set(report OFF)
    set(OPTRPT_FLAGS "")
endif()



### On rajoutes les flags des options  ###

message( STATUS "")
message( STATUS "=======================================================================================")
message( STATUS "")
message( STATUS "Compilation de CodeNS sur : " ${hostname})
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} " - " ${COMPILER_VERSION})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DOUBLE_FLAGS}")

if( openmp )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENMP_FLAGS}")
  message( STATUS "With  OpenMP")
else()
  message( STATUS "Without OpenMP")
endif()
if( profile )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PROFILE_FLAGS}")
  message( STATUS "With  Profiling informations")
endif()
if( report )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPTRPT_FLAGS}")
  message( STATUS "With  Optimization reports")
endif()
message( STATUS "")
if (${CMAKE_BUILD_TYPE} MATCHES "Release")
  message( STATUS "Optimized build")
  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  message( STATUS "Debugging build")
  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()
message( STATUS "")



### Définitions des dépendances  ###

set(BUILD_SHARED_LIBS OFF)
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(LINK_SEARCH_START_STATIC ON)
set(LINK_SEARCH_END_STATIC ON)
set(POSITION_INDEPENDENT_CODE ON)

FILE(GLOB SRC_FILES "*.f90")
#find_library(M_LIB m)
#find_library(OMP_LIB iomp5)
#find_library(PTHREAD_LIB pthread)



add_executable(solver_air ${SRC_FILES})
#target_link_libraries(solver_air ${M_LIB})
#target_link_libraries(solver_air ${OMP_LIB})
#target_link_libraries(solver_air ${PTHREAD_LIB})

### Autre (LTO et STATIC)  ###

if( static )
  set_target_properties(solver_air PROPERTIES LINK_SEARCH_START_STATIC ON)
  set_target_properties(solver_air PROPERTIES LINK_SEARCH_END_STATIC ON)
  set_property(TARGET solver_air PROPERTY LINK_SEARCH_START_STATIC ON)
  set_property(TARGET solver_air PROPERTY LINK_SEARCH_END_STATIC ON)
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
  set(CMAKE_SHARED_LINK_DYNAMIC_Fortran_FLAGS "")
  set(CMAKE_EXE_EXPORTS_Fortran_FLAG "")
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "")
  message( STATUS "Static build")
else()
  message( STATUS "Dynamic build")
endif()
if( lto AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
    set_target_properties(solver_air PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    set_property(TARGET solver_air PROPERTY LINK_INTERPROCEDURAL_OPTIMIZATION ON)
endif()

message( STATUS "")
message( STATUS "=======================================================================================")
message( STATUS "")
