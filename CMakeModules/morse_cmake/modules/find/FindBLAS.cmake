###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2019 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find BLAS library
# This module finds an installed fortran library that implements the BLAS
# linear-algebra interface (see http://www.netlib.org/blas/).
# The list of libraries searched for is taken
# from the autoconf macro file, acx_blas.m4 (distributed at
# http://ac-archive.sourceforge.net/ac-archive/acx_blas.html).
#
# This module sets the following variables:
#  BLAS_FOUND - set to true if a library implementing the BLAS interface
#    is found
#  BLAS_LDFLAGS_OTHER - list of required linker flags (excluding -l
#    and -L).
#  BLAS_CFLAGS_OTHER - list of required compiler flags (excluding -I).
#  BLAS_INCLUDE_DIRS - include directories
#  BLAS_LIBRARIES    - list of libraries (using full path name) to
#    link against to use BLAS, sequential or multithreaded version depending on BLA_VENDOR
#  BLAS_LIBRARIES_SEQ - list of libraries (using full path name) to
#    link against to use BLAS, sequential version (Intel MKL)
#  BLAS_LIBRARIES_PAR - list of libraries (using full path name) to
#    link against to use BLAS, multithreaded version (Intel MKL)
#  BLAS95_LIBRARIES - list of libraries (using full path name)
#    to link against to use BLAS95 interface
#  BLAS95_FOUND - set to true if a library implementing the BLAS f95 interface
#    is found
#  BLA_STATIC  if set on this determines what kind of linkage we do (static)
#  BLA_VENDOR  if set checks only the specified vendor, if not set checks
#     all the possibilities
#  BLAS_VENDOR_FOUND stores the BLAS vendor found
#  BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
#
#  BLA_PREFER_PKGCONFIG - try using pkg-config (ON by default)
#  BLAS_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = BLAS
#  <XPREFIX> = <PREFIX>        for common case
#  <XPREFIX> = <PREFIX>_STATIC for static linking
#  <XPREFIX>_FOUND          ... set to 1 if module(s) exist
#  <XPREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
#  <XPREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
#  <XPREFIX>_LDFLAGS        ... all required linker flags
#  <XPREFIX>_LDFLAGS_OTHER  ... all other linker flags
#  <XPREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
#  <XPREFIX>_CFLAGS         ... all required cflags
#  <XPREFIX>_CFLAGS_OTHER   ... the other compiler flags
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DBLAS_DIR=path/to/blas):
#  BLAS_DIR            - Where to find the base directory of blas
#  BLAS_INCDIR         - Where to find the header files
#  BLAS_LIBDIR         - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: BLAS_DIR, BLAS_INCDIR, BLAS_LIBDIR
# For MKL case and if no paths are given as hints, we will try to use the MKLROOT
# environment variable
# The user can also give directly the BLAS libraries to be used with
# BLAS_LIBRARIES_USER. If BLA_STATIC is ON users may also give some additional
# libraries to resolve dependencies with BLAS_LDFLAGS_USER
#  BLAS_VERBOSE Print some additional information during BLAS libraries detection
##########
### List of vendors (BLA_VENDOR) valid in this module
########## List of vendors (BLA_VENDOR) valid in this module
##  Open (for OpenBlas), Eigen (for EigenBlas), Goto, ATLAS PhiPACK,
##  CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, IBMESSLMT
##  Intel10_32 (intel mkl v10 32 bit), Intel10_64lp (intel mkl v10 64 bit,lp thread model, lp64 model),
##  Intel10_64lp_seq (intel mkl v10 64 bit,sequential code, lp64 model),
##  Intel( older versions of mkl 32 and 64 bit),
##  ARMPL, ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic, User (see BLAS_LIBRARIES_USER)
# C/CXX should be enabled to use Intel mkl
###
# We handle different modes to find the dependency
#
# - Detection if already installed on the system
#   - BLAS libraries can be detected from different ways
#     Here is the order of precedence:
#     1) we look in cmake variable BLAS_LIBDIR or BLAS_DIR (we guess the libdirs) if defined
#     2) we look in environment variable BLAS_LIBDIR or BLAS_DIR (we guess the libdirs) if defined
#     3) we look in common environnment variables depending on the system (INCLUDE, C_INCLUDE_PATH, CPATH - LIB, DYLD_LIBRARY_PATH, LD_LIBRARY_PATH, LIBRARY_PATH)
#     4) we look in common system paths depending on the system, see for example paths contained in the following cmake variables:
#       - CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES, CMAKE_C_IMPLICIT_LINK_DIRECTORIES
#

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# Common macros to use in finds
include(FindInit)

if (NOT BLAS_FOUND)
  set(BLAS_DIR "" CACHE PATH "Installation directory of BLAS library")
  if (NOT BLAS_FIND_QUIETLY)
    message(STATUS "A cache variable, namely BLAS_DIR, has been set to specify the install directory of BLAS")
  endif()
endif()

option(BLAS_VERBOSE "Print some additional information during BLAS libraries detection" OFF)
mark_as_advanced(BLAS_VERBOSE)
if(NOT DEFINED BLA_PREFER_PKGCONFIG)
  set(BLA_PREFER_PKGCONFIG ON CACHE BOOL "Try to find BLAS using pkg-config")
endif()
set(BLAS_LIBRARIES_USER "" CACHE STRING "BLAS libraries given by users")
set(BLAS_LIBRARIES_USER "" CACHE STRING "BLAS libraries given by users")
set(BLAS_LDFLAGS_USER   "" CACHE STRING "BLAS link flags given by users")
if(BLAS_LIBRARIES_USER)
  set(BLA_VENDOR "User")
endif()

set(_blas_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Check the language being used
get_property( _LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES )
if( _LANGUAGES_ MATCHES Fortran )
  set( _CHECK_FORTRAN TRUE )
  include(CheckFortranFunctionExists)
elseif( (_LANGUAGES_ MATCHES C) OR (_LANGUAGES_ MATCHES CXX) )
  set( _CHECK_FORTRAN FALSE )
  include(CheckFunctionExists)
else()
  if(BLAS_FIND_REQUIRED)
    message(FATAL_ERROR "FindBLAS requires Fortran, C, or C++ to be enabled.")
  else()
    message(STATUS "Looking for BLAS... - NOT found (Unsupported languages)")
    return()
  endif()
endif()

macro(Check_Fortran_Libraries LIBRARIES _prefix _name _flags _list _thread)
  # This macro checks for the existence of the combination of fortran libraries
  # given by _list.  If the combination is found, this macro checks (using the
  # Check_Fortran_Function_Exists macro) whether can link against that library
  # combination using the name of a routine given by _name using the linker
  # flags given by _flags.  If the combination of libraries is found and passes
  # the link test, LIBRARIES is set to the list of complete library paths that
  # have been found.  Otherwise, LIBRARIES is set to FALSE.

  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.

  set(_libdir ${ARGN})

  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)
  set(ENV_MKLROOT "$ENV{MKLROOT}")
  set(ENV_ARMPL "$ENV{ARMPL_DIR}")
  set(ENV_BLAS_DIR "$ENV{BLAS_DIR}")
  set(ENV_BLAS_LIBDIR "$ENV{BLAS_LIBDIR}")
  if (NOT _libdir)
    if (BLAS_LIBDIR)
      list(APPEND _libdir "${BLAS_LIBDIR}")
    elseif (BLAS_DIR)
      list(APPEND _libdir "${BLAS_DIR}")
      list(APPEND _libdir "${BLAS_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	list(APPEND _libdir "${BLAS_DIR}/lib64")
	list(APPEND _libdir "${BLAS_DIR}/lib/intel64")
      else()
	list(APPEND _libdir "${BLAS_DIR}/lib32")
	list(APPEND _libdir "${BLAS_DIR}/lib/ia32")
      endif()
    elseif(ENV_BLAS_LIBDIR)
      list(APPEND _libdir "${ENV_BLAS_LIBDIR}")
    elseif(ENV_BLAS_DIR)
      list(APPEND _libdir "${ENV_BLAS_DIR}")
      list(APPEND _libdir "${ENV_BLAS_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib64")
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib/intel64")
      else()
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib32")
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib/ia32")
      endif()
    else()
      if (ENV_MKLROOT)
	list(APPEND _libdir "${ENV_MKLROOT}/lib")
	if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	  list(APPEND _libdir "${ENV_MKLROOT}/lib64")
	  list(APPEND _libdir "${ENV_MKLROOT}/lib/intel64")
	else()
	  list(APPEND _libdir "${ENV_MKLROOT}/lib32")
	  list(APPEND _libdir "${ENV_MKLROOT}/lib/ia32")
	endif()
      endif()
      if (ENV_ARMPL)
	list(APPEND _libdir "${ENV_ARMPL}/lib")
      endif()
      list(APPEND _libdir "$ENV{LIBRARY_PATH}")
      if (WIN32)
	string(REPLACE ":" ";" _libdir2 "$ENV{LIB}")
      elseif (APPLE)
	string(REPLACE ":" ";" _libdir2 "$ENV{DYLD_LIBRARY_PATH}")
      else ()
	string(REPLACE ":" ";" _libdir2 "$ENV{LD_LIBRARY_PATH}")
      endif ()
      list(APPEND _libdir "${_libdir2}")
      list(APPEND _libdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    endif()
  endif ()

  if (BLAS_VERBOSE)
    message("${Cyan}Try to find BLAS libraries: ${_list}")
  endif ()

  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_work)
      if (BLA_STATIC)
	if (WIN32)
	  set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif ()
	if (APPLE)
	  set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
	else ()
	  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif ()
      else ()
	if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	  # for ubuntu's libblas3gf and liblapack3gf packages
	  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
	endif ()
      endif ()
      if (EXISTS ${_library})
	set(${_prefix}_${_library}_LIBRARY ${_library})
      else()
	find_library(${_prefix}_${_library}_LIBRARY
	  NAMES ${_library}
	  HINTS ${_libdir}
	  )
      endif()
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      # Print status if not found
      # -------------------------
      if (NOT ${_prefix}_${_library}_LIBRARY AND NOT BLAS_FIND_QUIETLY AND BLAS_VERBOSE)
	Print_Find_Library_Blas_Status(blas ${_library} ${_libdir})
      endif ()
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif(_libraries_work)
  endforeach(_library ${_list})

  if(_libraries_work)
    # Test this combination of libraries.
    if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND BLA_STATIC)
      list(INSERT ${LIBRARIES} 0 "-Wl,--start-group")
      list(APPEND ${LIBRARIES} "-Wl,--end-group")
    endif()
    set(CMAKE_REQUIRED_LIBRARIES "${_flags};${${LIBRARIES}};${_thread}")
    set(CMAKE_REQUIRED_INCLUDES "${BLAS_INCLUDE_DIRS}")
    set(CMAKE_REQUIRED_FLAGS "${BLAS_CFLAGS_OTHER}")
    if (BLAS_VERBOSE)
      message("${Cyan}BLAS libs found for BLA_VENDOR ${BLA_VENDOR}."
	"Try to compile symbol ${_name} with following libraries:"
	"${CMAKE_REQUIRED_LIBRARIES}")
    endif ()
    if(NOT BLAS_FOUND)
      unset(${_prefix}${_combined_name}_WORKS CACHE)
    endif()
    if (_CHECK_FORTRAN)
      check_fortran_function_exists("${_name}" ${_prefix}${_combined_name}_WORKS)
    else()
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    endif()
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
    # Print status if not found
    # -------------------------
    if (NOT _libraries_work AND NOT BLAS_FIND_QUIETLY AND BLAS_VERBOSE)
      Print_Find_Library_Blas_CheckFunc_Status(${_name} ${CMAKE_REQUIRED_LIBRARIES})
    endif ()
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()

  if(_libraries_work)
    set(${LIBRARIES} ${${LIBRARIES}} ${_thread})
  else(_libraries_work)
    set(${LIBRARIES} FALSE)
  endif(_libraries_work)

endmacro(Check_Fortran_Libraries)


set(BLAS_LDFLAGS_OTHER)
set(BLAS_LIBRARIES)
set(BLAS95_LIBRARIES)
if ($ENV{BLA_VENDOR} MATCHES ".+")
  set(BLA_VENDOR $ENV{BLA_VENDOR})
else ()
  if(NOT BLA_VENDOR)
    set(BLA_VENDOR "All")
  endif()
endif ()

set(ENV_BLAS_DIR "$ENV{BLAS_DIR}")
set(ENV_MKL_DIR "$ENV{MKLROOT}")
set(ENV_BLAS_INCDIR "$ENV{BLAS_INCDIR}")
set(ENV_BLAS_LIBDIR "$ENV{BLAS_LIBDIR}")
set(BLAS_GIVEN_BY_USER "FALSE")
if ( BLAS_LIBRARIES_USER OR BLAS_DIR OR ( BLAS_INCDIR AND BLAS_LIBDIR) OR ENV_BLAS_DIR OR (ENV_BLAS_INCDIR AND ENV_BLAS_LIBDIR) )
  set(BLAS_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if( PKG_CONFIG_EXECUTABLE AND NOT BLAS_GIVEN_BY_USER AND BLA_PREFER_PKGCONFIG)

  if (BLA_STATIC)
    set(MKL_STR_BLA_STATIC "static")
  else()
    set(MKL_STR_BLA_STATIC "dynamic")
  endif()
  # try different blas
  if ( (BLA_VENDOR MATCHES "Intel*") AND
      (NOT BLA_F95) AND NOT
      ((_LANGUAGES_ MATCHES Fortran) AND (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
      AND NOT (CMAKE_C_COMPILER_ID STREQUAL "GNU"))
  if (BLA_VENDOR STREQUAL "Intel10_64lp")
    pkg_search_module(BLAS mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
    if (BLAS_FOUND)
      set (BLAS_VENDOR_FOUND "Intel MKL Sequential")
      set (BLAS_LIBRARIES_SEQ "${BLAS_LIBRARIES}")
    endif()
    pkg_search_module(BLAS mkl-${MKL_STR_BLA_STATIC}-lp64-iomp)
    if (BLAS_FOUND)
      set (BLAS_VENDOR_FOUND "Intel MKL Multithreaded")
      set (BLAS_LIBRARIES_PAR "${BLAS_LIBRARIES}")
    endif()
  elseif(BLA_VENDOR STREQUAL "Intel10_64lp_seq")
    pkg_search_module(BLAS mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
    if (BLAS_FOUND)
      set (BLAS_VENDOR_FOUND "Intel MKL")
      set (BLAS_LIBRARIES_SEQ "${BLAS_LIBRARIES}")
    endif()
    endif()
  elseif(BLA_VENDOR STREQUAL "Open")
    pkg_search_module(BLAS openblas)
    if (BLAS_FOUND)
      set (BLAS_VENDOR_FOUND "Openblas")
    endif()
  elseif(BLA_VENDOR STREQUAL "Generic")
    pkg_search_module(BLAS blas)
    if (BLAS_FOUND)
      set (BLAS_VENDOR_FOUND "Netlib or other Generic libblas")
    endif()
  elseif(BLA_VENDOR STREQUAL "All")
    if ((NOT BLAS_FOUND) AND
	(NOT BLA_F95) AND NOT
	((_LANGUAGES_ MATCHES Fortran) AND (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
	AND NOT (CMAKE_C_COMPILER_ID STREQUAL "GNU"))
      pkg_search_module(BLAS mkl-${MKL_STR_BLA_STATIC}-lp64-iomp)
      if (BLAS_FOUND)
	set (BLAS_VENDOR_FOUND "Intel MKL Multithreaded")
	set (BLAS_LIBRARIES_PAR "${BLAS_LIBRARIES}")
      endif()
      pkg_search_module(BLAS mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      if (BLAS_FOUND)
	set (BLAS_VENDOR_FOUND "Intel MKL Sequential")
	set (BLAS_LIBRARIES_SEQ "${BLAS_LIBRARIES}")
      endif()
    endif()
    if (NOT BLAS_FOUND)
      pkg_search_module(BLAS openblas)
      if (BLAS_FOUND)
	set (BLAS_VENDOR_FOUND "Openblas")
      endif()
    endif()
    if (NOT BLAS_FOUND)
      pkg_search_module(BLAS blas)
      if (BLAS_FOUND)
	set (BLAS_VENDOR_FOUND "Netlib or other Generic libblas")
      endif()
    endif()
  endif()

  if (NOT BLAS_FIND_QUIETLY)
    if (BLAS_FOUND AND BLAS_LIBRARIES)
      message(STATUS "Looking for BLAS - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for BLAS - not found using PkgConfig."
	"\n   Perhaps you should add the directory containing blas.pc to"
	"\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (BLAS_FOUND AND BLAS_LIBRARIES)
    set(BLAS_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(BLAS)
  else()
    set(BLAS_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif()

if( (NOT BLAS_FOUND_WITH_PKGCONFIG) OR BLAS_GIVEN_BY_USER )
  #BLAS in intel mkl 10 library? (em64t 64bit)
  if (BLA_VENDOR MATCHES "Intel*" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES OR BLA_VENDOR MATCHES "Intel*")
      # Looking for include
      # -------------------

      # Add system include paths to search include
      # ------------------------------------------
      unset(_inc_env)
      set(ENV_MKLROOT "$ENV{MKLROOT}")
      set(ENV_BLAS_DIR "$ENV{BLAS_DIR}")
      set(ENV_BLAS_INCDIR "$ENV{BLAS_INCDIR}")
      if(ENV_BLAS_INCDIR)
	list(APPEND _inc_env "${ENV_BLAS_INCDIR}")
      elseif(ENV_BLAS_DIR)
	list(APPEND _inc_env "${ENV_BLAS_DIR}")
	list(APPEND _inc_env "${ENV_BLAS_DIR}/include")
      else()
	if (ENV_MKLROOT)
	  list(APPEND _inc_env "${ENV_MKLROOT}/include")
	endif()
	# system variables
	if(WIN32)
	  string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
	  list(APPEND _inc_env "${_path_env}")
	else()
	  string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
	  list(APPEND _inc_env "${_path_env}")
	  string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
	  list(APPEND _inc_env "${_path_env}")
	  string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
	  list(APPEND _inc_env "${_path_env}")
	  string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
	  list(APPEND _inc_env "${_path_env}")
	endif()
      endif()
      list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
      list(REMOVE_DUPLICATES _inc_env)

      # set paths where to look for
      set(PATH_TO_LOOK_FOR "${_inc_env}")

      # Try to find the mkl header in the given paths
      # -------------------------------------------------
      # call cmake macro to find the header path
      if(BLAS_INCDIR)
	set(BLAS_mkl.h_DIRS "BLAS_mkl.h_DIRS-NOTFOUND")
	find_path(BLAS_mkl.h_DIRS
	  NAMES mkl.h
	  HINTS ${BLAS_INCDIR})
      else()
	if(BLAS_DIR)
	  set(BLAS_mkl.h_DIRS "BLAS_mkl.h_DIRS-NOTFOUND")
	  find_path(BLAS_mkl.h_DIRS
	    NAMES mkl.h
	    HINTS ${BLAS_DIR}
	    PATH_SUFFIXES "include")
	else()
	  set(BLAS_mkl.h_DIRS "BLAS_mkl.h_DIRS-NOTFOUND")
	  find_path(BLAS_mkl.h_DIRS
	    NAMES mkl.h
	    HINTS ${PATH_TO_LOOK_FOR})
	endif()
      endif()
      mark_as_advanced(BLAS_mkl.h_DIRS)

      # If found, add path to cmake variable
      # ------------------------------------
      if (BLAS_mkl.h_DIRS)
	set(BLAS_INCLUDE_DIRS "${BLAS_mkl.h_DIRS}")
      else ()
	set(BLAS_INCLUDE_DIRS "BLAS_INCLUDE_DIRS-NOTFOUND")
	if(NOT BLAS_FIND_QUIETLY)
	  message(STATUS "Looking for BLAS -- mkl.h not found")
	endif()
      endif()

      list(APPEND _libdir "$ENV{LIBRARY_PATH}")
      if (WIN32)
	string(REPLACE ":" ";" _libdir2 "$ENV{LIB}")
      elseif (APPLE)
	string(REPLACE ":" ";" _libdir2 "$ENV{DYLD_LIBRARY_PATH}")
      else ()
	string(REPLACE ":" ";" _libdir2 "$ENV{LD_LIBRARY_PATH}")
      endif ()
      list(APPEND _libdir "${_libdir2}")
      list(APPEND _libdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
      # libiomp5
      # --------
      set(OMP_iomp5_LIBRARY "OMP_iomp5_LIBRARY-NOTFOUND")
      find_library(OMP_iomp5_LIBRARY
	NAMES iomp5
	HINTS ${_libdir}
	)
      mark_as_advanced(OMP_iomp5_LIBRARY)
      set(OMP_LIB "")
      # libgomp
      # -------
      set(OMP_gomp_LIBRARY "OMP_gomp_LIBRARY-NOTFOUND")
      find_library(OMP_gomp_LIBRARY
	NAMES gomp
	HINTS ${_libdir}
	)
      mark_as_advanced(OMP_gomp_LIBRARY)
      # choose one or another depending on the compilo
      if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
	if (OMP_gomp_LIBRARY)
	  set(OMP_LIB "${OMP_gomp_LIBRARY}")
	endif()
      else(CMAKE_C_COMPILER_ID STREQUAL "Intel")
	if (OMP_iomp5_LIBRARY)
	  set(OMP_LIB "${OMP_iomp5_LIBRARY}")
	endif()
      endif()

      if (UNIX AND NOT WIN32)
	# m
	find_library(M_LIBRARY
	  NAMES m
	  HINTS ${_libdir})
	mark_as_advanced(M_LIBRARY)
	if(M_LIBRARY)
	  set(LM "${M_LIBRARY}")
	else()
	  set(LM "")
	endif()
	# dl
	find_library(
	  DL_LIBRARY
	  NAMES dl
	  HINTS ${_libdir}
	  )
	mark_as_advanced(DL_LIBRARY)
	if(DL_LIBRARY)
	  set(LDL "${DL_LIBRARY}")
	else()
	  set(LDL "")
	endif()
	# Fortran
	set(LGFORTRAN "")
	if (CMAKE_C_COMPILER_ID MATCHES "GNU")
	  find_library(
	    FORTRAN_gfortran_LIBRARY
	    NAMES gfortran
	    HINTS ${_libdir}
	    )
	  mark_as_advanced(FORTRAN_gfortran_LIBRARY)
	  if (FORTRAN_gfortran_LIBRARY)
	    set(LGFORTRAN "${FORTRAN_gfortran_LIBRARY}")
	  endif()
	elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
	  find_library(
	    FORTRAN_ifcore_LIBRARY
	    NAMES ifcore
	    HINTS ${_libdir}
	    )
	  mark_as_advanced(FORTRAN_ifcore_LIBRARY)
	  if (FORTRAN_ifcore_LIBRARY)
	    set(LGFORTRAN "{FORTRAN_ifcore_LIBRARY}")
	  endif()
	endif()
	set(BLAS_CFLAGS_OTHER "")
	# OpenMP
	if (BLA_VENDOR STREQUAL "Intel10_64lp")
	  find_package(OpenMP REQUIRED)
	endif()
	if(OPENMP_C_FLAGS)
	  list(APPEND BLAS_CFLAGS_OTHER "${OPENMP_C_FLAGS}")
	endif()
	if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
	  if (BLA_VENDOR STREQUAL "Intel10_32")
	    list(APPEND BLAS_CFLAGS_OTHER "-m32")
	  else()
	    list(APPEND BLAS_CFLAGS_OTHER "-m64")
	  endif()
	endif()
	set(additional_flags "")
	if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
	  set(additional_flags "-Wl,--no-as-needed")
	endif()
      endif ()

      if (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
	find_package(Threads)
	if( THREADS_FOUND )
	  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
	endif ()

	set(BLAS_SEARCH_LIBS "")

	if(BLA_F95)

	  set(BLAS_mkl_SEARCH_SYMBOL sgemm_f95)
	  set(_LIBRARIES BLAS95_LIBRARIES)
	  if (WIN32)
	    if (BLA_STATIC)
	      set(BLAS_mkl_DLL_SUFFIX "")
	    else()
	      set(BLAS_mkl_DLL_SUFFIX "_dll")
	    endif()

	    # Find the main file (32-bit or 64-bit)
	    set(BLAS_SEARCH_LIBS_WIN_MAIN "")
	    if (BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS_WIN_MAIN
		"mkl_blas95${BLAS_mkl_DLL_SUFFIX} mkl_intel_c${BLAS_mkl_DLL_SUFFIX}")
	    endif()
	    if (BLA_VENDOR STREQUAL "Intel10_64lp*" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS_WIN_MAIN
		"mkl_blas95_lp64${BLAS_mkl_DLL_SUFFIX} mkl_intel_lp64${BLAS_mkl_DLL_SUFFIX}")
	    endif ()

	    # Add threading/sequential libs
	    set(BLAS_SEARCH_LIBS_WIN_THREAD "")
	    if (BLA_VENDOR STREQUAL "*_seq" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"mkl_sequential${BLAS_mkl_DLL_SUFFIX}")
	    endif()
	    if (NOT BLA_VENDOR STREQUAL "*_seq" OR BLA_VENDOR STREQUAL "All")
	      # old version
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"libguide40 mkl_intel_thread${BLAS_mkl_DLL_SUFFIX}")
	      # mkl >= 10.3
	      # sequential version
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"mkl_sequential${BLAS_mkl_DLL_SUFFIX}")
	      # multithreaded version
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"libiomp5md mkl_intel_thread${BLAS_mkl_DLL_SUFFIX}")
	    endif()

	    # Cartesian product of the above
	    foreach (MAIN ${BLAS_SEARCH_LIBS_WIN_MAIN})
	      foreach (THREAD ${BLAS_SEARCH_LIBS_WIN_THREAD})
		list(APPEND BLAS_SEARCH_LIBS
		  "${MAIN} ${THREAD} mkl_core${BLAS_mkl_DLL_SUFFIX}")
	      endforeach()
	    endforeach()
	  else (WIN32)
	    if( _LANGUAGES_ MATCHES Fortran )
	      if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		set(BLAS_mkl_intface "gf")
	      else()
		set(BLAS_mkl_intface "intel")
	      endif()
	    else()
	      set(BLAS_mkl_intface "intel")
	    endif()

	    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
	      set(BLAS_mkl_threading "gnu")
	    else()
	      set(BLAS_mkl_threading "intel")
	    endif()

	    if (BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_blas95 mkl_${BLAS_mkl_intface} mkl_${BLAS_mkl_threading}_thread mkl_core")
	    endif ()
	    if (BLA_VENDOR STREQUAL "Intel10_64lp" OR BLA_VENDOR STREQUAL "All")
	      # old version
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_blas95 mkl_${BLAS_mkl_intface}_lp64 mkl_${BLAS_mkl_threading}_thread mkl_core guide")
	      # mkl >= 10.3
	      # sequential version
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_blas95_lp64 mkl_${BLAS_mkl_intface}_lp64 mkl_sequential mkl_core")
	      # multithreaded version
		list(APPEND BLAS_SEARCH_LIBS
		"mkl_blas95_lp64 mkl_${BLAS_mkl_intface}_lp64 mkl_${BLAS_mkl_threading}_thread mkl_core")
	    endif ()
	    if (BLA_VENDOR STREQUAL "Intel10_64lp_seq" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_blas95_lp64 mkl_${BLAS_mkl_intface}_lp64 mkl_sequential mkl_core")
	      if (BLA_VENDOR STREQUAL "Intel10_64lp_seq")
		set(OMP_LIB "")
	      endif()
	    endif ()
	  endif (WIN32)

	else (BLA_F95)

	  set(BLAS_mkl_SEARCH_SYMBOL sgemm)
	  set(_LIBRARIES BLAS_LIBRARIES)
	  if (WIN32)
	    if (BLA_STATIC)
	      set(BLAS_mkl_DLL_SUFFIX "")
	    else()
	      set(BLAS_mkl_DLL_SUFFIX "_dll")
	    endif()

	    # Find the main file (32-bit or 64-bit)
	    set(BLAS_SEARCH_LIBS_WIN_MAIN "")
	    if (BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS_WIN_MAIN
		"mkl_intel_c${BLAS_mkl_DLL_SUFFIX}")
	    endif()
	    if (BLA_VENDOR STREQUAL "Intel10_64lp*" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS_WIN_MAIN
		"mkl_intel_lp64${BLAS_mkl_DLL_SUFFIX}")
	    endif ()

	    # Add threading/sequential libs
	    set(BLAS_SEARCH_LIBS_WIN_THREAD "")
	    if (NOT BLA_VENDOR STREQUAL "*_seq" OR BLA_VENDOR STREQUAL "All")
	      # old version
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"libguide40 mkl_intel_thread${BLAS_mkl_DLL_SUFFIX}")
	      # mkl >= 10.3
	      # sequential version
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"mkl_sequential${BLAS_mkl_DLL_SUFFIX}")
	      # multithreaded version
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"libiomp5md mkl_intel_thread${BLAS_mkl_DLL_SUFFIX}")
	    endif()
	    if (BLA_VENDOR STREQUAL "*_seq" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS_WIN_THREAD
		"mkl_sequential${BLAS_mkl_DLL_SUFFIX}")
	    endif()

	    # Cartesian product of the above
	    foreach (MAIN ${BLAS_SEARCH_LIBS_WIN_MAIN})
	      foreach (THREAD ${BLAS_SEARCH_LIBS_WIN_THREAD})
		list(APPEND BLAS_SEARCH_LIBS
		  "${MAIN} ${THREAD} mkl_core${BLAS_mkl_DLL_SUFFIX}")
	      endforeach()
	    endforeach()
	  else (WIN32)
	    if( _LANGUAGES_ MATCHES Fortran )
	      if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
		set(BLAS_mkl_intface "gf")
	      else()
		set(BLAS_mkl_intface "intel")
	      endif()
	    else()
	      set(BLAS_mkl_intface "intel")
	    endif()

	    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
	      set(BLAS_mkl_threading "gnu")
	    else()
	      set(BLAS_mkl_threading "intel")
	    endif()

	    if (BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_${BLAS_mkl_intface} mkl_${BLAS_mkl_threading}_thread mkl_core")
	    endif ()
	    if (BLA_VENDOR STREQUAL "Intel10_64lp" OR BLA_VENDOR STREQUAL "All")
	      # old version
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_${BLAS_mkl_intface}_lp64 mkl_${BLAS_mkl_threading}_thread mkl_core guide")
	      # mkl >= 10.3
	      # sequential version
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_${BLAS_mkl_intface}_lp64 mkl_sequential mkl_core")
	      # multithreaded version
		list(APPEND BLAS_SEARCH_LIBS
		"mkl_${BLAS_mkl_intface}_lp64 mkl_${BLAS_mkl_threading}_thread mkl_core")
	    endif ()
	    if (BLA_VENDOR STREQUAL "Intel10_64lp_seq" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_${BLAS_mkl_intface}_lp64 mkl_sequential mkl_core")
	      if (BLA_VENDOR STREQUAL "Intel10_64lp_seq")
		set(OMP_LIB "")
	      endif()
	    endif ()
	    #older versions of intel mkl libs
	    if (BLA_VENDOR STREQUAL "Intel" OR BLA_VENDOR STREQUAL "All")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_ia32")
	      list(APPEND BLAS_SEARCH_LIBS
		"mkl_em64t")
	    endif ()
	  endif (WIN32)

	endif (BLA_F95)

	foreach (IT ${BLAS_SEARCH_LIBS})
	  string(REPLACE " " ";" SEARCH_LIBS ${IT})
	  # determine if we look at a sequential or multithreaded version
	  set(MKL_MT OFF)
	  foreach(_lib ${SEARCH_LIBS})
	    if (${_lib} MATCHES "thread")
	      set (MKL_MT ON)
	    endif()
	  endforeach()
	  if (MKL_MT)
	    check_fortran_libraries(
	      ${_LIBRARIES}
	      BLAS
	      ${BLAS_mkl_SEARCH_SYMBOL}
	      "${additional_flags}"
	      "${SEARCH_LIBS}"
	      "${OMP_LIB};${CMAKE_THREAD_LIBS_INIT};${LM};${LDL}"
	      )
	  else()
	    check_fortran_libraries(
	      ${_LIBRARIES}
	      BLAS
	      ${BLAS_mkl_SEARCH_SYMBOL}
	      "${additional_flags}"
	      "${SEARCH_LIBS}"
	      "${CMAKE_THREAD_LIBS_INIT};${LM};${LDL}"
	      )
	  endif()
	  if(${_LIBRARIES})
	    set(BLAS_LDFLAGS_OTHER "${additional_flags}")
	    if (MKL_MT)
	      set (BLAS_LIBRARIES_PAR "${${_LIBRARIES}}")
	    else()
	      set (BLAS_LIBRARIES_SEQ "${${_LIBRARIES}}")
	    endif()
	    set(${_LIBRARIES}_COPY "${${_LIBRARIES}}")
	  endif()
	endforeach ()
	if (${_LIBRARIES}_COPY)
	  set(${_LIBRARIES} "${${_LIBRARIES}_COPY}")
	endif()
	if(NOT BLAS_FIND_QUIETLY)
	  if(${_LIBRARIES})
	    message(STATUS "Looking for MKL BLAS: found")
	  else()
	    message(STATUS "Looking for MKL BLAS: not found")
	  endif()
	endif()
	if (${_LIBRARIES})
	    set (BLAS_VENDOR_FOUND "Intel MKL")
	endif()
      endif (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
    endif(NOT BLAS_LIBRARIES OR BLA_VENDOR MATCHES "Intel*")
    if (NOT BLAS_LIBRARIES)
      set(BLAS_INCLUDE_DIRS "")
      set(BLAS_CFLAGS_OTHER "")
      set(BLAS_LDFLAGS_OTHER "")
    endif()
    unset(BLAS_mkl_intface)
    unset(BLAS_mkl_threading)
  endif (BLA_VENDOR MATCHES "Intel*" OR BLA_VENDOR STREQUAL "All")

  if (BLA_VENDOR MATCHES "ARMPL*" OR BLA_VENDOR STREQUAL "All")
    if(NOT BLAS_LIBRARIES OR BLA_VENDOR MATCHES "ARMPL*")
      # Looking for include
      # -------------------

      # Add system include paths to search include
      # ------------------------------------------
      unset(_inc_env)
      set(ENV_ARMPL "$ENV{ARMPL_DIR}")
      set(ENV_BLAS_DIR "$ENV{BLAS_DIR}")
      set(ENV_BLAS_INCDIR "$ENV{BLAS_INCDIR}")
      if(ENV_BLAS_INCDIR)
	list(APPEND _inc_env "${ENV_BLAS_INCDIR}")
      elseif(ENV_BLAS_DIR)
	list(APPEND _inc_env "${ENV_BLAS_DIR}")
	list(APPEND _inc_env "${ENV_BLAS_DIR}/include")
      else()
	if (ENV_ARMPL)
	  list(APPEND _inc_env "${ENV_ARMPL}/include")
	endif()
	# system variables
	if(WIN32)
	  string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
	  list(APPEND _inc_env "${_path_env}")
	else()
	  string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
	  list(APPEND _inc_env "${_path_env}")
	  string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
	  list(APPEND _inc_env "${_path_env}")
	  string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
	  list(APPEND _inc_env "${_path_env}")
	  string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
	  list(APPEND _inc_env "${_path_env}")
	endif()
      endif()
      list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
      list(REMOVE_DUPLICATES _inc_env)

      # set paths where to look for
      set(PATH_TO_LOOK_FOR "${_inc_env}")

      # Try to find the armpl header in the given paths
      # -------------------------------------------------
      # call cmake macro to find the header path
      if(BLAS_INCDIR)
	set(BLAS_armpl.h_DIRS "BLAS_armpl.h_DIRS-NOTFOUND")
	find_path(BLAS_armpl.h_DIRS
	  NAMES armpl.h
	  HINTS ${BLAS_INCDIR})
      else()
	if(BLAS_DIR)
	  set(BLAS_armpl.h_DIRS "BLAS_armpl.h_DIRS-NOTFOUND")
	  find_path(BLAS_armpl.h_DIRS
	    NAMES armpl.h
	    HINTS ${BLAS_DIR}
	    PATH_SUFFIXES "include")
	else()
	  set(BLAS_armpl.h_DIRS "BLAS_armpl.h_DIRS-NOTFOUND")
	  find_path(BLAS_armpl.h_DIRS
	    NAMES armpl.h
	    HINTS ${PATH_TO_LOOK_FOR})
	endif()
      endif()
      mark_as_advanced(BLAS_armpl.h_DIRS)

      # If found, add path to cmake variable
      # ------------------------------------
      if (BLAS_armpl.h_DIRS)
	set(BLAS_INCLUDE_DIRS "${BLAS_armpl.h_DIRS}")
      else ()
	set(BLAS_INCLUDE_DIRS "BLAS_INCLUDE_DIRS-NOTFOUND")
	if(NOT BLAS_FIND_QUIETLY)
	  message(STATUS "Looking for BLAS -- armpl.h not found")
	endif()
      endif()

      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"armpl"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for ARMPL BLAS: found")
	else()
	  message(STATUS "Looking for ARMPL BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "ARMPL")
      endif()

    endif (NOT BLAS_LIBRARIES OR BLA_VENDOR MATCHES "ARMPL*")
  endif (BLA_VENDOR MATCHES "ARMPL*" OR BLA_VENDOR STREQUAL "All")

  if (BLA_VENDOR STREQUAL "Goto" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      # gotoblas (http://www.tacc.utexas.edu/tacc-projects/gotoblas2)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"goto2"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for Goto BLAS: found")
	else()
	  message(STATUS "Looking for Goto BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "Goto")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "Goto" OR BLA_VENDOR STREQUAL "All")


  # OpenBlas
  if (BLA_VENDOR STREQUAL "Open" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      # openblas (http://www.openblas.net/)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"openblas"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for Open BLAS: found")
	else()
	  message(STATUS "Looking for Open BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "Openblas")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "Open" OR BLA_VENDOR STREQUAL "All")


  # EigenBlas
  if (BLA_VENDOR STREQUAL "Eigen" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      # eigenblas (http://eigen.tuxfamily.org/index.php?title=Main_Page)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"eigen_blas"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for Eigen BLAS: found")
	else()
	  message(STATUS "Looking for Eigen BLAS: not found")
	endif()
      endif()
    endif()

    if(NOT BLAS_LIBRARIES)
      # eigenblas (http://eigen.tuxfamily.org/index.php?title=Main_Page)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"eigen_blas_static"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for Eigen BLAS: found")
	else()
	  message(STATUS "Looking for Eigen BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "Eigen")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "Eigen" OR BLA_VENDOR STREQUAL "All")


  if (BLA_VENDOR STREQUAL "ATLAS" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      # BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	dgemm
	""
	"f77blas;atlas"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for Atlas BLAS: found")
	else()
	  message(STATUS "Looking for Atlas BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "Atlas")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "ATLAS" OR BLA_VENDOR STREQUAL "All")


  # BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
  if (BLA_VENDOR STREQUAL "PhiPACK" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"sgemm;dgemm;blas"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for PhiPACK BLAS: found")
	else()
	  message(STATUS "Looking for PhiPACK BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "PhiPACK")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "PhiPACK" OR BLA_VENDOR STREQUAL "All")


  # BLAS in Alpha CXML library?
  if (BLA_VENDOR STREQUAL "CXML" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"cxml"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for CXML BLAS: found")
	else()
	  message(STATUS "Looking for CXML BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "CXML")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "CXML" OR BLA_VENDOR STREQUAL "All")


  # BLAS in Alpha DXML library? (now called CXML, see above)
  if (BLA_VENDOR STREQUAL "DXML" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"dxml"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for DXML BLAS: found")
	else()
	  message(STATUS "Looking for DXML BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "DXML")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "DXML" OR BLA_VENDOR STREQUAL "All")


  # BLAS in Sun Performance library?
  if (BLA_VENDOR STREQUAL "SunPerf" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	"-xlic_lib=sunperf"
	"sunperf;sunmath"
	""
	)
      if(BLAS_LIBRARIES)
	set(BLAS_LDFLAGS_OTHER "-xlic_lib=sunperf")
      endif()
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for SunPerf BLAS: found")
	else()
	  message(STATUS "Looking for SunPerf BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "SunPerf")
      endif()
    endif()

  endif ()


  # BLAS in SCSL library?  (SGI/Cray Scientific Library)
  if (BLA_VENDOR STREQUAL "SCSL" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"scsl"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for SCSL BLAS: found")
	else()
	  message(STATUS "Looking for SCSL BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "SunPerf")
      endif()
    endif()

  endif ()


  # BLAS in SGIMATH library?
  if (BLA_VENDOR STREQUAL "SGIMATH" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"complib.sgimath"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for SGIMATH BLAS: found")
	else()
	  message(STATUS "Looking for SGIMATH BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "SGIMATH")
      endif()
    endif()

  endif ()


  # BLAS in IBM ESSL library (requires generic BLAS lib, too)
  if (BLA_VENDOR STREQUAL "IBMESSL" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"essl"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for IBM ESSL BLAS: found")
	else()
	  message(STATUS "Looking for IBM ESSL BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "IBM ESSL")
      endif()
    endif()

  endif ()

  # BLAS in IBM ESSL_MT library (requires generic BLAS lib, too)
  if (BLA_VENDOR STREQUAL "IBMESSLMT" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"esslsmp"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for IBM ESSL MT BLAS: found")
	else()
	  message(STATUS "Looking for IBM ESSL MT BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "IBM ESSL MT")
      endif()
    endif()

  endif ()


  #BLAS in acml library?
  if (BLA_VENDOR MATCHES "ACML.*" OR BLA_VENDOR STREQUAL "All")

    if( ((BLA_VENDOR STREQUAL "ACML") AND (NOT BLAS_ACML_LIB_DIRS)) OR
	((BLA_VENDOR STREQUAL "ACML_MP") AND (NOT BLAS_ACML_MP_LIB_DIRS)) OR
	((BLA_VENDOR STREQUAL "ACML_GPU") AND (NOT BLAS_ACML_GPU_LIB_DIRS)))

      # try to find acml in "standard" paths
      if( WIN32 )
	file( GLOB _ACML_ROOT "C:/AMD/acml*/ACML-EULA.txt" )
      else()
	file( GLOB _ACML_ROOT "/opt/acml*/ACML-EULA.txt" )
      endif()
      if( WIN32 )
	file( GLOB _ACML_GPU_ROOT "C:/AMD/acml*/GPGPUexamples" )
      else()
	file( GLOB _ACML_GPU_ROOT "/opt/acml*/GPGPUexamples" )
      endif()
      list(GET _ACML_ROOT 0 _ACML_ROOT)
      list(GET _ACML_GPU_ROOT 0 _ACML_GPU_ROOT)

      if( _ACML_ROOT )

	get_filename_component( _ACML_ROOT ${_ACML_ROOT} PATH )
	if( SIZEOF_INTEGER EQUAL 8 )
	  set( _ACML_PATH_SUFFIX "_int64" )
	else()
	  set( _ACML_PATH_SUFFIX "" )
	endif()
	if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
	  set( _ACML_COMPILER32 "ifort32" )
	  set( _ACML_COMPILER64 "ifort64" )
	elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro" )
	  set( _ACML_COMPILER32 "sun32" )
	  set( _ACML_COMPILER64 "sun64" )
	elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
	  set( _ACML_COMPILER32 "pgi32" )
	  if( WIN32 )
	    set( _ACML_COMPILER64 "win64" )
	  else()
	    set( _ACML_COMPILER64 "pgi64" )
	  endif()
	elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Open64" )
	  # 32 bit builds not supported on Open64 but for code simplicity
	  # We'll just use the same directory twice
	  set( _ACML_COMPILER32 "open64_64" )
	  set( _ACML_COMPILER64 "open64_64" )
	elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" )
	  set( _ACML_COMPILER32 "nag32" )
	  set( _ACML_COMPILER64 "nag64" )
	else()
	  set( _ACML_COMPILER32 "gfortran32" )
	  set( _ACML_COMPILER64 "gfortran64" )
	endif()

	if( BLA_VENDOR STREQUAL "ACML_MP" )
	  set(_ACML_MP_LIB_DIRS
	    "${_ACML_ROOT}/${_ACML_COMPILER32}_mp${_ACML_PATH_SUFFIX}/lib"
	    "${_ACML_ROOT}/${_ACML_COMPILER64}_mp${_ACML_PATH_SUFFIX}/lib" )
	else()
	  set(_ACML_LIB_DIRS
	    "${_ACML_ROOT}/${_ACML_COMPILER32}${_ACML_PATH_SUFFIX}/lib"
	    "${_ACML_ROOT}/${_ACML_COMPILER64}${_ACML_PATH_SUFFIX}/lib" )
	endif()

      endif(_ACML_ROOT)

    elseif(BLAS_${BLA_VENDOR}_LIB_DIRS)

      set(_${BLA_VENDOR}_LIB_DIRS ${BLAS_${BLA_VENDOR}_LIB_DIRS})

    endif()

    if( BLA_VENDOR STREQUAL "ACML_MP" )
      foreach( BLAS_ACML_MP_LIB_DIRS ${_ACML_MP_LIB_DIRS})
	check_fortran_libraries (
	  BLAS_LIBRARIES
	  BLAS
	  sgemm
	  "" "acml_mp;acml_mv" "" ${BLAS_ACML_MP_LIB_DIRS}
	  )
	if( BLAS_LIBRARIES )
	  set (BLAS_VENDOR_FOUND "ACML_MP")
	  break()
	endif()
      endforeach()
    elseif( BLA_VENDOR STREQUAL "ACML_GPU" )
      foreach( BLAS_ACML_GPU_LIB_DIRS ${_ACML_GPU_LIB_DIRS})
	check_fortran_libraries (
	  BLAS_LIBRARIES
	  BLAS
	  sgemm
	  "" "acml;acml_mv;CALBLAS" "" ${BLAS_ACML_GPU_LIB_DIRS}
	  )
	if( BLAS_LIBRARIES )
	  set (BLAS_VENDOR_FOUND "ACML_GPU")
	  break()
	endif()
      endforeach()
    else()
      foreach( BLAS_ACML_LIB_DIRS ${_ACML_LIB_DIRS} )
	check_fortran_libraries (
	  BLAS_LIBRARIES
	  BLAS
	  sgemm
	  "" "acml;acml_mv" "" ${BLAS_ACML_LIB_DIRS}
	  )
	if( BLAS_LIBRARIES )
	  set (BLAS_VENDOR_FOUND "ACML")
	  break()
	endif()
      endforeach()
    endif()

    # Either acml or acml_mp should be in LD_LIBRARY_PATH but not both
    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"acml;acml_mv"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for ACML BLAS: found")
	else()
	  message(STATUS "Looking for ACML BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "ACML")
      endif()
    endif()

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"acml_mp;acml_mv"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for ACML BLAS: found")
	else()
	  message(STATUS "Looking for ACML BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "ACML")
      endif()
    endif()

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	sgemm
	""
	"acml;acml_mv;CALBLAS"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for ACML BLAS: found")
	else()
	  message(STATUS "Looking for ACML BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "ACML")
      endif()
    endif()

  endif (BLA_VENDOR MATCHES "ACML.*" OR BLA_VENDOR STREQUAL "All") # ACML


  # Apple BLAS library?
  if (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")

    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	dgemm
	""
	"Accelerate"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for Apple BLAS: found")
	else()
	  message(STATUS "Looking for Apple BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "Apple Accelerate")
      endif()
    endif()

  endif (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")


  if (BLA_VENDOR STREQUAL "NAS" OR BLA_VENDOR STREQUAL "All")

    if ( NOT BLAS_LIBRARIES )
      check_fortran_libraries(
	BLAS_LIBRARIES
	BLAS
	dgemm
	""
	"vecLib"
	""
	)
      if(NOT BLAS_FIND_QUIETLY)
	if(BLAS_LIBRARIES)
	  message(STATUS "Looking for NAS BLAS: found")
	else()
	  message(STATUS "Looking for NAS BLAS: not found")
	endif()
      endif()
      if (BLAS_LIBRARIES)
	set (BLAS_VENDOR_FOUND "NAS")
      endif()
    endif ()

  endif (BLA_VENDOR STREQUAL "NAS" OR BLA_VENDOR STREQUAL "All")


  # Generic BLAS library?
  if (BLA_VENDOR STREQUAL "Generic" OR BLA_VENDOR STREQUAL "All")

    set(BLAS_SEARCH_LIBS "blas;blas_LINUX;blas_MAC;blas_WINDOWS;refblas")
    foreach (SEARCH_LIB ${BLAS_SEARCH_LIBS})
      if (BLAS_LIBRARIES)
      else ()
	check_fortran_libraries(
	  BLAS_LIBRARIES
	  BLAS
	  sgemm
	  ""
	  "${SEARCH_LIB}"
	  "${LGFORTRAN}"
	  )
	if(NOT BLAS_FIND_QUIETLY)
	  if(BLAS_LIBRARIES)
	    message(STATUS "Looking for Generic BLAS: found")
	  else()
	    message(STATUS "Looking for Generic BLAS: not found")
	  endif()
	endif()
	if (BLAS_LIBRARIES)
	  set (BLAS_VENDOR_FOUND "Netlib or other Generic libblas")
	endif()
      endif()
    endforeach ()

  endif (BLA_VENDOR STREQUAL "Generic" OR BLA_VENDOR STREQUAL "All")

  # BLAS given user, we need to test it
  if (BLA_VENDOR STREQUAL "User" AND BLAS_LIBRARIES_USER)
    check_fortran_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "${BLAS_LIBRARIES_USER}"
      "${BLAS_LDFLAGS_USER}"
      )
    if(NOT BLAS_FIND_QUIETLY)
      if(BLAS_LIBRARIES)
	message(STATUS "Looking for User's BLAS: found")
      else()
	message(STATUS "Looking for User's BLAS: not found")
      endif()
    endif()
    if (BLAS_LIBRARIES)
      set (BLAS_VENDOR_FOUND "User")
    endif()
  endif()

  if(BLA_F95)

    if(BLAS95_LIBRARIES)
      set(BLAS95_FOUND TRUE)
    else()
      set(BLAS95_FOUND FALSE)
    endif()

    if(NOT BLAS_FIND_QUIETLY)
      if(BLAS95_FOUND)
	message(STATUS "A library with BLAS95 API found.")
	message(STATUS "BLAS_LIBRARIES ${BLAS_LIBRARIES}")
      else(BLAS95_FOUND)
	message(WARNING "BLA_VENDOR has been set to ${BLA_VENDOR} but blas 95 libraries could not be found or check of symbols failed."
	  "\nPlease indicate where to find blas libraries. You have three options:\n"
	  "- Option 1: Provide the installation directory of BLAS library with cmake option: -DBLAS_DIR=your/path/to/blas\n"
	  "- Option 2: Provide the directory where to find BLAS libraries with cmake option: -DBLAS_LIBDIR=your/path/to/blas/libs\n"
	  "- Option 3: Update your environment variable (LIBRARY_PATH or Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
	  "\nTo follow libraries detection more precisely you can activate a verbose mode with -DBLAS_VERBOSE=ON at cmake configure."
	  "\nYou could also specify a BLAS vendor to look for by setting -DBLA_VENDOR=blas_vendor_name."
	  "\nList of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
	  "Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model), Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
	  "Intel( older versions of mkl 32 and 64 bit), ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
	if(BLAS_FIND_REQUIRED)
	  message(FATAL_ERROR
	    "A required library with BLAS95 API not found. Please specify library location.")
	else()
	  message(STATUS
	    "A library with BLAS95 API not found. Please specify library location.")
	endif()
      endif(BLAS95_FOUND)
    endif(NOT BLAS_FIND_QUIETLY)

    set(BLAS_FOUND TRUE)
    set(BLAS_LIBRARIES "${BLAS95_LIBRARIES}")
    if (NOT BLAS_LIBRARIES_DEP)
      set(BLAS_LIBRARIES_DEP "${BLAS95_LIBRARIES}")
    endif()

  else(BLA_F95)

    if(BLAS_LIBRARIES)
      set(BLAS_FOUND TRUE)
      if (NOT BLAS_LIBRARIES_DEP)
	set(BLAS_LIBRARIES_DEP "${BLAS_LIBRARIES}")
      endif()
    else()
      set(BLAS_FOUND FALSE)
    endif()

    if(NOT BLAS_FIND_QUIETLY)
      if(BLAS_FOUND)
	message(STATUS "A library with BLAS API found.")
	message(STATUS "BLAS_LIBRARIES ${BLAS_LIBRARIES}")
      else(BLAS_FOUND)
	message(WARNING "BLA_VENDOR has been set to ${BLA_VENDOR} but blas libraries could not be found or check of symbols failed."
	  "\nPlease indicate where to find blas libraries. You have three options:\n"
	  "- Option 1: Provide the installation directory of BLAS library with cmake option: -DBLAS_DIR=your/path/to/blas\n"
	  "- Option 2: Provide the directory where to find BLAS libraries with cmake option: -DBLAS_LIBDIR=your/path/to/blas/libs\n"
	  "- Option 3: Update your environment variable (LIBRARY_PATH or Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
	  "\nTo follow libraries detection more precisely you can activate a verbose mode with -DBLAS_VERBOSE=ON at cmake configure."
	  "\nYou could also specify a BLAS vendor to look for by setting -DBLA_VENDOR=blas_vendor_name."
	  "\nList of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
	  "Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model), Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
	  "Intel( older versions of mkl 32 and 64 bit), ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
	if(BLAS_FIND_REQUIRED)
	  message(FATAL_ERROR
	    "A required library with BLAS API not found. Please specify library location.")
	else()
	  message(STATUS
	    "A library with BLAS API not found. Please specify library location.")
	endif()
      endif(BLAS_FOUND)
    endif(NOT BLAS_FIND_QUIETLY)

  endif(BLA_F95)
endif( (NOT BLAS_FOUND_WITH_PKGCONFIG) OR BLAS_GIVEN_BY_USER )

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_blas_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

if (BLAS_FOUND)
  list(GET BLAS_LIBRARIES 0 first_lib)
  # first lib may be -Wl,--start-group (MKL) which is not a lib
  if (NOT EXISTS ${first_lib})
    list(GET BLAS_LIBRARIES 1 first_lib)
  endif()
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT BLAS_LIBRARY_DIRS)
    set(BLAS_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(BLAS_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of BLAS library" FORCE)
  else()
    set(BLAS_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of BLAS library" FORCE)
  endif()
endif()

if(BLAS_FOUND)
  if(NOT TARGET blas::blas AND BLAS_LIBRARIES)
    add_library(blas::blas INTERFACE IMPORTED)
    set_property(TARGET blas::blas APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LDFLAGS_OTHER}")
    set_property(TARGET blas::blas APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
    set_property(TARGET blas::blas APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "${BLAS_CFLAGS_OTHER}")
    set_property(TARGET blas::blas APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}")
    message(STATUS "Morse : target blas::blas available")
  endif()
  if(NOT TARGET blas::blas_seq AND BLAS_LIBRARIES_SEQ)
    add_library(blas::blas_seq INTERFACE IMPORTED)
    set_property(TARGET blas::blas_seq APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LDFLAGS_OTHER}")
    set_property(TARGET blas::blas_seq APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES_SEQ}")
    set_property(TARGET blas::blas_seq APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "${BLAS_CFLAGS_OTHER}")
    set_property(TARGET blas::blas_seq APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}")
    message(STATUS "Morse : target blas::blas_seq available")
  endif()
  if(NOT TARGET blas::blas_par AND BLAS_LIBRARIES_PAR)
    add_library(blas::blas_par INTERFACE IMPORTED)
    set_property(TARGET blas::blas_par APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LDFLAGS_OTHER}")
    set_property(TARGET blas::blas_par APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES_PAR}")
    set_property(TARGET blas::blas_par APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "${BLAS_CFLAGS_OTHER}")
    set_property(TARGET blas::blas_par APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}")
    message(STATUS "Morse : target blas::blas_par available")
  endif()
  if(NOT TARGET blas::blas_95 AND BLAS95_LIBRARIES AND BLAS95_FOUND)
    add_library(blas::blas_95 INTERFACE IMPORTED)
    set_property(TARGET blas::blas_95 APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS_LDFLAGS_OTHER}")
    set_property(TARGET blas::blas_95 APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${BLAS95_LIBRARIES}")
    set_property(TARGET blas::blas_95 APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "${BLAS_CFLAGS_OTHER}")
    set_property(TARGET blas::blas_95 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${BLAS_INCLUDE_DIRS}")
    message(STATUS "Morse : target blas::blas_95 available")
  endif()
endif()
mark_as_advanced(BLAS_DIR)
mark_as_advanced(BLAS_DIR_FOUND)
