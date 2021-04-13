#-----------------------------------------------------------------------------
#
# SCALFMMConfig.cmake - SCALFMM CMake configuration file for external projects.
#
# This file is configured by SCALFMM and used by the SCALFMM.cmake module
# to load SCALFMM's settings for an external project.
#


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was scalfmmConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

if(NOT TARGET scalfmm::scalfmm)
  list(APPEND CMAKE_MODULE_PATH "/mnt/c/scalfmm_2/scalfmm/CMakeModules/morse_cmake/modules/find")
  set(IMPORTED_LIBS OpenMP;MPI;BLAS;LAPACK;FFTW;OpenMP;MPI;BLAS;LAPACK;FFTW;STARPU)
  set(BLA_VENDOR All)
  include(CMakeFindDependencyMacro)
  foreach(lib IN LISTS IMPORTED_LIBS)
    find_dependency(${lib})
    if(NOT ${lib}_FOUND)
      message(FATAL_ERROR "MISSING ${lib} DEPENDENCY !")
    else()
      message(STATUS "Found ${lib} dependency.")
    endif()
  endforeach()
  include("${CMAKE_CURRENT_LIST_DIR}/scalfmm-targets.cmake")
endif()

