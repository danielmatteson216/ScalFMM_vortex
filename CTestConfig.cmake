## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "ScalFMM")
set(CTEST_NIGHTLY_START_TIME "00:00:00 GMT")
set(CTEST_SUBMIT_URL "http://cdash.inria.fr/CDash/submit.php?project=scalfmm")

#--------------------------------------------------------------------
# BUILDNAME variable construction
# This variable will be used to set the build name which will appear
# on the Chameleon dashboard http://cdash.inria.fr/CDash/
#--------------------------------------------------------------------
# Start with the short system name, e.g. "Linux", "FreeBSD" or "Windows"
if(NOT BUILDNAME)

  set(BUILDNAME "${CMAKE_SYSTEM_NAME}")

  # Add i386 or amd64
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(BUILDNAME "${BUILDNAME}-amd64")
  else()
    set(BUILDNAME "${BUILDNAME}-i386")
  endif()

  # Add compiler name
  get_filename_component(CMAKE_CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
  set(BUILDNAME "${BUILDNAME}-${CMAKE_CXX_COMPILER_NAME}")

  # Add the build type, e.g. "Debug, Release..."
  if(CMAKE_BUILD_TYPE)
    set(BUILDNAME "${BUILDNAME}-${CMAKE_BUILD_TYPE}")
  endif(CMAKE_BUILD_TYPE)

  # Specific options of Scalfmm
  if(SCALFMM_USE_SSE)
    set(BUILDNAME "${BUILDNAME}-sse")
  endif()

  if(SCALFMM_USE_BLAS)
    set(BUILDNAME "${BUILDNAME}-blas")
  endif()

  if(SCALFMM_USE_MPI)
    set(BUILDNAME "${BUILDNAME}-mpi")
  endif()

endif()
