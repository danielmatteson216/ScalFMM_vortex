Downloading, Building and Installing ScalFMM   {#install}
============================================

On this page, you will find all the requirements for building and installing ScalFMM.

[TOC]

# Download ScalFMM   {#download}

ScalFMM can be downloaded from Inria's GitLab repository

  - https://gitlab.inria.fr/solverstack/ScalFMM

It can also be installed using Spack, which will automatically manage its dependencies. You will find more information about Spack and ScalFMM here: http://morse.gforge.inria.fr/spack/spack.html#sec-2-6 .

    git clone https://github.com/fpruvost/spack.git
    cd spack
    bin/spack install scalfmm

To receive updates about new versions, you can register to the users mailing list at http://lists.gforge.inria.fr/cgi-bin/mailman/listinfo/scalfmm-public-users . This list has a very low traffic.

# Dependencies   {#deps}

  - CMake v2.8.12 or later
  - A C++ compiler that supports
    - C++14
    - OpenMP

**Optional dependencies**

  - Doxygen 1.8.8 or later to build the documentation
  - MPI, for distributed support
  - BLAS
  - FFTW

# Build   {#build}

## Setup   {#setup}

Building ScalFMM requires the standard CMake workflow.


    cd scalfmm/Build
    cmake .. [-DSCALFMM_USE_MPI=ON] # if MPI is needed


The build may be configured using `ccmake .`.

To build an executable, run:

    make exec_name

The documentation can be built using `make doc` if you have Doxygen available.

**Example**

    cd scalfmm/Build
    cmake -DSCALFMM_USE_BLAS=ON -DSCALFMM_USE_MKL_AS_BLAS=ON \
        -DSCALFMM_USE_SSE=OFF -DSCALFMM_USE_AVX=ON ..


## Configuration   {#conf}

The following options are available.

 * `CMAKE_INSTALL_PREFIX`: where to install ScalFmm
 * `SCALFMM_USE_MPI`: to use and enable MPI. Warning, you need to use this parameter at the first cmake command you write.
 * `SCALFMM_ATTACHE_SOURCE`: build with -g (which enables debugging with release binaries)
 * `SCALFMM_BUILD_EXAMPLES`: build the examples
 * `SCALFMM_BUILD_TESTS`: build the tests
 * `SCALFMM_BUILD_UTESTS`: build the unit tests
 * `SCALFMM_BUILD_DOC`: enable make doc generate the documentation
 * `SCALFMM_USE_ADDONS`: activate add ons
 * `SCALFMM_ADDON_FMMAPI`: build Fmm Api
 * `ScalFMMUSE_MEM_STATS`: use memory stats (which count any new/delete done during a simulation)
 * `SCALFMM_USE_BLAS`: enable BLAS (needed by Chebyshev interpolation kernel)
 * `SCALFMM_USE_MKL_AS_BLAS`: use MKL as blas
 * `SCALFMM_USE_FFT`: Use FFTW needed by the uniform interpolation kernel
     * `SCALFMM_USE_MKL_AS_FFTW`: use MKL as FFTW
 * `SCALFMM_USE_LOG`: print output debug information during the execution
 * `SCALFMM_USE_ASSERT`: enable safe tests during execution
 * `SCALFMM_USE_SSE`: compile with SEE support
 * `SCALFMM_USE_AVX`: compile with AVX support


# Installation {#installation}

To install ScalFMM, use the `make install` command.

To run small tests on using ScalFMM, you can create a source file in the `Tests/Utils` folder. Run `cmake .` to update the available targets and compile your file as any other ScalFMM test.
