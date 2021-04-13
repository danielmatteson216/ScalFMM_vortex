# *** THIS REPO AND ITS CONTENTS ARE A MODIFIED VERSION OF SCALFMM 2.0, USED FOR VORTEX INTERACTIONS ***
The installation instructions below are from the Inria ScalFMM gitlab/github and detail how to setup ScalFMM


# ScalFMM: Fast Multipole Method

[![pipeline status](https://gitlab.inria.fr/solverstack/ScalFMM/badges/develop/pipeline.svg)](https://gitlab.inria.fr/solverstack/ScalFMM/commits/develop)
[![coverage report](https://gitlab.inria.fr/solverstack/ScalFMM/badges/develop/coverage.svg)](https://gitlab.inria.fr/solverstack/ScalFMM/commits/develop)

----

:warning: ScalFMM has moved to Inria's GitLab: https://gitlab.inria.fr/solverstack/ScalFMM

----

**ScalFMM** is a C++ library that implements a kernel independent Fast Multipole Method.


Copyright Inria, please read the licence.

## Requirements

  - CMake v3.10.0 or later
  - C++ compiler that supports
    - C++14 [compiler support list](http://en.cppreference.com/w/cpp/compiler_support)
    - [OpenMP](http://www.openmp.org/resources/openmp-compilers/)

The following are optional:

  - [Doxygen](http://www.stack.nl/~dimitri/doxygen/) to build the documentation.
  - An MPI implementation to build the distributed files.
  - Custom BLAS, FFT implementations.
  - [StarPU](http://starpu.gforge.inria.fr/) for the relevant FMM implementations.

## Get and Build ScalFMM

### Cloning

To use last development states of ScalFMM, please clone the master
  branch. Note that ScalFMM contains two git submodules `morse_cmake` and `inastemp`.
  To get sources please use these commands:
``` bash
git clone --recursive git@gitlab.inria.fr:solverstack/ScalFMM.git -b requested_branch
```
or
```bash
git clone git@gitlab.inria.fr:solverstack/ScalFMM.git
cd ScalFMM
git submodule init
git submodule update

``` 
### Building
You can do an out-of-source build by creating a `build` folder out of your clone or you can use the `Build`
folder inside your clone.

``` bash
cd scalfmm/Build
# Use cmake, with relevant options
cmake .. # -DSCALFMM_USE_MPI=ON
```

The build may be configured after the first CMake invocation using, for instance, `ccmake` or `cmake-gui`.

```bash
# Still in the Build folder
ccmake ../
# Or
cmake-gui ../
```

The binaries are then compiled calling `make` (or `ninja` if you specified it at the configure step).
They can be found in `scalfmm/Build/Tests/{Release,Debug}/...`

Invoke `make help` to see the available targets.
Gloabal targets are available :
* `scalfmm_examples` builds the examples in the `Example` folder,
* `scalfmm_utests` builds the unit tests in the `Utests` folder.

An example build using StarPU:

```bash
cmake .. -DSCALFMM_USE_BLAS=ON -DSCALFMM_USE_MKL_AS_BLAS=ON -DSCALFMM_USE_FFT=ON -DSCALFMM_USE_STARPU=ON
make all
```

You can also specify your install directory with `-DCMAKE_INSTALL_PREFIX=/path/to/your/install` and then
call `make install`

## Using ScalFMM in your project

To find ScalFMM, `pkgconfig` can be used within your CMake and all ScalFMM dependencies will be found automatically.
Here is an example :

```cmake
find_package(scalfmm CONFIG REQUIRED)
if(scalfmm_FOUND)
  message(STATUS "ScalFMM Found")
  add_executable(my_exe program.cpp )
  target_link_libraries(my_exe scalfmm::scalfmm)
else()
  message(FATAL_ERROR "ScalFMM NOT FOUND")
endif()
```

## Documentation
The doc can be found [here](https://solverstack.gitlabpages.inria.fr/ScalFMM/) or you can build it locally.

```bash
cd scalfmm/Build
cmake .. -DSCALFMM_BUILD_DOC=ON # or if cmake has already been called, ccmake .
make doc
```

This will generate the documentation in HTML format in the `Build/Doc/html` folder.

```bash
# From the Build folder
cd Doc/html
firefox index.html
```
## Contributing and development guidelines

### Gitlab flow

Please, read the Gitlab flow article available [here](https://docs.gitlab.com/ee/workflow/gitlab_flow.html).

To make it simple, if you want to contribute to the library, create a branch from `master` with a meaningful name and develop
your feature in that branch. Keep your branch up to date by regularly rebasing your branch from the `master` branch to be up
to date. Once your are done, send a merge request.

## Help and News

You can subscribe to the scalfmm-public-users@lists.gforge.inria.fr mailing list (http://lists.gforge.inria.fr/cgi-bin/mailman/listinfo/scalfmm-public-users). The list is very low trafic (~ 2 mails per year), we will let you know of improvements and releases.

Contact the developers at : scalfmm-public-support@lists.gforge.inria.fr

## Folder structure
  - include : library core.
  - Data : particle distribution examples.
  - Examples : common usage examples.
  - Doc : documentation configuration.
  - UTests : unit tests.
  - Tests : examples to know how to use scalfmm/put particles in the tree/iterate on the tree...
  - Utils : some scripts and binaries to handle data files.
