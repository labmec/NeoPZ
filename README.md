# NeoPZ

NeoPZ is an open-source, cross-platform, finite element library
written in C++ with support for many advanced finite element techniques.

[![Run Unit Tests](https://github.com/labmec/neopz/actions/workflows/run_unittests.yml/badge.svg)](https://github.com/labmec/neopz/actions/workflows/run_unittests.yml)
[![Publish Docs](https://github.com/labmec/neopz/actions/workflows/publish_docs.yml/badge.svg)](https://github.com/labmec/neopz/actions/workflows/publish_docs.yml)
[![Compile External Projects](https://github.com/labmec/neopz/actions/workflows/compile_externalprojects.yml/badge.svg)](https://github.com/labmec/neopz/actions/workflows/compile_externalprojects.yml)

The documentation is available [here](https://labmec.github.io/neopz).


## Features
- Discontinuous, *H1, H(div)* and *H(curl)* conforming approximation spaces
- Multiphysics support
- *hp*-adaptivity
- Non-uniform meshes with support for hanging nodes
- Non-linear geometrical mappings (curved elements with exact representation)
- Run-time defined geometrical refinement patterns
- Forward automatic differentiation

And much more!

## Requirements
- A C++ 17 compiler
- [CMake](https://cmake.org/download/) 3.13.0+

### Optional external libraries
- [Paraview](https://www.paraview.org/download/) is recommended for display of VTK output files.

The usage of NeoPZ can be improved by linking against the following libraries:

- [log4cxx](https://logging.apache.org/log4cxx/latest_stable/), for logging information for debugging purposes.
- [Intel MKL](https://software.intel.com/en-us/mkl), for enabling sparse matrices solvers (in-house algorithms for skyline matrices are available, among other matrix storage formats).
- [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for (experimental) support of the [BDDC  technique](https://epubs.siam.org/doi/abs/10.1137/S1064827502412887?journalCode=sjoce3).
- [TBB](https://github.com/oneapi-src/oneTBB), also used for experimental techniques such as BDDC
- [blaze](https://bitbucket.org/blaze-lib/blaze), needed for projects using [SBFEM](https://www.cies.unsw.edu.au/scaled-boundary-finite-element-method-2a)
- [LAPACK](http://www.netlib.org/lapack/), for eigenvalues computations in full or banded matrices. If enabled, it is also internally used replacing in-house linear algebra algorithms with BLAS functions.
- [Boost](https://www.boost.org/), for experimental techniques.
- [Catch2](https://www.github.com/catchorg/Catch2), for Unit Testing.
- [Eigen](https://gitlab.com/libeigen/eigen), for integrating with Apple Accelerate on ARM-based CPUs and providing an alternative sparse solver (limited functionality, still in testing)

## Configuration and Install
The NeoPZ library uses CMake for configuring and installing the library. As a simple example, on UNIX systems, this could be done as:
```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DUSING_LOG4CXX=ON -DUSING_MKL=ON ..
make -j && sudo make install
```
where some CMake options were given as an illustration. The following CMake options are available:

- `REAL_TYPE`: which floating point type to use for real data types (geometry). It can be set as `float`, `double`, `long double`, `pzfpcounter` (last option is for internal usage). Default value: `double`.
- `STATE_TYPE`: which floating point type to use for the Finite Element matrices. It can be set as `float`, `double`, `long double`. Complex matrices will be created using `std::complex<STATE>`, which is aliased to `CSTATE`. Default value: `double`.
- `CMAKE_INSTALL_PREFIX`: where to install the NeoPZ library. Defaulted to `/opt/neopz` on UNIX systems.
- `USING_BLAZE`: Enable blaze library support
- `USING_BOOST`: Enable Boost libraries (`Boost::date_time`, `Boost::graph` and `Boost::unit_test_framework`) support
- `USING_LAPACK`: Enable LAPACK library support
- `USING_MKL`: Enable MKL library support. It will also enable usage of LAPACK functions.
- `USING_LOG4CXX`: Enable logging information through log4cxx support.
- `USING_METIS`: Enable Metis library support.
- `CMAKE_BUILD_TYPE`: `Release`, `Debug` and `RelWithDebugInfo`
- `BUILD_UNITTESTING`: this option will include the Unit Test suite. It automatically downloads `Catch2` if needed.
- `BUILD_PROJECTS`: this option will include the examples from [NeoPZExamples](https://github.com/labmec/NeoPZExamples) in the build tree.
- `BUILD_PLASTICITY_MATERIALS`: to include support for classes modelling plastic materials.
- `BUILD_DOCS`: build documentation using [Doxygen](https://www.doxygen.nl/index.html). It requires Doxygen and [Graphviz](https://graphviz.org/).
- `BUILD_SPHINX_DOCS`: use [Breathe](https://breathe.readthedocs.io/) for generating [Sphinx](https://www.sphinx-doc.org/en/master/) from the Doxygen output.

*Note:* MKL has recently changed how they distribute their libraries, using the [oneAPI](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html) programming model. See [here](https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-base-linux/top/before-you-begin.html) how to set up oneAPI MKL for usage with NeoPZ. In UNIX systems, adding the 
```sh
source /opt/intel/oneapi/setvars.sh intel64
```
command to your shell startup file (`~/.bashrc_profile`, `~/.zshenv`, *etc*) should suffice.

For older MKL distributions, the script is located in:
```sh
source MKL_DIR/bin/compilervars.sh intel64
```

### Internal options/deprecated options
The remaining `USING_XXX` CMake options are for internal usage only.
The following options are listed for completeness:

- `BUILD_PERF_TESTS`: for building a few performance tests. The tests are in need of a revision.
- `BUILD_PUBLICATIONS`: build legacy projects with code used to generate results that were published.


Should any problems arise during the installation of NeoPZ library, please contact us on <neopz@googlegroups.com>.

## Using NeoPZ in your C++ project
The installation of NeoPZ will provide the `add_pz_target` function to be used in the CMake files of your project. A minimal example of the `CMakeLists.txt` files for a project using the NeoPZ library follows:

`proj_root/CMakelists.txt`

```cmake
cmake_minimum_required(VERSION 3.13)

project (MySuperProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding the pz package
find_package(NeoPZ REQUIRED)
# Finding other packages
# Let us assume that the following call will define
# OtherDep::OtherDep
# OTHERDEP_INCLUDE_DIRS
find_package(OtherDep REQUIRED)
add_subdirectory(MySubDirectory)
```

`proj_root/MySubDirectory/CMakeLists.txt`

```cmake
add_pz_target(
  NAME MyTarget
  SOURCES MySources.cpp
  FILES file_to_be_copied_to_binary_dir.abc
  REQUIRED PZ_USING_MKL)
  
  target_link_libraries(MyTarget PRIVATE OtherDep::OtherDep)
  target_include_directories(MyTarget PRIVATE ${OTHERDEP_INCLUDE_DIRS})
```

*Note:* If you don't want to use the function `add_pz_target()`, the PZ library target is available under the alias `NeoPZ::pz`. It can still be linked through CMake's `target_link_libraries(my_target <PRIVATE|PUBLIC> NeoPZ::pz)`

The CMake options defined when configuring NeoPZ will be available to the external targets with the `PZ` prefix (`USING_MKL` becomes `PZ_USING_MKL`). So, if a certain option is needed in an external project, the `REQUIRED` field in `add_pz_targets` can be used. *Note:* `PZ_LOG` is used both internally and externally for identifying if there is support for logging.

*Note:* NeoPZ will be installed in UNIX systems in `/opt/neopz` by default. If `/opt` or `/opt/bin` is not on your `PATH`, the CMake variable `NEOPZ_DIR` can be used in your project using NeoPZ:

```sh
mkdir build && cd build
cmake -DNEOPZ_DIR=pz_install_dir ..
```

or, if you prefer, the following line

```sh
export PATH=pz_install_dir:$PATH
```

can be added to a startup file of your shell. In both examples, `pz_install_dir` is `/opt/neopz`.

## NeoPZ documentation

A Doxygen documentation can be found 
[here](http://www.labmec.org.br/pz/arquivos-html/html/index.html).

## How to cite NeoPZ

Devloo, P.R. B., 1997. PZ: An object oriented environment
for scientific programming. Computer Methods in Applied Mechanics
and Engineering 150, 133–153.
https://doi.org/10.1016/s0045-7825(97)00097-2

```bibtex
@article{neopz97,
author = {Devloo, Philippe Remy Bernard},
year = {1997},
month = {12},
pages = {},
title = {PZ: An object oriented environment for scientific programming},
volume = {150},
booktitle = {Computer Methods in Applied Mechanics and Engineering},
doi = {https://doi.org/10.1016/s0045-7825(97)00097-2}
}
```
