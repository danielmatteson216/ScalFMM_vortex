#!/bin/bash

# Performs an analysis of ScalFMM source code
# We consider:
# 1) to be in ScalFMM's source code root
# 2) ScalFMM built into the ./Build directory
# 3) ctests with coverage results should have been performed

# capture coverage
lcov --directory Build --capture --output-file scalfmm.lcov
lcov_cobertura.py scalfmm.lcov --output scalfmm-coverage.xml
# to get it displayed and captured by gitlab to expose the badge on the main page
gcovr -r $PWD

## filter sources:
export SOURCES_TO_ANALYZE="include src Examples Tests UTests"

# run cppcheck analysis
# commented for now because too long
export CPPCHECK_INCLUDES="-IBuild/include -Iinclude"
export SOURCES_TO_ANALYZE="src Examples Tests UTests"
export CPPCHECK_UNDEF="-USCALFMM_USE_AVX -USCALFMM_USE_AVX2 -USCALFMM_USE_SSE -USCALFMM_USE_EZTRACE -USCALFMM_BLAS_NOCHANGE -USCALFMM_BLAS_UPCASE -U__AVXPE_INTEL_COMPILER -U__INTEL_COMPILER -U__MIC__ -U__SSEPE_INTEL_COMPILER -U__SSSE3__ -U__SSSE4_1__ -U_WIN32 -Uming -U_MSC_VER -U__IBMCPP__ -U__PGI -U__SUNPRO_CC -U__clang__ -U__APPLE__ -U__HP_aCC -U__HP_cc -U__ICC"
cppcheck -v -f --language=c++ --std=c++14 --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem ${CPPCHECK_UNDEF} ${CPPCHECK_INCLUDES} ${SOURCES_TO_ANALYZE} 2> scalfmm-cppcheck.xml

# run rats analysis
rats -w 3 --xml ${SOURCES_TO_ANALYZE} > scalfmm-rats.xml

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.bordeaux.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN
sonar.links.homepage=https://gitlab.inria.fr/solverstack/ScalFMM
sonar.links.scm=https://fpruvost@gitlab.inria.fr/solverstack/ScalFMM.git
sonar.links.issue=https://gitlab.inria.fr/solverstack/ScalFMM/issues
sonar.projectKey=hiepacs:scalfmm:gitlab:develop
sonar.projectDescription=C++ library that implements a kernel independent Fast Multipole Method (LGPL+CeCILL-C)
sonar.projectVersion=2.0
sonar.language=c++
sonar.sourceEncoding=UTF-8
sonar.sources=Examples, src, Tests, UTests
sonar.exclusions=include/Kernels/FKernelConcepts.hpp, Tests/noDist/PerfTest/TestDriver.hpp
sonar.cxx.includeDirectories=$(echo | gcc -E -Wp,-v - 2>&1 | grep "^ " | tr '\n' ',')Build/include, include
sonar.cxx.compiler.charset=UTF-8
sonar.cxx.compiler.regex=^(.*):(\\\d+):\\\d+: warning: (.*)\\\[(.*)\\\]$
sonar.cxx.compiler.reportPath=scalfmm-build.log
sonar.cxx.coverage.reportPath=scalfmm-coverage.xml
sonar.cxx.rats.reportPath=scalfmm-rats.xml
sonar.cxx.cppcheck.reportPath=scalfmm-cppcheck.xml
#sonar.cxx.clangsa.reportPath=build/analyzer_reports/*/*.plist
EOF
