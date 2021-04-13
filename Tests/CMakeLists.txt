# Find all code files
set(source_tests_files
  GroupTree/compareDAGmapping.cpp 
  GroupTree/testBlockedAlgorithm.cpp
  GroupTree/testBlockedAlgorithmDyn.cpp
  GroupTree/testBlockedChebyshev.cpp 
  GroupTree/testBlockedMpiAlgorithm.cpp
  GroupTree/testBlockedMpiChebyshev.cpp
  GroupTree/testBlockedMpiInterpolation.cpp
  GroupTree/testBlockedMpiUniform.cpp
  GroupTree/testBlockedRotationCompare.cpp
  GroupTree/testBlockedRotation.cpp
  GroupTree/testBlockedRotationCuda.cpp
  GroupTree/testBlockedRotationMpi.cpp
  GroupTree/testBlockedRotationMpiCuda.cpp
  GroupTree/testBlockedSpeedOperators.cpp
  GroupTree/testBlockedTaylor.cpp
  GroupTree/testBlockedTree.cpp
  GroupTree/testBlockedUnifCudaBench.cpp 
  GroupTree/testBlockedUnifCuda.cpp
  GroupTree/testBlockedUniformBench.cpp
  GroupTree/testBlockedUniformCompare.cpp
  GroupTree/testBlockedUniform.cpp
  GroupTree/testBlockedWithCudaAlgorithm.cpp
  GroupTree/testBlockedWithOpenCLAlgorithm.cpp
  Kernels/DirectAlgorithm.cpp
  Kernels/testChebAlgorithm.cpp
  Kernels/testChebAlgorithmProc.cpp
  Kernels/testChebTensorialAlgorithm.cpp
  Kernels/testFlopsChebAlgorithm.cpp
  Kernels/testOmniPath.cpp
  Kernels/testP2PEfficency.cpp
  Kernels/testRotationAlgorithm.cpp
  Kernels/testRotationAlgorithmProc.cpp
  Kernels/testRotationPeriodicBench.cpp
  Kernels/testSmoothUnifAlgorithm.cpp
  Kernels/testSphericalAlgorithm.cpp
  Kernels/testSphericalBench.cpp
  Kernels/testSphericalBlasAlgorithm.cpp
  Kernels/testSphericalBlockBlasAlgorithm.cpp
  Kernels/testSphericalDlpolyAlgorithm.cpp
  Kernels/testSphericalGalaxyCsv.cpp
  Kernels/testSphericalProcAlgorithm.cpp
  Kernels/testSphericalRotationAlgorithm.cpp
  Kernels/testTaylor.cpp
  Kernels/testTaylorKernel.cpp
  Kernels/testTaylorSimple.cpp
  Kernels/testTsmAlgorithm.cpp
  Kernels/testTuneSphericalBlockBlas.cpp
  Kernels/testUnifAlgorithm.cpp
  Kernels/testUnifAlgorithmProc.cpp
  Kernels/testUnifTensorialAlgorithm.cpp
  Utils/testBlas.cpp
  Utils/testChebBinaryM2L.cpp
  Utils/testChebBinarySymM2L.cpp
  Utils/testChebInterpolator.cpp
  Utils/testChebInterpolatorExtendedBBox.cpp
  Utils/testChebM2Lprecomputation.cpp
  Utils/testChebOctree.cpp
  Utils/testChebSxUCBSy.cpp
  Utils/testChebSymmetries.cpp
  Utils/testChebTensorProduct.cpp
  Utils/testCompareIOTree.cpp
  Utils/testDebug.cpp
  Utils/testFastDiscreteConvolution.cpp
  Utils/testFFT.cpp
  Utils/testFmmAlgorithm.cpp
  Utils/testFmmAlgorithmOmp4.cpp
  Utils/testFmmAlgorithmPeriodic.cpp
  Utils/testFmmAlgorithmProc.cpp
  Utils/testFmmAlgorithmProcPeriodic.cpp
  Utils/testFmmAlgorithmProcRandom.cpp
  Utils/testFmmAlgorithmProcRotation.cpp
  Utils/testFmmAlgorithmProcTsm.cpp
  Utils/testFmmAlgorithmTsm.cpp
  Utils/testLapack.cpp
  Utils/testLapackQR.cpp
  Utils/testLoader.cpp
  Utils/testLoaderFMATsm.cpp
  Utils/testLoaderFMATsmCreate.cpp
  Utils/testMemStats.cpp
  Utils/testMpiCsvSaver.cpp
  Utils/testOctree.cpp
  Utils/testOctreeFuncteur.cpp
  Utils/testOctreeIter.cpp
  Utils/testOctreeParallelFuncteur.cpp
  Utils/testOctreePrintMorton.cpp
  Utils/testOctreeRearrange.cpp
  Utils/testOctreeRearrangePeriodic.cpp
  Utils/testOctreeRearrangeProc.cpp
  Utils/testOctreeRearrangeTsm.cpp
  Utils/testParameterNames.cpp
  Utils/testPartitionsMapping.cpp
  Utils/testStatsTree.cpp
  Utils/testTemplateExample.cpp
  Utils/testTic.cpp
  Utils/testTreeBuilderThread.cpp
  Utils/testTreeIO.cpp
  Utils/testUnifInterpolator.cpp
  Utils/testUnifTensorialInterpolator.cpp    
  )

# Add execs - 1 cpp = 1 exec
foreach(exec ${source_tests_files}) 
  set(compile_exec TRUE)
  get_filename_component( execname ${exec} NAME_WE )
	
  foreach(fuse_key ${FUSE_DEP_AVAILABLE})
    file(STRINGS "${exec}" lines_fuse REGEX "@FUSE_${fuse_key}")
    if(lines_fuse AND NOT ${fuse_key} IN_LIST FUSE_LIST)
      message( STATUS "This needs ${fuse_key} = ${exec}" )
      set(compile_exec FALSE)
    endif()
  endforeach()

  # Dependency are OK
  if( compile_exec )
    add_executable( ${execname} ${exec})
    list(APPEND SCALFMM_TESTS_TARGETS ${execname})
    set_target_properties(${execname} PROPERTIES ENABLE_EXPORTS TRUE
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BUILD_TYPE}
      )
    target_link_libraries( ${execname} ${CMAKE_PROJECT_NAME})
  endif()
endforeach(exec)

add_custom_target(scalfmm_tests ALL DEPENDS ${SCALFMM_TESTS_TARGETS})

