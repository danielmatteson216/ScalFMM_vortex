# Utesting

# Find all code files
set(source_tests_files
  utestAlgorithmBuilder.cpp
  utestAlignedMemory.cpp
  utestBoolArray.cpp
  utestBuffer.cpp
  utestChebyshev.cpp
  utestChebyshevDirectPeriodic.cpp
  utestChebyshevDirectTsm.cpp
  utestChebyshevMpi.cpp
  utestChebyshevThread.cpp
  utestFBasicParticleContainer.cpp
  utestFBasicParticle.cpp
  utestFmmAlgorithmProc.cpp
  utestInterpolationMultiRhs.cpp
  utestLagrange.cpp
  utestLagrangeMpi.cpp
  utestLagrangeThread.cpp
  utestLetTree.cpp
  utestListBlockAllocator.cpp
  utestList.cpp
  utestMorton.cpp
  utestMpiBitonic.cpp
  utestMPILoader.cpp
  utestMpiQs.cpp
  utestMpiTreeBuilder.cpp
  utestNeighborIndexes.cpp
  utestOctree.cpp
  utestP2PExclusion.cpp
  utestQuicksort.cpp
  utestRotation.cpp
  utestRotationDirectSeveralTime.cpp
  utestRotationDirectTsm.cpp
  utestRotationThread.cpp
  utestSphericalBlasAlgorithm.cpp
  utestSpherical.cpp
  utestSphericalDirect.cpp
  utestSphericalWithPrevious.cpp
  utestStaticMpiTreeBuilder.cpp
  utestTest.cpp
  utestVector.cpp
  Utils/variadic_vector/utest_variadic_vector.cpp
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

    target_link_libraries( ${execname} PRIVATE ${CMAKE_PROJECT_NAME} )
    list(APPEND SCALFMM_UTESTS_TARGETS ${execname})

    if(lines_mpi)
      add_test(NAME ${execname} COMMAND $<IF:${MPIEXEC_EXECUTABLE}, ${MPIEXEC_EXECUTABLE} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${execname}>)
    else()
      add_test(NAME ${execname} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${execname})
    endif()
  endif()
endforeach(exec)

add_custom_target(scalfmm_utests ALL DEPENDS ${SCALFMM_UTESTS_TARGETS})

