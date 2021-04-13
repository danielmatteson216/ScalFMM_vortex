# Examples
# --------

# List of source files
set(source_tests_files
  changeFmaFormat.cpp
  ChebyshevHybridFMM.cpp
  ChebyshevOpenMPAdaptiveFMM.cpp
  ChebyshevOpenMPFMM.cpp
  ChebyshevStarpuImplicit.cpp
  compare2Files.cpp
  compareAllPoissonKernels.cpp
  CutOffAlgorithm.cpp
  DirectComputation.cpp
  generateDistributions.cpp
  LagrangeHybridFMM.cpp
  LagrangeInterpolationAdaptiveFMM.cpp
  LagrangeOpenMPFMM.cpp
  LagrangeStarpuImplicit.cpp
  RotationFMM.cpp
  RotationMPIFMM.cpp
  statisticsOnOctree.cpp
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
    list(APPEND SCALFMM_EXAMPLES_TARGETS ${execname})
    list(APPEND SCALFMM_TESTS_TARGETS ${execname})
    set_target_properties(${execname} PROPERTIES ENABLE_EXPORTS TRUE
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BUILD_TYPE}
      )
    target_link_libraries( ${execname} PRIVATE ${CMAKE_PROJECT_NAME})
    install(TARGETS ${execname} RUNTIME DESTINATION bin)
  endif()
endforeach(exec)

add_custom_target(scalfmm_examples ALL DEPENDS ${SCALFMM_EXAMPLES_TARGETS})


