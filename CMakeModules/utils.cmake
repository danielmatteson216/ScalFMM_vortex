# helper function for generating targets

# This function creates an executable target with a list
# of sources and all requirements passed as arguments
# if the expression passed as boolean is true.
# ----------------------------------------------------
function(scalfmm_create_target target_name sources_list dependencies options )
  # Here, we generate a vector of boolean to check if dependencies (targets)
  # exists. 1 is yes, 0 otherwise
  foreach(dep IN LISTS dependencies)
    if(TARGET ${dep})
      list(APPEND VECTOR_BOOL 1)
    else()
      list(APPEND VECTOR_BOOL 0)
    endif()
  endforeach()
  # If all dependencies are met, we generate the target, otherwise not.
  if(NOT 0 IN_LIST VECTOR_BOOL)
    add_executable(${target_name} ${sources_list})
    target_link_libraries(${target_name} PRIVATE ${dependencies})
    target_compile_options(${target_name} PRIVATE ${options})
    set_target_properties(${target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BUILD_TYPE})
    install(TARGETS ${target_name} RUNTIME DESTINATION bin)
  endif()
endfunction(scalfmm_create_target)

