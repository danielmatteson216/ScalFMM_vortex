#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "scalfmm::scalfmm" for configuration "RELEASE"
set_property(TARGET scalfmm::scalfmm APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(scalfmm::scalfmm PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libscalfmm.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS scalfmm::scalfmm )
list(APPEND _IMPORT_CHECK_FILES_FOR_scalfmm::scalfmm "${_IMPORT_PREFIX}/lib/libscalfmm.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
