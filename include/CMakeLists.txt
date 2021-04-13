# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

# Configuration du projet
project(Lib_scalfmm)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

ADD_DEFINITIONS(${ScaLFMM_CXX_FLAGS})

# Searching all cpp file
file(
	GLOB_RECURSE
	source_lib_files
	./*.cpp
)

set(source_lib_files_enabled)
foreach(exec ${source_lib_files}) 
	set(compile_exec "TRUE")
	
	foreach(fuse_key ${FUSE_LIST})
	    file(STRINGS "${exec}" lines_fuse REGEX "@FUSE_${fuse_key}")
	    if(lines_fuse)
		    if( NOT SCALFMM_USE_${fuse_key} )
			    MESSAGE( STATUS "This needs ${fuse_key} = ${exec}" )
			    set(compile_exec "FALSE")
		    endif()
	    endif()
	endforeach()

	# Dependency are OK
	if( compile_exec )
        list(APPEND source_lib_files_enabled ${exec})
	endif()
endforeach(exec)


# Add CUDA files once they are compiled from cu to .o
if(SCALFMM_USE_CUDA)
    # Find all the CU files in my project
    file(GLOB_RECURSE source_cu_files ./*.cu)
 
    # Iterate and add builind command for each file
    set( SCALFMM_CUDA_SOURCES "" )
    FOREACH (_file ${source_cu_files})
 
        GET_FILENAME_COMPONENT (_filewe ${_file} NAME_WE)
        SET (_filehpp_output ${CMAKE_CURRENT_BINARY_DIR}/${_filewe}.o)
 
        ADD_CUSTOM_COMMAND(OUTPUT  ${_filehpp_output}
                  DEPENDS ${_file}
                  COMMAND echo ARGS Compiling ${_filewe}
                  COMMAND nvcc ARGS -I${CUDA_NEEDED_INCLUDE_DIRS} -c ${_file} -o ${_filehpp_output} ${CUSTOM_CUDA_FLAGS})
 
        LIST (APPEND SCALFMM_CUDA_SOURCES ${_filehpp_output})
    ENDFOREACH ()
 
    MESSAGE( STATUS "SCALFMM_CUDA_SOURCES          = ${SCALFMM_CUDA_SOURCES}" )
endif()

# Adding cpp files to project
add_library(
	scalfmm
	${source_lib_files_enabled}
	${SCALFMM_CUDA_SOURCES}
)

# Add blas library (even if it is set to off)
target_link_libraries(
	scalfmm
	${scalfmm_lib}
	${SCALFMM_LIBRARIES}
)

# Adding the entire project dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_CURRENT_SOURCE_DIR}
     ${CMAKE_CURRENT_BINARY_DIR}
     ${SCALFMM_INCLUDES}
)

#########################################################################
#
# Install lib
install( TARGETS scalfmm DESTINATION lib )
#
# Install headers
SET(my_include_dirs "."  "Adaptive"
  "Arranger"       "Components"  "Containers"  "Core"    "Extensions"
  "Files" "Utils" )
INSTALL( DIRECTORY "Kernels" DESTINATION include FILES_MATCHING PATTERN   "*.hpp")
INSTALL( DIRECTORY "GroupTree" DESTINATION include )
#    FILES_MATCHING PATTERN "*.hpp"  PATTERN "*.h" PATTERN "*.cu" PATTERN "*.cu")

FOREACH(my_dir ${my_include_dirs})	
   INSTALL( DIRECTORY ${my_dir} DESTINATION include
           FILES_MATCHING PATTERN   "*.hpp"   PATTERN "*.h")
ENDFOREACH()
INSTALL( FILES "${SCALFMM_BINARY_DIR}/Src/ScalFmmConfig.h" DESTINATION include )
#
#########################################################################
#

