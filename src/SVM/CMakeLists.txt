# ----------------------------------------------------------------------
# CMAKE file for the Utils lib
# ----------------------------------------------------------------------
set(lib_name "SVM")


# ----------------------------------------------------------------------
# Get the sources for the library
set(lib_srcs 
	"svm.cpp"
	)



# ----------------------------------------------------------------------
# Get the headers for the install include folder
set(lib_headers 
	"svm.h"
	)

set(lib_include_dir "${include_path}/${lib_name}")



# ----------------------------------------------------------------------
# Do the compilation, or not
if(NOT "${lib_srcs}" STREQUAL "")

	# ----------------------------------------------------------------------
	# Compile the library and install it
	add_library(${lib_name} SHARED ${lib_srcs} ${lib_headers})
	target_link_libraries(${lib_name} ${ITK_LIBRARIES} ${Qt4_LIBRARIES} hdf5_hl)

	# ----------------------------------------------------------------------
	# install command for the library
	install(TARGETS ${lib_name} 
		RUNTIME DESTINATION ${bin_path}
		LIBRARY DESTINATION ${lib_path}
		ARCHIVE DESTINATION ${lib_path}) 

endif()


# ----------------------------------------------------------------------
# No library to compile, just install the headers
install(FILES ${lib_headers} DESTINATION ${lib_include_dir})

