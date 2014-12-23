# ----------------------------------------------------------------------
# CMAKE file for the Utils lib
# ----------------------------------------------------------------------
set(lib_name "Utils")


# ----------------------------------------------------------------------
# Get the sources for the library
set(lib_srcs 
	"FilenamesReader.cpp"
	"Directory.cpp"
   	"MatrixIO.cpp"
	"VonMisesDistribution.cpp"
	"CardiacDicomExtractor.cpp"
	"PatientData.cpp"
	"MatrixWriter.cpp"
	"MatrixReader.cpp"
	"DicomParser.cpp"
	"SeriesExtractor.cpp"
	"RegistrationHelpers.cpp")


# ----------------------------------------------------------------------
# Get the headers for the install include folder
set(lib_headers 
	"Helpers.h"
	"Directory.h"
	"Helpers.hpp"
	"ImageIO.h"
	"CardiacDicomExtractor.h"
	"VonMisesDistribution.h"
	"ImageIO.hpp"
	"FilenamesReader.h"
	"Knn.h"
	"Knn.hpp"
	"SeriesExtractor.h"
	"PatientData.h"
	"MatrixIO.h"
	"MatrixWriter.h"
	"DicomParser.h"
	"MatrixReader.h"
	"MatrixCommon.h"
	"RegistrationHelpers.h")

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



