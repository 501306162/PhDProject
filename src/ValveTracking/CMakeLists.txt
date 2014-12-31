set(lib_name "ValveTracking")


set(lib_srcs
	"JsonConverter.cpp"
	"ValveLine.cpp"
	"CMRFileExtractor.cpp"
	"BoundingBox.cpp"
	"ValveIO.cpp"
	"ValveHelpers.cpp"
	"ValveAligner.cpp"
	"OpenCVValve.cpp"
	"ValveOriginFinder.cpp"
	"JsonReader.cpp")


set(lib_headers
	"JsonConverter.h"
	"CMRFileExtractor.h"
	"BoundingBox.h"
	"ValveOriginFinder.h"
	"ValveAligner.h"
	"OpenCVValve.h"
	"ValveLine.h"
	"ValveHelpers.h"
	"ValveIO.h"
	"JsonReader.h")




set(lib_include_dir "${include_path}/${lib_name}")

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_DIRS})

# ----------------------------------------------------------------------
# Compile the library and install it
add_library(${lib_name} SHARED ${lib_srcs} ${lib_headers})
target_link_libraries(${lib_name} ${ITK_LIBRARIES} ${Qt4_LIBRARIES} hdf5_hl ${OpenCV_LIBS} qjson Utils)

# ----------------------------------------------------------------------
# install command for the library
install(TARGETS ${lib_name} 
	RUNTIME DESTINATION ${bin_path}
	LIBRARY DESTINATION ${lib_path}
	ARCHIVE DESTINATION ${lib_path}) 



# ----------------------------------------------------------------------
# install the headers
install(FILES ${lib_headers} DESTINATION ${lib_include_dir})
