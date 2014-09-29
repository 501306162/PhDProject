# ----------------------------------------------------------------------
# CMAKE file for the manifold lib
# ----------------------------------------------------------------------
set(lib_name "Filters")


# ----------------------------------------------------------------------
# Get the sources for the library
set(lib_srcs
	"")


# ----------------------------------------------------------------------
# Get the headers for the install include folder
set(lib_headers 
	"ImageToDoGImageFilter.h"
	"ImageToDoGImageFilter.hpp"
	"DoGKeyPointExtractor.h"
	"DoGKeyPointExtractor.hpp"
	"KeyPointOrientator.h"
	"KeyPointOrientator.hpp"
	"FeaturePointGrouper.h"
	"FeaturePointGrouper.hpp"
	"ImagePatchExtractor.h"
	"ImagePatchExtractor.hpp"
	"FeatureClusterLearner.h"
	"FeatureClusterLearner.hpp"
	"FeatureCommon.h"
	"HistogramOfGradientsFeatureExtractor.h"
	"HistogramOfGradientsFeatureExtractor.hpp"
	"StructureTensorKeyPointValidator.h"
	"StructureTensorKeyPointValidator.hpp"
	"ClusterDistributionGenerator.h"
	"ClusterDistributionGenerator.hpp"
	)

set(lib_include_dir "${include_path}/${lib_name}")



# ----------------------------------------------------------------------
# Do the compilation, or not
if(NOT "${lib_srcs}" STREQUAL "")

	# ----------------------------------------------------------------------
	# Compile the library and install it
	add_library(${lib_name} SHARED ${lib_srcs} ${lib_headers})
	target_link_libraries(${lib_name} ${ITK_LIBRARIES} ${Qt4_LIBRARIES})

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


# ----------------------------------------------------------------------
# 
