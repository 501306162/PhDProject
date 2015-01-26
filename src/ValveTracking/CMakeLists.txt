set(lib_name "ValveTracking")


set(lib_srcs
	"JsonConverter.cpp"
	"ValveLine.cpp"
	"CMRFileExtractor.cpp"
	"PointData.cpp"
	"lbp.cpp"
	"PatchTrainingData.cpp"
	"lbp_hist.cpp"
	"FlipChecker.cpp"
	"BoundingBox.cpp"
	"LineIntersectionFinder.cpp"
	"ValveIO.cpp"
	"ResultViewer.cpp"
	"PatchExtractor.cpp"
	"ComputeProbImage.cpp"
	"FeatureExtractor.cpp"
	"LBPFeatureExtractor.cpp"
	"ValvePointCandidateFinder.cpp"
	"BinaryPatchFeatureExtractor.cpp"
	"ValvePlane.cpp"
	"LengthData.cpp"
	"TestData.cpp"
	"TrainingData.cpp"
	"CandidateLineFinder.cpp"
	"ValveNormaliser.cpp"
	"ParameterHelpers.cpp"
	"ValveHelpers.cpp"
	"SVMClassifier.cpp"
	"ValveAligner.cpp"
	"SimpleMRFSegmenter.cpp"
	"LBPFeatureExtractor.cpp"
	"OpenCVValve.cpp"
	"ValveOriginFinder.cpp"
	"JsonReader.cpp")


set(lib_headers
	"JsonConverter.h"
	"PointData.h"
	"ComputeProbImage.h"
	"ResultViewer.h"
	"PatchExtractor2.h"
	"PatchExtractor2.hpp"
	"LengthData.h"
	"CandidateLineFinder.h"
	"ValvePointCandidateFinder.h"
	"SVMClassifier.h"
	"TestData.h"
	"lbp.h"
	"PatchTrainingData.h"
	"ParameterHelpers.h"
	"lbp_hist.h"
	"PatchExtractor.h"
	"LBPFeatureExtractor.h"
	"LineIntersectionFinder.h"
	"BinaryPatchFeatureExtractor.h"
	"LBPFeatureExtractor.h"
	"TrainingData.h"
	"ValvePlane.h"
	"FeatureExtractor.h"
	"FlipChecker.h"
	"ValveNormaliser.h"
	"CMRFileExtractor.h"
	"BoundingBox.h"
	"ValveOriginFinder.h"
	"ValveAligner.h"
	"SimpleMRFSegmenter.h"
	"OpenCVValve.h"
	"ValveLine.h"
	"ValveHelpers.h"
	"ValveIO.h"
	"JsonReader.h")




set(lib_include_dir "${include_path}/${lib_name}")


include_directories(${OpenCV_DIRS})

# ----------------------------------------------------------------------
# Compile the library and install it
add_library(${lib_name} SHARED ${lib_srcs} ${lib_headers})
target_link_libraries(${lib_name} ${ITK_LIBRARIES} ${Qt4_LIBRARIES} ${OpenCV_LIBS} hdf5_hl MRFLib SVM qjson Utils)

# ----------------------------------------------------------------------
# install command for the library
install(TARGETS ${lib_name} 
	RUNTIME DESTINATION ${bin_path}
	LIBRARY DESTINATION ${lib_path}
	ARCHIVE DESTINATION ${lib_path}) 



# ----------------------------------------------------------------------
# install the headers
install(FILES ${lib_headers} DESTINATION ${lib_include_dir})
