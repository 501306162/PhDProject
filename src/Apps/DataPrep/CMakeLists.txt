# ----------------------------------------------------------------------
# Set of programs for testing
set(APP_NAMES 
	"AlignShapes"
	"AlignImages"
	"ResampleLabelImages"
	"BuildManifoldShapeModel"
	"LearnGMMs"
	"ComputeDistanceMaps"
	"VolumePointPicker"
	"FeatureLearning"
	"NonRigidShapeAlignment"
	"ThresholdFunction"
	#"ContourExtraction"
	)

include_directories(${ITK_INCLUDE_DIRS})


foreach(APP ${APP_NAMES})
	add_executable(${APP} "${APP}.cpp")
	add_dependencies(${APP} Utils Manifold)
	target_link_libraries(${APP} ${ITK_LIBRARIES} nlopt MRFLib Utils Manifold)
	install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${APP}" DESTINATION ${bin_path})
endforeach()
