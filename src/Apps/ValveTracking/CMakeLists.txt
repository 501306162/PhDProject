set(APP_NAMES 
	"ParseResults"
	"ValveAlignment"
	#"FeatureTest"
	"ROITest"
	"SegmentROIs"
	"ExtractPatch"
	"TestFitting"
	"ValveROIFinder"
	"ExtractTrainingFeatures"
	"SortOutDicom"
	"ValveDataAppender"
	"PlaneCalculator"
	"ComputeAlignedPlanes"
	#"ModelTest"
	"ValveNameReplacer"
	"Testing"
	"SAStackBuilder")

set(DataSorter_SRC "DataSorter.cpp")
set(DataSorter_HDRS "main_window.h")



include_directories(${ITK_INCLUDE_DIRS})

find_package("qjson")
include(${VTK_USE_FILE})
include_directories(${qjson_INCLUDE_DIRS})




foreach(APP ${APP_NAMES})
	add_executable(${APP} "${APP}.cpp")
	target_link_libraries(${APP} ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
	install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${APP}" DESTINATION "${bin_path}/ValveTracking")
endforeach()



qt4_wrap_cpp(headers_moc "main_window.h")
add_executable(DataSorter "DataSorter.cpp" "main_window.cpp" ${headers_moc})
target_link_libraries(DataSorter ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/DataSorter" DESTINATION "${bin_path}/ValveTracking")
