set(APP_NAMES 
	"ParseResults"
	"ValveAlignment"
	"PatchTest"
	#"FeatureTest"
	"ROITest"
	"TMCheck"
	#"CreatePointLocations"
	"TrainClassifiers"
	#"HaarPatches"
	"ComputeNewPlanes"
	"SegmentROIs"
	"ExtractPatch"
	#"TestFitting"
	"ValveROIFinder"
	"SortOutDicom"
	"ValveDataAppender"
	"FlipPoints"
	"PlaneCalculator"
	"CreatePatches"
	"ComputeAlignedPlanes"
	"ComputeLineLengths"
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

set(Extractor_SRCS "ExtractTrainingFeatures.cpp" "ExtractTrainingFeaturesFunctions.cpp")

add_executable(ExtractTrainingFeatures ${Extractor_SRCS})
target_link_libraries(ExtractTrainingFeatures ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/ExtractTrainingFeatures" DESTINATION "${bin_path}/ValveTracking")


set(BigProgs 
	"Fitting2"
	"SMOTEer"
	"BadTracking"
	"MeanShiftTest"
	"PF")


foreach(Prog ${BigProgs})
	add_executable(${Prog} "${Prog}.cpp"  "PFFunctions.cpp" "FittingFunctions.cpp" "ExtractTrainingFeaturesFunctions.cpp")
	target_link_libraries(${Prog} ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
	install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${Prog}" DESTINATION "${bin_path}/ValveTracking")
endforeach()

#set(Fitting_SRCS "Fitting.cpp" "FittingFunctions.cpp" "ExtractTrainingFeaturesFunctions.cpp")
#add_executable(Fitting ${Fitting_SRCS})
#target_link_libraries(Fitting ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
#install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/Fitting" DESTINATION "${bin_path}/ValveTracking")

#set(Fitting2_SRCS "Fitting2.cpp" )
#add_executable(Fitting2 ${Fitting2_SRCS})


#set(PF_SRCS "PF.cpp" "PFFunctions.cpp" "FittingFunctions.cpp" "ExtractTrainingFeaturesFunctions.cpp")
#add_executable(PF ${PF_SRCS})
#target_link_libraries(PF ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
#install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/PF" DESTINATION "${bin_path}/ValveTracking")


#set(MeanShiftTest_SRCS "MeanShiftTest.cpp" "FittingFunctions.cpp" "ExtractTrainingFeaturesFunctions.cpp")
#add_executable(MeanShiftTest ${MeanShiftTest_SRCS})
#target_link_libraries(MeanShiftTest ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
#install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/MeanShiftTest" DESTINATION "${bin_path}/ValveTracking")






qt4_wrap_cpp(headers_moc "main_window.h")
add_executable(DataSorter "DataSorter.cpp" "main_window.cpp" ${headers_moc})
target_link_libraries(DataSorter ${ITK_LIBRARIES} ${QT_LIBRARIES} qjson nlopt MRFLib Utils ValveTracking Manifold)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/DataSorter" DESTINATION "${bin_path}/ValveTracking")
