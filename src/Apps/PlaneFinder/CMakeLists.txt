set( SRCS 
	"main.cpp"
   	"main_window.cpp"
	"main_window_setup.cpp"
	"data_loader.cpp"
	"image_list.cpp"
	"image_viewer.cpp"
	"point_picker.cpp"
	"line.cpp"
	"line_list.cpp"
	"containers.cpp"
	"key_controls.cpp"
	"io.cpp"
	"patch_extractor.cpp"
	"dicom_parser.cpp"
	"image_exporter.cpp"
	)


find_package("qjson")
include_directories(${qjson_INCLUDE_DIRS})

set(qt_headers "main_window.h" "key_controls.h" "dicom_parser.h")
include(${VTK_USE_FILE})

qt4_wrap_cpp( headers_moc ${qt_headers} )


add_executable( PlaneFinder ${SRCS} ${headers_moc} )
target_link_libraries( PlaneFinder ${ITK_LIBRARIES} ${QT_LIBRARIES} Utils qjson )
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/PlaneFinder" DESTINATION ${bin_path})
