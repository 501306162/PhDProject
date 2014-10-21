set( SRCS 
	"main.cpp"
   	"main_window.cpp"
	"data_loader.cpp"
	"image_list.cpp"
	"image_viewer.cpp"
	"point_picker.cpp"
	"line.cpp"
	"containers.cpp"
	)

set(qt_headers "main_window.h")
include(${VTK_USE_FILE})

qt4_wrap_cpp( headers_moc ${qt_headers} )


add_executable( PlaneFinder ${SRCS} ${headers_moc} )
target_link_libraries( PlaneFinder ${ITK_LIBRARIES} ${QT_LIBRARIES} )
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/PlaneFinder" DESTINATION ${bin_path})
