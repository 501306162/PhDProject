set( sources "main.cpp"
	"SeriesExtractor.cpp"
	"MainWindow.cpp" )

set(qt_headers "MainWindow.h")


qt4_wrap_cpp( headers_moc ${qt_headers} )


add_executable( DicomSorter ${sources} ${headers_moc} )
target_link_libraries( DicomSorter ${ITK_LIBRARIES} ${QT_LIBRARIES} )
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/DicomSorter" DESTINATION ${bin_path})
