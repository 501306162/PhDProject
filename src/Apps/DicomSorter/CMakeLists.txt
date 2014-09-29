set( sources "main.cpp"
	"SeriesExtractor.cpp"	)


add_executable( DicomSorter ${sources} )
target_link_libraries( DicomSorter ${ITK_LIBRARIES} )
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/DicomSorter" DESTINATION ${bin_path})
