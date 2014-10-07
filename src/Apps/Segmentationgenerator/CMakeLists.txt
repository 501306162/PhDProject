# ----------------------------------------------------------------------
# Gather all the sources
set( sources "main.cpp" "functions.cpp" "InitialTransformExtractor.cpp" )

add_executable( SegmentationGenerator ${sources} )
add_dependencies( SegmentationGenerator Utils Manifold)
target_link_libraries( SegmentationGenerator ${ITK_LIBRARIES} nlopt MRFLib Utils Manifold )
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/SegmentationGenerator" DESTINATION ${bin_path})
