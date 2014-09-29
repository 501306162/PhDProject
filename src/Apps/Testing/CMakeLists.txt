# ----------------------------------------------------------------------
# Set of programs for testing
set(APP_NAMES 
	"InitTesting"
	"MRFTest"
	"DOGTest"
	"PatientDataTest"
	"DiffusionMapTest")


include_directories(${ITK_INCLUDE_DIRS})


foreach(APP ${APP_NAMES})
	add_executable(${APP} "${APP}.cpp")
	add_dependencies(${APP} Utils Manifold)
	target_link_libraries(${APP} ${ITK_LIBRARIES} Utils Manifold MRFLib)
	install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${APP}" DESTINATION ${bin_path})
endforeach()


