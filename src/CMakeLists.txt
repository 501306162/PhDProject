# ----------------------------------------------------------------------
# CMAKE file for the source directory
# ----------------------------------------------------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# ----------------------------------------------------------------------
# Get the includes together that will be used by all the different parts
find_package(ITK REQUIRED)
find_package(VTK REQUIRED)
find_package(Qt4 REQUIRED)
find_package(Eigen3 REQUIRED)


message(${QT_INCLUDE_DIR})

include(${QT_USE_FILE})
include(${ITK_USE_FILE})

if(APPLE)
	include_directories("/usr/local/include/eigen3")
else()
	include(${EIGEN_USE_FILE})
endif()

enable_testing()

# ----------------------------------------------------------------------
# Compile and include gtest and gmock
add_subdirectory("External")

# ----------------------------------------------------------------------
# Set up the tests 
include_directories("${gmock_SOURCE_DIR}/include")
include_directories("${gtest_SOURCE_DIR}/include")


# ----------------------------------------------------------------------
# Gather together all the source subdirs and add them to the cmake process
set(SRC_SUBDIRS 
	"MRFLib"
	"Utils"
	"Filters"
	"Manifold"
	"Apps"
	"Registration"
	"Segmentation"
	)

include_directories(${CMAKE_CURRENT_DIR} ${SRC_SUBDIRS} ${EIGEN_INCLUDE_DIRS})

foreach(SUBDIR ${SRC_SUBDIRS})
	add_subdirectory(${SUBDIR})
endforeach()

