# ----------------------------------------------------------------------
# Standard definitions
cmake_minimum_required( VERSION 2.8 ) 
project( iresisd )

set( iresisd_VERSION_MAJOR 1 )
set( iresisd_VERSION_MINOR 1 )


# ----------------------------------------------------------------------
# Set up the dependencies
find_package( OpenCV REQUIRED )
#find_package( DCMTK REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
#include_directories( ${DCMT_INCLUDE_DIRS} )

include_directories( "/usr/local/include/c++")
include_directories( "/usr/local/include/nifti" )
include_directories( "/usr/local/include/dcmtk" )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	add_definitions(-DNO_OMP)
endif()

set(DCMTK_LIBRARIES charls
dcmimage
dcmjpls
dcmqrdb
dcmtls
ijg12
oflog
dcmdata
dcmimgle
dcmnet
dcmwlm
ijg16
ofstd
dcmdsig
dcmjpeg
dcmpstat
dcmsr
ijg8
bz2
)

# ----------------------------------------------------------------------
# Set up the source code to build the executable
set( SRC_DIR "src" )
set( iresisd_SRCS 
	"${SRC_DIR}/Sdata.cpp"
	"${SRC_DIR}/Phi.cpp"
	"${SRC_DIR}/XMLreader.cpp"
	"${SRC_DIR}/ISISD.cpp"
	"${SRC_DIR}/SpeedFunction.cpp"
	"${SRC_DIR}/VolumeData.cpp"
	"${SRC_DIR}/IReSD.cpp"
	"${SRC_DIR}/main.cpp"
	"${SRC_DIR}/Interface.cpp"
	"${SRC_DIR}/RBF.cpp"
	"${SRC_DIR}/Sgeom.cpp"
	"${SRC_DIR}/AuxClasses.cpp"
	"${SRC_DIR}/tinyxml2.cpp" )

# ----------------------------------------------------------------------
# Build the executable
add_executable( iresisd ${iresisd_SRCS} )
target_link_libraries( iresisd fftw3 fftw3_threads  ${OpenCV_LIBRARIES} ${DCMTK_LIBRARIES} niftiio teem )

install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/iresisd" DESTINATION bin)

