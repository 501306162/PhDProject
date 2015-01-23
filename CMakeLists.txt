cmake_minimum_required(VERSION 2.8)
set(project_name "PhDProject")
project(${project_name})

enable_testing()

set(QT_QMAKE_EXECUTABLE /usr/bin/qmake-qt4)


# ----------------------------------------------------------------------
# includes the cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
include(AddApplication)


# ----------------------------------------------------------------------
# Set the paths for the install output
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})

set(install_path ${CMAKE_INSTALL_PREFIX})
set(include_path "${install_path}/include/${project_name}")
set(lib_path "${install_path}/lib/${project_name}")
set(bin_path "${install_path}/bin/${project_name}")
set(doc_path "${install_path}/doc/${project_name}")


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------------------------------------------------
# add the src directory
add_subdirectory(src)



