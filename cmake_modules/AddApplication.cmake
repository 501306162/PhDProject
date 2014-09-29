# ----------------------------------------------------------------------
# Functions and macros for adding an application
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# Main add application function
function(AddApplication folder_path)

	# get the folder name
	get_filename_component(folder_name ${folder_path} NAME)

	#message("Adding Application in ${folder_name}")
	

	# check if a cmakelists file exists
	set(lists_file "${folder_path}/CMakeLists.txt")
	if(EXISTS ${lists_file})

		# do the normal add directory
		#message("Cmake File Exists: Adding Sub Directory")
		add_subdirectory(${folder_path})
	else()
		#message("Cmake File Doesn't exist: Doing defualt Adding")
		#AddDefualt(${folder_path})

		

	endif()
	#message("-------------------------------------------------------")		
endfunction()


# ----------------------------------------------------------------------
# Function to add an application in the defualt manner
function(AddDefualt folder_path)
	set(src_files "")
	GetCppFiles(${folder_path} ${src_files})
endfunction()


# ----------------------------------------------------------------------
# Function to get all the source files from a folder
function(GetCppFiles folder cpp_files)
	file(GLOB ${cpp_files} RELATIVE ${folder} *.cpp)	
	#message(${cpp_files})
endfunction()

