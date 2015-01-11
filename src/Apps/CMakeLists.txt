# ----------------------------------------------------------------------
# CMAKE file for the Apps directory
set(APPS 
	"Testing"
	"Testing2"
	"DataPrep"
	"SegmentationGenerator"
	"PlaneFinder"
	"ValveDataAppender"
	#"iresisd-1.1"
	"ValveTracking"
	"DicomSorter")

# ----------------------------------------------------------------------
# Add all the apps
#message("---------------------------------------------------------------")
foreach(APP ${APPS})
	set(APP_FULL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${APP}")
	AddApplication(${APP_FULL_PATH})	
endforeach()
