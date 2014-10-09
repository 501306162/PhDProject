bin_PROGRAMS=$(top_builddir)/bin/IReSISD

__top_builddir__bin_IReSISD_SOURCES=AuxClasses.cpp tinyxml2.cpp XMLreader.cpp VolumeData.cpp RBF.cpp Phi.cpp SpeedFunction.cpp Sdata.cpp Sgeom.cpp ISISD.cpp IReSD.cpp Interface.cpp main.cpp

AM_CXXFLAGS=-Wall -fopenmp

__top_builddir__bin_IReSISD_LDADD=-ldcmdata -loflog -lofstd -L/usr/lib/ -lopencv_core -lopencv_imgproc -lopencv_features2d -lopencv_flann -lopencv_contrib -lopencv_calib3d -lopencv_gpu -lopencv_legacy -lopencv_ml -lopencv_objdetect -lopencv_video -lopencv_highgui -lfftw3_threads -lfftw3 -lniftiio -lnifticdf -lgomp -lteem
