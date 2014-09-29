import dicom
import os
import os.path
import sys
import gdcm

def read_directory(dir_path):
    files = [f for f in os.listdir(dir_path)]
    print files





# First we need to load up all the dicom file in the directory
read_directory(sys.argv[1])




