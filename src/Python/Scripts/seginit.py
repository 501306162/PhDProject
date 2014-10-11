""" 
class to run the initialisation for the segmentation.
What this does, is it runs the iresisd executable with only a single iteration.
This will build the output volume. The output volume is then loaded into the 
MITK executable where we can select the circles for the iresisd to be run 
from
"""
import sys
import subprocess
import os
import os.path

class SegInit(object):

    """Docstring for SegInit. """

    def __init__(self):
        """@todo: to be defined1. """

    def data(self, data):
        self._data = data
    
    def mitk_exe(self, exe):
        self._mitk_exe = exe

    def iresisd_exe(self, exe):
        self._iresisd_exe = exe

    def simple_xml(self, xml):
        self._simple_xml = xml

    def cwd(self, cwd):
        self._cwd = cwd

    def run(self):
        # setup the subprocess
        proc = subprocess.Popen([ self._iresisd_exe, self._simple_xml ], cwd=self._cwd)
        success = proc.wait()

        if not success == 0:
            print "There was an error running iresisd on this data"
            sys.exit()

        # now run mitk with the output volume as the argument
        output_volume = os.path.join(self._cwd, "VolumeData.nii")


        proc2 = subprocess.Popen([ self._mitk_exe, output_volume ], cwd=self._cwd)
        success = proc2.wait()

        



