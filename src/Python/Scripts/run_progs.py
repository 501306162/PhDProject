"""
Set of classes that run the MITK and iresisd programs
"""

import subprocess


class IresisdSorter(object):

    """Docstring for IresisdSorter. """

    def __init__(self):
        """@todo: to be defined1. """
        self._input_dir = ""
        self._output_dir = ""
        self._exe_path = ""

    def exe(self, exe):
        self._exe_path = exe

    def input_dir(self, input_dir):
        self._input_dir = input_dir

    def output_dir(self, output_dir):
        self._output_dir = output_dir

    def run(self):

        proc = subprocess.Popen([ self._exe_path, self._input_dir, self._output_dir ])
        success = proc.wait()

        if not success == 0:
            print "There was an error running iresisd on this data"
            sys.exit()



class Iresisd(object):

    """Docstring for Iresisd. """

    def __init__(self):
        """@todo: to be defined1. """
        self._cwd = ""
        self._simple_xml = ""
        self._exe_path = ""

    def exe(self, exe):
        self._exe_path = exe

    def simple_xml(self, xml):
        self._simple_xml = xml

    def cwd(self, cwd):
        self._cwd = cwd

    def run(self):
        # setup the subprocess
        proc = subprocess.Popen([ self._exe_path, self._simple_xml ], cwd=self._cwd)
        success = proc.wait()

        if not success == 0:
            print "There was an error running iresisd on this data"
            sys.exit()


class MITK(object):

    """Docstring for MITK. """

    def __init__(self):
        """@todo: to be defined1. """
        self._exe_path = ""
        self._images = []

    def exe(self, exe):
        self._exe_path = exe

    def images(self, images):
        self._images = images

    def run(self):

        # create the args
        args = [self._exe_path] + self._images

        proc2 = subprocess.Popen(args)
        success = proc2.wait()

        
