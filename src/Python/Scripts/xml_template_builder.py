"""
class that creates the xml files for the iresisd program
"""

class XmlTemplateBuilder(object):

    """Docstring for XmlTemplateBuilder. """

    def __init__(self, template):
        """@todo: to be defined1.

        :options: @todo

        """
        self._template = template
        self._iterations = 1
        self._r = 5
        self._x = 100
        self._y = 100
        self._z = 100
        self._filename = "simple_example.xml"
        self._dicom_dir = ""
        self._display = 0
        self._first_slice = True
    
    def display(self, display):
        self._display = display

    def dicom_dir(self, dir):
        self._dicom_dir = dir

    def iterations(self, iter):
        self._iterations = iter;

    def circle(self, r, x, y, z):
        self._r = r
        self._x = x
        self._y = y
        self._z = z

    def first_slice(self, val):
        self._first_slice = val

    def filename(self, filename):
        self._filename = filename

    def generate(self):
        # open the template
        input_file = open(self._template, "r")
        lines = input_file.read()
        input_file.close()

        # do the string substitutions
        lines = lines.replace("{NUM_ITERATIONS}", str(self._iterations))
        lines = lines.replace("{DATA_FOLDER}", str(self._dicom_dir))
        lines = lines.replace("{DISPLAY}", str(self._display))
        lines = lines.replace("{X}", str(self._x))
        lines = lines.replace("{Y}", str(self._y))
        lines = lines.replace("{Z}", str(self._z))
        lines = lines.replace("{R}", str(self._r))
        
        if self._first_slice:
            lines = lines.replace("{FIRST_SLICE}", "true")
        else:
            lines = lines.replace("{FIRST_SLICE}", "false")

        # write the output file
        output_file = open(self._filename, "w")
        output_file.write(lines)
        output_file.close()


