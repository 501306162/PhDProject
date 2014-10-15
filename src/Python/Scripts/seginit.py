"""
Set of functions to run the initialisation
"""
import os
import os.path

from xml_template_builder import *
from run_progs import *


# run the initialisation
def run_initialisation(input_data, data_item):

    xml_file = create_xml_file(input_data, data_item)
    working_dir = input_data["initialisation_dump"]
    
    iresisd = Iresisd()
    iresisd.exe(input_data["iresisd_path"])
    iresisd.simple_xml(xml_file)
    iresisd.cwd(working_dir)
    iresisd.run()


    # get the path for the viewer 
    image_path = os.path.join(working_dir, "VolumeData.nii")

    mitk = MITK()
    mitk.exe(input_data["mitk_path"])
    mitk.images([image_path])
    mitk.run()



# create the simple xml file that will be used by the iresisd program
def create_xml_file(input_data, data_item):

    xml_output = os.path.join(input_data["initialisation_dump"],"simple_example.xml")
    data_path = os.path.join(input_data["dicom_base_path"], data_item["name"])

    builder = XmlTemplateBuilder(input_data["xml_template"])
    builder.iterations(1)
    builder.dicom_dir(data_path)
    builder.filename(xml_output)
    builder.generate()

    return xml_output

