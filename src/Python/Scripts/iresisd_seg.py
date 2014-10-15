"""
Functions for running the iresisd program
"""
import os 
import os.path
import sys

from xml_template_builder import *
from run_progs import *

def run_iresisd(input_data, data_item, choices=["rv", "lv"]):

    # are the items initialised ready for the iresisd
    if not check_initialistations(data_item):
        print "Data Item: " + data_item["name"] + " not fully initialised"
        sys.exit(1)
        
    # make sure the output directories exist
    prepare_directories(input_data, data_item)

    
    # iterate through the parts
    for part in data_item["parts"]:
        if part["name"] not in choices:
            continue

        xml_file = create_xml_file(input_data, data_item, part)
        cwd = os.path.dirname(xml_file)

        iresisd = Iresisd()
        iresisd.exe(input_data["iresisd_path"])
        iresisd.simple_xml(xml_file)
        iresisd.cwd(cwd)
        iresisd.run()



def create_xml_file(input_data, data_item, part):

    # create the xml file 
    xml_output = os.path.join(
            input_data["iresisd_output_base_path"], 
            data_item["name"], 
            part["name"], 
            "simple_example.xml")
    
    data_path = os.path.join(input_data["dicom_base_path"], data_item["name"])



    xml_builder = XmlTemplateBuilder(input_data["xml_template"])
    xml_builder.iterations(part["iterations"])
    xml_builder.circle(
            part["circle"]["r"],
            part["circle"]["c"][0],
            part["circle"]["c"][1],
            part["circle"]["c"][2])
    xml_builder.filename(xml_output)
    xml_builder.dicom_dir(data_path)
    xml_builder.first_slice(False)
    xml_builder.display(1)
    xml_builder.generate()

    return xml_output



def check_initialistations(data_item):
    for part in data_item["parts"]:
        if not part["initialised"]:
            return False

    return True


def prepare_directories(input_data, data_item):
    
    # check the base path exists
    base_path = os.path.join(input_data["sorted_iresisd_base_path"], data_item["name"])

    if not os.path.isdir(base_path):
        os.makedirs(base_path)

    # create the sub paths
    for part in data_item["parts"]:
        sub_path = os.path.join(base_path, part["name"])

        if not os.path.isdir(sub_path):
            os.makedirs(sub_path)
    





