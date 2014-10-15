"""
Functions for running the program for the iresisd cleanup
"""
import sys
import os
import os.path

from run_progs import *


def iresisd_cleanup(input_data, data_item):
    # are the items initialised ready for the iresisd
    if not check_initialistations(data_item):
        print "Data Item: " + data_item["name"] + " not segmented"
        sys.exit(1)


    prepare_directories(input_data, data_item)


    for part in data_item["parts"]:
        input_dir = os.path.join(
                input_data["iresisd_output_base_path"], 
                data_item["name"], 
                part["name"])

        output_dir = os.path.join(
                input_data["sorted_iresisd_base_path"], 
                data_item["name"], 
                part["name"]) 

        sorter = IresisdSorter()
        sorter.exe(input_data["iresisd_sorter_path"])
        sorter.input_dir(input_dir)
        sorter.output_dir(output_dir)
        sorter.run()


    



def check_initialistations(data_item):
    for part in data_item["parts"]:
        if not part["initial_seg"]:
            return False

    return True


def prepare_directories(input_data, data_item):
    
    # check the base path exsorted_iresisd_base_pathists
    base_path = os.path.join(input_data["sorted_iresisd_base_path"], data_item["name"])

    if not os.path.isdir(base_path):
        os.makedirs(base_path)

    # create the sub paths
    for part in data_item["parts"]:
        sub_path = os.path.join(base_path, part["name"])

        if not os.path.isdir(sub_path):
            os.makedirs(sub_path)
    



