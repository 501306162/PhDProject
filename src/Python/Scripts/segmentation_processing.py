import json
import argparse
import sys
import os
import os.path
from xml_template_builder import *
from seginit import *


# decoding functions
def _decode_list(data):
    rv = []
    for item in data:
        if isinstance(item, unicode):
            item = item.encode('utf-8')
        elif isinstance(item, list):
            item = _decode_list(item)
        elif isinstance(item, dict):
            item = _decode_dict(item)
        rv.append(item)
    return rv

# decoding functions 
def _decode_dict(data):
    rv = {}
    for key, value in data.iteritems():
        if isinstance(key, unicode):
            key = key.encode('utf-8')
        if isinstance(value, unicode):
            value = value.encode('utf-8')
        elif isinstance(value, list):
            value = _decode_list(value)
        elif isinstance(value, dict):
            value = _decode_dict(value)
        rv[key] = value
    return rv


# function to read the json data
def read_json(fname):
    json_data = open(fname);
    data = json.load(json_data, object_hook=_decode_dict);
    return data



# function to parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file", help="The json file containing the data for labelling process")
    parser.add_argument("-d", "--data_set", type=int, help="Choose a specific data set to work on or use all")
    parser.add_argument("-m", "--mode", type=int, help="Select the mode to run with [initialise (0), iresisd (1), sort_iresisd (2), amend (3), interpolate (4)]", default=0)
    return parser.parse_args()


# function to extract the data that we are going to be working on
def extract_data(json_data, data_set):
    data = []
    if not data_set:
        data = [p for p in json_data['instances']]
    else:
        #check the request patient is valid
        data = [p for p in json_data['instances'] if p["number"] == data_set]

    if len(data) == 0:
        print "No data sets found"
        sys.exit(1)
    else:
        return data





def main():
    args = parse_args()
    input_data = read_json(args.json_file)

    #iterate through the set of data sets we are working with
    for data_item in extract_data(input_data, args.data_set):

        # get some variables
        data_path = os.path.join(input_data["dicom_base_path"], data_item["name"])

        # initialisation mode
        if args.mode == 0:

            xml_output = os.path.join(input_data["initialisation_dump"],"simple_example.xml")

            builder = XmlTemplateBuilder(input_data["xml_template"])
            builder.iterations(1)
            builder.dicom_dir(data_path)
            builder.filename(xml_output)
            builder.generate()


            seg_initialiser = SegInit()
            seg_initialiser.simple_xml(xml_output)
            seg_initialiser.iresisd_exe(input_data["iresisd_path"])
            seg_initialiser.mitk_exe(input_data["mitk_path"])
            seg_initialiser.cwd(input_data["initialisation_dump"])
            seg_initialiser.run()
            

    

    





if __name__ == '__main__':
    main()

