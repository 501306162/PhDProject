import sys
import json
import nrrd
import os
import os.path
import re
import numpy as np

class ValveIO:

    def __init__(self):
        self.filename = "";
        self.image_ending = "_image.nrrd"

    # Read the input file and return the data 
    def load(self, filename, time_step=0):
        part = self.parse_json(filename)[time_step];
        return self.build_entry(part)

    def load_all(self, filename):
        return [self.build_entry(x) for x in self.parse_json(filename)]

    def parse_json(self, filename):
        f = open(filename)
        string_data = f.read().encode('utf-8') 
        return json.loads(string_data)

    def construct_filenames(self, base):
        image_fname = base + self.image_ending
        return str(image_fname)

    def parse_points(self, point):
        return np.array([point["x"], point["y"], point["z"]])

    def parse_index(self, index):
        rows = index["y"]
        cols = index["x"]
        return np.array([rows,cols])

    def build_entry(self, part):
        output = {}

        # read the image
        image_filename = self.construct_filenames(part["filename"])
        output["image_data"] = nrrd.read(image_filename)
        output["image"] = output["image_data"][0][:,:,0]

        # parse the points
        output["p1"] = self.parse_points(part["p1"])
        output["p2"] = self.parse_points(part["p2"])        
        output["i1"] = self.parse_index(part["i1"])
        output["i2"] = self.parse_index(part["i2"])

        return output

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]  
    

# procedural loading of a valve
def read_valve(filename, time_step=0):
    io = ValveIO()
    return io.load(filename, time_step)


# read all the valves from a given directory
def read_valves_from_directory(directory, time_step=0):
    filenames = [os.path.join(directory,f) for f in os.listdir(directory) if ".txt" in f]

    #sort the filenames
    filenames.sort(key=natural_sort_key)
    return [read_valve(x, time_step) for x in filenames]
    


# test the class
def main():
    io = ValveIO()
    entry = io.load(sys.argv[1])
    print entry["p1"]
    print entry["p2"]
    print entry["i1"]
    print entry["i2"]

    out = read_valves_from_directory(sys.argv[2])
    print out[0]["image"]
    

if __name__ == '__main__':
    main()

