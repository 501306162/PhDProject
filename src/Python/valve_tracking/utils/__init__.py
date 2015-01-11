import numpy as np
import sys
import skimage.util


def extract_image_patch(image, location, size=(3,3)):

    # pad the input by the size
    padded = skimage.util.pad(image,size,mode='minimum')

    startr = (location[1]+size[0]) - (size[0]/2)
    startc = (location[0]+size[1]) - (size[1]/2)
    endr = startr+size[0]
    endc = startc+size[1]
    
    #print "Starts: ", startr, startc


    patch = padded[startr:endr,startc:endc]

    return patch
    


def neighbors(arr,x,y,size=(3,3)):
    ''' Given a 2D-array, returns an nxn array whose "center" element is arr[x,y]'''
    arr=np.roll(np.roll(arr,shift=-x+1,axis=0),shift=-y+1,axis=1)
    return arr[:size[0],:size[1]]

def extract_image_patch_float(image, location, size=(3,3)):
    loc = []
    loc.append(int(np.round(location[0])))
    loc.append(int(np.round(location[1])))
    
    return extract_image_patch(image, loc, size)


def extract_valve_patch(image, point, size=(3,3)):
    location = []
    location.append(point[1])
    location.append(point[0])
    return extract_image_patch_float(image, location, size)

    


    




