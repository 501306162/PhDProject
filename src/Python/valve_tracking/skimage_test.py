import sys
import utils
import utils.valve_io
import geometry.alignment
import matplotlib.pyplot as plt
from skimage import exposure
import math
import random
from sklearn import svm
import skimage.feature
import skimage.filter
from sklearn.neighbors import KNeighborsRegressor
import filters.segmentation
import skimage.feature
from skimage.filter import canny
import skimage.measure

import sklearn.linear_model
import numpy as np

from sklearn.feature_extraction import image as imext
from  sklearn.decomposition import PCA
from sklearn.tree import DecisionTreeRegressor


def flip_valves(input_valves):
    to_flip = [2,7,10,12,19,22,25,26,27,28,29,33,35,36,37,40,42]
    for n in to_flip:
        input_valves[n] = geometry.alignment.flip_valve(input_valves[n])

    return input_valves

def align_valves(input_valves):
    for n, v in enumerate(input_valves):
        v = geometry.alignment.align_valves(input_valves[0],v)
    
    return input_valves

def segment_image(image, location, patch_size, sigma):
    # get the segmented image
    binary, smooth = filters.segmentation.segment_by_valve_patch(
            image, location, sigma, 
            (patch_size, patch_size))
    
    return binary


def extract_test_training(input_valves, test_number):
    valves = []
    test_valve = None
    for n, v in enumerate(input_valves):
        if n == test_number:
            test_valve = v
        else:
            valves.append(v)
        
    return valves, test_valve

def extract_feature(patch):
    lbp = skimage.feature.local_binary_pattern(patch, 3, 3*8)

    #plt.imshow(lbp)
    #plt.gray()
    #plt.show()

    return lbp.ravel()

def extract_contours(binary):
    contours  = skimage.measure.find_contours(binary,0)
    contour_indexs = []
    for contour in contours:
        for index in contour:
            x = int(index[1])
            y = int(index[0])
            contour_indexs.append([x,y])

    return contour_indexs


def main():

    input_folder = sys.argv[1]
    test_image_number = int(sys.argv[2])
    sigma = float(sys.argv[3])
    patch_size = 30
    n_features = patch_size*patch_size
    n_components = 20
    n_negative_patches = 100
    segmentation_patch_size = 50
    point_choice = "i1"

    # read the valves
    input_valves = utils.valve_io.read_valves_from_directory(input_folder, time_step=5)

    # flip and align the valves
    input_valves = flip_valves(input_valves)
    input_valves = align_valves(input_valves)

    # extract traingin and test set
    valves, test_valve = extract_test_training(input_valves, test_image_number)
    n_training_valves = len(valves)
    
            
    positive_patches = np.zeros((n_training_valves, n_features))
    negative_patch_list = []

    for n, v in enumerate(valves):
        training_image = v["image"]
        training_point = v[point_choice]
        im_shape = training_image.shape

        binary = segment_image(training_image, training_point, segmentation_patch_size, sigma)
        
        location = []
        location.append(int(training_point[0]))
        location.append(int(training_point[1]))
        patch = utils.extract_image_patch(binary, location, (patch_size, patch_size))
        
        #plt.imshow(patch)
        #plt.plot(location[0],location[1], 'o')
        #plt.gray()
        #plt.show()

        positive_patches[n,:] = extract_feature(patch)

        contours = extract_contours(binary)
        

        total_contours = len(contours)
        step = int(total_contours / n_negative_patches)
        for n in range(n_negative_patches):
            contour = contours[n*step]
            patch = utils.extract_image_patch(binary, contour, (patch_size, patch_size))
            negative_patch_list.append(extract_feature(patch))


        '''        
        for i in range(n_negative_patches):
            neg_index = n*n_negative_patches+i
            loc_rand = []
            
            start_range = 0 
            end_range_x = im_shape[0]-1
            end_range_y = im_shape[1]-1
            
            ind1 = random.randint(start_range, end_range_y)
            ind2 = random.randint(start_range, end_range_x)
            loc_rand.append(ind1)
            loc_rand.append(ind2)
        '''


    # construct the labels
    n_negative_patches = len(negative_patch_list)
    labels = np.zeros((positive_patches.shape[0]+n_negative_patches,1))
    labels[0:positive_patches.shape[0]-1] = 1
    labels[positive_patches.shape[0]:positive_patches.shape[0]+n_negative_patches] = 0

    negative_patches = np.zeros((n_negative_patches,n_features))
    for n, p in enumerate(negative_patch_list):
        negative_patches[n,:] = p


    X = np.concatenate((positive_patches, negative_patches))
    

    
    #clf = svm.SVR()
    clf = svm.SVR()
    clf.fit(X, labels.ravel())

    
    # extract the test patches 
    test_image = test_valve["image"]
    
    # get a set of indices 
    binary = segment_image(test_image, test_valve[point_choice], segmentation_patch_size, sigma)
    contours  = extract_contours(binary)

    
    test_patches = np.zeros((len(contours), n_features))
    for n, contour in enumerate(contours):
        patch = utils.extract_image_patch(binary, contour, (patch_size, patch_size))
        test_patches[n,:] = extract_feature(patch)
        

    y = clf.predict(test_patches)
    #y = clf.predict(test_patches[3])
    #print test_patches[3]

    #print y
    #y = clf.predict(test_patches[70])
    #print test_patches[70]
    #print y
    #y = clf.predict(X[10])
    #print X[10]
    #print y

    out = np.zeros(test_image.shape)
    for n, ind in enumerate(contours):
        out[ind[1],ind[0]] = y[n];



    coords = skimage.feature.peak_local_max(out, min_distance=5)


    #print negative_patches

    fig, (ax1,ax2,ax3) = plt.subplots(1,3)
    plt.gray()
    ax1.imshow(test_image)
    ax1.plot(coords[:,1], coords[:,0], 'o')
    ax1.plot(test_valve[point_choice][0], test_valve[point_choice][1], 'o')
    ax2.imshow(out)
    ax1.contour(binary, [0.5], linewidths=1.2, colors='y')
    ax3.imshow(binary)
    plt.show()



if __name__ == '__main__':
    np.set_printoptions(threshold=np.nan)
    main()
