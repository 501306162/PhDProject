import skimage.transform as tf
import numpy as np
import matplotlib.pylab as plt

# ----------------------------------------------------------------------
def compute_valve_alignment(valve1, valve2, flatten=False):
    v1p1 = valve1["i1"]
    v1p2 = valve1["i2"]
    v2p1 = valve2["i1"]
    v2p2 = valve2["i2"]

    # estimate the transform
    src = np.array(( (v2p1[0],v2p1[1]), (v2p2[0],v2p2[1]) ));
    dst = np.array(( (v1p1[0],v1p1[1]), (v1p2[0],v1p2[1]) ));

    tform = tf.SimilarityTransform()
    tform.estimate(dst,src)

    return tform


# ----------------------------------------------------------------------
def align_valves(valve1, valve2):

    tform = compute_valve_alignment(valve1, valve2)

    valve2["image"] = tf.warp(valve2["image"], tform)
    valve2["i1"] = valve1["i1"]
    valve2["i2"] = valve1["i2"]

    return valve2


# ----------------------------------------------------------------------
def transform_index(index, im_shape, transform):
    return transform(index)[0]
    point = []
    point.append(index[0])
    point.append(im_shape[1]-index[1])

    tpoint = ptransform(point)[0]
    outindex = []
    outindex.append(tpoint[0])


# ----------------------------------------------------------------------
def flip_valve(valve):
    valve["image"] = np.swapaxes(valve["image"],0,1)
    i1x = valve["i1"][0]
    i1y = valve["i1"][1]
    i2x = valve["i2"][0]
    i2y = valve["i2"][1]
    valve["i1"][0] = i1y
    valve["i1"][1] = i1x
    valve["i2"][0] = i2y
    valve["i2"][1] = i2x

    return valve

