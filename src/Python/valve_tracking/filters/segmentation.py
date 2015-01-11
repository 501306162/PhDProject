import skimage.exposure
import utils

def segment_by_valve_patch(image, location, sigma=1.0, size=(3,3)):

    print location, image.shape
    # smooth the input
    image_smooth = skimage.filter.gaussian_filter(image, sigma)

    patch = utils.extract_valve_patch(image_smooth, location, size)
    seg = segment_by_patch(image_smooth, patch)

    #return image_smooth, image_smooth
    return seg, image_smooth



def segment_by_patch(image, patch):
    threshold = skimage.filter.threshold_otsu(patch)
    return image < threshold



