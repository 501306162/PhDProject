import utils
import utils.valve_io
import filters.segmentation
import sys
import matplotlib.pylab as plt
from skimage.filter import threshold_otsu


def main():
    valve = utils.valve_io.read_valve(sys.argv[1])
    image = valve["image"]


    # extract the patch
    patch_size = 30
    seg, smooth = filters.segmentation.segment_by_valve_patch(
            image, valve["i1"], 
            float(sys.argv[2]), (patch_size,patch_size))


       
    # display the image
    fig, (ax1,ax2,ax3) = plt.subplots(1,3)

    plt.gray()
    ax1.imshow(image)
    ax2.imshow(smooth)
    ax3.imshow(seg)


    plt.show()


if __name__ == '__main__':
    main()
