#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from skimage.feature import peak_local_max

def getCenterPoints(cframe, rframe, M, oX, oY, threshold, SaturationOffset, R, isBackgroundDark):
    # start_init = time.time()

    #threshold = 40
    #R = 7          #Center area Radius
    #SaturationOffset = 15

    FrameSize = cframe.shape
    X = FrameSize[1]
    Y = FrameSize[0]

    neighborhood_size = 10
    data_max = filters.maximum_filter(cframe, neighborhood_size)
    maxima = (cframe == data_max)
    data_min = filters.minimum_filter(cframe, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)

    #Center Intensity threashold
    if(isBackgroundDark):
        LOCALTHREASHOLD = SaturationOffset
    else:
        LOCALTHREASHOLD = 255 - SaturationOffset


    x, y = [], []
    for dy,dx in slices:
        # x_center = (dx.start + dx.stop - 1)/2
        # y_center = (dy.start + dy.stop - 1)/2
        x_center = dx.start
        y_center = dy.start

        #Keep center if it is within original frame and meets intensity threashold.
        if( (x_center>M ) & (x_center<(M+oX)) & (y_center>M ) & (y_center<(M+oY)) ):
            LOCALmean = np.mean(rframe[int(y_center-R):int(y_center+R+1), int(x_center-R):int(x_center+R+1)])
            if(isBackgroundDark):
                if(LOCALmean > LOCALTHREASHOLD):
                    x.append(x_center)
                    y.append(y_center)
            else:
                if(LOCALmean < LOCALTHREASHOLD):
                    x.append(x_center)
                    y.append(y_center)

    sframe = np.zeros((Y, X), dtype='uint8')
    sframe[maxima == 1] = 255


    return [np.int64(x), np.int64(y), sframe]

def getCenterPoints_TopHat(cframe, rframe, M, oX, oY, threshold, SaturationOffset, R, isBackgroundDark):
    # start_init = time.time()

    #threshold = 40
    #R = 7          #Center area Radius
    #SaturationOffset = 15

    FrameSize = cframe.shape
    X = FrameSize[1]
    Y = FrameSize[0]

    # neighborhood_size = 10
    # data_max = filters.maximum_filter(cframe, neighborhood_size)
    # maxima = (cframe == data_max)
    # data_min = filters.minimum_filter(cframe, neighborhood_size)
    # diff = ((data_max - data_min) > threshold)
    # maxima[diff == 0] = 0

    #struct = ndimage.generate_binary_structure(2, 1)
    struct = np.ones((5,5),dtype='bool')
    struct[0,0] = False
    struct[0,4] = False
    struct[4,0] = False
    struct[4,4] = False
    tmp = ndimage.grey_erosion(cframe, structure=struct)
    tmp = ndimage.grey_dilation(tmp, structure=struct)
    maxima = np.uint8( np.clip( np.double(cframe) - np.double(tmp) , 0, 255 ) )
    #maxima = ndimage.white_tophat(cframe, structure=struct)
    maxima[maxima <= 40] = 0

    coordinates = peak_local_max(maxima, min_distance=3)
    N = len(coordinates)

    #labeled, num_objects = ndimage.label(maxima)
    #slices = ndimage.find_objects(labeled)


    #Center Intensity threashold
    if(isBackgroundDark):
        LOCALTHREASHOLD = SaturationOffset
    else:
        LOCALTHREASHOLD = 255 - SaturationOffset

    sframe = np.zeros((Y, X), dtype='uint8')
    #sframe[maxima == 1] = 255

    x, y = [], []
    for n in range(0, N):
        # x_center = (dx.start + dx.stop - 1)/2
        # y_center = (dy.start + dy.stop - 1)/2
        x_center = coordinates[n,1]
        y_center = coordinates[n,0]
        sframe[y_center, x_center] = 255

        #Keep center if it is within original frame and meets intensity threashold.
        if( (x_center>M ) & (x_center<(M+oX)) & (y_center>M ) & (y_center<(M+oY)) ):
            LOCALmean = np.mean(rframe[int(y_center-R):int(y_center+R+1), int(x_center-R):int(x_center+R+1)])
            if(isBackgroundDark):
                if(LOCALmean > LOCALTHREASHOLD):
                    x.append(x_center)
                    y.append(y_center)
            else:
                if(LOCALmean < LOCALTHREASHOLD):
                    x.append(x_center)
                    y.append(y_center)



    return [np.int64(x), np.int64(y), sframe]