#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np
import scipy.ndimage as ndimage

def generateOverlapImage(frames, tsize, overlap, numX, numY, mX, mY, isBackgroundDark):

    FrameSize = frames.shape
    if(isBackgroundDark==0):
        for n in range(0, FrameSize[0]):
            frames[n,:,:,0] = np.invert(frames[n,:,:,0])

    dsize = tsize - overlap

    #X = tsize + (numX-1)*dsize
    #Y = tsize + (numY-1)*dsize
    rframe = np.zeros( (mY , mX) , dtype='uint8')


    for xx in range(0, (numX)):
        for yy in range(0, (numY)):

            if( (xx != (numX-1)) and (yy != (numY-1)) ):
                rframe[(yy*dsize):(yy*dsize+tsize) , (xx*dsize):(xx*dsize+tsize)] = frames[ (xx*numY)+yy, :, :, 0 ]

            elif( (xx == (numX-1)) and (yy == (numY-1)) ):
                rframe[(mY-tsize-1):(mY-1) , (mX-tsize-1):(mX-1)] = frames[ (xx*numY)+yy, :, :, 0 ]

            elif( xx == (numX-1) ):
                rframe[(yy*dsize):(yy*dsize+tsize) , (mX-tsize-1):(mX-1)] = frames[ (xx*numY)+yy, :, :, 0 ]

            elif( yy == (numY-1) ):
                rframe[(mY-tsize-1):(mY-1) , (xx*dsize):(xx*dsize+tsize)] = frames[ (xx*numY)+yy, :, :, 0 ]

                #rframe[(yy*dsize):(yy*dsize+tsize) , (xx*dsize):(xx*dsize+tsize)] = frames[ (xx*numY)+yy, :, :, 0 ]

    return rframe


    # dsize = tsize - overlap

    # X = tsize + (numX-1)*dsize
    # Y = tsize + (numY-1)*dsize
    # rframe = np.zeros( (Y , X) , dtype='uint8')

    # if(isBackgroundDark):
    #     for xx in range(0, (numX)):
    #         for yy in range(0, (numY)):
    #             rframe[(yy*dsize):(yy*dsize+tsize) , (xx*dsize):(xx*dsize+tsize)] = frames[ (xx*numY)+yy, :, :, 0 ]
    # else:
    #     for xx in range(0, (numX)):
    #         for yy in range(0, (numY)):
    #             rframe[(yy*dsize):(yy*dsize+tsize) , (xx*dsize):(xx*dsize+tsize)] = np.invert(frames[ (xx*numY)+yy, :, :, 0 ])

    # return rframe

def generateOverlapImage_WithCenter(frames, tsize, overlap, numX, numY, mX, mY):

    dsize = tsize - overlap
    #qsize = int(tsize / 4)
    qsize = int((tsize/2) - (dsize/2))
    hsize = int(dsize/2)
    h2size = int(tsize/2)
    #X = tsize + (numX-1)*dsize
    #Y = tsize + (numY-1)*dsize
    cframe = np.zeros( (mY , mX) , dtype='uint8')
    #NOT USING IMAGE CENTERS
    for xx in range(0, (numX)):
        for yy in range(0, (numY)):

            if( (xx != (numX-1)) and (yy != (numY-1)) ):
                cframe[(yy*dsize):(yy*dsize+tsize) , (xx*dsize):(xx*dsize+tsize)] = frames[ (xx*numY)+yy, :, :, 0 ]

            elif( (xx == (numX-1)) and (yy == (numY-1)) ):
                cframe[(mY-tsize-1):(mY-1) , (mX-tsize-1):(mX-1)] = frames[ (xx*numY)+yy, :, :, 0 ]

            elif( xx == (numX-1) ):
                cframe[(yy*dsize):(yy*dsize+tsize) , (mX-tsize-1):(mX-1)] = frames[ (xx*numY)+yy, :, :, 0 ]

            elif( yy == (numY-1) ):
                cframe[(mY-tsize-1):(mY-1) , (xx*dsize):(xx*dsize+tsize)] = frames[ (xx*numY)+yy, :, :, 0 ]



            #cframe[((yy*dsize+1)-hsize+h2size):((yy*dsize+1)+hsize+h2size) , ((xx*dsize+1)-hsize+h2size):((xx*dsize+1)+hsize+h2size)] = frames[ (xx*numY)+yy, qsize:(qsize+dsize), qsize:(qsize+dsize), 0 ]

    cframe = ndimage.gaussian_filter(cframe, sigma=1);

    return cframe