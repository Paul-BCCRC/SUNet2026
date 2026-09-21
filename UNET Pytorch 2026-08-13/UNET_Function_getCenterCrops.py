#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np

def getCenterCrops(rframe, sframe, xp, yp, B_X, B_Y, isBackgroundDark):
    FrameSize = rframe.shape
    X = FrameSize[1]
    Y = FrameSize[0]
    try:
        Z = FrameSize[2]
    except:
        Z = 1

    B_X2 = B_X/2
    B_Y2 = B_Y/2

    numcells = len(xp)
    Xpredict = np.zeros((numcells, B_Y, B_X, Z), dtype='uint8')
    Centerframes = np.zeros((numcells, B_Y, B_X, 1), dtype='uint8')

    for n in range(0, numcells):

        L = int(xp[n]-B_X2)
        R = int(xp[n]+B_X2)
        T = int(yp[n]-B_Y2)
        B = int(yp[n]+B_Y2)

        if(  (L>=0) & (R<(X-1)) & (T>=0) & (B<(Y-1))  ):
            Xpredict[n, :, :, :] = rframe[T:B, L:R, :]
            if(isBackgroundDark==0):
                Xpredict[n, :, :, 0] = np.invert(rframe[T:B, L:R, 0])

            Centerframes[n, :, :, 0] = sframe[T:B, L:R]

    return [Xpredict, Centerframes]