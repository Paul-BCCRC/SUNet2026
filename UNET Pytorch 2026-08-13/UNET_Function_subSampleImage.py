#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np

def subSampleImage(Mframe, tsize, overlap, isBackgroundDark):

    if(isBackgroundDark==0):
        Mframe = np.invert(Mframe)
    dsize = tsize - overlap
    FrameSize = Mframe.shape
    X = FrameSize[1]
    Y = FrameSize[0]
    numX=int(np.floor( (X-tsize) / dsize ) + 2)
    numY=int(np.floor( (Y-tsize) / dsize ) + 2)

    Xpredict = np.zeros(( (numX)*(numY), tsize, tsize, 1), dtype='uint8')

    for xx in range(0, (numX)):
        for yy in range(0, (numY)):

            if( (xx != (numX-1)) and (yy != (numY-1)) ):
                temp = Mframe[ (yy*dsize):(yy*dsize+tsize) , (xx*dsize):(xx*dsize+tsize) ]

            elif( (xx == (numX-1)) and (yy == (numY-1)) ):
                temp = Mframe[ (Y-tsize-1):(Y-1) , (X-tsize-1):(X-1) ]

            elif( xx == (numX-1) ):
                temp = Mframe[ (yy*dsize):(yy*dsize+tsize) , (X-tsize-1):(X-1) ]

            elif( yy == (numY-1) ):
                temp = Mframe[ (Y-tsize-1):(Y-1) , (xx*dsize):(xx*dsize+tsize) ]

            Xpredict[(xx*numY)+yy, :, :, 0] = temp

    return [Xpredict, numX, numY, int(np.floor( (X-tsize) / dsize )*dsize +tsize), int(np.floor( (Y-tsize) / dsize )*dsize +tsize)]