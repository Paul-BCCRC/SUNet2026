#Writen by Paul Gallagher (BC Cancer Research Centre)
import numpy as np

def cropMirrored_getScreenPos(xp, yp, XPredict, YPredict, CenterCrops, FrameSize, border_mirror):

    numImage = len(XPredict)
    for n in range(0, numImage):
        xp[n] = xp[n] - border_mirror
        yp[n] = yp[n] - border_mirror

    if(len(FrameSize)==2):
        Xfull = FrameSize[1]
        Yfull = FrameSize[0]
    else:
        if(FrameSize[0] == min(FrameSize)):
            Xfull = FrameSize[2]
            Yfull = FrameSize[1]
        else:
            Xfull = FrameSize[1]
            Yfull = FrameSize[0]

    N = len(XPredict)
    Screenx = np.zeros( (N), dtype='int' )
    Screeny = np.zeros( (N), dtype='int' )

    nXPredict = [0]*N
    nYPredict = [0]*N
    nCenterCrops = [0]*N
    for n in range(0, N):
        nXPredict[n] = XPredict[n,:,:,:]
        nYPredict[n] = YPredict[n,:,:,:]
        nCenterCrops[n] = CenterCrops[n,:,:,:]

    for n in range(0, N):

        [N, Y, X, K] = XPredict.shape
        Screenx[n] = xp[n] - (X/2)
        Screeny[n] = yp[n] - (Y/2)

        #Horizontal
        if(Screenx[n] < 0):
            offsetx = int(np.abs(Screenx[n]))
            nXPredict[n] = XPredict[n, :, offsetx:(X-1), :]
            nYPredict[n] = YPredict[n, :, offsetx:(X-1), :]
            nCenterCrops[n] = CenterCrops[n, :, offsetx:(X-1), :]
            Screenx[n] = 0;
        elif(Screenx[n]+X > Xfull):
            offsetx = int(Screenx[n]+X - Xfull)
            nXPredict[n] = XPredict[n, :, 0:(X-offsetx), :]
            nYPredict[n] = YPredict[n, :, 0:(X-offsetx), :]
            nCenterCrops[n] = CenterCrops[n, :, 0:(X-offsetx), :]
            Screenx[n] = Xfull - (X - offsetx);
        else:
            nXPredict[n] = XPredict[n, :, :, :]
            nYPredict[n] = YPredict[n, :, :, :]
            nCenterCrops[n] = CenterCrops[n, :, :, :]

        #Vertical
        if(Screeny[n] < 0):
            offsety = int(np.abs(Screeny[n]))
            nXPredict[n] = nXPredict[n][offsety:(Y-1), :, :]
            nYPredict[n] = nYPredict[n][offsety:(Y-1), :, :]
            nCenterCrops[n] = nCenterCrops[n][offsety:(Y-1), :, :]
            Screeny[n] = 0;
        elif(Screeny[n]+Y > Yfull):
            offsety = int(Screeny[n]+Y - Yfull)
            nXPredict[n] = nXPredict[n][0:(Y-offsety), :, :]
            nYPredict[n] = nYPredict[n][0:(Y-offsety), :, :]
            nCenterCrops[n] = nCenterCrops[n][0:(Y-offsety), :, :]
            Screeny[n] = Yfull - (Y - offsety);

    Screenx = np.uint32(Screenx)
    Screeny = np.uint32(Screeny)

    return[ xp, yp, Screenx, Screeny, nXPredict, nYPredict, nCenterCrops ]