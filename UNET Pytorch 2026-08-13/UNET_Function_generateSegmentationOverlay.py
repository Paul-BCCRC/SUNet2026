#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np
import skimage.morphology as morph
import skimage.io as skio

def generateSegmentationOverlay(rframe, cframe, sframe, Ypredict, xp, yp, border_mirror, oX, oY):

    # numImage = len(xp)
    # for n in range(0, numImage):
    #     xp[n] = xp[n] + border_mirror
    #     yp[n] = yp[n] + border_mirror

    invframe = np.invert(rframe)

    FrameSize = rframe.shape
    tileSize = Ypredict.shape
    X = FrameSize[1]
    Y = FrameSize[0]
    # oX = OFrameSize[1]
    # oY = OFrameSize[0]
    xx = tileSize[1]
    yy = tileSize[2]
    Overlay = np.zeros((Y, X, 3), dtype='uint8');
    OverlayALL = np.zeros((Y, X), dtype='uint8');
    tframe = np.zeros((Y, X, 3), dtype='float32');
    iframe = np.zeros((Y, X, 3), dtype='float32');
#    frame2 = np.zeros((Y, X, 3), dtype='uint8');
#    frame3 = np.zeros((Y, X, 3), dtype='uint8');
    numcells = len(xp)

    colors = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 0.5, 0],
            [1, 0, 0.5],
            [0.5, 1, 0],
            [0, 1, 0.5],
            [0, 0.5, 1],
            [0.5, 0, 1],
            [1, 0.5, 0.5],
            [0.5, 1, 0.5],
            [0.5, 0.5, 1],
            [1, 1, 0.5],
            [1, 0.5, 1],
            [0.5, 1, 1]
        ]

    for n in range(0, numcells):
        L = int(xp[n] - 64)
        R = int(xp[n] + 64)
        U = int(yp[n] - 64)
        D = int(yp[n] + 64)
        colorindex = n % 18
        mag4 = np.ones((yy, xx, 3), dtype='float32');
        for m in range(0, 3):
            mag4[:,:,m] = mag4[:,:,m] * colors[colorindex][m]

        mag = morph.binary_dilation(( Ypredict[n,:,:,0]*255 ) );
        #mag = morph.dilation(( Ypredict[n,:,:,0]*255 ) );
        mag2 = np.multiply(np.ones((xx, yy), dtype='uint8')*255, mag)
        mag3 = mag2 - Ypredict[n,:,:,0]*255

        OverlayALL[ U:D , L:R ] = OverlayALL[ U:D , L:R ] + mag3

        for m in range(0, 3):
            mag4[:,:,m] = np.multiply(mag4[:,:,m], mag3)
            #Overlay[ U:D , L:R , m] = Overlay[ U:D , L:R , m] - mag3
            Overlay[ U:D , L:R , m] = Overlay[ U:D , L:R , m] + mag4[:,:,m]

        # Overlay[ U:D , L:R , 0] = Overlay[ U:D , L:R , 0] + mag3
        # Overlay[ U:D , L:R , 1] = Overlay[ U:D , L:R , 1] + mag3
        # Overlay[ U:D , L:R , 2] = Overlay[ U:D , L:R , 2] + mag3

        #mask = ndimage.binary_fill_holes(Ypredict[n,:,:,0])

    # for m in range(0, 3):
    #     frame2[:,:,m] = rframe
    # for m in range(0, 3):
    #     frame3[:,:,m] = invframe
    OverlayALL = np.invert(OverlayALL)
    OverlayALL = OverlayALL/255

    for m in range(0, 3):
        tframe[:,:,m] = np.multiply(np.float32(rframe), OverlayALL);
        iframe[:,:,m] = np.multiply(np.float32(invframe), OverlayALL);
        #tframe[:,:,m] = np.multiply(tframe[:,:,m], np.invert(tframe )/255);
        #iframe[:,:,m] = np.multiply(iframe[:,:,m], np.invert(tframe )/255);
        #tframe[:,:,m] = tframe[:,:,m] - cframe
        #iframe[:,:,m] = iframe[:,:,m] - cframe


    #Overlay = mag;
    tframe[:,:,0] = tframe[:,:,0] + Overlay[:,:,0]
    tframe[:,:,1] = tframe[:,:,1] + Overlay[:,:,1] + sframe
    tframe[:,:,2] = tframe[:,:,2] + Overlay[:,:,2] + cframe
    tframe = np.uint8(np.clip(tframe, 0, 255));

    iframe[:,:,0] = iframe[:,:,0] + Overlay[:,:,0]
    iframe[:,:,1] = iframe[:,:,1] + Overlay[:,:,1] + sframe
    iframe[:,:,2] = iframe[:,:,2] + Overlay[:,:,2] + cframe
    iframe = np.uint8(np.clip(iframe, 0, 255));

    Overlay = np.float32(Overlay)

    Overlay[:,:,1] = Overlay[:,:,1] + np.float32(sframe)
    Overlay[:,:,2] = Overlay[:,:,2] + np.float32(cframe)
    Overlay = np.uint8(np.clip(Overlay, 0, 255));

    tframe = tframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror,:]
    iframe = iframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror,:]
    rframe = rframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror]
    invframe = invframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror]
    Overlay = Overlay[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror,:]

    return [tframe, iframe, rframe, invframe, Overlay]

def generateSegmentationOverlay_dashline(rframe, cframe, sframe, Ypredict, xp, yp, border_mirror, oX, oY):

    xp = xp + border_mirror
    yp = yp + border_mirror
    starburst = skio.imread(r'A:\UNET_NucleiOverlap\UNET_Python_Code\UNET_5x5_02-10-2021\black-and-white-sunburst-pattern.tif', plugin="tifffile")
    starburst = starburst[:,:,0]
    starburst2 = skio.imread(r'A:\UNET_NucleiOverlap\UNET_Python_Code\UNET_5x5_02-10-2021\black-and-white-sunburst-pattern2.tif', plugin="tifffile")
    starburst2 = starburst2[:,:,0]
    invframe = np.invert(rframe)
    FrameSize = rframe.shape
    tileSize = Ypredict.shape
    X = FrameSize[1]
    Y = FrameSize[0]
    # oX = OFrameSize[1]
    # oY = OFrameSize[0]
    xx = tileSize[1]
    yy = tileSize[2]
    Overlay = np.zeros((Y, X, 3), dtype='uint8');
    OverlayALL = np.zeros((Y, X), dtype='uint8');
    tframe = np.zeros((Y, X, 3), dtype='float32');
    iframe = np.zeros((Y, X, 3), dtype='float32');
#    frame2 = np.zeros((Y, X, 3), dtype='uint8');
#    frame3 = np.zeros((Y, X, 3), dtype='uint8');
    numcells = len(xp)

    colors = [
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1]
        ]

    for n in range(0, numcells):
        L = int(xp[n] - 64)
        R = int(xp[n] + 64)
        U = int(yp[n] - 64)
        D = int(yp[n] + 64)
        colorindex = n % 18
        mag4 = np.ones((yy, xx, 3), dtype='float32');
        for m in range(0, 3):
            mag4[:,:,m] = mag4[:,:,m] * colors[colorindex][m]

        mag = morph.binary_dilation(( Ypredict[n,:,:,0]*255 ) );
        mag2 = np.multiply(np.ones((xx, yy), dtype='uint8')*255, mag)
        mag3 = mag2 - Ypredict[n,:,:,0]*255
        linetype = n % 3
        if(linetype == 0):
            mag3 = np.multiply(mag3, starburst)*255
        elif(linetype == 1):
            mag3 = np.multiply(mag3, starburst2)*255
        else:
            mag3 = mag3

        OverlayALL[ U:D , L:R ] = OverlayALL[ U:D , L:R ] + mag3

        for m in range(0, 3):
            mag4[:,:,m] = np.multiply(mag4[:,:,m], mag3)
            #Overlay[ U:D , L:R , m] = Overlay[ U:D , L:R , m] - mag3
            Overlay[ U:D , L:R , m] = Overlay[ U:D , L:R , m] + mag4[:,:,m]

        # Overlay[ U:D , L:R , 0] = Overlay[ U:D , L:R , 0] + mag3
        # Overlay[ U:D , L:R , 1] = Overlay[ U:D , L:R , 1] + mag3
        # Overlay[ U:D , L:R , 2] = Overlay[ U:D , L:R , 2] + mag3

        #mask = ndimage.binary_fill_holes(Ypredict[n,:,:,0])

    # for m in range(0, 3):
    #     frame2[:,:,m] = rframe
    # for m in range(0, 3):
    #     frame3[:,:,m] = invframe
    OverlayALL = np.invert(OverlayALL)
    OverlayALL = OverlayALL/255

    for m in range(0, 3):
        tframe[:,:,m] = np.multiply(np.float32(rframe), OverlayALL);
        iframe[:,:,m] = np.multiply(np.float32(invframe), OverlayALL);


    #Overlay = mag;
    tframe[:,:,0] = tframe[:,:,0]
    tframe[:,:,1] = tframe[:,:,1]
    tframe[:,:,2] = tframe[:,:,2]
    tframe = np.uint8(np.clip(tframe, 0, 255));

    iframe[:,:,0] = iframe[:,:,0]
    iframe[:,:,1] = iframe[:,:,1]
    iframe[:,:,2] = iframe[:,:,2]
    iframe = np.uint8(np.clip(iframe, 0, 255));

    Overlay = np.float32(Overlay)

    Overlay[:,:,1] = Overlay[:,:,1]
    Overlay[:,:,2] = Overlay[:,:,2]
    Overlay = np.uint8(np.clip(Overlay, 0, 255));

    tframe = tframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror,:]
    iframe = iframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror,:]
    rframe = rframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror]
    invframe = invframe[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror]
    Overlay = Overlay[border_mirror:oY+border_mirror,border_mirror:oX+border_mirror,:]

    return [tframe, iframe, rframe, invframe, Overlay]