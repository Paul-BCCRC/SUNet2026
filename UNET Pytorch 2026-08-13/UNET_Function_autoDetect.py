#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np

def autodetect_isBackgroundDark(frame):

    # borderpixels = np.ones((1,) , 'uint8' )*128;
    # borderpixels = np.append(borderpixels, frame[0,:], 0);
    # borderpixels = np.append(borderpixels, frame[:,0], 0);
    # borderpixels = np.append(borderpixels, frame[(frame.shape[0]-1),:], 0);
    # borderpixels = np.append(borderpixels, frame[:,(frame.shape[1]-1)], 0);

    # AvgBP = np.mean(borderpixels);

    hist = np.histogram(frame, bins=256);
    B = np.argmax(hist[0]);

    if(B > 128):
        isB=0;
    else:
        isB=1;

    return [isB, B]

def autodetect_isOME(imagefilename):
    isOME = 0
    UP = imagefilename.upper()
    if(len(UP) >= 8):
        if(UP[-9:].find(".OME.TIF") != -1):
            isOME = 1;
    return isOME