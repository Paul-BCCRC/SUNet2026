#Writen by Paul Gallagher (BC Cancer Research Centre)
"""
Created on Tue Oct 21 14:35:37 2025

Notes:
https://docs.python.org/3/library/struct.html

@author: x-paul
"""
import numpy as np

import struct
#import multiprocessing as mp
#from threading import Thread, current_thread
#from multiprocessing import Process, current_process


def readCBX(path, filename, ext):

    slash = '/'
    cbpath = path + slash + filename + '.' + ext


    cbfile = open(cbpath, "rb")
    cbdata = cbfile.read()

    BYTE_TOTAL = len(cbdata)
    CURRENT_BYTE = np.uint64(0);

    dummy = int.from_bytes(cbdata[CURRENT_BYTE:CURRENT_BYTE+1], byteorder='little')
    numFeats = int.from_bytes(cbdata[CURRENT_BYTE+2:CURRENT_BYTE+4], byteorder='little')
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(4);

    BYTE_FEATS = BYTE_TOTAL - (numFeats*16) - 4;
    numCells = np.uint64(np.floor((BYTE_FEATS/4)/numFeats))

    #Read feature labels
    flabels = [''] * numFeats
    for p in range(0, numFeats):
        flabels[p] = cbdata[CURRENT_BYTE:CURRENT_BYTE+15].decode("utf-8")
        flabels[p] = flabels[p].replace(" ", "")
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(16);

    fvalues = np.zeros((numCells, numFeats), dtype='float')
    for m in range(0, numCells):
        for n in range(0, numFeats):
            fvalues[m, n] = struct.unpack('<f', cbdata[CURRENT_BYTE : CURRENT_BYTE+4])[0]
            CURRENT_BYTE = CURRENT_BYTE + np.uint64(4);

    cbfile.close()

    return [flabels, fvalues]

def writeCBX(flabels, fvalues, path, filename, ext):

    numlabels = len(flabels)
    [numCell, numfeat] = fvalues.shape

    if(numlabels == numfeat):

        slash = '/'
        cbpath = path + slash + filename + '.' + ext
        cbfile = open(cbpath, "wb")

        cb_bytes = bytearray(struct.pack( '<H', np.uint16(0) ) );
        cb_bytes.extend(np.uint16(numfeat).astype('uint16').tobytes())

        for n in range(0, numfeat):
            for m in range(0, 15):
                try:
                    cb_bytes.extend(np.array([flabels[n][m]], 'S1').tobytes())
                except:
                    cb_bytes.extend(np.array([r' '], 'S1').tobytes())
            cb_bytes.extend(np.array([r','], 'S1').tobytes())

        for m in range(0, numCell):
            for n in range(0, numfeat):
                cb_bytes.extend(fvalues[m][n].astype('float32').tobytes())

        cbfile.write(cb_bytes)
        cbfile.close()
    else:
        print(r'Number of Labels and Features mismatch!')