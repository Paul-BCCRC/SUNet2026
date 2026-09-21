#Writen by Paul Gallagher (BC Cancer Research Centre)

import time
#import os
import numpy as np
import copy

#import scipy.ndimage as ndimage
import skimage.morphology
#import torch

from CytoFeatures import findOverlap, CytoFeatures, CytoLabels, getNonOverlap, makeHeaderFeatures#, ShapeFeatures_batch, FFTFeatures_batch, ODFeatures_batch
#from read_write_CCG_MultiProcessing import readCCG
#from read_write_CMG import readCMG
#from read_write_CBX import readCBX, writeCBX



# imagefilepath = r'G:\DATA\SpilloverCompensationTEST\Out'
# #imagefilepath = r'G:\DATA\SpilloverCompensationTEST\Out\SPEEDTEST'
# #imagefilename = r'FeatureTEST.cmg'
# #imagefilename = r'SpeedTest_CMG.cmg'
# imagefilename = r'Scaled_Top_Left_C23-01C1-AbII-Subtracted_crop3.cmg'
# outpath = r'G:\DATA\SpilloverCompensationTEST\Out'



# #-----Get CMG or CCG Data------------------------------------------------------
# print(r'Loading CMG/CCG Data...')
# start_load = time.time()
# imagefilename = '\\' + imagefilename

# imagefilenameSplit = imagefilename.split(".");
# isCMG = imagefilenameSplit[1][1]=='m'

# if(isCMG):
#     [Images, Masks, Header] = readCMG(imagefilepath, imagefilenameSplit[0])
# else:
#     [Images, Masks, Header, ImageLabel, MaskLabel, ilist] = readCCG(imagefilepath, imagefilenameSplit[0])
#     del ilist
def SSC(Images, Masks, Header, isCMG, ImageLabel, MaskLabel, isCBfile, in_flabels, in_fvalues, SSCweight):
    start_ = time.time()
    if(isCMG):
        ImageLabel = [r''] * (Header[1][0])
        for i in range(0, Header[1][0]):
            ImageLabel[i] = str(i+1)
        MaskLabel = [r''] * (Header[26][0])
        for i in range(0, Header[26][0]):
            MaskLabel[i] = str(i+1)


    #------------------------------------------------------------------------------

    N = len(Masks)              #Number of Cells
    X = Images[0].shape[1]      #Width
    Y = Images[0].shape[0]      #Height
    P = Images[0].shape[2]      #Number of Image Planes
    M = Masks[0].shape[2]       #Number of Mask Planes

    #Number of batches. Increase if not enough VRAM for all the cells at the same time
    B = int(np.ceil(N/10000))

    #Reshape Image and Mask data from lists to arrays
    ImageTensor = np.zeros((N, Y, X, P), dtype='uint8')
    MaskTensor = np.zeros((N, Y, X, M), dtype='uint8')
    for n in range(0, N):
        ImageTensor[n, :, :, :] = Images[n]
        MaskTensor[n, :, :, :] = Masks[n]

    #Calculate which cells interact (ilist) and by how much (slist)
    print(r'Finding cell interactions...')
    start_itract = time.time()
    ilist = [[]]*3
    slist = [[]]*3

    #Nuclei interations
    [ilist[0], slist[0]] = findOverlap(Masks, Header, 0)

    #struct1 = ndimage.generate_binary_structure(2, 1)
    struct1 = skimage.morphology.disk(3, dtype=np.uint8)

    dMask = np.zeros((N,Y,X,1),dtype='uint8')
    cMask = np.zeros((N,Y,X,1),dtype='uint8')

    for n in range(0, N):
        #dMask[n,:,:,0] = ndimage.binary_dilation(MaskTensor[n,:,:,0], structure=struct1,iterations=3).astype(MaskTensor[n,:,:,0].dtype)
        dMask[n,:,:,0] = skimage.morphology.dilation(MaskTensor[n,:,:,0], footprint=struct1)
        cMask[n,:,:,0] = MaskTensor[n,:,:,1] - MaskTensor[n,:,:,0]
        cMask[cMask > 1] = 0

    #Annulus interations
    [ilist[1], slist[1]] = findOverlap(dMask, Header, 0)
    #Cyto-Nuclei interations
    [ilist[2], slist[2]] = findOverlap(cMask, Header, 0)

    del dMask
    del cMask
    del Images
    del Masks
    end_itract = time.time()
    print(r'Time to find cell interactions: ' + str(end_itract - start_itract))

    # #-----Get existing CB Feature Data if it exists--------------------------------
    # print(r'Checking for existing CB file...')
    # CB_file_path = imagefilepath + '\\' + imagefilenameSplit[0] + r'.cb4'
    # isCBfile = os.path.exists(CB_file_path)
    # if(isCBfile):
    #     [in_flabels, in_fvalues] = readCBX(imagefilepath, imagefilenameSplit[0], r'cb4')
    # #------------------------------------------------------------------------------



    #----Use GPU to get features---------------------------------------------------
    print(r'Calculating features...')
    start_calcfeat = time.time()


    #[testreturn] = ShapeFeatures_batch(MaskTensor, 1)
    #[testreturn] = FFTFeatures_batch(MaskTensor, 64, 1)
    #[testreturn] = ODFeatures_batch(ImageTensor, MaskTensor, 1)

    [OMasks, nOMasks] = getNonOverlap(MaskTensor, Header, 1, ilist)
    #nOMasks=-1
    out_flabels = CytoLabels(ImageLabel)
    [iMean, iTotal, mArea, out_fvalues] = CytoFeatures(ImageTensor, MaskTensor, OMasks, nOMasks, B)

    #------------------------------------------------------------------------------

    print(r'Calculating spillover adjustments...')
    I_adjusted_m = copy.deepcopy(iMean)
    I_adjusted_t = copy.deepcopy(iTotal)
    numInteractions = np.zeros((N, 3), dtype='float32')
    sumFractions = np.zeros((N, 3), dtype='float32')

    for n in range(0, N):
        for m in range(0, 3):
            mm = 2*m

            NN = len(ilist[m][n])
            numInteractions[n][m] = NN

            if(NN > 0):
                frac_ThisCell = np.zeros(NN, dtype='float32')
                frac_OtherCell = np.zeros(NN, dtype='float32')

                #ThisCell_ilist = copy.deepcopy(np.uint64(ilist[m][n]))

                for nn in range(0, NN):
                    frac_ThisCell[nn] = slist[m][n][nn] / mArea[n][mm]
                    frac_OtherCell[nn] = slist[m][n][nn] / mArea[ np.uint64(ilist[m][n][nn]-1) ][ mm ]

                    top_m = I_adjusted_m[n, :, mm] - ( SSCweight * (frac_OtherCell[nn]) * iMean[ np.uint64(ilist[m][n][nn]-1), :, mm ] )
                    top_t = I_adjusted_t[n, :, mm] - ( SSCweight * (frac_OtherCell[nn]) * iTotal[ np.uint64(ilist[m][n][nn]-1), :, mm ] )
                    bottom = 1 - ( frac_OtherCell[nn] * frac_ThisCell[nn] )

                    adj_m = top_m / bottom
                    adj_t = top_t / bottom

                    #Set any negative values to zero.
                    adj_m[adj_m < 0] = 0
                    adj_t[adj_t < 0] = 0

                    I_adjusted_m[n, :, mm] = adj_m
                    I_adjusted_t[n, :, mm] = adj_t

                    sumFractions[n,m] = sumFractions[n,m] + frac_ThisCell[nn] #+ frac_OtherCell[nn]



    end_calcfeat = time.time()
    print(r'Time to calculate all features: ' + str(end_calcfeat - start_calcfeat))

    #DIFF = np.int64(I_adjusted) - np.int64(I_initSum)

    out_fvalues[:, 8] = numInteractions[:, 0]
    out_fvalues[:, 9] = numInteractions[:, 2]
    out_fvalues[:, 10] = sumFractions[:, 0]
    out_fvalues[:, 11] = sumFractions[:, 2]
    count= (7*3)*P + 66
    for p in range(0, P):
        out_fvalues[:, count+p] = I_adjusted_m[:, p, 0]     #SPNm_m
    count+=P
    for p in range(0, P):
        out_fvalues[:, count+p] = I_adjusted_m[:, p, 2]     #SPAnd_m
    count+=P
    for p in range(0, P):
        out_fvalues[:, count+p] = I_adjusted_m[:, p, 4]     #SPCmNm_m
    count+=P
    for p in range(0, P):
        out_fvalues[:, count+p] = I_adjusted_t[:, p, 0]     #SPNm_t
    count+=P
    for p in range(0, P):
        out_fvalues[:, count+p] = I_adjusted_t[:, p, 2]     #SPAnd_t
    count+=P
    for p in range(0, P):
        out_fvalues[:, count+p] = I_adjusted_t[:, p, 4]     #SPCmNm_t
    count+=P

    [header_flabels, header_fvalues] = makeHeaderFeatures(Header)
    out_flabels = header_flabels + out_flabels
    out_fvalues = np.concatenate( (header_fvalues, out_fvalues), axis=1 )

    # if(isCBfile):
    #     print(r'Appending CB file...')
    #     #append to existing features
    #     flabels = in_flabels + out_flabels
    #     fvalues = np.concatenate( (in_fvalues, out_fvalues), axis=1 )
    # else:
    #     print(r'Creating CB file...')
    #     #create all new label/values
    #     flabels = out_flabels
    #     fvalues = out_fvalues

    # writeCBX(flabels, fvalues, outpath, imagefilenameSplit[0], r'cb6')

    end_ = time.time()
    print('SSC Done: ' + str(end_ - start_))
    print('---------------------------------')

    return [out_flabels, out_fvalues]

def makeCMGHeaderStructure(Images,Masks,Screenx,Screeny,Resolution,Objective,CassettePosition,vorx,vory,BackgroundFloat):
    Header = -1

    NumberOfImage = len(Images)

    Width = np.zeros(NumberOfImage,dtype='uint32')
    Height = np.zeros(NumberOfImage,dtype='uint32')
    NbColorMap = np.zeros(NumberOfImage,dtype='uint8')
    NbBitMap = np.zeros(NumberOfImage,dtype='uint8')

    for n in range(0, NumberOfImage):
        Height[n] = Images[n].shape[0]
        Width[n] = Images[n].shape[1]

        ImagesLen = len(Images[n].shape)
        MasksLen = len(Masks[n].shape)
        if(ImagesLen==2):
            numImage = 1;
        elif(ImagesLen==3):
            numImage = Images[n].shape[2];

        if(MasksLen==2):
            numMask = 1;
        elif(MasksLen==3):
            numMask = Masks[n].shape[2];

        NbColorMap[n] = numImage
        NbBitMap[n] = numMask


    Mode = np.ones(NumberOfImage,dtype='uint8')*99
    Class = np.zeros(NumberOfImage,dtype='uint32')

    # if(len(Screenx)!=NumberOfImage):
    #     Screenx = np.zeros(NumberOfImage,dtype='uint32')

    # if(len(Screeny)!=NumberOfImage):
    #     Screeny = np.zeros(NumberOfImage,dtype='uint32')

    Stagex = np.zeros(NumberOfImage,dtype='uint64')
    Stagey = np.zeros(NumberOfImage,dtype='uint64')
    Stagez = np.zeros(NumberOfImage,dtype='uint64')

    # if(len(Resolution)!=NumberOfImage):
    #     Resolution = np.ones(NumberOfImage,dtype='float32') * 0.18

    LowThreshold = np.ones(NumberOfImage,dtype='uint16') * 9728
    MidThreshold = np.ones(NumberOfImage,dtype='uint16') * 4864
    Group = np.zeros(NumberOfImage,dtype='uint8')
    Accession = np.zeros(NumberOfImage,dtype='uint32')
    for x in range(1, NumberOfImage):
        Accession[x:] = Accession[x:] + np.ones((NumberOfImage-x),dtype='uint32')

    Iod = np.ones(NumberOfImage,dtype='float32') * 22500
    Fluor = np.zeros(NumberOfImage,dtype='uint8')
    Diagnosis = np.zeros(NumberOfImage,dtype='uint16')
    RedFaction = np.ones(NumberOfImage,dtype='float32') *0#* 0.3
    GreenFaction = np.ones(NumberOfImage,dtype='float32') *0#* 0.7
    BlueFaction = np.zeros(NumberOfImage,dtype='float32') *0#
    Index = np.zeros(NumberOfImage,dtype='uint32')

    # if(len(Objective)!=NumberOfImage):
    #     Objective = np.ones(NumberOfImage,dtype='uint32') * 538980404

    Calibrated = np.zeros(NumberOfImage,dtype='uint8')
    StackX_int = np.zeros(NumberOfImage,dtype='uint32')
    StackY_int = np.zeros(NumberOfImage,dtype='uint32')

    # if(len(CassettePosition)!=NumberOfImage):
    #     CassettePosition = np.zeros(NumberOfImage,dtype='uint8')

    # if(len(vorx)!=NumberOfImage):
    #     vorx = np.zeros(NumberOfImage,dtype='uint32')

    # if(len(vory)!=NumberOfImage):
    #     vory = np.zeros(NumberOfImage,dtype='uint32')


    BestFocusFrame = np.zeros(NumberOfImage,dtype='uint8')

    # if(len(BackgroundFloat)!=NumberOfImage):
    #     BackgroundFloat = np.ones(NumberOfImage,dtype='float32') * 254

    PrimaryColourChannel = np.zeros(NumberOfImage,dtype='uint8')
    Layer = np.zeros(NumberOfImage,dtype='uint16')
    Points = np.zeros((NumberOfImage, 10),dtype='uint8')
    NumFeature = np.ones(NumberOfImage,dtype='uint8') * 3
    RGB_Order = np.zeros((NumberOfImage, 3),dtype='uint8')
    RGB_Order[:,0] = np.zeros((NumberOfImage),dtype='uint8')
    RGB_Order[:,1] = np.ones((NumberOfImage),dtype='uint8') * 1
    RGB_Order[:,2] = np.ones((NumberOfImage),dtype='uint8') * 2

    Header = [Mode, NbColorMap, Class, Screenx, Screeny, Stagex, Stagey, Stagez, Resolution, LowThreshold, MidThreshold, Group, Width, Height, Accession, Iod, Fluor, Diagnosis, RedFaction, GreenFaction, BlueFaction, Index, Objective, Calibrated, StackX_int, StackY_int, NbBitMap, CassettePosition, vorx, vory, BestFocusFrame, BackgroundFloat, PrimaryColourChannel, Layer, Points, NumFeature, RGB_Order]

    return Header