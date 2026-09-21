# PHSA NON-COMMERCIAL SOURCE-AVAILABLE LICENCE
# Draft for PHSA review | Prepared September 9, 2026

# Software: SUnet version 0.9
# Copyright (c) 2026 Provincial Health Services Authority.
# All rights reserved.
# Licensing contact: PHSA Technology Development Office (TDO)
# tdoadmin@phsa.ca

# 1. SCOPE AND ACCEPTANCE

# This licence applies to the PHSA-owned source code, executable code, and
# documentation identified above (the "Software"). Provincial Health Services
# Authority ("PHSA") grants only the permissions expressly stated here. By using,
# copying, or modifying the Software, you accept these terms. If acting for an
# organization, you must have authority to accept them on its behalf.

# 2. PERMITTED NON-COMMERCIAL USE

# Subject to this licence, PHSA grants a non-exclusive, royalty-free,
# non-transferable copyright licence to use, copy, and modify the Software:

# (a) for internal non-commercial research and education by non-profit
#     institutions and Canadian federal, provincial, territorial, and local
#     government bodies; or
# (b) for internal, non-commercial evaluation by other organizations, solely
#     to determine whether to seek a separate licence from PHSA.

# "Non-commercial" excludes use for commercial advantage, monetary compensation,
# product or service development for commercial exploitation, paid services,
# consulting, or work performed for or on behalf of a commercial entity. An
# organization's non-profit status does not itself make a use non-commercial.

# Any other use requires a separate written licence from PHSA through TDO.

# 3. COPYING, MODIFICATION, AND REDISTRIBUTION

# Copies and modifications must remain within the licensed organization and be
# used only as permitted above. Preserve this licence, copyright notices, and
# attributions in every copy; clearly identify modifications and their dates.
# Do not sell, sublicense, publish, redistribute, or provide third-party access
# to the Software or modifications without PHSA's prior written permission.

# 4. AUTOMATED INGESTION AND AI

# Without PHSA's prior written permission, do not submit the Software or its
# modifications to AI or machine-learning training systems, code-ingestion
# services, crawlers, indexing services, or text or data mining systems. This
# restriction concerns ingestion of the Software itself; it does not prohibit
# ordinary compilation, testing, or execution for a use permitted by section 2.

# 5. RESERVED RIGHTS AND THIRD-PARTY MATERIALS

# PHSA retains ownership of the Software. No patent or trademark licence is
# granted. If a proposed use requires patent rights, obtain a separate written
# licence from PHSA before that use. Separately identified third-party materials
# remain governed by their own licences; this licence does not replace them.
# Do not imply PHSA or BC Cancer endorsement or use their logos without written
# permission.

# 6. RESEARCH SOFTWARE; WARRANTY AND LIABILITY

# The Software is provided "AS IS", without warranties, including merchantability,
# fitness for a particular purpose, accuracy, or non-infringement. PHSA has no
# obligation to support, maintain, or update it. This licence does not authorize
# clinical diagnosis, treatment, or patient-care decisions.

# To the extent permitted by applicable law, PHSA and its personnel are not
# liable for losses or damages arising from the Software or its use, including
# lost data, lost profits, or direct, indirect, incidental, or consequential
# damages, even if advised of their possibility.

# 7. TERMINATION AND GOVERNING LAW

# Your permissions terminate upon breach. You must then stop using the Software
# and delete your copies and modifications. Sections concerning reserved rights,
# warranties, liability, and governing law survive termination. This licence is
# governed by the laws of British Columbia and the applicable laws of Canada.
# Revised terms apply only to releases supplied under those revised terms.

# END OF LICENCE


#Writen by Paul Gallagher (BC Cancer Research Centre)

#-----Initialize/Import--------------------------------------------------------
#from UNET_Function_mirrorBorder import MirrorBorder
from UNET_Function_subSampleImage import subSampleImage
#from UNET_Function_generateOverlapImage import generateOverlapImage
from UNET_Function_generateOverlapImage import generateOverlapImage_WithCenter
from UNET_Function_getCenterPoints import getCenterPoints#, getCenterPoints_TopHat
from UNET_Function_getCenterCrops import getCenterCrops
from UNET_Function_generateSegmentationOverlay import generateSegmentationOverlay
from UNET_Function_cropMirrored_getScreenPos import cropMirrored_getScreenPos
from UNET_Function_cleanMasks_MP import cleanMasks_MP
from UNET_Function_autoDetect import autodetect_isOME, autodetect_isBackgroundDark
from UNET_Function_proximityFilter import  proximityFilter_MP#, proximityFilter
from UNET_Function_readOMEmetadata import readOMEmetadata

import os
import numpy as np
import time
import copy
import sys
#import math

from UNET_Structures import UNet_A, UNet_B

from read_write_CMG import writeCMG, CMG_header
#from read_write_CCG import writeCCG, CCG_header
#from read_write_CCG_bytearray import writeCCG, CCG_header
from read_write_CCG_MultiProcessing import writeCCG#, CCG_header

from PIL import Image
from ome_types import from_tiff
#from libtiff import TIFF
#import matplotlib.pyplot as plt
#import skimage.io as skio
#import scipy.ndimage as ndimage
import skimage.morphology

import pyvips

import torch
#from torch import nn
#from torch import optim#, nn
from tqdm import tqdm



#------------------------------------------------------------------------------
def UNET_Segmentation_Function(imagefilepath, imagefilename):
    #start_ = time.time()
    start_init = time.time()


    #-----Run Parameters-------------------------------------------------------

    #General Parameters
    display_times = 1;                          #Print out display times in seconds
    append_centers = 0;                         #Append cell centers masks to prediction masks for write out.
    TIFF_page_number = 1;                       #Page number of input image to segment, should be Nuclei image.
    #Include_Other_Pages = 1;                   #Save other page image data in cmg/ccg

    #Input Image Parameters
    maxImageLength = 5000;                      #Max length for subimage
    isOME = autodetect_isOME(imagefilename);
    is12bit = 0                                 #Correct for 12 bit data within 16 bit range
    overlap = 64;                               #Number of pixels to overlap tiles for netowrk A.
    border_mirror = 0;                          #Number of pixels to extend/mirror input frame for all edges.

    #Output Parameters
    outpath = imagefilepath;                    #Set path to create folder containing output files.
    CMG_out = 1                                 #Save total result as a CMG
    CCG_out = 0                                 #Save total result as a CCG
    CCG_out_I = 0                               #Calculate mask interactions for CCG_out
    CMG_subimg_out = 1                          #Save subimage result as CMG
    CCG_subimg_out = 0                          #Save subimage result as CCG
    CCG_subimg_out_I = 0                        #Calculate mask interactions for CCG__subimg_out

    #Background Parameters
    isAutoDetectBackground = 0;                 #Attempt to determine background brightness.
    isNucleiBackgroundDark = 1;
    isBackgroundDark =      [1, 1, 1, 1, 1, 1, 1, 1]              #Background flag. Set 1 if the background of input image is dark, 0 if bright. Must be same length as channel_map

    #Center Finding Parameters
    threshold = 6                               #
    SaturationOffset = 15                       #
    R = 7                                       #

    #Mask Clean Up Parameters
    nucleiMaskThresh = 0.5;                     #Pixel prediction values above this are considered positive.
    cytoMaskThresh = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
    num_CleanProccess = 0                       #Number of Processes for mask clean. (0="auto" Scales based on number of cells)

    #Proximity Filter Parameters
    doProxFilt = 1;                             #Perform proximity filter
    prox_distance = 20;                         #Nuclei below this distance in pixels are concidered for removal
    iThresh = 0.85;                             #Removal Threshold
    num_ProxProccess = 0                        #Number of Processes for prox. filt. (0="auto" Scales based on number of cells)

    #Cytoplasm Segmentation Parameters
    doCytoSegmention = 1;
    channel_map =           [5, 6, 7, 9, 11, 14, 15, 17]           #Channels/pages to use
    OrPriorityGrouping =    [1, 1, 1, 1, 1, 1, 1, 1]
    dilation_radius = 3;                        #Dilates this many pixels around mask. Nuclei for minimum cytoplasm prediction (Dilated nuclei Or with cyto predictions)

    #Spatial Spillover Compensation Parameters
    doSSC = 1;                                  #Call Spatial Spillover Compensation
    SSCweight = 0.75;
    #--------------------------------------------------------------------------

    if(doSSC):
        from SpatialSpilloverCompensation import SSC, makeCMGHeaderStructure
        #from CytoFeatures import makeHeaderFeatures
        from read_write_CBX import writeCBX

    numChannelMap = len(channel_map)
    maxPriority = np.max(OrPriorityGrouping)

    #-----Network Dimensions----------------------------------------------
    A_X = 128                                 #Input size of network A.
    A_Y = 128
    B_X = 128                                 #Input size of network B.
    B_Y = 128
    A_predictbatch = 100                      #Set how many samples per prediction iteration for both networks.
    B_predictbatch = 100
    Nfeatures = 32                            #Number of feature maps for first level of network. This number is multiplyed by 2, 4, 8, 16 deeper into UNET.
    #--------------------------------------------------------------------------


    #-----Weight Paths---------------------------------------------------------
    #Nuclei
    wpathA = r'D:\FromCalum\SegmentationWeights\UNET_ARound4Tweaked.pth'                    #Nuclei centre weights
    #wpathA = r'D:\FromCalum\SegmentationWeights\UNET_A_Thionin.pth'
    #wpathB = r'D:\FromCalum\SegmentationWeights\UNET_B_Thionin.pth'                         #Nuclei segmentation weights
    wpathB = r'D:\FromCalum\SegmentationWeights\UNET_BRound4Tweaked.pth'

    #Cytoplasm
    wpathC = [r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',  #Cyto weights. Weights in this list should align with channels in "channel_map"
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth',
              r'D:\FromCalum\SegmentationWeights\Segmentation-Cytoplasm_model_weights.pth'
              ]
    #--------------------------------------------------------------------------


    #-----CCG Labels-----------------------------------------------------------
    iLabels = [r'DAPI',                   #1
               r'TRITC',                  #2
               r'Cy5',                    #3
               r'CD3 - TRITC',            #4
               r'FOXP3 - CY5',            #5
               r'CD4 - TRITC',            #6
               r'CD8 - Cy5',              #7
               r'CD20 - TRITC',           #8
               r'CD56 - Cy5',             #9
               r'CD11c - TRITC',          #10
               r'CD68 - Cy5',             #11
               r'SMA - TRITC',            #12
               r'PD-L1 - Cy5',            #13
               r'CD45 - TRITC',           #14
               r'PD-1 - Cy5',             #15
               r'CK_TRITC - TRITC',       #16
               r'Ki-67 - Cy5']            #17

    #iLabels = [r'Thionin']

    mLabels = [r'Nuclei',
               r'Cyto',
               r'Center']
    #--------------------------------------------------------------------------

    device = (
        "cuda"
        if torch.cuda.is_available()
        else "mps"
        if torch.backends.mps.is_available()
        else "cpu"
    )
    torch.cuda.empty_cache()                  #Clear GPU memory from possible previous runs
    end_init = time.time()
    if(display_times):
        print('Time to Initialize: ' + str(end_init - start_init))



    numargs = len(sys.argv)
    print('-----------------------')
    #print('Number of arguments:', numargs, 'arguments.')
    #print( 'Argument List:', str(sys.argv))
    print( 'File Path:', str(imagefilepath))
    print( 'File Name:', str(imagefilename))
    print('-----------------------')
    print('Segmentating...', str(imagefilepath), str(imagefilename))
    #imagefilepath=r'A:\UNET_NucleiOverlap\UNET_SCANTESTS\LCR_Thionin_Unmixed'
    #imagefilename=r'\Scaled_InvertedConcentrations_row_5_col_9'
    #imagefilename=r'\HE1'

    if(numargs > 1):
        imagefilepath=sys.argv[1]
        imagefilename=sys.argv[2]



    #Configure networks structure, compile, and load wieghts.
    start_compile = time.time()
    #-----Compile network------------------------------------------------------
    Amodel = UNet_A(Nfeatures).to(device)
    Bmodel = UNet_B(Nfeatures).to(device)
    #--------------------------------------------------------------------------
    end_compile = time.time()
    if(display_times):
        print(r'Time to Construct UNETS: ' + str(end_compile - start_compile))




    #-----Load Input Image-----------------------------------------------------
    #Read in entire input images
    #frame = skio.imread(os.path.join(imagefilepath+imagefilename))

    image = pyvips.Image.new_from_file(os.path.join(imagefilepath+imagefilename), access="random")
    page_num = image.get("n-pages")
    oX = image.width
    oY = image.height

    pages = [pyvips.Image.new_from_file(os.path.join(imagefilepath+imagefilename), page=i, access="random")
             for i in range(0, page_num)]
    image = pages[0].bandjoin(pages[1:])

    #Extract OME metadata
    if(isOME):
        ome = from_tiff(os.path.join(imagefilepath+imagefilename))
        [Objective, Resolution, ImageType, ScaleF, iLabels] = readOMEmetadata(ome, TIFF_page_number)
    #--------------------------------------------------------------------------

    max_A_X = int( np.floor( maxImageLength / A_X ))
    max_A_Y = int( np.floor( maxImageLength / A_Y ))

    oAX = int( np.ceil( oX / A_X ))
    oAY = int( np.ceil( oY / A_Y ))

    num_sX = int( np.ceil( oAX / (max_A_X) ) )
    num_sY = int( np.ceil( oAY / (max_A_Y) ) )

    num_XX = int( np.ceil( oAX / num_sX ))
    num_YY = int( np.ceil( oAY / num_sY ))

    sX = A_X*num_XX
    sY = A_Y*num_YY

    numSubImage = num_sX * num_sY
    count = 0;

    isFirstTile = 1
    pY = 0
    for yy in range(0, num_sY):
        pX = 0
        if(yy!=0):
            pY = pY + (sY - B_Y)
        for xx in range(0, num_sX):
            print(r'-------' + str(int(count*100/numSubImage)) + r'%----------------------------------------')
            print(r'Segmenting...(' + str(xx) + r'/' + str(num_sX) + r', ' + str(yy) + r'/' + str(num_sY) + r')' )

            start_read = time.time()
            if(xx!=0):
                pX = pX + (sX - B_X)


            tX=sX
            tY=sY
            if(xx == (num_sX-1)):
                tX=(oX-(num_sX-1)*(sX-B_X))
            if(yy == (num_sY-1)):
                tY=(oY-(num_sY-1)*(sY-B_Y))


            cropROI = image.crop(pX, pY, tX, tY)



            frame = cropROI.numpy()


            # if(page_num > 1):
            #     frame = frame[:,:,TIFF_page_number-1]
            del cropROI
            FrameSize = frame.shape

            if(len(FrameSize)==2):
                frame = frame[:,:,np.newaxis]
                # tempframe = np.zeros((FrameSize[0], FrameSize[1], 1), frame.dtype)
                # tempframe[:,:,0] = frame
                # frame = tempframe
                # del tempframe

            #If image is 16 bit convert to 8 bit.
            if(frame.dtype == 'uint16'):
                if(is12bit):
                    frame = frame*16
                frame = np.uint8(frame/256);

            # #Get dimentions of input image.
            # FrameSize = frame.shape

            if(isAutoDetectBackground):
                [isNucleiBackgroundDark, background_value] = autodetect_isBackgroundDark( frame[:,:,TIFF_page_number-1] );
                if(doCytoSegmention):
                    for n in range(0, len(channel_map)):
                        [isBackgroundDark[n], background_value] = autodetect_isBackgroundDark( frame[:,:,(channel_map[n]-1)] );


            end_read = time.time()
            if(display_times):
                print(r'Time to Read TIFF: ' + str(end_read - start_read))
            #--------------------------------------------------------------------------



            #-----Subsample Input Image------------------------------------------------
            start_sub = time.time()
            #Tiles mirrored input image for cell center location.
            [XpredictA, numX, numY, resX, resY] = subSampleImage(frame[:,:,TIFF_page_number-1], A_X, overlap, isNucleiBackgroundDark)

            end_sub = time.time()
            if(display_times):
                print(r'Time to Subsample: ' + str(end_sub - start_sub))
            mFrameSize = frame.shape
            mX = mFrameSize[1]
            mY = mFrameSize[0]
            #del Mframe
            #--------------------------------------------------------------------------



            #-----Predict Nuclei Centers ('A')-----------------------------------------
            YpredictA, YpredictATime = makePrediction(XpredictA, device, Amodel, wpathA, A_predictbatch)
            if(display_times):
                print(r'Prediction Time(Cell Centers): ' + YpredictATime)
            #--------------------------------------------------------------------------



            #-----Reconstruct input image with centers---------------------------------
            start_genimage = time.time()
            #Reconstruct input image with prediction tiles. rframe = full size raw image.
            #rframe = GenerateOverlapImage(XpredictA, A_X, overlap, numX, numY, mX, mY, isBackgroundDark)

            #Construct full size cell center prediction image. cframe = full size centres prediction image.
            cframe = generateOverlapImage_WithCenter(YpredictA*255, A_X, overlap, numX, numY, mX, mY)
            end_genimage = time.time()
            if(display_times):
                print(r'Time to recreate input and center images: ' + str(end_genimage - start_genimage))
            del XpredictA, YpredictA
            #--------------------------------------------------------------------------



            #-----Locate centerpoints--------------------------------------------------
            start_centerpoints = time.time()
            #Find center pixel cordinates of each prediction using "Top hat". sframe = full size pixel centres image.
            [xp, yp, sframe] = getCenterPoints(cframe, frame[:,:,TIFF_page_number-1], border_mirror, tX, tY, threshold, SaturationOffset, R, isBackgroundDark[0])
            end_centerpoints = time.time()
            if(display_times):
                print(r'Time to locate center points: ' + str(end_centerpoints - start_centerpoints))
            #--------------------------------------------------------------------------

            if(len(xp)>=1):

                #-----Crop centerpoints----------------------------------------------------
                start_centercrops = time.time()
                #Crop rframe at each center point in sframe for segmentation prediction
                [XpredictB, CenterCrops] = getCenterCrops(frame, sframe, xp, yp, B_X, B_Y, isNucleiBackgroundDark)
                end_centercrops = time.time()
                if(display_times):
                    print(r'Time to locate center crops: ' + str(end_centercrops - start_centercrops))
                #--------------------------------------------------------------------------

                #img = Image.fromarray(RGB)
                #img.save(outpath + imagefilename + '_SelectedCenters.TIFF')


                #-----Predict Nuclei Segmentation ('B')------------------------------------
                YpredictB, YpredictBTime = makePrediction(XpredictB[:,:,:,[TIFF_page_number-1]], device, Bmodel, wpathB, B_predictbatch)
                if(display_times):
                    print(r'Prediction Time(Cell Segmentation): ' + YpredictBTime)
                #--------------------------------------------------------------------------



                #-----Clean Predictions----------------------------------------------------
                start_clean = time.time()
                #Modify mask predictions.
                #-Set 0 or 1 threshold
                #-Remove non connected pixels
                #-Fill any holes
                #-Remove nuclei with borders within 2 pixels of border.
                #-Remove nuclei within mirrored area.
                #cleanMaskTimes = np.zeros((5), 'float');
                #XpredictB, YpredictB, CenterCrops, xp, yp, cleanMaskTimes = cleanMasks(XpredictB, YpredictB, CenterCrops, xp, yp, frame, mask_Thresh, border_mirror, cleanMaskTimes)

                # start_clean = time.time()
                # XpredictB, YpredictB, CenterCrops, xp, yp = cleanMasks(1, mX, mY, mask_Thresh, border_mirror, XpredictB, YpredictB, CenterCrops, xp, yp)
                # end_clean = time.time()
                # print(r'clean Time: ' + str(end_clean - start_clean))

                start_cleanMP = time.time()
                print(r'Fixing Predictions...')
                XpredictB, YpredictB, CenterCrops, xp, yp = cleanMasks_MP(1, mX, mY, nucleiMaskThresh, border_mirror, XpredictB, YpredictB, CenterCrops, xp, yp, num_CleanProccess)
                end_cleanMP = time.time()
                print(r'cleanMP Time: ' + str(end_cleanMP - start_cleanMP))


                end_clean = time.time()
                N, H, W, C = XpredictB.shape

                if(display_times):
                    print(r'Time to Clean Masks: ' + str(end_clean - start_clean))
                #     print(r'    Threashold: ' + str(cleanMaskTimes[0]))
                #     print(r'    Remove non connected areas: ' + str(cleanMaskTimes[1]))
                #     print(r'    Fill holes: ' + str(cleanMaskTimes[2]))
                #     print(r'    Cells borders : ' + str(cleanMaskTimes[3]))
                #     print(r'    Remove cells exceeding input image: ' + str(cleanMaskTimes[4]))
                #--------------------------------------------------------------------------

                if(YpredictB.shape[0] >= 1):

                    #-----Proximity Filter-----------------------------------------------------

                    if(doProxFilt):
                        start_proxMP = time.time()
                        print(r'Proximity...')

                        [XpredictB, YpredictB, CenterCrops, xp, yp] = proximityFilter_MP(XpredictB, YpredictB, CenterCrops, xp, yp, B_X, B_Y, prox_distance, iThresh, num_ProxProccess)
                        #[XpredictB, YpredictB, CenterCrops, xp, yp] = proximityFilter(XpredictB, YpredictB, CenterCrops, xp, yp, B_X, B_Y, prox_distance, iThresh)
                        N, H, W, C = XpredictB.shape
                        end_proxMP = time.time()
                        print(r'proxMP Time: ' + str(end_proxMP - start_proxMP))
                    #--------------------------------------------------------------------------


                    #-----Output Segmentation Overlay Images-----------------------------------
                    start_gensegimage = time.time()
                    CC = imagefilename

                    outfilenameSplit = CC.split(".");
                    outfilename = outfilenameSplit[0];

                    # for i in range(0, len(outfilenameSplit)-2):
                    #     outfilename = outfilename + r'.' + outfilenameSplit[i+1]
                    if(numSubImage>1):
                        outfilenametile = outfilenameSplit[0] + r'(' + str(xx) + r', ' + str(yy) + r')'
                    else:
                        outfilenametile = outfilenameSplit[0] + r'.' + outfilenameSplit[1]

                    outpath = imagefilepath + outfilename + r'\\'
                    if not os.path.exists(outpath):
                        os.mkdir(outpath)

                    if(CMG_subimg_out):
                        if not os.path.exists(outpath+r'CMG\\'):
                            os.mkdir(outpath+r'CMG\\')
                    if(CCG_subimg_out):
                        if not os.path.exists(outpath+r'CCG\\'):
                            os.mkdir(outpath+r'CCG\\')

                    #Create overlay image with segmentation and center predictions.
                    [segOverlay, isegOverlay, segRaw, isegRaw, OverlayMask] = generateSegmentationOverlay(frame[:,:,TIFF_page_number-1], cframe, sframe, YpredictB, xp, yp, border_mirror, tX, tY);
                    del frame, cframe, sframe

                    # if(isBackgroundDark[TIFF_page_number-1]):
                    savepath = outpath + outfilenametile + r'_SegmentationMap.TIFF'
                    seg = []
                    seg.append(Image.fromarray(segRaw))
                    seg.append(Image.fromarray(segOverlay[:,:,:],'RGB'))
                    seg.append(Image.fromarray(OverlayMask[:,:,:],'RGB'))
                    seg[0].save(savepath, save_all=True, append_images=seg[1:])
                    # else:
                    #     isavepath = outpath + outfilenametile + r'_InvertedSegmentationMap.TIFF'
                    #     iseg = []
                    #     iseg.append(Image.fromarray(isegRaw))
                    #     iseg.append(Image.fromarray(isegOverlay[:,:,:],'RGB'))
                    #     iseg.append(Image.fromarray(OverlayMask[:,:,:],'RGB'))
                    #     iseg[0].save(isavepath, save_all=True, append_images=iseg[1:])

                    del isegRaw, isegOverlay, OverlayMask, segRaw, segOverlay
                    #write_tiffs2(X, Y, XpredictRED, YpredictRED, folder=outpath)
                    end_gensegimage = time.time()
                    if(display_times):
                        print(r'Time to Generate Segmentation Image: ' + str(end_gensegimage - start_gensegimage))
                    #--------------------------------------------------------------------------


                    #-----Cytoplasm segmentation-----------------------------------------------

                    if(doCytoSegmention):

                        if(append_centers):
                            tMasks = np.zeros((N, H, W, 3), dtype='uint8')
                            tMasks[:,:,:,-1:] = CenterCrops[:,:,:,0]
                        else:
                            tMasks = np.zeros((N, H, W, 2), dtype='uint8')
                        tMasks[:,:,:,0] = YpredictB[:,:,:,0]

                        for ii in range(1, numChannelMap):
                            if(isBackgroundDark[ii]==0):
                                for n in range(0, N):
                                    XpredictB[n,:,:,(channel_map[ii]-1) ] = np.invert(XpredictB[n,:,:,(channel_map[ii]-1) ]);

                        activeList = np.ones(N, dtype='uint8')
                        for p in range(0, maxPriority):
                            cPriority = p+1

                            for ii in range(1, numChannelMap):

                                if(OrPriorityGrouping[ii] == cPriority):
                                    del YpredictB
                                    #-----Predict Segmentation ('B')-------------------------------------------
                                    YpredictB, YpredictBTime = makePrediction(XpredictB[:,:,:,[channel_map[ii]-1]], device, Bmodel, wpathC[ii], B_predictbatch)
                                    if(display_times):
                                        print(r'Prediction Time(Cyto Segmentation): ' + YpredictBTime)
                                    #--------------------------------------------------------------------------

                                    #-----Clean Predictions----------------------------------------------------
                                    [XpredictB, YpredictB, CenterCrops, xp, yp] = cleanMasks_MP(0, mX, mY, cytoMaskThresh[ii], border_mirror, XpredictB, YpredictB, CenterCrops, xp, yp, 2)
                                    #[XpredictB, YpredictB] = cleanMasks_CMG(XpredictB, YpredictB, 0.3)
                                    #[XpredictB, YpredictB] = cleanMasks_CMG_weights(XpredictB, YpredictB)
                                    #--------------------------------------------------------------------------

                                    tMasks[:,:,:,1] = OrMasks_GPU( tMasks[:,:,:,[1]], YpredictB, activeList, device )

                            for n in range(0, N):
                                if(np.sum(tMasks[n,:,:,1]) < 1):
                                    activeList[n] = 1
                                else:
                                    activeList[n] = 0

                            tMasks[:,:,:,1] = OrMasks_GPU( tMasks[:,:,:,[1]], tMasks[:,:,:,[0]], np.invert(activeList*255)/255, device )

                        for ii in range(1, numChannelMap):
                            if(isBackgroundDark[ii]==0):
                                for n in range(0, N):
                                    XpredictB[n,:,:,(channel_map[ii]-1) ] = np.invert(XpredictB[n,:,:,(channel_map[ii]-1) ]);


                        #"or" all cyto predictions with dilated nuclei mask
                        dMasks = np.zeros((N, H, W, 1), dtype='uint8')
                        #struct1 = ndimage.generate_binary_structure(2, 1)
                        struct1 = skimage.morphology.disk(dilation_radius, dtype=np.uint8)
                        #struct1 = ndimage.generate_binary_structure(2, 1)
                        for n in range(0, N):
                            #dMasks[n,:,:,0] = ndimage.binary_dilation(tMasks[n,:,:,0], structure=struct1,iterations=1).astype(tMasks[n,:,:,0].dtype)
                            dMasks[n,:,:,0] = skimage.morphology.dilation(tMasks[n,:,:,0], footprint=struct1)
                            #dMasks[n,:,:,0] = ndimage.binary_dilation(tMasks[n,:,:,0], structure=struct1,iterations=dilation_radius).astype(tMasks[n,:,:,0].dtype)
                        tMasks[:,:,:,1] = OrMasks_GPU(tMasks[:,:,:,[1]], dMasks[:,:,:,[0]], activeList, device )
                        del dMasks
                        del struct1

                    #--------------------------------------------------------------------------

                        del YpredictB
                        YpredictB = copy.deepcopy(tMasks)
                        del tMasks

                    if(isBackgroundDark[0]==0):
                        for n in range(0, N):
                            XpredictB[n,:,:,0 ] = np.invert(XpredictB[n,:,:,0 ]);

                    #-----Crop Mirrored Areas--------------------------------------------------
                    start_mirrorcrop = time.time()
                    [ xp, yp, Screenx, Screeny, XpredictB, YpredictB, CenterCrops ] = cropMirrored_getScreenPos(xp, yp, XpredictB, YpredictB, CenterCrops, FrameSize, border_mirror);
                    end_mirrorcrop = time.time()
                    if(display_times):
                        print(r'Time to crop nuclei in mirrored areas: ' + str(end_mirrorcrop - start_mirrorcrop))
                    #--------------------------------------------------------------------------


                    # #-----Append other pages of input image to segmented image-----------------
                    # if(not(Include_Other_Pages)):
                    # #     XpredictB = XpredictB
                    # #     #XpredictB = includePages(imagefilepath, imagefilename, XpredictB, xp, yp, Screenx, Screeny);
                    # # else:
                    #     for n in range(0, len(XpredictB)):
                    #         tempframe = XpredictB[n][:,:,TIFF_page_number-1]
                    #         XpredictB[n] = tempframe[:,:,np.newaxis]
                    #     #XpredictB = XpredictB[:,:,:,TIFF_page_number-1]
                    # #--------------------------------------------------------------------------





                    Z = np.array([0]);
                    # if(isAutoDetectBackground):
                    #     h_Background = np.ones(N,dtype='float32') * background_value
                    # else:
                    h_Background = np.array([0]);
                    if(isOME):
                        if Objective is not None:
                            h_Objective = np.ones(N,dtype='uint32') * np.uint32(Objective)
                        else:
                            h_Objective = np.array([0]);
                        if Resolution is not None:
                            h_Resolution = np.ones(N,dtype='float32') * Resolution
                        else:
                            h_Resolution = np.array([0]);
                        if ImageType is not None:
                            if(ImageType == r'Absorption'):
                                h_ImageType = np.zeros(N,dtype='uint8')
                            elif(ImageType == r'Concentration'):
                                h_ImageType = np.ones(N,dtype='uint8')
                            elif(ImageType == r'Fluorescence'):
                                h_ImageType = np.ones(N,dtype='uint8') * 2
                        else:
                            h_ImageType = np.array([0]);

                        if ScaleF is not None:
                            h_ScaleF = np.uint16( np.ones(N,dtype='uint16') * (ScaleF*100))
                        else:
                            h_ScaleF = np.array([0]);
                    else:
                        h_Objective = np.array([0]);
                        h_Resolution = np.array([0]);
                        h_ImageType = np.array([0]);
                        h_ScaleF = np.array([0]);

                    if(CMG_subimg_out or CCG_subimg_out):
                        cmg_h = CMG_header(N, Z,Z,np.uint32(Screenx+pX),np.uint32(Screeny+pY),Z,Z,Z,h_Resolution,Z,Z, Z,Z,Z,Z,Z,Z,Z,Z,Z,h_Objective, Z,Z,Z,h_ImageType,np.uint32(xp+pX),np.uint32(yp+pY),Z,h_Background,Z,h_ScaleF, Z,Z,Z)
                        if(CMG_subimg_out):
                            writeCMG(cmg_h, XpredictB, YpredictB, outpath+r'CMG\\', outfilename + r'(' + str(xx) + r', ' + str(yy) + r')')
                        if(CCG_subimg_out):
                            writeCCG(cmg_h, iLabels, mLabels, XpredictB, YpredictB, outpath+r'CCG\\', outfilename + r'(' + str(xx) + r', ' + str(yy) + r')', CCG_subimg_out_I, 0)

                    Z = np.array([0]);
                    if(doSSC and doCytoSegmention):
                        cb4_path = outpath + r'\CB4\\'
                        cb4_name = outfilename + r'(' + str(xx) + r', ' + str(yy) + r')'
                        if not os.path.exists(cb4_path):
                            os.mkdir(cb4_path)
                        Header_struct = makeCMGHeaderStructure(XpredictB,YpredictB,np.uint32(Screenx+pX),np.uint32(Screeny+pY),h_Resolution,h_Objective,h_ImageType,np.uint32(xp+pX),np.uint32(yp+pY),h_Background)
                        #Header = CMG_header(N, Z,Z,np.uint32(Screenx),np.uint32(Screeny),Z,Z,Z,h_Resolution,Z,Z, Z,Z,Z,Z,Z,Z,Z,Z,Z,h_Objective, Z,Z,Z,h_ImageType,np.uint32(xp),np.uint32(yp),Z,h_Background,Z,h_ScaleF, Z,Z,Z);

                        if(isFirstTile):

                            [flabels, fvalues] = SSC(XpredictB, YpredictB, Header_struct, 0, iLabels, mLabels, 0, -1, -1, SSCweight)
                            writeCBX(flabels, fvalues, cb4_path, cb4_name, r'cb4')
                        else:
                            [flabels, tile_out_fvalues] = SSC(XpredictB, YpredictB, Header_struct, 0, iLabels, mLabels, 0, -1, -1, SSCweight)
                            fvalues = np.concatenate( (fvalues, tile_out_fvalues), axis=0 )
                            writeCBX(flabels, tile_out_fvalues, cb4_path, cb4_name, r'cb4')

                    #Append iteration data
                    if(CMG_out or CCG_out):
                        if(isFirstTile):
                            T_N = N
                            T_Screenx = Screenx+pX
                            T_Screeny = Screeny+pY
                            T_h_Resolution = h_Resolution
                            T_h_Objective = h_Objective
                            T_h_ImageType = h_ImageType
                            T_xp = xp+pX
                            T_yp = yp+pY
                            T_h_Background = h_Background
                            T_h_ScaleF = h_ScaleF

                            T_XpredictB = XpredictB
                            T_YpredictB = YpredictB
                        else:
                            T_N = T_N + N
                            T_Screenx = np.append(T_Screenx, Screenx+pX)
                            T_Screeny = np.append(T_Screeny, Screeny+pY)
                            T_h_Resolution = np.append(T_h_Resolution, h_Resolution)
                            T_h_Objective = np.append(T_h_Objective, h_Objective)
                            T_h_ImageType = np.append(T_h_ImageType, h_ImageType)
                            T_xp = np.append(T_xp, xp+pX)
                            T_yp = np.append(T_yp, yp+pY)
                            T_h_Background = np.append(T_h_Background, h_Background)
                            T_h_ScaleF = np.append(T_h_ScaleF, h_ScaleF)

                            T_XpredictB = T_XpredictB + XpredictB
                            T_YpredictB = T_YpredictB + YpredictB

                    #End of iteration
                    isFirstTile = 0;
            count = count+1

    #####ITERATIONS OVER##########
    #-----WRITE TOTAL CMG/CCG--------------------------------------------------
    if(CMG_out or CCG_out):
        Header = CMG_header(T_N, Z,Z,np.uint32(T_Screenx),np.uint32(T_Screeny),Z,Z,Z,T_h_Resolution,Z,Z, Z,Z,Z,Z,Z,Z,Z,Z,Z,T_h_Objective, Z,Z,Z,T_h_ImageType,np.uint32(T_xp),np.uint32(T_yp),Z,T_h_Background,Z,T_h_ScaleF, Z,Z,Z);
        if(CMG_out):
            start_cmg_out = time.time()
            #cmg_h = CMG_header(T_N, Z,Z,np.uint32(T_Screenx),np.uint32(T_Screeny),Z,Z,Z,T_h_Resolution,Z,Z, Z,Z,Z,Z,Z,Z,Z,Z,Z,T_h_Objective, Z,Z,Z,T_h_ImageType,np.uint32(T_xp),np.uint32(T_yp),Z,T_h_Background,Z,T_h_ScaleF, Z,Z,Z);
            #writeCMG(cmg_h, T_XpredictB, T_YpredictB, outpath, outfilename)
            writeCMG(Header, T_XpredictB, T_YpredictB, outpath, outfilename)
            end_cmg_out = time.time()
            if(display_times):
                print(r'Time to Write CMG: ' + str(end_cmg_out - start_cmg_out))
        if(CCG_out):
            start_ccg_out = time.time()
            #ccg_h = CCG_header(T_N, Z,Z,np.uint32(T_Screenx),np.uint32(T_Screeny),Z,Z,Z,T_h_Resolution,Z,Z, Z,Z,Z,Z,Z,Z,Z,Z,Z,T_h_Objective, Z,Z,Z,T_h_ImageType,np.uint32(T_xp),np.uint32(T_yp),Z,T_h_Background,Z,T_h_ScaleF, Z,Z,Z);
            #writeCCG(ccg_h, iLabels, mLabels, T_XpredictB, T_YpredictB, outpath, outfilename, 0)
            writeCCG(Header, iLabels, mLabels, T_XpredictB, T_YpredictB, outpath, outfilename, CCG_out_I, 0)
            end_ccg_out = time.time()
            if(display_times):
                print(r'Time to Write CCG: ' + str(end_ccg_out - start_ccg_out))
    #--------------------------------------------------------------------------

    if(doSSC and doCytoSegmention):
        #[out_flabels, out_fvalues] = SSC(T_XpredictB, T_YpredictB, Header, 0, iLabels, mLabels, 0, -1, -1, SSCweight)
        #[in_flabels, in_fvalues] = makeHeaderFeatures(Header)
        #flabels = in_flabels + out_flabels
        #fvalues = np.concatenate( (in_fvalues, out_fvalues), axis=1 )
        writeCBX(flabels, fvalues, outpath, outfilename, r'cb4')

    print('Full Segmentation Complete.')

    return []


def makePrediction(Xpredict, device, model, wPath, predictbatch):
    start_predict = time.time()
    model.load_state_dict(torch.load(wPath, map_location=torch.device(device), weights_only=True))

    N, Y, X, C = Xpredict.shape
    Xpredict = Xpredict.reshape((N, C, Y, X))
    N_batches = int(np.ceil(N / predictbatch))
    Ypredict = np.ones((N, 1, X, Y), dtype='float32');
    pbi=0
    for i in tqdm(range(0, N_batches), position=0, leave=True):
        if(i == (N_batches-1)):
            img = Xpredict[pbi:,[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict[pbi:,[0],:,:] = y_pred.cpu().detach().numpy()
        else:
            img = Xpredict[pbi:(pbi+predictbatch),[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict[pbi:(pbi+predictbatch),[0],:,:] = y_pred.cpu().detach().numpy()
            pbi=pbi+predictbatch

    Xpredict = Xpredict.reshape((N, Y, X, C))
    Ypredict = Ypredict.reshape((N, Y, X, C))
    end_predict = time.time()
    predictTime = str(end_predict - start_predict)
    return Ypredict, predictTime

def OrMasks_GPU(A, B, activeList, device):
    C = A
    gputA = torch.from_numpy(A).float().to(device)
    gputB = torch.from_numpy(B).float().to(device)

    gputC = torch.logical_or(gputA, gputB)

    C = np.uint8( gputC.cpu().detach().numpy() )

    #Undo OR for cells in previous priorities
    N = len(activeList)
    for n in range(0, N):
        if(activeList[n]==0):
            C[n,:,:,:] = A[n,:,:,:]

    del [gputA, gputB, gputC]
    torch.cuda.empty_cache()
    return C[:,:,:,0]