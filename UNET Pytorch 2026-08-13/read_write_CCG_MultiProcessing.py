#Writen by Paul Gallagher (BC Cancer Research Centre)

#import os
import time
import numpy as np
import copy
#import io
import struct
import multiprocessing as mp
#from threading import Thread, current_thread
#from multiprocessing import Process, current_process
#from PIL import Image, ImageOps
import skimage.io as skio
from skimage import measure
import scipy.ndimage as ndimage
#from itertools import chain

#[Images, Masks, Header] = readCMG(r'\\crcfile11\CI\Paul_Gallagher\MATLAB Applications\Utility Scripts&Functions\CMG_Creator\TestCMGTiffs', 'TestCMG')

def CCG_header(NumberOfImage,Mode,Class,Screenx,Screeny,Stagex,Stagey,Stagez,Resolution,LowThreshold,MidThreshold,Group,Accession,Iod,Fluor,Diagnosis,RedFaction,GreenFaction,BlueFaction,Index,Objective,Calibrated,StackX_int,StackY_int,CassettePosition,vorx,vory,BestFocusFrame,BackgroundFloat,PrimaryColourChannel,Layer,Points,NumFeature,RGB_Order):
    ###example: [Header] = header(3,[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0])

    Header = -1

    if(len(Mode)!=NumberOfImage):
        Mode = np.ones(NumberOfImage,dtype='uint8')*99

    if(len(Class)!=NumberOfImage):
        Class = np.zeros(NumberOfImage,dtype='uint32')

    if(len(Screenx)!=NumberOfImage):
        Screenx = np.zeros(NumberOfImage,dtype='uint32')

    if(len(Screeny)!=NumberOfImage):
        Screeny = np.zeros(NumberOfImage,dtype='uint32')

    if(len(Stagex)!=NumberOfImage):
        Stagex = np.zeros(NumberOfImage,dtype='uint64')

    if(len(Stagey)!=NumberOfImage):
        Stagey = np.zeros(NumberOfImage,dtype='uint64')

    if(len(Stagez)!=NumberOfImage):
        Stagez = np.zeros(NumberOfImage,dtype='uint64')

    if(len(Resolution)!=NumberOfImage):
        Resolution = np.ones(NumberOfImage,dtype='float') * 0.18

    if(len(LowThreshold)!=NumberOfImage):
        LowThreshold = np.ones(NumberOfImage,dtype='uint16') * 9728

    if(len(MidThreshold)!=NumberOfImage):
        MidThreshold = np.ones(NumberOfImage,dtype='uint16') * 4864

    if(len(Group)!=NumberOfImage):
        Group = np.zeros(NumberOfImage,dtype='uint8')

    if(len(Accession)!=NumberOfImage):
        Accession = np.zeros(NumberOfImage,dtype='uint32')
        for x in range(1, NumberOfImage):
            Accession[x:] = Accession[x:] + np.ones((NumberOfImage-x),dtype='uint32')

    if(len(Iod)!=NumberOfImage):
        Iod = np.ones(NumberOfImage,dtype='float') * 22500

    if(len(Fluor)!=NumberOfImage):
        Fluor = np.zeros(NumberOfImage,dtype='uint8')

    if(len(Diagnosis)!=NumberOfImage):
        Diagnosis = np.zeros(NumberOfImage,dtype='uint16')

    if(len(RedFaction)!=NumberOfImage):
        RedFaction = np.ones(NumberOfImage,dtype='float') *0#* 0.3

    if(len(GreenFaction)!=NumberOfImage):
        GreenFaction = np.ones(NumberOfImage,dtype='float') *0#* 0.7

    if(len(BlueFaction)!=NumberOfImage):
        BlueFaction = np.zeros(NumberOfImage,dtype='float') *0#

    if(len(Index)!=NumberOfImage):
        Index = np.zeros(NumberOfImage,dtype='uint32')

    if(len(Objective)!=NumberOfImage):
        Objective = np.ones(NumberOfImage,dtype='uint32') * 538980404

    if(len(Calibrated)!=NumberOfImage):
        Calibrated = np.zeros(NumberOfImage,dtype='uint8')

    if(len(StackX_int)!=NumberOfImage):
        StackX_int = np.zeros(NumberOfImage,dtype='uint32')

    if(len(StackY_int)!=NumberOfImage):
        StackY_int = np.zeros(NumberOfImage,dtype='uint32')

    if(len(CassettePosition)!=NumberOfImage):
        CassettePosition = np.zeros(NumberOfImage,dtype='uint8')

    if(len(vorx)!=NumberOfImage):
        vorx = np.zeros(NumberOfImage,dtype='uint32')

    if(len(vory)!=NumberOfImage):
        vory = np.zeros(NumberOfImage,dtype='uint32')

    if(len(BestFocusFrame)!=NumberOfImage):
        BestFocusFrame = np.zeros(NumberOfImage,dtype='uint8')

    if(len(BackgroundFloat)!=NumberOfImage):
        BackgroundFloat = np.ones(NumberOfImage,dtype='float') * 254

    if(len(PrimaryColourChannel)!=NumberOfImage):
        PrimaryColourChannel = np.zeros(NumberOfImage,dtype='uint8')

    if(len(Layer)!=NumberOfImage):
        Layer = np.zeros(NumberOfImage,dtype='uint16')
    # sLayer = Layer.shape
    # if(sLayer[0]!=NumberOfImage or 1 == NumberOfImage):
    #     Layer = np.zeros((NumberOfImage, 2),dtype='uint8')

    sPoints = Points.shape
    if(sPoints[0]!=NumberOfImage or 1 == NumberOfImage):
        Points = np.zeros((NumberOfImage, 10),dtype='uint8')

    if(len(NumFeature)!=NumberOfImage):
        NumFeature = np.ones(NumberOfImage,dtype='uint8') * 3

    if(len(RGB_Order)!=NumberOfImage or 1 == NumberOfImage):
        RGB_Order = np.zeros((NumberOfImage, 3),dtype='uint8')
        RGB_Order[:,0] = np.zeros((NumberOfImage),dtype='uint8')
        RGB_Order[:,1] = np.ones((NumberOfImage),dtype='uint8') * 1
        RGB_Order[:,2] = np.ones((NumberOfImage),dtype='uint8') * 2

    Width = 0
    Height = 0
    NbColorMap = 0
    NbBitMap = 0

    Header = [Mode, NbColorMap, Class, Screenx, Screeny, Stagex, Stagey, Stagez, Resolution, LowThreshold, MidThreshold, Group, Width, Height, Accession, Iod, Fluor, Diagnosis, RedFaction, GreenFaction, BlueFaction, Index, Objective, Calibrated, StackX_int, StackY_int, NbBitMap, CassettePosition, vorx, vory, BestFocusFrame, BackgroundFloat, PrimaryColourChannel, Layer, Points, NumFeature, RGB_Order]

    return Header



def readCCG(path, filename):

    slash = '/'
    tiffpath = path + slash + filename + '.TIFF'
    ccgpath = path + slash + filename + '.ccg'

    #----------------------------------Read Multi page TIFF-------------------------------------

    frame = skio.imread(tiffpath)

    #-------------------------------------------------------------------------------------------

    #with open(cmg_path, "rb") as cmg_file:
    ccgfile = open(ccgpath, "rb")
    ccgdata = ccgfile.read()

    BYTE_TOTAL = len(ccgdata)
    CURRENT_BYTE = np.uint64(0);
    icount = 0

    #Get number of cell objects.
    N = int.from_bytes(ccgdata[CURRENT_BYTE:CURRENT_BYTE+np.uint64(7)], byteorder='little')
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(8);

    #Get number of image planes.
    P = int.from_bytes(ccgdata[CURRENT_BYTE:CURRENT_BYTE+np.uint64(1)], byteorder='little')
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(2);

    #Read plane labels
    ImageLabel = [''] * P
    for p in range(0, P):
        ImageLabel[p] = ccgdata[CURRENT_BYTE:CURRENT_BYTE+np.uint64(15)].decode("utf-8")
        ImageLabel[p] = ImageLabel[p].replace(" ", "")
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(16);

    #Read vorx for each object
    h_vorx = np.zeros(N, dtype='uint32')
    for n in range(0, N):
        h_vorx[n] = int.from_bytes(ccgdata[CURRENT_BYTE:CURRENT_BYTE+np.uint64(3)], byteorder='little')
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(4);

    #Read vory for each object
    h_vory = np.zeros(N, dtype='uint32')
    for n in range(0, N):
        h_vory[n] = int.from_bytes(ccgdata[CURRENT_BYTE:CURRENT_BYTE+np.uint64(3)], byteorder='little')
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(4);

    #Read memory locations for each object
    memoryloc = np.zeros(N, dtype='uint64')
    for n in range(0, N):
        memoryloc[n] = int.from_bytes(ccgdata[CURRENT_BYTE:CURRENT_BYTE+np.uint64(7)], byteorder='little')
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(8);

    Mode = np.zeros(N, dtype='uint8')
    NbColorMap = np.zeros(N, dtype='uint8')
    Class = np.zeros(N, dtype='uint32')
    Screenx = np.zeros(N, dtype='uint32')
    Screeny = np.zeros(N, dtype='uint32')
    Stagex = np.zeros(N, dtype='uint64')
    Stagey = np.zeros(N, dtype='uint64')
    Stagez = np.zeros(N, dtype='uint64')
    Resolution = np.zeros(N, dtype='float')
    LowThreshold = np.zeros(N, dtype='uint16')
    MidThreshold = np.zeros(N, dtype='uint16')
    Group = np.zeros(N, dtype='uint8')
    Width = np.zeros(N, dtype='uint32')
    Height = np.zeros(N, dtype='uint32')
    Accession = np.zeros(N, dtype='uint32')
    Iod = np.zeros(N, dtype='float')
    Fluor = np.zeros(N, dtype='uint8')
    Diagnosis = np.zeros(N, dtype='uint16')
    RedFaction = np.zeros(N, dtype='float')
    GreenFaction = np.zeros(N, dtype='float')
    BlueFaction = np.zeros(N, dtype='float')
    Index = np.zeros(N, dtype='uint32')
    Objective = np.zeros(N, dtype='uint32')
    Calibrated = np.zeros(N, dtype='uint8')
    StackX_int = np.zeros(N, dtype='uint32')
    StackY_int = np.zeros(N, dtype='uint32')
    NbBitMap = np.zeros(N, dtype='uint8')
    CassettePosition = np.zeros(N, dtype='uint8')
    vorx = np.zeros(N, dtype='uint32')
    vory = np.zeros(N, dtype='uint32')
    BestFocusFrame = np.zeros(N, dtype='uint8')
    BackgroundFloat = np.zeros(N, dtype='float')
    PrimaryColourChannel = np.zeros(N, dtype='uint8')
    Layer = np.zeros(N, dtype='uint16')
    Points = np.zeros((N, 9), dtype='uint8')
    NumFeature = np.zeros(N, dtype='uint8')
    RGB_Order = np.zeros((N, 3), dtype='uint8')
    Images = []
    Masks = []
    MaskLabel = []

    ImageCount = np.zeros(N, dtype='uint16')
    MaskCount = np.zeros(N, dtype='uint16')
    ilist = [r''] * N


    while(BYTE_TOTAL > CURRENT_BYTE):
        Mode[icount] = ccgdata[CURRENT_BYTE+np.uint64(1)]
        NbColorMap[icount] = ccgdata[CURRENT_BYTE+np.uint64(2)]
        Class[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(3) : CURRENT_BYTE+np.uint64(7)], byteorder='little')
        Screenx[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(7) : CURRENT_BYTE+np.uint64(11)], byteorder='little')
        Screeny[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(11) : CURRENT_BYTE+np.uint64(15)], byteorder='little')
        Stagex[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(15) : CURRENT_BYTE+np.uint64(23)], byteorder='little')
        Stagey[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(23) : CURRENT_BYTE+np.uint64(31)], byteorder='little')
        Stagez[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(31) : CURRENT_BYTE+np.uint64(39)], byteorder='little')
        [Resolution[icount]] = struct.unpack('<f', ccgdata[CURRENT_BYTE+np.uint64(39) : CURRENT_BYTE+np.uint64(43)])
        LowThreshold[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(43) : CURRENT_BYTE+np.uint64(45)], byteorder='little')
        MidThreshold[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(45) : CURRENT_BYTE+np.uint64(47)], byteorder='little')
        Group[icount] = ccgdata[CURRENT_BYTE+np.uint64(47)]
        Width[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(48) : CURRENT_BYTE+np.uint64(52)], byteorder='little')
        Height[icount]= int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(52) : CURRENT_BYTE+np.uint64(56)], byteorder='little')
        Accession[icount]= int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(56) : CURRENT_BYTE+np.uint64(60)], byteorder='little')
        [Iod[icount]] = struct.unpack('<f', ccgdata[CURRENT_BYTE+np.uint64(60) : CURRENT_BYTE+np.uint64(64)])
        Fluor[icount] = ccgdata[CURRENT_BYTE+np.uint64(64)]
        Diagnosis[icount]= int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(65) : CURRENT_BYTE+np.uint64(67)], byteorder='little')
        [RedFaction[icount]] = struct.unpack('<f', ccgdata[CURRENT_BYTE+np.uint64(67) : CURRENT_BYTE+np.uint64(71)])
        [GreenFaction[icount]] = struct.unpack('<f', ccgdata[CURRENT_BYTE+np.uint64(71) : CURRENT_BYTE+np.uint64(75)])
        [BlueFaction[icount]] = struct.unpack('<f', ccgdata[CURRENT_BYTE+np.uint64(75) : CURRENT_BYTE+np.uint64(79)])
        Index[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(79) : CURRENT_BYTE+np.uint64(83)], byteorder='little')
        Objective[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(83) : CURRENT_BYTE+np.uint64(87)], byteorder='little')
        Calibrated[icount] = ccgdata[CURRENT_BYTE+np.uint64(87)]
        StackX_int[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(88) : CURRENT_BYTE+np.uint64(92)], byteorder='little')
        StackY_int[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(92) : CURRENT_BYTE+np.uint64(96)], byteorder='little')
        NbBitMap[icount] = ccgdata[CURRENT_BYTE+np.uint64(96)]
        CassettePosition[icount] = ccgdata[CURRENT_BYTE+np.uint64(97)]
        vorx[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(98) : CURRENT_BYTE+np.uint64(102)], byteorder='little')
        vory[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(102) : CURRENT_BYTE+np.uint64(106)], byteorder='little')
        BestFocusFrame[icount] = ccgdata[CURRENT_BYTE+np.uint64(106)]
        [BackgroundFloat[icount]] = struct.unpack('<f', ccgdata[CURRENT_BYTE+np.uint64(107) : CURRENT_BYTE+np.uint64(111)])
        PrimaryColourChannel[icount] = ccgdata[CURRENT_BYTE+np.uint64(111)]
        # for i in range(0, 2):
        #     Layer[icount, i] = ccgdata[CURRENT_BYTE+112+i]
        Layer[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(112) : CURRENT_BYTE+np.uint64(113)], byteorder='little')
        for i in range(0, 9):
            Points[icount, i] = ccgdata[CURRENT_BYTE+np.uint64(113)+np.uint64(i)]
        NumFeature[icount] = ccgdata[CURRENT_BYTE+np.uint64(123)]
        for i in range(0, 3):
            RGB_Order[icount, i] = ccgdata[CURRENT_BYTE+np.uint64(123)+np.uint64(i)]

        CURRENT_BYTE = CURRENT_BYTE + np.uint64(128);


        #---------------------------Parse Image Data-----------------------------------
        frames = np.zeros((Height[icount], Width[icount], P))

        for p in range(0, P):
            frames[:,:,p] = frame[p, Screeny[icount]:(Screeny[icount]+Height[icount]), Screenx[icount]:(Screenx[icount]+Width[icount])]

        Images.append(frames)
        #------------------------------------------------------------------------------

        #--------------------------Parse Mask Data-------------------------------------

        #Get number of images for this cell
        ImageCount[icount] = int.from_bytes(ccgdata[CURRENT_BYTE : CURRENT_BYTE+np.uint64(1)], byteorder='little')

        #Get number of masks for this cell
        MaskCount[icount] = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(2) : CURRENT_BYTE+np.uint64(3)], byteorder='little')
        mask = np.zeros((Height[icount], Width[icount], MaskCount[icount]), 'uint8')

        CURRENT_BYTE = CURRENT_BYTE + np.uint64(4)

        #Get cell interaction list
        ilen = ccgdata[CURRENT_BYTE]
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(1)
        ilist[icount] = np.zeros(ilen, 'float')

        for i in range(0, ilen):
            ilist[icount][i] = struct.unpack('<f', ccgdata[CURRENT_BYTE : CURRENT_BYTE+np.uint64(4)])[0]
            CURRENT_BYTE = CURRENT_BYTE + np.uint64(4)

        #Temp space to store mask labels
        #readStr = np.ones((MaskCount[icount], 16), dtype='uint8')*32
        readStr = [''] * MaskCount[icount]



        for n in range(0, MaskCount[icount]):

            #Read mask label
            readStr[n] = ccgdata[CURRENT_BYTE : CURRENT_BYTE+np.uint64(15)].decode("utf-8")
            readStr[n] = readStr[n].replace(" ", "")

            #Read chaincount
            ChainCount = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(16) : CURRENT_BYTE+np.uint64(17)], byteorder='little')

            CURRENT_BYTE = CURRENT_BYTE + np.uint64(18)



            for m in range(0, ChainCount):

                #Read starting X
                StartX = int.from_bytes(ccgdata[CURRENT_BYTE : CURRENT_BYTE+np.uint64(3)], byteorder='little')

                #Read starting Y
                StartY = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(4) : CURRENT_BYTE+np.uint64(7)], byteorder='little')

                #Read chain length
                CLength = int.from_bytes(ccgdata[CURRENT_BYTE+np.uint64(8) : CURRENT_BYTE+np.uint64(11)], byteorder='little')
                CCode = np.zeros(CLength, dtype='uint8')

                #If no chaincode
                if(CLength==0):

                    #Set starting pixel
                    mask[StartY-1, StartX-1, n] = 1

                    #Update current byte
                    CURRENT_BYTE = CURRENT_BYTE + np.uint64(12)

                else:

                    #Read in chaincode bytes
                    for i in range(0, CLength):
                        CCode[i] = ccgdata[CURRENT_BYTE+np.uint64(12)+np.uint64(i)]

                    #Convert chaincode to mask
                    mask[:, :, n] = generateMaskfromCCode( mask[:, :, n], CCode, StartX, StartY )

                    #Update current byte
                    CURRENT_BYTE = CURRENT_BYTE + np.uint64(12 + CLength)

        Masks.append(mask)

        MaskLabel.append(readStr)
        icount = icount+1
        #------------------------------------------------------------------------------



    Header = [Mode, NbColorMap, Class, Screenx, Screeny, Stagex, Stagey, Stagez, Resolution, LowThreshold, MidThreshold, Group, Width, Height, Accession, Iod, Fluor, Diagnosis, RedFaction, GreenFaction, BlueFaction, Index, Objective, Calibrated, StackX_int, StackY_int, NbBitMap, CassettePosition, vorx, vory, BestFocusFrame, BackgroundFloat, PrimaryColourChannel, Layer, Points, NumFeature, RGB_Order]

    ccgfile.close()

    return [Images, Masks, Header, ImageLabel, MaskLabel, ilist]











def writeCCG(Header, ImageLabel, MaskLabel, Images, Masks, path, filename, calcOverlap, MP):

    print('-------------------------------------------------')
    print(path + filename + '.ccg')
    #print('Writing CCG...')

    numRuns = len(Masks)
    Header[12] = np.zeros(numRuns, dtype='uint32')
    Header[13] = np.zeros(numRuns, dtype='uint32')
    Header[1] = np.zeros(numRuns, dtype='uint8')
    Header[26] = np.zeros(numRuns, dtype='uint8')
    for n in range(0, numRuns):
        Header[12][n] = Images[n].shape[1]
        Header[13][n] = Images[n].shape[0]
        Header[1][n] = Images[n].shape[2]
        Header[26][n] = Masks[n].shape[2]


    minX = np.min(Header[3]);
    minY = np.min(Header[4]);

    Header[3] = Header[3] - minX;
    Header[4] = Header[4] - minY;

    Header[28] = Header[28] - minX;
    Header[29] = Header[29] - minY;




    #--------------------------Saving TIFF-------------------------------------
    start_maketiff = time.time()

    print('Constructing TIFF...')

    maxX = np.max(Header[3]);
    I_X = np.argmax(Header[3]);
    maxY = np.max(Header[4]);
    I_Y = np.argmax(Header[4]);

    if(numRuns > 100):
        Nbackimg = 100;
    else:
        Nbackimg = numRuns;

    borderpixels = np.ones((1,) , 'uint8' )*128;
    planes = np.zeros(numRuns, 'uint16');
    for n in range(0, numRuns):
        planes[n] = np.uint16( Images[n].shape[2] );
        # if(Images[n].shape[2] > maxZ):
        #     maxZ = Images[n].shape[2];
    for n in range(0, Nbackimg):
        borderpixels = np.append(borderpixels, Images[n][0,:,0], 0);
        borderpixels = np.append(borderpixels, Images[n][:,0,0], 0);
        borderpixels = np.append(borderpixels, Images[n][(Images[n].shape[0]-1),:,0], 0);
        borderpixels = np.append(borderpixels, Images[n][:,(Images[n].shape[1]-1),0], 0);

    maxZ = np.max(planes);
    AvgBP = np.mean(borderpixels);
    #I_Z = np.argmax(Header[1]);
    #[ maxY, I_Y] = np.max(Header[4]);
    #[ maxZ, I_Z] = np.max(Header[1]);

    #XX = maxX - minX + Images[I_X].shape[1];
    #YY = maxY - minY + Images[I_Y].shape[0];

    XX = maxX + Images[I_X].shape[1];
    YY = maxY + Images[I_Y].shape[0];

    if(AvgBP > 128):
        DATA_frame = np.ones((YY, XX, maxZ) , 'uint8')*255;
    else:
        DATA_frame = np.zeros((YY, XX, maxZ) , 'uint8');

    for nn in range(0, numRuns):
        # posx = Header[3][nn] - minX;
        # posy = Header[4][nn] - minY;
        posx = Header[3][nn];
        posy = Header[4][nn];
        tY = Header[12][nn]
        tX = Header[13][nn]
        try:
            tZ = Images[nn].shape[2]
        except:
            tZ = 1;

        #[tY, tX, tZ] = Images[nn].shape;
        for zz in range(0, tZ):
            DATA_frame[posy:(posy+tY), posx:(posx+tX), zz] = Images[nn][:,:,zz];

        #Header[3][nn] = posx;
        #Header[4][nn] = posy;

    end_maketiff = time.time()
    maketifftime = (end_maketiff - start_maketiff)
    print('make Tiff time: ' + str(maketifftime))


    start_savetiff = time.time()

    writepath = path + '\\' + filename + r'.TIFF';

    print('Writing TIFF...')
    saveThread = Thread(target=saveCCGTIFF, args =(DATA_frame, maxZ, writepath))
    saveThread.start()

    end_savetiff = time.time()
    savetifftime = (end_savetiff - start_savetiff)
    print('save Tiff time: ' + str(savetifftime))
    del DATA_frame
    #--------------------------------------------------------------------------



    slash = '/'
    ccgpath = path + slash + filename + '.ccg'
    #ccgfile = open(ccgpath, "wb")
    N = len(Images);

    freemantime = 0;

    CURRENT_BYTE = np.uint64(0);

    #Find Interacting Cells
    if(calcOverlap):
        print('Finding Interacting Cells...')
        I_Obj = findOverlap(Masks, Header, 0);
    else:
        I_Obj = [[]]*numRuns


    #Write number of objects
    ccg_bytes = bytearray(struct.pack( '<Q', np.uint64(N) ) );
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(8);

    #Write number of image planes
    ccg_bytes.extend(np.uint16(maxZ).astype('uint16').tobytes())
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(2);

    #Write image labels
    for z in range(0, maxZ):
        for m in range(0, 16):
            try:
                ccg_bytes.extend(np.array([ImageLabel[z][m]], 'S1').tobytes())
            except:
                ccg_bytes.extend(np.array([r' '], 'S1').tobytes())
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(16*maxZ);

    #Write vorx
    ccg_bytes.extend(Header[28][:].astype('uint32').tobytes())
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(4*N);

    #Write vory
    ccg_bytes.extend(Header[29][:].astype('uint32').tobytes())
    CURRENT_BYTE = CURRENT_BYTE + np.uint64(4*N);

    MEMSTART = CURRENT_BYTE;
    #Allocate space for pointers
    memlist = np.zeros((N), 'uint64');
    ccg_bytes.extend(memlist.astype('uint64').tobytes())

    CURRENT_BYTE = CURRENT_BYTE + np.uint64(8*N);

    ccg_header_len = len(ccg_bytes)

    out_ccg_bytes = mp.Queue()
    out_memlist = mp.Queue()

    doMP = 0
    maskThresh = 1000


    print('Writing Chain Code...')
    if(MP <= 0):
        [ccg_bytes2, memlist] = processCellsObjects(Images, Masks, Header, MaskLabel, I_Obj, out_ccg_bytes, out_memlist, 0);
        ccg_bytes.extend(ccg_bytes2)
    else:
        if(MP == 1):
            if(sum(Header[26]) > maskThresh):
                numProc = np.int32(mp.cpu_count() / 2)
                doMP = 1
            else:
                [ccg_bytes2, memlist] = processCellsObjects(Images, Masks, Header, MaskLabel, I_Obj, out_ccg_bytes, out_memlist, 0);
                ccg_bytes.extend(ccg_bytes2)
                doMP = 0


        elif(MP > 1):
            doMP = 1
            numProc = MP


        if(doMP):

            procs = []
            B = np.int64( np.ceil(N/numProc) );
            cCell = 0;
            for i in range(0, numProc):
                if(i+1 != numProc):
                    Img = Images[cCell:(cCell+B)]
                    Msk = Masks[cCell:(cCell+B)]
                    MskLb = MaskLabel#[cCell:(cCell+B)]
                    iObj = I_Obj[cCell:(cCell+B)]
                    Hdr = copy.deepcopy(Header)
                    for j in range(0, 37):
                        Hdr[j] = Header[j][cCell:(cCell+B)]
                    cCell = cCell+B;

                    p = mp.Process(
                            target=processCellsObjects,
                            args=( Img, Msk, Hdr, MskLb, iObj, out_ccg_bytes, out_memlist, i) )

                    # p = mp.Process(
                    #         target=testMultiProcess,
                    #         args=( i, out_ccg_bytes, out_memlist) )

                    procs.append(p)
                    p.start()

                else:
                    Img = Images[cCell:]
                    Msk = Masks[cCell:]
                    MskLb = MaskLabel#[cCell:]
                    iObj = I_Obj[cCell:]
                    Hdr = copy.deepcopy(Header)
                    for j in range(0, 37):
                        Hdr[j] = Header[j][cCell:]

                    p = mp.Process(
                            target=processCellsObjects,
                            args=(Img, Msk, Hdr, MskLb, iObj, out_ccg_bytes, out_memlist, i) )

                    procs.append(p)
                    p.start()
                    i=i

            resultdict_bytes = {}
            resultdict_mem = {}
            for i in range(0, numProc):
                resultdict_bytes.update(out_ccg_bytes.get())
                resultdict_mem.update(out_memlist.get())

            for p in procs:
                p.join()

            numbytes = 0
            curcell = 0
            for i in range(0, numProc):
                ccg_bytes.extend(resultdict_bytes[i])
                if(i+1 != numProc):
                    memlist[curcell:(curcell+B)] = [x+numbytes for x in resultdict_mem[i] ]
                else:
                    memlist[curcell:] = [x+numbytes for x in resultdict_mem[i] ]
                numbytes = numbytes + len(resultdict_bytes[i])
                curcell = curcell + B

    memlist = [x+ccg_header_len for x in memlist]







    ccg_bytes[  MEMSTART  : (MEMSTART + np.uint64(N*8))  ] = np.uint64(memlist).astype('uint64').tobytes();
         #ccgfile.write(np.uint64(memlist[n]))
    start_saveccg = time.time()

    ccgfile = open(ccgpath, "wb")
    ccgfile.write(ccg_bytes)
    ccgfile.close()

    end_saveccg = time.time()
    saveccgtime = (end_saveccg - start_saveccg)
    print('save ccg time: ' + str(saveccgtime))

    #print('Freeman time: ' + str(freemantime))
    print('Done.')
    print('-------------------------------------------------')
    return []



def saveCCGTIFF(DATA_frame, maxZ, writepath):
    im = []
    for z in range(0, maxZ):
        im.append(Image.fromarray(DATA_frame[:,:,z]))

    im[0].save(writepath, save_all=True, append_images=im[1:], compression="tiff_lzw")
    return []


def testMultiProcess(INPUT, out_ccg_bytes, out_memlist):
    outdict = {}
    outdict[INPUT] = INPUT
    out_ccg_bytes.put(outdict)
    out_memlist.put(outdict)


def processCellsObjects(Images, Masks, Header, MaskLabel, I_Obj, out_ccg_bytes, out_memlist, procID):
    N = len(Masks)
    memlist = np.zeros((N), 'uint64');
    CURRENT_BYTE = np.uint64(0);
    ccg_bytes = []#bytearray(struct.pack( '<Q', np.uint64(0) ) );
    ccg_bytes_dict = {}
    memlist_dict = {}

    for n in range(0, N):

        memlist[n] = CURRENT_BYTE;
        #dictkey = str(procID) + '_' + str(n)
        #memlist_dict[dictkey] = CURRENT_BYTE

        H = Images[n].shape[0]
        W = Images[n].shape[1]
        numImage = Images[n].shape
        numMask = Masks[n].shape

        #Header[3][n] = Header[28][n] - (W/2);
        #Header[4][n] = Header[29][n] - (H/2);

        #I_RGB = Images[n]
        #I_Bitmap = Masks[n]


        ccg_bytes.extend(b'c')
        ccg_bytes.extend(Header[0][n].astype('uint8').tobytes())
        ccg_bytes.extend(np.uint8(numImage[2]).astype('uint8').tobytes())
        ccg_bytes.extend(Header[2][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[3][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[4][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[5][n].astype('uint64').tobytes())
        ccg_bytes.extend(Header[6][n].astype('uint64').tobytes())
        ccg_bytes.extend(Header[7][n].astype('uint64').tobytes())
        ccg_bytes.extend(Header[8][n].astype('float32').tobytes())
        ccg_bytes.extend(Header[9][n].astype('uint16').tobytes())
        ccg_bytes.extend(Header[10][n].astype('uint16').tobytes())
        ccg_bytes.extend(Header[11][n].astype('uint8').tobytes())
        ccg_bytes.extend(np.uint32(W).astype('uint32').tobytes())
        ccg_bytes.extend(np.uint32(H).astype('uint32').tobytes())
        ccg_bytes.extend(Header[14][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[15][n].astype('float32').tobytes())
        ccg_bytes.extend(Header[16][n].astype('uint8').tobytes())
        ccg_bytes.extend(Header[17][n].astype('uint16').tobytes())
        ccg_bytes.extend(Header[18][n].astype('float32').tobytes())
        ccg_bytes.extend(Header[19][n].astype('float32').tobytes())
        ccg_bytes.extend(Header[20][n].astype('float32').tobytes())
        ccg_bytes.extend(Header[21][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[22][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[23][n].astype('uint8').tobytes())
        ccg_bytes.extend(Header[24][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[25][n].astype('uint32').tobytes())
        ccg_bytes.extend(np.uint8(numMask[2]).astype('uint8').tobytes())
        ccg_bytes.extend(Header[27][n].astype('uint8').tobytes())
        ccg_bytes.extend(Header[28][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[29][n].astype('uint32').tobytes())
        ccg_bytes.extend(Header[30][n].astype('uint8').tobytes())
        ccg_bytes.extend(Header[31][n].astype('float32').tobytes())
        ccg_bytes.extend(Header[32][n].astype('uint8').tobytes())
        ccg_bytes.extend(Header[33][n].astype('uint16').tobytes())
        #ccg_bytes.extend(Header[33][n][0].astype('uint8').tobytes())
        #ccg_bytes.extend(Header[33][n][1].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][0].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][1].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][2].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][3].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][4].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][5].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][6].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][7].astype('uint8').tobytes())
        ccg_bytes.extend(Header[34][n][8].astype('uint8').tobytes())
        ccg_bytes.extend(Header[35][n].astype('uint8').tobytes())
        ccg_bytes.extend(Header[36][n][0].astype('uint8').tobytes())
        ccg_bytes.extend(Header[36][n][1].astype('uint8').tobytes())
        ccg_bytes.extend(Header[36][n][2].astype('uint8').tobytes())
        ccg_bytes.extend(b'$')

        CURRENT_BYTE = CURRENT_BYTE + np.uint64(128);

        #write number of images
        ccg_bytes.extend(np.uint16(numImage[2]).astype('uint16').tobytes())
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(2);

        #write number of masks
        ccg_bytes.extend(np.uint16(numMask[2]).astype('uint16').tobytes())
        CURRENT_BYTE = CURRENT_BYTE + np.uint64(2);

        #Write number of cell interactions
        I_Obj_num = len(I_Obj[n])
        ccg_bytes.extend(np.uint8(I_Obj_num).astype('uint8').tobytes())
        for i in range(0, I_Obj_num):
            ccg_bytes.extend(I_Obj[n][i].astype('float32').tobytes())

        CURRENT_BYTE = CURRENT_BYTE + np.uint64(4*I_Obj_num + 1);

        for z in range(0, numMask[2]):

            mask = Masks[n][:,:,z]

            #Write mask label
            for m in range(0, 16):
                try:
                    ccg_bytes.extend(np.array([MaskLabel[z][m]], 'S1').tobytes())
                except:
                    ccg_bytes.extend(np.array([r' '], 'S1').tobytes())
            CURRENT_BYTE = CURRENT_BYTE + np.uint64(16);

            masksum = np.sum(mask)

            if(masksum==0):
                #Mask is blank
                #Write chain count = zero
                ccg_bytes.extend(np.uint16(0).astype('uint16').tobytes())

                CURRENT_BYTE = CURRENT_BYTE + np.uint64(2);
            else:
                [labelmask, M] = measure.label(mask, return_num=True, connectivity=1 );

                #Write chain count = M
                ccg_bytes.extend(np.uint16(M).astype('uint16').tobytes())
                #ccgfile.write(np.uint16(M))
                CURRENT_BYTE = CURRENT_BYTE + np.uint64(2);

                for m in range(0, M):
                    submask = labelmask == (m+1);

                    submasksum = np.sum(submask);

                    if(submasksum==1):
                        [I, J] = np.where(submask==1)
                        #Write StartX
                        ccg_bytes.extend(np.uint32(J+1).astype('uint32').tobytes())

                        #Write StartY
                        ccg_bytes.extend(np.uint32(I+1).astype('uint32').tobytes())

                        #Write Chain Length=0
                        ccg_bytes.extend(np.uint32(0).astype('uint32').tobytes())


                        CURRENT_BYTE = CURRENT_BYTE + np.uint64(12);
                    else:
                        #start_fman = time.time()
                        [X0, Code, C] = Freeman_Chain_Code(submask)
                        #end_fman = time.time()
                        #freemantime = freemantime + (end_fman - start_fman)

                        #Write StartX
                        ccg_bytes.extend(np.uint32(X0[0]).astype('uint32').tobytes())

                        #Write StartY
                        ccg_bytes.extend(np.uint32(X0[1]).astype('uint32').tobytes())

                        #Write Chain Length
                        ccg_bytes.extend(np.uint32(C).astype('uint32').tobytes())

                        #Write Chain Data
                        ccg_bytes.extend(np.uint8(Code).astype('uint8').tobytes())


                        CURRENT_BYTE = CURRENT_BYTE + np.uint64(12+C);

    ccg_bytes_dict[procID] = ccg_bytes
    memlist_dict[procID] = memlist

    out_ccg_bytes.put(ccg_bytes_dict)
    out_memlist.put(memlist_dict)

    return [ccg_bytes, memlist]



def Freeman_Chain_Code(mask):
    #Code translated here from Freeman_chain_code.m
    #start_ = time.time()

    X0 = np.zeros(2, 'uint32');
    Code = np.zeros(1, 'uint8');
    [Y, X] = mask.shape

    #start_bpoint = time.time()

    colsum = np.sum(mask, axis=0);
    col_idx0 = np.nonzero(colsum)[0][0];
    row_idx0 = np.nonzero(mask[:, col_idx0])[0][0];

    X0[0] = row_idx0 + 1;
    X0[1] = col_idx0 + 1;

    #end_bpoint = time.time()
    #print('start time: ' + str(end_bpoint - start_bpoint))

    move_index = np.array([[1, 0, -1, -1, -1, 0, 1, 1], [-1, -1, -1, 0, 1, 1, 1, 0]]);

    #while

    row_idx1 = row_idx0;
    col_idx1 = col_idx0;

    dirr = 8;
    dir0 = opp_dir(dirr);

    mask[ row_idx0, col_idx0  ] = 0;

    while ( row_idx1 < 0 or row_idx1 >= Y or col_idx1 < 0 or col_idx1 >= X or mask[row_idx1,col_idx1] == 0 ):
        dir0 = prv_dir(dir0);
        row_idx1 = row_idx0 + move_index[1,dir0-1];
        col_idx1 = col_idx0 + move_index[0,dir0-1];

    Code[0] = np.uint8(dir0);
    mask[row_idx0, col_idx0] = 1;
    dirr = dir0;

    row_idx2 = row_idx1;
    col_idx2 = col_idx1;
    C = 1;

    while( row_idx2 != row_idx0 or col_idx2 != col_idx0 ):
        dir0 = opp_dir(dirr);
        mask[row_idx2, col_idx2] = 0;
        while ( row_idx1 < 0 or row_idx1 >= Y or col_idx1 < 0 or col_idx1 >= X or mask[row_idx1,col_idx1] == 0 ):
            dir0 = prv_dir(dir0);

            row_idx1 = row_idx2 + move_index[1,dir0-1]; # next point coordinates
            col_idx1 = col_idx2 + move_index[0,dir0-1]; # from index matrix

        Code = np.append(Code, np.uint8(dir0));
        dirr = dir0;

        mask[row_idx2,col_idx2] = 1; # reset pixel initial value
        row_idx2 = row_idx1;
        col_idx2 = col_idx1;

        C = C + 1;

    #end_ = time.time()
    #print('total time: ' + str(end_ - start_))

    return[X0, Code, C]

def generateMaskfromCCode(mask, CCode, starty, startx):
    move_index = np.array( [ [1, 0, -1, -1, -1, 0, 1, 1], [-1, -1, -1, 0, 1, 1, 1, 0] ] )

    Cx = startx-1
    Cy = starty-1
    mask[Cy, Cx] = 1

    Clen = len(CCode)

    for c in range(0, Clen):
        Cx = Cx + move_index[0, (CCode[c]-1)]
        Cy = Cy + move_index[1, (CCode[c]-1)]
        mask[Cy, Cx] = 1

    mask = ndimage.binary_fill_holes(mask).astype(int)


    return mask

def findOverlap(Masks, header, mask_plane):
    NN = len(Masks);
    total_list = [[]]*NN
    for nn in range(0, NN):
        SX = np.double( header[3][nn] );        #Screenx
        SY = np.double( header[4][nn] );        #Screeny
        XX = np.uint32( header[12][nn] );       #Width
        YY = np.uint32( header[13][nn] );       #Height

        top = np.floor(SY);
        bottom = np.floor(SY) + np.ceil(YY);
        left = np.floor(SX);
        right = np.floor(SX) + np.ceil(XX);

        xp = header[28]+1;                        #vorx;
        yp = header[29]+1;                        #vory;

        xrange = np.array(np.where(np.logical_and(xp>=left, xp<=right)))+1
        yrange = np.array(np.where(np.logical_and(yp>=top, yp<=bottom)))+1

        L = np.isin(xrange, yrange)
        L = np.multiply(L, xrange)
        L = L[L != 0]
        List = copy.deepcopy(L);
        LL = L.shape[0]

        Xmax=0;
        Ymax=0;
        for l in range(0, LL):
            if(Xmax < header[12][ L[l]-1 ]):
                Xmax = header[12][ L[l]-1 ];

            if(Ymax < header[13][ L[l]-1 ]):
                Ymax = header[13][ L[l]-1 ];

        for l in range(0, LL):
            ll = LL - (l+1);
            Onum = L[ ll ]-1;
            nSX = np.double( header[3][ Onum ] );
            nSY = np.double( header[4][ Onum ] );
            nXX = np.uint32( header[12][ Onum ] );
            nYY = np.uint32( header[13][ Onum ] );

            OX = np.int64(Xmax + (nSX - SX));
            OY = np.int64(Ymax + (nSY - SY));

            Omask = np.zeros((YY+2*Ymax, XX+2*Xmax), 'uint8');
            Omask[ (Ymax):(YY+Ymax) , (Xmax):(XX+Xmax) ] = Masks[nn][:,:,mask_plane];
            Omask[ (OY):(nYY+OY), (OX):(nXX+OX) ] = Omask[ (OY):(nYY+OY), (OX):(nXX+OX) ] + Masks[ Onum ][:,:,mask_plane];

            if( np.max(Omask) <= 1 or (List[ll]-1)==nn ):
                List = np.delete(List, ll)

        total_list[nn] = np.float32(List)

    return total_list

def prv_dir(d):

   dd = np.mod(d-2, 8) + 1;

   return dd

def opp_dir(d):

   dd = np.mod(d+2, 8) + 1;

   return dd