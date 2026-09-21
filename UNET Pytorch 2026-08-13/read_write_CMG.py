#Writen by Paul Gallagher (BC Cancer Research Centre)

#import os
import numpy as np
#import io
import struct

#[Images, Masks, Header] = readCMG(r'\\crcfile11\CI\Paul_Gallagher\MATLAB Applications\Utility Scripts&Functions\CMG_Creator\TestCMGTiffs', 'TestCMG')

def readCMG(path, filename):

#    path=r'\\crcfile11\CI\Paul_Gallagher\MATLAB Applications\Utility Scripts&Functions\CMG_Creator'
#    filename='cancer3'
    slash = '/'
    cmgpath = path + slash + filename + '.cmg'

    #with open(cmg_path, "rb") as cmg_file:
    cmgfile = open(cmgpath, "rb")
    cmgdata = cmgfile.read()

    btotal = len(cmgdata)
    bcurrent = 0
    icount = 0


    while (btotal > bcurrent):###Find out how many objects
        NbColorMap = cmgdata[bcurrent+2]
        NbBitMap = cmgdata[bcurrent+96]
        width = int.from_bytes(cmgdata[bcurrent+48 : bcurrent+51], byteorder='little')
        height= int.from_bytes(cmgdata[bcurrent+52 : bcurrent+55], byteorder='little')
        imagesize = width * height

        metapadding = 0;
        isMeta = True;
        while(isMeta):
            if(cmgdata[bcurrent+127+metapadding] != 36):
                metapadding = metapadding+1;
            else:
                isMeta = False;

        bcurrent = bcurrent + 128 + metapadding
        bcurrent = bcurrent + int(imagesize)*int(NbColorMap)
        bcurrent = bcurrent + int(imagesize)*int(NbBitMap)
        icount = icount + 1

    ctotal = icount

    Images = []
    Masks = []

    Mode = np.zeros(ctotal, dtype='uint8')
    NbColorMap = np.zeros(ctotal, dtype='uint8')
    Class = np.zeros(ctotal, dtype='uint32')
    Screenx = np.zeros(ctotal, dtype='uint32')
    Screeny = np.zeros(ctotal, dtype='uint32')
    Stagex = np.zeros(ctotal, dtype='uint64')
    Stagey = np.zeros(ctotal, dtype='uint64')
    Stagez = np.zeros(ctotal, dtype='uint64')
    Resolution = np.zeros(ctotal, dtype='float')
    LowThreshold = np.zeros(ctotal, dtype='uint16')
    MidThreshold = np.zeros(ctotal, dtype='uint16')
    Group = np.zeros(ctotal, dtype='uint8')
    Width = np.zeros(ctotal, dtype='uint32')
    Height = np.zeros(ctotal, dtype='uint32')
    Accession = np.zeros(ctotal, dtype='uint32')
    Iod = np.zeros(ctotal, dtype='float')
    Fluor = np.zeros(ctotal, dtype='uint8')
    Diagnosis = np.zeros(ctotal, dtype='uint16')
    RedFaction = np.zeros(ctotal, dtype='float')
    GreenFaction = np.zeros(ctotal, dtype='float')
    BlueFaction = np.zeros(ctotal, dtype='float')
    Index = np.zeros(ctotal, dtype='uint32')
    Objective = np.zeros(ctotal, dtype='uint32')
    Calibrated = np.zeros(ctotal, dtype='uint8')
    StackX_int = np.zeros(ctotal, dtype='uint32')
    StackY_int = np.zeros(ctotal, dtype='uint32')
    NbBitMap = np.zeros(ctotal, dtype='uint8')
    CassettePosition = np.zeros(ctotal, dtype='uint8')
    vorx = np.zeros(ctotal, dtype='uint32')
    vory = np.zeros(ctotal, dtype='uint32')
    BestFocusFrame = np.zeros(ctotal, dtype='uint8')
    BackgroundFloat = np.zeros(ctotal, dtype='float')
    PrimaryColourChannel = np.zeros(ctotal, dtype='uint8')
    #Layer = np.zeros((ctotal, 2), dtype='uint8')
    Layer = np.zeros(ctotal, dtype='uint16')
    Points = np.zeros((ctotal, 9), dtype='uint8')
    NumFeature = np.zeros(ctotal, dtype='uint8')
    RGB_Order = np.zeros((ctotal, 3), dtype='uint8')

    bcurrent = 0
    icount = 0

    #while (BYTE_TOTAL > BYTE_CURRENT):###
    for n in range(0, ctotal):
    #----------------------------Parse Header Data---------------------------------
        Mode[icount] = cmgdata[bcurrent+1]
        NbColorMap[icount] = cmgdata[bcurrent+2]
        Class[icount] = int.from_bytes(cmgdata[bcurrent+3 : bcurrent+7], byteorder='little')
        Screenx[icount] = int.from_bytes(cmgdata[bcurrent+7 : bcurrent+11], byteorder='little')
        Screeny[icount] = int.from_bytes(cmgdata[bcurrent+11 : bcurrent+15], byteorder='little')
        Stagex[icount] = int.from_bytes(cmgdata[bcurrent+15 : bcurrent+23], byteorder='little')
        Stagey[icount] = int.from_bytes(cmgdata[bcurrent+23 : bcurrent+31], byteorder='little')
        Stagez[icount] = int.from_bytes(cmgdata[bcurrent+31 : bcurrent+39], byteorder='little')
        [Resolution[icount]] = struct.unpack('<f', cmgdata[bcurrent+39 : bcurrent+43])
        LowThreshold[icount] = int.from_bytes(cmgdata[bcurrent+43 : bcurrent+45], byteorder='little')
        MidThreshold[icount] = int.from_bytes(cmgdata[bcurrent+45 : bcurrent+47], byteorder='little')
        Group[icount] = cmgdata[bcurrent+47]
        Width[icount] = int.from_bytes(cmgdata[bcurrent+48 : bcurrent+52], byteorder='little')
        Height[icount]= int.from_bytes(cmgdata[bcurrent+52 : bcurrent+56], byteorder='little')
        Accession[icount]= int.from_bytes(cmgdata[bcurrent+56 : bcurrent+60], byteorder='little')
        [Iod[icount]] = struct.unpack('<f', cmgdata[bcurrent+60 : bcurrent+64])
        Fluor[icount] = cmgdata[bcurrent+64]
        Diagnosis[icount]= int.from_bytes(cmgdata[bcurrent+65 : bcurrent+67], byteorder='little')
        [RedFaction[icount]] = struct.unpack('<f', cmgdata[bcurrent+67 : bcurrent+71])
        [GreenFaction[icount]] = struct.unpack('<f', cmgdata[bcurrent+71 : bcurrent+75])
        [BlueFaction[icount]] = struct.unpack('<f', cmgdata[bcurrent+75 : bcurrent+79])
        Index[icount] = int.from_bytes(cmgdata[bcurrent+79 : bcurrent+83], byteorder='little')
        Objective[icount] = int.from_bytes(cmgdata[bcurrent+83 : bcurrent+87], byteorder='little')
        Calibrated[icount] = cmgdata[bcurrent+87]
        StackX_int[icount] = int.from_bytes(cmgdata[bcurrent+88 : bcurrent+92], byteorder='little')
        StackY_int[icount] = int.from_bytes(cmgdata[bcurrent+92 : bcurrent+96], byteorder='little')
        NbBitMap[icount] = cmgdata[bcurrent+96]
        CassettePosition[icount] = cmgdata[bcurrent+97]
        vorx[icount] = int.from_bytes(cmgdata[bcurrent+98 : bcurrent+102], byteorder='little')
        vory[icount] = int.from_bytes(cmgdata[bcurrent+102 : bcurrent+106], byteorder='little')
        BestFocusFrame[icount] = cmgdata[bcurrent+106]
        [BackgroundFloat[icount]] = struct.unpack('<f', cmgdata[bcurrent+107 : bcurrent+111])
        PrimaryColourChannel[icount] = cmgdata[bcurrent+111]
        Layer[icount] = int.from_bytes(cmgdata[bcurrent+112 : bcurrent+114], byteorder='little')
        # for i in range(0, 2):
        #     Layer[icount, i] = cmgdata[bcurrent+112+i]
        for i in range(0, 9):
            Points[icount, i] = cmgdata[bcurrent+114+i]
        NumFeature[icount] = cmgdata[bcurrent+123]
        for i in range(0, 3):
            RGB_Order[icount, i] = cmgdata[bcurrent+124+i]

        imagesize = Width[icount] * Height[icount]

        metapadding = 0;
        isMeta = True;
        while(isMeta):
            if(cmgdata[bcurrent+127+metapadding] != 36):
                metapadding = metapadding+1;
            else:
                isMeta = False;
        #print(metapadding)
        bcurrent = bcurrent + 128 + metapadding
    #------------------------------------------------------------------------------


    #---------------------------Parse Image Data-----------------------------------
        frames = np.zeros((Height[icount], Width[icount], NbColorMap[icount]))
        for n in range(0, NbColorMap[icount]):
            bytevector = np.frombuffer(cmgdata[  bcurrent + int(imagesize)*n  :  bcurrent + int(imagesize)*(n+1) ], dtype=np.uint8)
            frames[:,:,n] = np.reshape(bytevector, (Height[icount], Width[icount]))
            #bytevector.reshape( Height[icount], Width[icount] )

        Images.append(frames)
        bcurrent = bcurrent + int(imagesize)*int(NbColorMap[icount])

    #------------------------------------------------------------------------------

    #--------------------------Parse Mask Data-------------------------------------
        frames = np.zeros((Height[icount], Width[icount], NbBitMap[icount]))
        for n in range(0, NbBitMap[icount]):
            bytevector = np.frombuffer(cmgdata[  bcurrent + int(imagesize)*n  :  bcurrent + int(imagesize)*(n+1) ], dtype=np.uint8)
            frames[:,:,n] = np.reshape(bytevector, (Height[icount], Width[icount]))
            #bytevector.reshape( Height[icount] , Width[icount] )

        Masks.append(frames)
        bcurrent = bcurrent + int(imagesize)*int(NbBitMap[icount])


    #------------------------------------------------------------------------------
        icount = icount + 1

    Header = [Mode, NbColorMap, Class, Screenx, Screeny, Stagex, Stagey, Stagez, Resolution, LowThreshold, MidThreshold, Group, Width, Height, Accession, Iod, Fluor, Diagnosis, RedFaction, GreenFaction, BlueFaction, Index, Objective, Calibrated, StackX_int, StackY_int, NbBitMap, CassettePosition, vorx, vory, BestFocusFrame, BackgroundFloat, PrimaryColourChannel, Layer, Points, NumFeature, RGB_Order]

    cmgfile.close()

    return [Images, Masks, Header]


# def writeCMG(Header, Images, Masks, path, filename):
#     numRuns = len(Images)

#     slash = '/'
#     cmgpath = path + slash + filename + '.cmg'

#     cmgfile = open(cmgpath, "wb")

#     for n in range(0, numRuns):

#         H = Images[n].shape[0]
#         W = Images[n].shape[1]
#         numImage = Images[n].shape[2]
#         numMask = Masks[n].shape[2]

#         #Header[3][n] = Header[28][n] - (W/2);
#         #Header[4][n] = Header[29][n] - (H/2);

#         I_RGB = Images[n]
#         I_Bitmap = Masks[n]

#         cmgfile.write(b'c')
#         cmgfile.write(Header[0][n])#Mode
#         cmgfile.write(np.uint8(numImage))#NBColorMap
#         cmgfile.write(Header[2][n])#Class
#         cmgfile.write(Header[3][n])#Screenx
#         cmgfile.write(Header[4][n])#Screeny
#         cmgfile.write(Header[5][n])#Stagex
#         cmgfile.write(Header[6][n])#Stagey
#         cmgfile.write(Header[7][n])#Stagez
#         cmgfile.write(bytearray(struct.pack("f",Header[8][n])))#Resolution
#         cmgfile.write(Header[9][n])#LowThreshold
#         cmgfile.write(Header[10][n])#MidThreshold
#         cmgfile.write(Header[11][n])#Group
#         cmgfile.write(np.uint32(W))#Width
#         cmgfile.write(np.uint32(H))#Height
#         cmgfile.write(Header[14][n])#Accession
#         cmgfile.write(bytearray(struct.pack("f",Header[15][n])))#Iod
#         cmgfile.write(Header[16][n])#Fluor
#         cmgfile.write(Header[17][n])#Diagnosis
#         cmgfile.write(bytearray(struct.pack("f",Header[18][n])))#RedFraction
#         cmgfile.write(bytearray(struct.pack("f",Header[19][n])))#GreenFraction
#         cmgfile.write(bytearray(struct.pack("f",Header[20][n])))#BlueFraction
#         cmgfile.write(Header[21][n])#Index
#         cmgfile.write(Header[22][n])#Objective
#         cmgfile.write(Header[23][n])#Calibrated
#         cmgfile.write(Header[24][n])#StackX_int
#         cmgfile.write(Header[25][n])#StackY_int
#         cmgfile.write(np.uint8(numMask))#NbBitMap
#         cmgfile.write(Header[27][n])#CassettePosition
#         cmgfile.write(Header[28][n])#vorx
#         cmgfile.write(Header[29][n])#vory
#         cmgfile.write(Header[30][n])#BestFocusFrame
#         cmgfile.write(bytearray(struct.pack("f",Header[31][n])))#BackgroundFloat
#         cmgfile.write(Header[32][n])#PrimaryColorChannel
#         cmgfile.write(Header[33][n][0])#Layer
#         cmgfile.write(Header[33][n][1])#Layers
#         cmgfile.write(Header[34][n][0])#Points
#         cmgfile.write(Header[34][n][1])#Points
#         cmgfile.write(Header[34][n][2])#Points
#         cmgfile.write(Header[34][n][3])#Points
#         cmgfile.write(Header[34][n][4])#Points
#         cmgfile.write(Header[34][n][5])#Points
#         cmgfile.write(Header[34][n][6])#Points
#         cmgfile.write(Header[34][n][7])#Points
#         cmgfile.write(Header[34][n][8])#Points
#         cmgfile.write(Header[35][n])#NumFeatures
#         cmgfile.write(Header[36][n][0])#RGB_Order
#         cmgfile.write(Header[36][n][1])#RGB_Order
#         cmgfile.write(Header[36][n][2])#RGB_Order
#         cmgfile.write(b'$')

#         for z in range(0, numImage):
#             for y in range(0, H):
#                 cmgfile.write(bytearray(np.uint8(I_RGB[y, :, z])))

#         for z in range(0, numMask):
#             for y in range(0, H):
#                 cmgfile.write(bytearray(np.uint8(I_Bitmap[y, :, z])))

#     cmgfile.close()
#     #return []


def writeCMG(Header, Images, Masks, path, filename):
    numRuns = len(Images)

    slash = '/'
    cmgpath = path + slash + filename + '.cmg'

    cmgfile = open(cmgpath, "wb")

    for n in range(0, numRuns):

        H = Images[n].shape[0]
        W = Images[n].shape[1]
        ImagesLen = len(Images[n].shape)
        MasksLen = len(Masks[n].shape)

        if(ImagesLen==2):
            numImage = 1;
            I_RGB = np.zeros((H, W, numImage), dtype='uint8')
            I_RGB[:,:,0] = Images[n]
        elif(ImagesLen==3):
            numImage = Images[n].shape[2];
            I_RGB = np.zeros((H, W, numImage), dtype='uint8')
            I_RGB[:,:,:] = Images[n]

        if(MasksLen==2):
            numMask = 1;
            I_Bitmap = np.zeros((H, W, numMask), dtype='uint8')
            I_Bitmap[:,:,0] = Masks[n]
        elif(MasksLen==3):
            numMask = Masks[n].shape[2];
            I_Bitmap = np.zeros((H, W, numMask), dtype='uint8')
            I_Bitmap[:,:,:] = Masks[n]

        #Header[3][n] = Header[28][n] - (W/2);
        #Header[4][n] = Header[29][n] - (H/2);

        #I_RGB = Images[n]
        #I_Bitmap = Masks[n]

        cmgfile.write(b'c')
        cmgfile.write(Header[0][n])#Mode
        cmgfile.write(np.uint8(numImage))#NBColorMap
        cmgfile.write(Header[2][n])#Class
        cmgfile.write(Header[3][n])#Screenx
        cmgfile.write(Header[4][n])#Screeny
        cmgfile.write(Header[5][n])#Stagex
        cmgfile.write(Header[6][n])#Stagey
        cmgfile.write(Header[7][n])#Stagez
        cmgfile.write(bytearray(struct.pack("f",Header[8][n])))#Resolution
        cmgfile.write(Header[9][n])#LowThreshold
        cmgfile.write(Header[10][n])#MidThreshold
        cmgfile.write(Header[11][n])#Group
        cmgfile.write(np.uint32(W))#Width
        cmgfile.write(np.uint32(H))#Height
        cmgfile.write(Header[14][n])#Accession
        cmgfile.write(bytearray(struct.pack("f",Header[15][n])))#Iod
        cmgfile.write(Header[16][n])#Fluor
        cmgfile.write(Header[17][n])#Diagnosis
        cmgfile.write(bytearray(struct.pack("f",Header[18][n])))#RedFraction
        cmgfile.write(bytearray(struct.pack("f",Header[19][n])))#GreenFraction
        cmgfile.write(bytearray(struct.pack("f",Header[20][n])))#BlueFraction
        cmgfile.write(Header[21][n])#Index
        cmgfile.write(Header[22][n])#Objective
        cmgfile.write(Header[23][n])#Calibrated
        cmgfile.write(Header[24][n])#StackX_int
        cmgfile.write(Header[25][n])#StackY_int
        cmgfile.write(np.uint8(numMask))#NbBitMap
        cmgfile.write(Header[27][n])#CassettePosition
        cmgfile.write(Header[28][n])#vorx
        cmgfile.write(Header[29][n])#vory
        cmgfile.write(Header[30][n])#BestFocusFrame
        cmgfile.write(bytearray(struct.pack("f",Header[31][n])))#BackgroundFloat
        cmgfile.write(Header[32][n])#PrimaryColorChannel
        cmgfile.write(Header[33][n])#Layer
        #cmgfile.write(Header[33][n][0])#Layers
        #cmgfile.write(Header[33][n][1])#Layers
        cmgfile.write(Header[34][n][0])#Points
        cmgfile.write(Header[34][n][1])#Points
        cmgfile.write(Header[34][n][2])#Points
        cmgfile.write(Header[34][n][3])#Points
        cmgfile.write(Header[34][n][4])#Points
        cmgfile.write(Header[34][n][5])#Points
        cmgfile.write(Header[34][n][6])#Points
        cmgfile.write(Header[34][n][7])#Points
        cmgfile.write(Header[34][n][8])#Points
        cmgfile.write(Header[35][n])#NumFeatures
        cmgfile.write(Header[36][n][0])#RGB_Order
        cmgfile.write(Header[36][n][1])#RGB_Order
        cmgfile.write(Header[36][n][2])#RGB_Order
        cmgfile.write(b'$')

        for z in range(0, numImage):
            for y in range(0, H):
                cmgfile.write(bytearray(np.uint8(I_RGB[y, :, z])))

        for z in range(0, numMask):
            for y in range(0, H):
                cmgfile.write(bytearray(np.uint8(I_Bitmap[y, :, z])))

    cmgfile.close()
    #return []

def CMG_header(NumberOfImage,Mode,Class,Screenx,Screeny,Stagex,Stagey,Stagez,Resolution,LowThreshold,MidThreshold,Group,Accession,Iod,Fluor,Diagnosis,RedFaction,GreenFaction,BlueFaction,Index,Objective,Calibrated,StackX_int,StackY_int,CassettePosition,vorx,vory,BestFocusFrame,BackgroundFloat,PrimaryColourChannel,Layer,Points,NumFeature,RGB_Order):
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
        Resolution = np.ones(NumberOfImage,dtype='float32') * 0.18

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
        Iod = np.ones(NumberOfImage,dtype='float32') * 22500

    if(len(Fluor)!=NumberOfImage):
        Fluor = np.zeros(NumberOfImage,dtype='uint8')

    if(len(Diagnosis)!=NumberOfImage):
        Diagnosis = np.zeros(NumberOfImage,dtype='uint16')

    if(len(RedFaction)!=NumberOfImage):
        RedFaction = np.ones(NumberOfImage,dtype='float32') *0#* 0.3

    if(len(GreenFaction)!=NumberOfImage):
        GreenFaction = np.ones(NumberOfImage,dtype='float32') *0#* 0.7

    if(len(BlueFaction)!=NumberOfImage):
        BlueFaction = np.zeros(NumberOfImage,dtype='float32') *0#

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
        BackgroundFloat = np.ones(NumberOfImage,dtype='float32') * 254

    if(len(PrimaryColourChannel)!=NumberOfImage):
        PrimaryColourChannel = np.zeros(NumberOfImage,dtype='uint8')

    # sLayer = Layer.shape
    # if(sLayer[0]!=NumberOfImage or 1 == NumberOfImage):
    #     Layer = np.zeros((NumberOfImage, 2),dtype='uint8')
    if(len(Layer)!=NumberOfImage):
        Layer = np.zeros(NumberOfImage,dtype='uint16')

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