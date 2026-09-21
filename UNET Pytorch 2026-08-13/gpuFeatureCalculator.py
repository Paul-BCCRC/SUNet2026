#Writen by Paul Gallagher (BC Cancer Research Centre)
"""
Created on Tue Mar 10 10:55:35 2026

@author: pgallagher
"""
#import os
import numpy as np
import math
#import time
#import copy
from scipy.interpolate import interp1d

import torch as t
#import scipy.ndimage as ndimage

def IntensityFeatures(ImageTensor, MaskTensor):
    #Calculates features below:
        #area           = mask sum
        #imean          = intensity average
        #itotal         = intensity sum
        #istd           = intensity standard deviation
    [N, Y, X, P] = ImageTensor.shape
    [N, Y, X, M] = MaskTensor.shape

    area = np.zeros((N, M), dtype='float32')
    imean = np.zeros((N, P, M), dtype='float32')
    itotal = np.zeros((N, P, M), dtype='float32')
    istd = np.zeros((N, P, M), dtype='float32')

    #Load Image and Mask data to GPU
    gpuImage = t.from_numpy(ImageTensor).float().to("cuda")
    gpuMask = t.from_numpy(MaskTensor).float().to("cuda")

    #Subtract nucleai mask from dilated and Cyto masks
    gpuMask[:,:,:,[2]] = t.sub(gpuMask[:,:,:,[1]], gpuMask[:,:,:,[0]], alpha=1)
    gpuMask[:,:,:,[4]] = t.sub(gpuMask[:,:,:,[3]], gpuMask[:,:,:,[0]], alpha=1)

    #Calculate Area of all masks
    gpuArea = t.sum( gpuMask, [1,2] )

    #Retrieve and free Areas from VRAM
    area = np.float32( gpuArea.cpu().detach().numpy() )
    #del gpuArea
    #torch.cuda.empty_cache()

    #For each mask
    for m in range(0, M):

        #Apply mask to all Image channel data
        gpuBitMasked = t.mul(gpuImage, gpuMask[:,:,:,[m]])

        #Calculate total intensities for each channel and retrieve from VRAM
        gpuTotal = t.sum(gpuBitMasked, [1, 2])
        itotal[:, :, m] = np.float32( gpuTotal.cpu().detach().numpy() )

        #Unfortuately the mean and sd functions are applied to the whole image including zero values outside the mask.
        #So I calculate them manually.
        #Calculate mean intensities for each channel and retrieve from VRAM
        gpuMean = t.div(gpuTotal, gpuArea[:,[m]])
        imean[:, :, m] = np.float32( gpuMean.cpu().detach().numpy() )

        #Calculate standard dev. intensities for each channel and retrieve from VRAM
        gpuSTD = t.sub(gpuBitMasked, gpuMean[:,np.newaxis,np.newaxis,:])
        gpuSTD = t.mul(gpuSTD, gpuMask[:,:,:,[m]])
        gpuSTD = t.square(gpuSTD)
        gpuSTD = t.sum(gpuSTD, [1, 2])
        gpuSTD = t.div(gpuSTD, gpuArea[:,[m]])
        gpuSTD = t.sqrt(gpuSTD)
        istd[:, :, m] = np.float32( gpuSTD.cpu().detach().numpy() )

        #Clear VRAM
        del [gpuBitMasked, gpuTotal, gpuMean, gpuSTD]
        t.cuda.empty_cache()

    #Cleanup to prevent VRAM memory leaks
    del [gpuImage, gpuMask, gpuArea]
    t.cuda.empty_cache()

    return [area, imean, itotal, istd]


def ShapeFeatures(MaskTensor):
    #start_shapef = time.time()
    #Calculates features below:
        #area           = mask sum
        #parea          = perimeter sum
        #cir            = circumference
        #centroid_x     = x center
        #centroid_y     = y center
        #min_r          = radial min
        #max_r          = radial max
        #mean_r         = radial mean
        #var_r          = radial variance
        #spher          = sphericity
        #eccentr        = eccentricity
        #cell_orient    = cell_orientation
        #inert          = inertia
        #comp           = compact
    gpuCuda0 = t.device('cuda:0')
    [N, Y, X, M] = MaskTensor.shape

    #Allocate memory for features
    area = np.zeros((N, M), dtype='float32')
    parea = np.zeros((N, M), dtype='float32')
    cir = np.zeros((N, M), dtype='float32')
    centroid_x = np.zeros((N, M), dtype='float32')
    centroid_y = np.zeros((N, M), dtype='float32')
    mean_r = np.zeros((N, M), dtype='float32')
    max_r = np.zeros((N, M), dtype='float32')
    min_r = np.zeros((N, M), dtype='float32')
    var_r = np.zeros((N, M), dtype='float32')
    spher = np.zeros((N, M), dtype='float32')
    eccentr = np.zeros((N, M), dtype='float32')
    cell_orient = np.zeros((N, M), dtype='float32')
    inert = np.zeros((N, M), dtype='float32')
    comp = np.zeros((N, M), dtype='float32')

    #Create X and Y position matrixes.
    Gx = np.zeros((Y,X),dtype='float')
    G_row = np.zeros((1,X),dtype='float')
    for x in range(0, X):# make row vector [0 1 2 3 4...]
        G_row[0, x] = x+0.5
    for y in range(0, Y):#set each row to the row vector
        Gx[y, :] = G_row
    Gy = Gx.T#Transpose matrix for y indexes

    #Pass Mask to the GPU
    gpuMask = t.from_numpy(MaskTensor).float().to("cuda")

    #Pass index matrixes to the GPU
    gpuGx = t.from_numpy(Gx).float().to("cuda")
    gpuGy = t.from_numpy(Gy).float().to("cuda")

    #Calculate Perimeter
    gpuT = t.zeros([N,X,Y,M,4], dtype=t.float32, device=gpuCuda0)
    gpuT[:,:-1,:,:,0] = gpuMask[:, 1:, :, :]
    gpuT[:,1:,:,:,1] = gpuMask[:, :-1, :, :]
    gpuT[:,:,:-1,:,2] = gpuMask[:, :, 1:, :]
    gpuT[:,:,1:,:,3] = gpuMask[:, :, :-1, :]
    gpuT = t.sub(gpuMask[:,:,:,:,np.newaxis], gpuT)
    gpuT = t.mul(gpuMask[:,:,:,:,np.newaxis], gpuT)
    gpuT = t.sum(gpuT, [4])
    gpuPerimeter = t.mul(gpuT, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0))
    gpuPerimeter[(gpuPerimeter>1)] = 1

    #Calculate Circumference
    gpuT[(gpuT==2)] = np.sqrt(2)
    gpuT[(gpuT==3)] = 2
    gpuT[(gpuT==4)] = 1
    gpuCir = t.sum(gpuT, [1, 2])
    cir = np.float32( gpuCir.cpu().detach().numpy() )

    #Calculate area of masks
    gpuArea = t.sum( gpuMask, [1,2] )
    gpuPArea = t.sum( gpuPerimeter, [1,2] )
    area = np.float32( gpuArea.cpu().detach().numpy() )
    parea = np.float32( gpuPArea.cpu().detach().numpy() )
    #Set any zeros to at least 1 (divide by zero)
    gpuArea[(gpuArea<=1)] = 2
    gpuPArea[(gpuPArea<=1)] = 2

    #Calculate index/positions of each pixel
    gpuXind = t.mul(gpuMask, gpuGx[np.newaxis, :, :, np.newaxis])
    gpuYind = t.mul(gpuMask, gpuGy[np.newaxis, :, :, np.newaxis])

    #Calculate x and y mean index/positions for each mask (MASK CENTERS)
    gpuMean_x = t.div( t.sum(gpuXind, [1, 2]), gpuArea )
    gpuMean_y = t.div( t.sum(gpuYind, [1, 2]), gpuArea )
    #gpuCentroid_x = t.add( gpuMean_x, t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 0.5 )
    #gpuCentroid_y = t.add( gpuMean_y, t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 0.5 )
    centroid_x = np.float32( gpuMean_x.cpu().detach().numpy() )
    centroid_y = np.float32( gpuMean_y.cpu().detach().numpy() )

    #Squared distance from mean position
    gpuX2 = t.mul( t.square( t.sub(gpuXind, gpuMean_x[:,np.newaxis,np.newaxis,:]) ), gpuMask )
    gpuY2 = t.mul( t.square( t.sub(gpuYind, gpuMean_y[:,np.newaxis,np.newaxis,:]) ), gpuMask )
    gpuX = t.mul( t.sub(gpuMean_x[:,np.newaxis,np.newaxis,:], gpuXind), gpuMask )
    gpuY = t.mul( t.sub(gpuMean_y[:,np.newaxis,np.newaxis,:], gpuYind), gpuMask )

    gpuX_moment2 = t.sum(gpuX2, [1, 2])
    gpuY_moment2 = t.sum(gpuY2, [1, 2])
    gpuXY_cross_moment2 = t.sum( t.mul( gpuX , gpuY ), [1, 2] )

    #Calculate Inertia
    gpuInert = t.div( t.mul( t.add(gpuX_moment2, gpuY_moment2), t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 2*math.pi), t.square(gpuArea) )
    inert = np.float32( gpuInert.cpu().detach().numpy() )

    #Calculate Eccentricity
    gpuCalcEcc = t.add( t.square(gpuX_moment2), t.square(gpuY_moment2)  )
    gpuCalcEcc = t.add(gpuCalcEcc, t.mul( t.square(gpuXY_cross_moment2), t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 4 )  )
    gpuCalcEcc = t.sub(gpuCalcEcc, t.mul( t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 2, t.mul( gpuX_moment2, gpuY_moment2 ) ) )

    gpuCalcEcc_a = t.div( t.sub( t.add(gpuX_moment2,gpuY_moment2), t.sqrt(gpuCalcEcc) ) , t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 2)
    gpuCalcEcc_b = t.div( t.add( t.add(gpuX_moment2,gpuY_moment2), t.sqrt(gpuCalcEcc) ) , t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 2)
    gpuCalcEcc_a[(gpuCalcEcc_a==0)] = 1

    gpuEccentricity = t.div(t.sqrt(gpuCalcEcc_b), t.sqrt(gpuCalcEcc_a) )
    eccentr = np.float32( gpuEccentricity.cpu().detach().numpy() )

    #Calculate Cell Orientation
    gpuXY_cross_moment2[(gpuXY_cross_moment2==0)] = 1
    gpuCell_orient_atan = t.atan( t.div( t.sub(gpuCalcEcc_b, gpuX_moment2), gpuXY_cross_moment2 ) )
    gpuCell_orient = t.mul(t.ones([N, M], dtype=t.float32, device=gpuCuda0) *(180/math.pi) , t.add(t.ones([N, M], dtype=t.float32, device=gpuCuda0) *(math.pi/2) , gpuCell_orient_atan) )
    cell_orient = np.float32( gpuCell_orient.cpu().detach().numpy() )

    #Calculate Compactness
    gpuCompact = t.div( t.square(gpuCir),  t.mul( t.ones([N, M], dtype=t.float32, device=gpuCuda0) *(4*math.pi) , gpuArea) )
    comp = np.float32( gpuCompact.cpu().detach().numpy() )

    #Calculate radial distances
    gpuD2 = t.add(gpuX2, gpuY2)
    gpuRad = t.mul(t.sqrt(gpuD2), gpuPerimeter)

    #Find radial min/max
    gpuMax_r = t.max(t.max(gpuRad, 1)[0], 1)[0]
    gpuminR = t.mul(t.sqrt(gpuD2), gpuPerimeter)
    gpuminR[(gpuminR==0)] = X*Y     #Set Zeros to a big number
    gpuMin_r = t.min(t.min(gpuminR, 1)[0], 1)[0]
    max_r = np.float32( gpuMax_r.cpu().detach().numpy() )
    min_r = np.float32( gpuMin_r.cpu().detach().numpy() )

    #Calculate radial mean
    gpuMean_r = t.div( t.sum(gpuRad, [1,2]), gpuPArea)
    mean_r = np.float32( gpuMean_r.cpu().detach().numpy() )

    #Calculate radial variance
    gpuVarDiff = t.mul( t.sub( gpuRad, gpuMean_r[:, np.newaxis, np.newaxis, :] ), gpuPerimeter )
    gpuNumi = t.sum( t.square( gpuVarDiff ) , [1, 2] )
    gpuDenomi = t.sub(gpuPArea, t.ones([N, M], dtype=t.float32, device=gpuCuda0))
    gpuVar_r =  t.div( gpuNumi  , gpuDenomi )
    var_r = np.float32( gpuVar_r.cpu().detach().numpy() )

    #Calculate sphericity
    gpuSpher = t.div(gpuMin_r, gpuMax_r)
    spher = np.float32( gpuSpher.cpu().detach().numpy() )

    #Cleanup to prevent VRAM memory leaks
    del [gpuMask, gpuPerimeter, gpuGx, gpuGy, gpuArea, gpuPArea, gpuXind, gpuYind, gpuMean_x, gpuMean_y]#, gpuCentroid_x, gpuCentroid_y]
    del [gpuX2, gpuY2, gpuX, gpuY, gpuX_moment2, gpuY_moment2, gpuXY_cross_moment2, gpuInert]
    del [gpuCalcEcc, gpuCalcEcc_a, gpuCalcEcc_b, gpuEccentricity, gpuCell_orient_atan, gpuCell_orient]
    del [gpuT, gpuCir, gpuCompact]
    del [gpuD2, gpuRad, gpuMean_r, gpuMax_r, gpuminR, gpuMin_r, gpuVarDiff, gpuNumi, gpuDenomi,  gpuSpher, gpuVar_r, gpuCuda0]
    t.cuda.empty_cache()

    # end_shapef = time.time()
    # print(r'Time to calculate shape features: ' + str(end_shapef - start_shapef))

    return [area, parea, centroid_x, centroid_y, mean_r, max_r, min_r, var_r, spher, eccentr, cell_orient, inert, cir, comp]


def FFTFeatures(MaskTensor, RADIAL_QUANTIZATION):
    #Calculates features below:
        #harmonics      = fft amplitudes
        #freq_low_fft   = amplitude sum 3 through 6
        #freq_high_fft  = amplitude sum 7 through 31
        #elongation     = elongation using 2nd harmonic

    gpuCuda0 = t.device('cuda:0')
    [N, Y, X, M] = MaskTensor.shape

    #Allocate memory for features
    #TEST = np.zeros((N, M, RADIAL_QUANTIZATION), dtype='float32')
    #TEST2 = np.zeros((N, M, RADIAL_QUANTIZATION), dtype='float32')
    #TEST3 = np.zeros((N, X, Y, M), dtype='float32')
    harmonics = np.zeros((N, M, RADIAL_QUANTIZATION), dtype='float32')
    freq_low_fft = np.zeros((N, M), dtype='float32')
    freq_high_fft = np.zeros((N, M), dtype='float32')
    elongation = np.zeros((N, M), dtype='float32')

    #Create X and Y position matrixes.
    Gx = np.zeros((Y,X),dtype='float32')
    G_row = np.zeros((1,X),dtype='float32')
    for x in range(0, X):# make row vector [0 1 2 3 4...]
        G_row[0, x] = x+0.5
    for y in range(0, Y):#set each row to the row vector
        Gx[y, :] = G_row
    Gy = Gx.T#Transpose matrix for y indexes

    #Pass Mask to the GPU
    gpuMask = t.from_numpy(MaskTensor).float().to("cuda")

    #Pass index matrixes to the GPU
    gpuGx = t.from_numpy(Gx).float().to("cuda")
    gpuGy = t.from_numpy(Gy).float().to("cuda")

    #Calculate Perimeter
    gpuT = t.zeros([N,X,Y,M,4], dtype=t.float32, device=gpuCuda0)
    gpuT[:,:-1,:,:,0] = gpuMask[:, 1:, :, :]
    gpuT[:,1:,:,:,1] = gpuMask[:, :-1, :, :]
    gpuT[:,:,:-1,:,2] = gpuMask[:, :, 1:, :]
    gpuT[:,:,1:,:,3] = gpuMask[:, :, :-1, :]
    gpuT = t.sub(gpuMask[:,:,:,:,np.newaxis], gpuT)
    gpuT = t.mul(gpuMask[:,:,:,:,np.newaxis], gpuT)
    gpuT = t.sum(gpuT, [4])
    gpuPerimeter = t.mul(gpuT, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0))
    gpuPerimeter[(gpuPerimeter>1)] = 1

    #Calculate area of masks
    gpuArea = t.sum( gpuMask, [1,2] )
    gpuPArea = t.sum( gpuPerimeter, [1,2] )
    #Set any zeros to at least 1 (divide by zero)
    gpuArea[(gpuArea<=1)] = 2
    gpuPArea[(gpuPArea<=1)] = 2

    #Calculate index/positions of each pixel
    gpuXind = t.mul(gpuMask, gpuGx[np.newaxis, :, :, np.newaxis])
    gpuYind = t.mul(gpuMask, gpuGy[np.newaxis, :, :, np.newaxis])

    #Calculate x and y mean index/positions for each mask (MASK CENTERS)
    gpuMean_x = t.div( t.sum(gpuXind, [1, 2]), gpuArea )
    gpuMean_y = t.div( t.sum(gpuYind, [1, 2]), gpuArea )

    #Distance from mean position
    gpuX = t.mul( t.sub(gpuXind, gpuMean_x[:,np.newaxis,np.newaxis,:]), gpuPerimeter )
    gpuY = t.mul( t.sub(gpuMean_y[:,np.newaxis,np.newaxis,:], gpuYind), gpuPerimeter )

    #Squared distance from mean position
    gpuX2 = t.square( gpuX )
    gpuY2 = t.square( gpuY )

    #Calculate radial distances
    gpuRad = t.sqrt( t.add(gpuX2, gpuY2) )
    gpuRad_NZ = t.sqrt( t.add(gpuX2, gpuY2) )
    gpuRad_NZ[(gpuRad==0)] = 1

    #Calculate radial mean
    #gpuMean_r = t.div( t.sum(gpuRad, [1,2]), gpuPArea)
    #mean_r = np.float32( gpuMean_r.cpu().detach().numpy() )

    #Calculate Cosine ratio
    gpuW = t.div( gpuX, gpuRad_NZ )
    gpuW[(gpuW > 1)] = 1
    gpuW[(gpuW < -1)] = -1
    gpuW = t.mul( gpuW, gpuPerimeter )

    #Radial bin each perimeter pixel
    gpuW = t.acos(gpuW)

    #Convert to degrees
    gpuW = t.mul( gpuW, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * (360/(2*math.pi)) )
    gpuW[(gpuY<0)] = 360 - gpuW[(gpuY<0)]
    #TEST = np.float32( gpuW.cpu().detach().numpy() )
    ASTEP = 360/RADIAL_QUANTIZATION
    gpuW = t.add( t.div( gpuW, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * ASTEP ) , t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * 0.5)
    gpuW = t.floor(gpuW)
    #TEST2 = np.float32( gpuW.cpu().detach().numpy() )
    gpuW2 = t.sub(gpuW, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * (RADIAL_QUANTIZATION-1))
    gpuW[(gpuW>(RADIAL_QUANTIZATION-1))] = gpuW2[(gpuW>(RADIAL_QUANTIZATION-1))]
    #TEST3 = np.float32( gpuW.cpu().detach().numpy() )
    gpuW = t.mul( t.add( gpuW , gpuPerimeter ), gpuPerimeter )




    #-------------------------------
    # gpuW = t.mul( gpuW, gpuPerimeter )
    # gpuW = t.div( gpuW, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * math.pi )
    # gpuW = t.mul( gpuW, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * (RADIAL_QUANTIZATION/2) )
    # # gpuW = t.add( gpuW, t.ones([N, X, Y, M], dtype=t.float32, device=gpuCuda0) * 0.5 )
    # gpuW = t.floor(gpuW)
    # gpuW[int(RADIAL_QUANTIZATION/2)] = (RADIAL_QUANTIZATION/2)-1

    # TEST = np.float32( gpuW.cpu().detach().numpy() )

    # gpuW[(gpuY<0)] = (RADIAL_QUANTIZATION-1) - gpuW[(gpuY<0)]
    # TEST2 = np.float32( gpuW.cpu().detach().numpy() )

    # gpuW = t.mul( t.add( gpuW , gpuPerimeter ), gpuPerimeter )
    # TEST3 = np.float32( gpuW.cpu().detach().numpy() )
    #--------------------

    gpuQRad = t.zeros([N, M, RADIAL_QUANTIZATION], dtype=t.float32, device=gpuCuda0)
    gpuQCount = t.zeros([N, M, RADIAL_QUANTIZATION], dtype=t.float32, device=gpuCuda0)


    for q in range(0, RADIAL_QUANTIZATION):
        Qnum = q+1
        gpuQMask = t.zeros([N, X, Y, M], dtype=t.float32, device=gpuCuda0)
        gpuQMask[(gpuW==Qnum)] = 1
        # if(q==0):
        #     gpuQMask[(gpuW==1+RADIAL_QUANTIZATION)] = 1

        gpuQCount[:, :, q] = t.sum(gpuQMask, [1, 2])
        gpuQRad[:, :, q] = t.sum(t.mul(gpuRad, gpuQMask), [1, 2])

    gpuQCount[(gpuQCount==0)] = 1
    gpuQRad = t.div(gpuQRad, gpuQCount)

    #TEST2 = np.float32( gpuQRad.cpu().detach().numpy() )
    QRad = np.float32( gpuQRad.cpu().detach().numpy() )
    #QRad2 = np.float32( gpuQRad.cpu().detach().numpy() )

    # for n in range(0, N):
    #     for m in range(0, M):
    #         radius = np.zeros((RADIAL_QUANTIZATION+1),dtype=float)
    #         radius[:-1] = QRad2[n, m, :]
    #         radius[RADIAL_QUANTIZATION] = radius[0]
    #         space = 0
    #         last = 0
    #         i = RADIAL_QUANTIZATION

    #         while((last==np.float32(0)) and (i >= 0)):
    #             if(radius[i] > 0):
    #                 last = radius[i]
    #             else:
    #                 space=space+1
    #                 i=i-1

    #         for i in range(0, RADIAL_QUANTIZATION+1):
    #             if(radius[i] > 0):
    #                 fnew = radius[i]
    #                 if(space > 0):
    #                     diff = (fnew - last) / (space + 1)
    #                     for j in range(1, space+1):
    #                         k = i-j
    #                         if(k<0):
    #                             k += RADIAL_QUANTIZATION
    #                         radius[k] = fnew - diff * j
    #                     space = 0
    #                 last = fnew
    #             if(radius[i] == 0):
    #                 space=space+1

    #         radius[RADIAL_QUANTIZATION] = radius[0]
    #         QRad2[n, m, :] = radius[:-1]

    for n in range(0, N):
        for m in range(0, M):
            y = np.zeros((3*RADIAL_QUANTIZATION),dtype=float)
            y[0:RADIAL_QUANTIZATION] = QRad[n,m,:]
            y[RADIAL_QUANTIZATION:(2*RADIAL_QUANTIZATION)] = QRad[n,m,:]
            y[(2*RADIAL_QUANTIZATION):(3*RADIAL_QUANTIZATION)] = QRad[n,m,:]
            x = np.arange(len(y))
            idx = np.nonzero(y)
            x_known = x[idx]
            y_known = y[idx]
            f = interp1d(x_known, y_known, kind='linear', fill_value="extrapolate")
            intp = f(x)
            QRad[n,m,:] = intp[RADIAL_QUANTIZATION:(2*RADIAL_QUANTIZATION)]

    gpuQRad = t.from_numpy(QRad).float().to("cuda")
    gpuQRadMean = t.div( t.sum(gpuQRad, [2]), t.ones([N, M], dtype=t.float32, device=gpuCuda0) * RADIAL_QUANTIZATION)
    gpuQRad = t.div(gpuQRad, gpuQRadMean[:, :, np.newaxis])
    gpuQRad = t.div(gpuQRad, t.ones([N, M, RADIAL_QUANTIZATION], dtype=t.float32, device=gpuCuda0) * RADIAL_QUANTIZATION)
    #TEST3 = np.float32( gpuQRad.cpu().detach().numpy() )

    # QRad = np.float32( gpuQRad.cpu().detach().numpy() )

    #Calculate fft
    gpuFFT = t.fft.fft(gpuQRad, dim=2)
    #FFT_real = np.float32( gpuFFT.real.cpu().detach().numpy() )
    #FFT_imag = np.float32( gpuFFT.imag.cpu().detach().numpy() )

    #Combine real and imag numbers
    gpuHarmonic = t.add( t.square(gpuFFT.real) , t.square(gpuFFT.imag) )
    gpuHarmonic = t.sqrt(gpuHarmonic)

    #Calculate low freq
    gpuFFT_low = t.sum(gpuHarmonic[:,:,3:7],[2])
    freq_low_fft = np.float32( gpuFFT_low.cpu().detach().numpy() )
    #Calculate high freq
    gpuFFT_high = t.sum(gpuHarmonic[:,:,7:32],[2])
    freq_high_fft = np.float32( gpuFFT_high.cpu().detach().numpy() )
    #Finish harmonic amplitude

    harmonics = np.float32( gpuHarmonic.cpu().detach().numpy() )

    # FFT_Test_in = (np.arange(32)+1)*3;
    # gpuFFTtest_in = t.from_numpy(FFT_Test_in).float().to("cuda")
    # gpuFFTtest_out = t.fft.fft(gpuFFTtest_in)
    # gpuFFTtest_out_real = gpuFFTtest_out.real
    # gpuFFTtest_out_imag = gpuFFTtest_out.imag
    # FFT_Test_out_real = np.float32( gpuFFTtest_out_real.cpu().detach().numpy() )
    # FFT_Test_out_imag = np.float32( gpuFFTtest_out_imag.cpu().detach().numpy() )

    #Calculate Elongation with second harmonic
    gpuAVER = t.div( gpuFFT.real[:,:,0], t.ones([N, M], dtype=t.float32, device=gpuCuda0) * 2)
    gpuC = gpuHarmonic[:,:,2]

    gpuNumi = t.add(gpuAVER, gpuC)
    gpuDenomi = t.sub(gpuAVER, gpuC)
    gpuDenomi[(gpuDenomi==0)] = 1
    gpuELONG = t.div(gpuNumi, gpuDenomi)
    elongation = np.float32( gpuELONG.cpu().detach().numpy() )

    #Cleanup to prevent VRAM memory leaks
    del [gpuMask, gpuPerimeter, gpuGx, gpuGy, gpuArea, gpuPArea, gpuXind, gpuYind, gpuMean_x, gpuMean_y ]
    del [gpuX2, gpuY2, gpuX, gpuY, gpuRad, gpuRad_NZ, gpuT]
    del [gpuW, gpuW2, gpuQRad, gpuQRadMean, gpuQCount, gpuQMask, gpuFFT]
    del [gpuHarmonic, gpuFFT_low, gpuFFT_high]
    del [gpuAVER, gpuC, gpuNumi, gpuDenomi, gpuELONG, gpuCuda0]
    t.cuda.empty_cache()

    return [harmonics, freq_low_fft, freq_high_fft, elongation]

def fft(TYPE, QRad, num):
    N = len(QRad)
    fft_real = np.zeros(N)
    fft_imag = np.zeros(N)
    data = np.zeros((2*N)+2)
    for n in range( 0, N):
        nn = n*2
        data[nn+2] = 0.0;
        data[nn+1] = QRad[n];

    fsize=int(0)
    i=int(0)
    m=int(0)
    mmax=int(0)
    step=int(0)

    theta = 0
    wi = 0
    wr = 0
    j = int(1)

    fsize = 2 * num
    for i in range(1, fsize+1, 2):
        if(i<j):
            tempr = data[j]
            tempi = data[j+1]
            data[j] = data[i]
            data[j + 1] = data[i + 1]
            data[i] = tempr
            data[i + 1] = tempi

        m = int(fsize / 2)

        while ((m > 1) and (j > m)):
            j -= int(m)
            m /= int(2)

        j += int(m);

    mmax = 2
    while (mmax < fsize):
        step = 2 * mmax
        for m in range(1, mmax+1, 2):
            if (mmax != 0):
                theta = math.pi * TYPE * (1 - m) / mmax

            wr = math.cos(theta)
            wi = math.sin(theta)
            for i in range(m, fsize+1, step):
                j = i + mmax
                tempr = (wr * data[j]) - (wi * data[j + 1])
                tempi = (wr * data[j + 1]) + (wi * data[j])
                data[j] = data[i] - tempr
                data[j + 1] = data[i + 1] - tempi
                data[i] = data[i] + tempr
                data[i + 1] = data[i + 1] + tempi


        mmax = step;

    for n in range( 0, N ):
        nn = n*2
        fft_real[n] = data[nn+1]
        fft_imag[n] = data[nn+2]

    return [fft_real, fft_imag]


def ODFeatures(ImageTensor, MaskTensor):
    #Calculates features below:
        #IOD            = total OD
        #ODmean         = mean OD
        #ODvar          = variance OD
        #ODskew         = skewness OD
        #ODkurt         = kurtosis OD

    [N, Y, X, P] = ImageTensor.shape
    [N, Y, X, M] = MaskTensor.shape
    gpuCuda0 = t.device('cuda:0')

    BKG = np.zeros((N, P), dtype='float32')
    IOD = np.zeros((N, P, M), dtype='float32')
    AvgIOD = np.zeros((P, M), dtype='float32')
    ODmean = np.zeros((N, P, M), dtype='float32')
    ODvar = np.zeros((N, P, M), dtype='float32')
    ODskew = np.zeros((N, P, M), dtype='float32')
    ODkurt = np.zeros((N, P, M), dtype='float32')

    #Load Image and Mask data to GPU
    gpuImage = t.from_numpy(ImageTensor).float().to("cuda")
    gpuMask = t.from_numpy(MaskTensor).float().to("cuda")

    gpuBKG = t.mode( t.reshape(gpuImage, (N, -1, P) ), 1)
    BKG = np.float32( gpuBKG.cpu().detach().numpy() )

    #Calculate pixel count (Area)
    gpuArea = t.sum( gpuMask, [1,2] )
    gpuArea[(gpuArea<=1)] = 2

    #Allocate space for OD data
    gpuODImage = t.zeros([N, X, Y, P, M], dtype=t.float32, device=gpuCuda0)

    gpuODImage = gpuImage[:,:,:,:,np.newaxis]
    gpuODImage[(gpuODImage==0)] = 1
    gpuODImage = t.mul( t.ones([N, X, Y, P, M], dtype=t.float32, device=gpuCuda0) * -1 , t.log10(t.div(gpuODImage, t.ones([N, X, Y, P, M], dtype=t.float32, device=gpuCuda0) * 255)) )
    gpuODImage = t.mul( gpuODImage, gpuMask[:,:,:,np.newaxis,:] )

    gpuIOD = t.sum(gpuODImage, [1, 2])
    IOD = np.float32( gpuIOD.cpu().detach().numpy() )

    gpuAvgIOD = t.div( t.sum(gpuIOD, [0]), t.ones([P, M], dtype=t.float32, device=gpuCuda0) * N)
    AvgIOD = np.float32( gpuAvgIOD.cpu().detach().numpy() )

    gpuODmean = t.div( gpuIOD, gpuArea[:, np.newaxis, :])
    ODmean = np.float32( gpuODmean.cpu().detach().numpy() )

    gpuODImage = t.mul(  t.sub(gpuODImage, gpuODmean[:, np.newaxis, np.newaxis, :, :]),  gpuMask[:,:,:,np.newaxis,:] )
    gpuODvar = t.sum(t.mul(gpuODImage, gpuODImage), [1,2])
    gpuODskew = t.sum(t.mul(gpuODImage, t.mul(gpuODImage, gpuODImage)), [1,2])
    gpuODkurt = t.sum(t.mul(gpuODImage, t.mul(gpuODImage, t.mul(gpuODImage, gpuODImage))), [1,2])

    gpuODvar2 = t.div(gpuODvar , gpuArea[:,np.newaxis,:] )
    gpuODskew = t.div(gpuODskew , gpuArea[:,np.newaxis,:] )
    gpuODkurt = t.div(gpuODkurt , gpuArea[:,np.newaxis,:] )

    gpuODvar = t.div(gpuODvar , t.square(gpuODmean) )
    ODvar = np.float32( gpuODvar.cpu().detach().numpy() )

    gpuODskew = t.div(gpuODskew , t.mul(t.sqrt(gpuODvar2), gpuODvar2) )
    ODskew = np.float32( gpuODskew.cpu().detach().numpy() )

    gpuODkurt = t.div(gpuODkurt , t.square(gpuODvar2) )
    ODkurt = np.float32( gpuODkurt.cpu().detach().numpy() )

    del [gpuCuda0, gpuImage, gpuMask, gpuArea, gpuBKG]
    del [gpuODImage, gpuIOD, gpuAvgIOD, gpuODmean]
    del [gpuODvar, gpuODvar2, gpuODskew, gpuODkurt]
    t.cuda.empty_cache()

    return[IOD, AvgIOD, ODmean, ODvar, ODskew, ODkurt]

def TextureFeatures():

    return[-1]