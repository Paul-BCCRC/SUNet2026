#Writen by Paul Gallagher (BC Cancer Research Centre)

#import os
import numpy as np
import copy
import time

#import torch as t
#import scipy.ndimage as ndimage
import skimage.morphology

from gpuFeatureCalculator import IntensityFeatures, ShapeFeatures, FFTFeatures, ODFeatures

def CytoLabels(ImageLabel):
    P = len(ImageLabel)
    out_flabels = [r''] * (12 + 2*27 + 3*(7+2)*P )
    count=0

    out_flabels[count+0] = r'Nm_a'
    out_flabels[count+1] = r'Nmd_a'
    out_flabels[count+2] = r'And_a'
    out_flabels[count+3] = r'Cm_a'
    out_flabels[count+4] = r'CmNm_a'
    out_flabels[count+5] = r'OCm_a'
    out_flabels[count+6] = r'nOCm_a'
    out_flabels[count+7] = r'isCmNm'
    out_flabels[count+8] = r'Nm_#ofItract'
    out_flabels[count+9] = r'CmNm_#ofItract'
    out_flabels[count+10] = r'Nm_+fracs'
    out_flabels[count+11] = r'CmNm_+fracs'

    out_flabels[count+12] = r'Nm_parea'
    out_flabels[count+13] = r'Nm_cir'
    out_flabels[count+14] = r'Nm_centroid_x'
    out_flabels[count+15] = r'Nm_centroid_y'
    out_flabels[count+16] = r'Nm_min_r'
    out_flabels[count+17] = r'Nm_max_r'
    out_flabels[count+18] = r'Nm_mean_r'
    out_flabels[count+19] = r'Nm_var_r'
    out_flabels[count+20] = r'Nm_spher'
    out_flabels[count+21] = r'Nm_eccentr'
    out_flabels[count+22] = r'Nm_cell_orient'
    out_flabels[count+23] = r'Nm_inert'
    out_flabels[count+24] = r'Nm_compac'

    out_flabels[count+25] = r'Nm_elongation'
    out_flabels[count+26] = r'Nm_bndvarlow'
    out_flabels[count+27] = r'Nm_bndvarhigh'
    out_flabels[count+28] = r'Nm_fft2'
    out_flabels[count+29] = r'Nm_fft3'
    out_flabels[count+30] = r'Nm_fft4'
    out_flabels[count+31] = r'Nm_fft5'
    out_flabels[count+32] = r'Nm_fft6'
    out_flabels[count+33] = r'Nm_fft7'
    out_flabels[count+34] = r'Nm_fft8'
    out_flabels[count+35] = r'Nm_fft9'
    out_flabels[count+36] = r'Nm_fft10'
    out_flabels[count+37] = r'Nm_fft11'
    out_flabels[count+38] = r'Nm_fft12'

    out_flabels[count+39] = r'Cm_parea'
    out_flabels[count+40] = r'Cm_cir'
    out_flabels[count+41] = r'Cm_centroid_x'
    out_flabels[count+42] = r'Cm_centroid_y'
    out_flabels[count+43] = r'Cm_min_r'
    out_flabels[count+44] = r'Cm_max_r'
    out_flabels[count+45] = r'Cm_mean_r'
    out_flabels[count+46] = r'Cm_var_r'
    out_flabels[count+47] = r'Cm_spher'
    out_flabels[count+48] = r'Cm_eccentr'
    out_flabels[count+49] = r'Cm_cell_orient'
    out_flabels[count+50] = r'Cm_inert'
    out_flabels[count+51] = r'Cm_compac'

    out_flabels[count+52] = r'Cm_elongation'
    out_flabels[count+53] = r'Cm_bndvarlow'
    out_flabels[count+54] = r'Cm_bndvarhigh'
    out_flabels[count+55] = r'Cm_fft2'
    out_flabels[count+56] = r'Cm_fft3'
    out_flabels[count+57] = r'Cm_fft4'
    out_flabels[count+58] = r'Cm_fft5'
    out_flabels[count+59] = r'Cm_fft6'
    out_flabels[count+60] = r'Cm_fft7'
    out_flabels[count+61] = r'Cm_fft8'
    out_flabels[count+62] = r'Cm_fft9'
    out_flabels[count+63] = r'Cm_fft10'
    out_flabels[count+64] = r'Cm_fft11'
    out_flabels[count+65] = r'Cm_fft12'


    count+=66

    for p in range(0, P):
        out_flabels[count+p] = r'Nm_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'Nm_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'Nm_' + ImageLabel[p] + r'_sd'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'Nmd_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'Nmd_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'Nmd_' + ImageLabel[p] + r'_sd'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'And_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'And_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'And_' + ImageLabel[p] + r'_sd'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'Cm_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'Cm_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'Cm_' + ImageLabel[p] + r'_sd'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'CmNm_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'CmNm_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'CmNm_' + ImageLabel[p] + r'_sd'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'OCm_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'OCm_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'OCm_' + ImageLabel[p] + r'_sd'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'nOCm_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'nOCm_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'nOCm_' + ImageLabel[p] + r'_sd'
    count+=P




    for p in range(0, P):
        out_flabels[count+p] = r'SPNm_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'SPAnd_' + ImageLabel[p] + r'_m'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'SPCmNm_' + ImageLabel[p] + r'_m'
    count+=P

    for p in range(0, P):
        out_flabels[count+p] = r'SPNm_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'SPAnd_' + ImageLabel[p] + r'_t'
    count+=P
    for p in range(0, P):
        out_flabels[count+p] = r'SPCmNm_' + ImageLabel[p] + r'_t'
    count+=P

    return out_flabels

def CytoFeatures(ImageTensor, MaskTensor, OMasks, nOMasks, B):
    #Get input data dimentions
    [N, Y, X, P] = ImageTensor.shape
    [N, Y, X, M] = MaskTensor.shape

    #Calculates features below:
        #mean_x         = xmean
        #mean_y         = ymean
        #mean_r         = radial mean
        #max_r          = radial max
        #min_r          = radial min
        #var_r          = radial variance
        #spher          = sphericity
        #eccentr        = eccentricity
        #cell_orient    = cell_orientation
        #inert          = inertia
        #comp           = compact

        #harmonics      = fft amplitudes
        #freq_low_fft   = amplitude sum 3 through 6
        #freq_high_fft  = amplitude sum 7 through 31
        #elongation     = elongation using 2nd harmonic

        #area           = mask sum
        #imean          = intensity average
        #itotal         = intensity sum
        #istd           = intensity standard deviation
    #cuda0 = t.device('cuda:0')

#-----Get Shape features-------------------------------------------------------

    Area = np.zeros((N, M), dtype='float32')
    PArea = np.zeros((N, M), dtype='float32')
    Cir = np.zeros((N, M), dtype='float32')
    Centroid_x = np.zeros((N, M), dtype='float32')
    Centroid_y = np.zeros((N, M), dtype='float32')
    Min_r = np.zeros((N, M), dtype='float32')
    Max_r = np.zeros((N, M), dtype='float32')
    Mean_r = np.zeros((N, M), dtype='float32')
    Var_r = np.zeros((N, M), dtype='float32')
    Spher = np.zeros((N, M), dtype='float32')
    Eccentr = np.zeros((N, M), dtype='float32')
    Cell_orient = np.zeros((N, M), dtype='float32')
    Inert = np.zeros((N, M), dtype='float32')
    Comp = np.zeros((N, M), dtype='float32')

    cCell = 0
    N_Batch = np.int64( np.ceil(N/B) )
    N_Remain = np.mod( N, N_Batch )
    if(N_Remain==0):
        N_Remain = N_Batch

    start_shapef = time.time()
    for b in range(0, B):
        if(b+1 != B):
            [Area[cCell:(cCell+N_Batch), :]
             , PArea[cCell:(cCell+N_Batch), :]
             , Centroid_x[cCell:(cCell+N_Batch), :]
             , Centroid_y[cCell:(cCell+N_Batch), :]
             , Mean_r[cCell:(cCell+N_Batch), :]
             , Max_r[cCell:(cCell+N_Batch), :]
             , Min_r[cCell:(cCell+N_Batch), :]
             , Var_r[cCell:(cCell+N_Batch), :]
             , Spher[cCell:(cCell+N_Batch), :]
             , Eccentr[cCell:(cCell+N_Batch), :]
             , Cell_orient[cCell:(cCell+N_Batch), :]
             , Inert[cCell:(cCell+N_Batch), :]
             , Cir[cCell:(cCell+N_Batch), :]
             , Comp[cCell:(cCell+N_Batch), :]] = ShapeFeatures(MaskTensor[cCell:(cCell+N_Batch), :,:,:])
            cCell = cCell + N_Batch
        else:
            [Area[cCell:, :]
             , PArea[cCell:, :]
             , Centroid_x[cCell:, :]
             , Centroid_y[cCell:, :]
             , Mean_r[cCell:, :]
             , Max_r[cCell:, :]
             , Min_r[cCell:, :]
             , Var_r[cCell:, :]
             , Spher[cCell:, :]
             , Eccentr[cCell:, :]
             , Cell_orient[cCell:, :]
             , Inert[cCell:, :]
             , Cir[cCell:, :]
             , Comp[cCell:, :]] = ShapeFeatures(MaskTensor[cCell:,:,:,:])

    end_shapef = time.time()
    print(r'Time to calculate shape features: ' + str(end_shapef - start_shapef))
#------------------------------------------------------------------------------

#-----Get Fourier features-----------------------------------------------------
    RADIAL_QUANTIZATION = 64
    #Allocate memory for features
    harmonics = np.zeros((N, M, RADIAL_QUANTIZATION), dtype='float32')
    freq_low_fft = np.zeros((N, M), dtype='float32')
    freq_high_fft = np.zeros((N, M), dtype='float32')
    elongation = np.zeros((N, M), dtype='float32')

    cCell = 0
    N_Batch = np.int64( np.ceil(N/B) )
    N_Remain = np.mod( N, N_Batch )
    if(N_Remain==0):
        N_Remain = N_Batch

    start_fftf = time.time()
    for b in range(0, B):
        if(b+1 != B):
            [harmonics[cCell:(cCell+N_Batch), :, :]
             , freq_low_fft[cCell:(cCell+N_Batch), :]
             , freq_high_fft[cCell:(cCell+N_Batch), :]
             , elongation[cCell:(cCell+N_Batch), :] ] = FFTFeatures(MaskTensor[cCell:(cCell+N_Batch), :,:,:], RADIAL_QUANTIZATION)
            cCell = cCell + N_Batch
        else:
            [harmonics[cCell:, :, :]
             , freq_low_fft[cCell:, :]
             , freq_high_fft[cCell:, :]
             , elongation[cCell:, :] ] = FFTFeatures(MaskTensor[cCell:,:,:,:], RADIAL_QUANTIZATION)

    end_fftf = time.time()
    print(r'Time to calculate fft features: ' + str(end_fftf - start_fftf))
#------------------------------------------------------------------------------

#-----Get Intensity features---------------------------------------------------

    #Create space for new mask channels with Nm, Nmd, And, Cm, CmNm
    NewMaskTensor = np.zeros((N,Y,X,7),dtype='uint8')
    NewMaskTensor[:,:,:,0] = MaskTensor[:,:,:,0]
    NewMaskTensor[:,:,:,3] = MaskTensor[:,:,:,1]
    NewMaskTensor[:,:,:,5] = OMasks
    NewMaskTensor[:,:,:,6] = nOMasks

    #As far as I know GPU cannot dilate so do it before hand
    #Dilate each Nuclei by 3


    struct1 = skimage.morphology.disk(3, dtype=np.uint8)
    #struct1 = ndimage.generate_binary_structure(2, 1)
    for n in range(0, N):
        NewMaskTensor[n,:,:,1] = skimage.morphology.dilation(NewMaskTensor[n,:,:,0], footprint=struct1)
        #NewMaskTensor[n,:,:,1] = ndimage.binary_dilation(NewMaskTensor[n,:,:,0], structure=struct1,iterations=3).astype(NewMaskTensor[n,:,:,0].dtype)

    MaskTensor = NewMaskTensor
    del NewMaskTensor

    [N, Y, X, M] = MaskTensor.shape

    #Allocate memory for features

    mArea = np.zeros((N, M), dtype='float32')
    isNmCm = np.zeros((N, 1), dtype='float32')
    iMean = np.zeros((N, P, M), dtype='float32')
    iTotal = np.zeros((N, P, M), dtype='float32')
    iSTD = np.zeros((N, P, M), dtype='float32')

    start_intenf = time.time()
    cCell = 0
    N_Batch = np.int64( np.ceil(N/B) )
    for b in range(0, B):
        if(b+1 != B):
            [mArea[cCell:(cCell+N_Batch), :]
             , iMean[cCell:(cCell+N_Batch), :, :]
             , iTotal[cCell:(cCell+N_Batch), :, :]
             , iSTD[cCell:(cCell+N_Batch), :, :]] = IntensityFeatures(ImageTensor[cCell:(cCell+N_Batch),:,:,:], MaskTensor[cCell:(cCell+N_Batch),:,:,:])
            cCell = cCell + N_Batch
        else:
            [mArea[cCell:, :]
             , iMean[cCell:, :, :]
             , iTotal[cCell:, :, :]
             , iSTD[cCell:, :, :]] = IntensityFeatures(ImageTensor[cCell:,:,:,:], MaskTensor[cCell:,:,:,:])

    isNmCm_diff = mArea[:,0] - mArea[:,3]
    isNmCm[(isNmCm_diff==0)] = 1

    end_intenf = time.time()
    print(r'Time to calculate intensity features: ' + str(end_intenf - start_intenf))
#------------------------------------------------------------------------------

    J = M+5
    fvalues = np.zeros( (N, (J + 2*27 + 3*(M+2)*P)), dtype='float32' )
    #Populate fvalues for CB write
    fvalues[:, 0:M] = mArea
    fvalues[:, M] = isNmCm[:, 0]
    fvalues[:, J+0] = PArea[:, 0]
    fvalues[:, J+1] = Cir[:, 0]
    fvalues[:, J+2] = Centroid_x[:, 0]
    fvalues[:, J+3] = Centroid_y[:, 0]
    fvalues[:, J+4] = Min_r[:, 0]
    fvalues[:, J+5] = Max_r[:, 0]
    fvalues[:, J+6] = Mean_r[:, 0]
    fvalues[:, J+7] = Var_r[:, 0]
    fvalues[:, J+8] = Spher[:, 0]
    fvalues[:, J+9] = Eccentr[:, 0]
    fvalues[:, J+10] = Cell_orient[:, 0]
    fvalues[:, J+11] = Inert[:, 0]
    fvalues[:, J+12] = Comp[:, 0]

    fvalues[:, J+13] = elongation[:,0]
    fvalues[:, J+14] = freq_low_fft[:,0]
    fvalues[:, J+15] = freq_high_fft[:,0]
    fvalues[:, J+16:J+27] = harmonics[:, 0, 2:13]

    fvalues[:, J+27] = PArea[:, 1]
    fvalues[:, J+28] = Cir[:, 1]
    fvalues[:, J+29] = Centroid_x[:, 1]
    fvalues[:, J+30] = Centroid_y[:, 1]
    fvalues[:, J+31] = Min_r[:, 1]
    fvalues[:, J+32] = Max_r[:, 1]
    fvalues[:, J+33] = Mean_r[:, 1]
    fvalues[:, J+34] = Var_r[:, 1]
    fvalues[:, J+35] = Spher[:, 1]
    fvalues[:, J+36] = Eccentr[:, 1]
    fvalues[:, J+37] = Cell_orient[:, 1]
    fvalues[:, J+38] = Inert[:, 1]
    fvalues[:, J+39] = Comp[:, 1]

    fvalues[:, J+40] = elongation[:,1]
    fvalues[:, J+41] = freq_low_fft[:,1]
    fvalues[:, J+42] = freq_high_fft[:,1]
    fvalues[:, J+43:J+54] = harmonics[:, 1, 2:13]

    count = J + 2*27
    for m in range(0, M):
        for p in range(0, P):
            fvalues[:, count+p] = iMean[:,p,m]
        count+=P
        for p in range(0, P):
            fvalues[:, count+p] = iTotal[:,p,m]
        count+=P
        for p in range(0, P):
            fvalues[:, count+p] = iSTD[:,p,m]
        count+=P

    return [iMean, iTotal, mArea, fvalues]




def ShapeFeatures_batch(MaskTensor, B):
    #start_shapef = time.time()
    #Calculates features below:
        #mean_x         = xmean
        #mean_y         = ymean
        #mean_r         = radial mean
        #max_r          = radial max
        #min_r          = radial min
        #var_r          = radial variance
        #spher          = sphericity
        #eccentr        = eccentricity
        #cell_orient    = cell_orientation
        #inert          = inertia
        #comp           = compact
    [N, Y, X, M] = MaskTensor.shape

    #Allocate memory for features
    #fvalues = np.zeros( (N, (24*P + 10)), dtype='float32' )
    Area = np.zeros((N, M), dtype='float32')
    PArea = np.zeros((N, M), dtype='float32')
    Centroid_x = np.zeros((N, M), dtype='float32')
    Centroid_y = np.zeros((N, M), dtype='float32')
    Mean_r = np.zeros((N, M), dtype='float32')
    Max_r = np.zeros((N, M), dtype='float32')
    Min_r = np.zeros((N, M), dtype='float32')
    Var_r = np.zeros((N, M), dtype='float32')
    Spher = np.zeros((N, M), dtype='float32')
    Eccentr = np.zeros((N, M), dtype='float32')
    Cell_orient = np.zeros((N, M), dtype='float32')
    Inert = np.zeros((N, M), dtype='float32')
    Cir = np.zeros((N, M), dtype='float32')
    Comp = np.zeros((N, M), dtype='float32')

    cCell = 0
    N_Batch = np.int64( np.ceil(N/B) )
    N_Remain = np.mod( N, N_Batch )
    if(N_Remain==0):
        N_Remain = N_Batch

    for b in range(0, B):
        if(b+1 != B):
            [Area[cCell:(cCell+N_Batch), :], PArea[cCell:(cCell+N_Batch), :], Centroid_x[cCell:(cCell+N_Batch), :], Centroid_y[cCell:(cCell+N_Batch), :], Mean_r[cCell:(cCell+N_Batch), :], Max_r[cCell:(cCell+N_Batch), :], Min_r[cCell:(cCell+N_Batch), :], Var_r[cCell:(cCell+N_Batch), :], Spher[cCell:(cCell+N_Batch), :], Eccentr[cCell:(cCell+N_Batch), :], Cell_orient[cCell:(cCell+N_Batch), :], Inert[cCell:(cCell+N_Batch), :], Cir[cCell:(cCell+N_Batch), :], Comp[cCell:(cCell+N_Batch), :]] = ShapeFeatures(MaskTensor[cCell:(cCell+N_Batch), :,:,:])
            cCell = cCell + N_Batch
        else:
            [Area[cCell:, :], PArea[cCell:, :], Centroid_x[cCell:, :], Centroid_y[cCell:, :], Mean_r[cCell:, :], Max_r[cCell:, :], Min_r[cCell:, :], Var_r[cCell:, :], Spher[cCell:, :], Eccentr[cCell:, :], Cell_orient[cCell:, :], Inert[cCell:, :], Cir[cCell:, :], Comp[cCell:, :]] = ShapeFeatures(MaskTensor[cCell:,:,:,:])

    #end_shapef = time.time()
    #print(r'Time to calculate shape features: ' + str(end_shapef - start_shapef))

    return [-1]

def FFTFeatures_batch(MaskTensor, QUANT, B):

   [N, Y, X, M] = MaskTensor.shape



   cCell = 0
   N_Batch = np.int64( np.ceil(N/B) )
   N_Remain = np.mod( N, N_Batch )
   if(N_Remain==0):
       N_Remain = N_Batch

   for b in range(0, B):
       if(b+1 != B):
           [R] = FFTFeatures(MaskTensor[cCell:(cCell+N_Batch), :,:,:], QUANT)
           cCell = cCell + N_Batch
       else:
           [R] = FFTFeatures(MaskTensor[cCell:,:,:,:], QUANT)


   return [-1]

def ODFeatures_batch(ImageTensor, MaskTensor, B):
   [N, Y, X, P] = ImageTensor.shape
   [N, Y, X, M] = MaskTensor.shape



   cCell = 0
   N_Batch = np.int64( np.ceil(N/B) )
   N_Remain = np.mod( N, N_Batch )
   if(N_Remain==0):
       N_Remain = N_Batch

   for b in range(0, B):
       if(b+1 != B):
           [R] = ODFeatures(ImageTensor[cCell:(cCell+N_Batch), :,:,:], MaskTensor[cCell:(cCell+N_Batch), :,:,:])
           cCell = cCell + N_Batch
       else:
           [R] = ODFeatures(ImageTensor[cCell:,:,:,:], MaskTensor[cCell:,:,:,:])

   return [-1]

#Function finds cells with overlap(ilist) and the sum of the overlap(slist)
def findOverlap(Masks, header, mask_plane):
    NN = len(Masks);
    ilist = [[]]*NN
    slist = [[]]*NN
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
        sList = copy.deepcopy(L);
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

            Omask = np.zeros((YY+2*Ymax, XX+2*Xmax), 'uint8')
            if(mask_plane>=1):
                cmask = Masks[nn][:,:,mask_plane] - Masks[nn][:,:,0]
            else:
                cmask = Masks[nn][:,:,mask_plane]
            Omask[ (Ymax):(YY+Ymax) , (Xmax):(XX+Xmax) ] = cmask;

            if(mask_plane>=1):
                omask = Masks[ Onum ][:,:,mask_plane] - Masks[ Onum ][:,:,0]
            else:
                omask = Masks[ Onum ][:,:,mask_plane]
            Omask[ (OY):(nYY+OY), (OX):(nXX+OX) ] = Omask[ (OY):(nYY+OY), (OX):(nXX+OX) ] + omask;

            if( np.max(Omask) <= 1 or (List[ll]-1)==nn ):
                List = np.delete(List, ll)
                sList = np.delete(sList, ll)
            else:
                sList[ll] = np.count_nonzero(Omask >= 2)

        ilist[nn] = np.float32(List)
        slist[nn] = np.float32(sList)

    return ilist, slist

def getNonOverlap(Masks, header, mask_plane, ilist):
    [N, Y, X, M] = Masks.shape
    Omasks = np.zeros( (N, Y, X), dtype='uint8' )
    nOmasks = np.zeros( (N, Y, X), dtype='uint8' )
    ilist = ilist[mask_plane]

    for n in range(0, N):
        LL = len(ilist[n])
        Omasks[n,:,:] = Masks[n,:,:,mask_plane]

        SX = np.double( header[3][n] );        #Screenx
        SY = np.double( header[4][n] );        #Screeny
        XX = np.uint32( header[12][n] );       #Width
        YY = np.uint32( header[13][n] );       #Height

        Xmax=0;
        Ymax=0;
        for l in range(0, LL):
            if(Xmax < header[12][ np.uint32( ilist[n][l]-1 )]):
                Xmax = header[12][ np.uint32( ilist[n][l]-1 )];

            if(Ymax < header[13][ np.uint32( ilist[n][l]-1 )]):
                Ymax = header[13][ np.uint32( ilist[n][l]-1 )];

        OOmask = np.zeros((YY+2*Ymax, XX+2*Xmax), 'float32')
        OOmask[ (Ymax):(YY+Ymax) , (Xmax):(XX+Xmax) ] = np.float32( Masks[n,:,:,mask_plane] );
        for l in range(0, LL):
            ll = LL - (l+1);
            Onum = np.uint32(ilist[n][ll]-1);
            nSX = np.double( header[3][ Onum ] );
            nSY = np.double( header[4][ Onum ] );
            nXX = np.uint32( header[12][ Onum ] );
            nYY = np.uint32( header[13][ Onum ] );

            OX = np.int64(Xmax + (nSX - SX));
            OY = np.int64(Ymax + (nSY - SY));

            OOmask[ (OY):(nYY+OY), (OX):(nXX+OX) ] = OOmask[ (OY):(nYY+OY), (OX):(nXX+OX) ] - np.float32( Masks[Onum,:,:,mask_plane] );

        OOmask[OOmask < 0] = 0
        nOmasks[n , :, :] = OOmask[(Ymax):(YY+Ymax) , (Xmax):(XX+Xmax)]
        Omasks[n, :, :] = Omasks[n, :, :] - nOmasks[n , :, :]
        Omasks[Omasks > 1] = 0

    return [Omasks, nOmasks]

def makeHeaderFeatures(Header):
    N = len(Header[0])
    F = 13
    in_flabels = [r''] * F
    count=0
    in_flabels[count] = r'stage_x'
    in_flabels[count+1] = r'stage_y'
    in_flabels[count+2] = r'stage_z'
    in_flabels[count+3] = r'screenx'
    in_flabels[count+4] = r'screeny'
    in_flabels[count+5] = r'stackx'
    in_flabels[count+6] = r'stacky'
    in_flabels[count+7] = r'GROUP'
    in_flabels[count+8] = r'prob'
    in_flabels[count+9] = r'locks'
    in_flabels[count+10] = r'diagnosis'
    in_flabels[count+11] = r'store'
    in_flabels[count+12] = r'num'


    in_fvalues = np.zeros( (N, F), dtype='float32' )

    in_fvalues[:, 0] = np.float32(Header[5])        #StageX
    in_fvalues[:, 1] = np.float32(Header[6])        #StageY
    in_fvalues[:, 2] = np.float32(Header[7])        #StageZ
    in_fvalues[:, 3] = np.float32(Header[3])        #ScreenX
    in_fvalues[:, 4] = np.float32(Header[4])        #ScreenY
    in_fvalues[:, 5] = np.float32(Header[24])       #StackX
    in_fvalues[:, 6] = np.float32(Header[25])       #StackY
    in_fvalues[:, 7] = np.float32(Header[11])       #Group
    in_fvalues[:, 10] = np.float32(Header[17])       #Diagnosis



    return [in_flabels, in_fvalues]