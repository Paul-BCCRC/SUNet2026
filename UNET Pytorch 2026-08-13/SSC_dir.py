#Writen by Paul Gallagher (BC Cancer Research Centre)
"""
Created on Wed Sep  2 14:58:36 2026

@author: pgallagher
"""
if __name__ == '__main__': #REQUIRED TO PREVENT INFINITE FORKING

    import os
    import time
    from read_write_CMG import readCMG#, writeCMG, CMG_header
    from read_write_CCG_MultiProcessing import readCCG#, writeCCG#, CCG_header

    from SpatialSpilloverCompensation import SSC#, makeCMGHeaderStructure
    #from CytoFeatures import makeHeaderFeatures
    from read_write_CBX import writeCBX


    #imagefilepath=r'A:\UNET_NucleiOverlap\UNET_SCANTESTS\LCR_Thionin_Unmixed\Test'
    #imagefilepath = r'D:\Test2'
    #imagefilepath = r'F:\Thionin_feb2026\UNMIXED\03_VS20-27501_F4\TileUnmixing(2-Spectra_(2024-2-1)(12;45.7.042).xlsx)_Scalefactor=1.7 1.2'
    imagefilepath = r'G:\DATA\SpilloverCompensationTEST\Out\2026_08_25TEST\CMG2CB4_TEST'

    #imagefilename = r'Scaled_Top_Left_C23-01C1-AbII-Subtracted_crop.tiff'
    #imagefilename = r'Thionin_InvertedConcentrationMap(^).ome.TIFF'
    #imagefilename = r'001_NORM.TIFF'
    #imagefilename = r'Scaled_Top_Left_C23-01C1-AbII-Subtracted_crop3'

    isCMG = 1
    SSCweight = 0.75

    outpath = imagefilepath
    #outpath = r'F:\Thionin_feb2026\UNMIXED\03_VS20-27501_F4\TileUnmixing(2-Spectra_(2024-2-1)(12;45.7.042).xlsx)_Scalefactor=1.7 1.2'

    #imagefilename = '\\' + imagefilename
    start_ = time.time()

    if(isCMG):
        cmg_files = [f for f in os.listdir(imagefilepath) if f.endswith(r'cmg')]
        N = len(cmg_files)
        print( str(N) + ' CMG(s) found.')
        print( r'Calculating features on CMGs')
        for n in range(0, N):
            imagefilename = cmg_files[n]
            imagefilename = imagefilename[:-4]
            print( r'Calculating features on ... ' + imagefilename)
            [Images, Masks, Header] = readCMG(imagefilepath, imagefilename)
            [flabels, fvalues] = SSC(Images, Masks, Header, isCMG, -1, -1, 0, -1, -1, SSCweight)
            writeCBX(flabels, fvalues, outpath, imagefilename, r'cb4')

    else:
        ccg_files = [f for f in os.listdir(imagefilepath) if f.endswith(r'ccg')]
        N = len(ccg_files)
        print( str(N) + ' CCG(s) found.')
        print( r'Calculating features on CCGs')
        for n in range(0, N):
            imagefilename = ccg_files[n]
            imagefilename = imagefilename[:-4]
            print( r'Calculating features on ... ' + imagefilename)
            [Images, Masks, Header, ImageLabel, MaskLabel, ilist] = readCCG(imagefilepath, imagefilename)
            [flabels, fvalues] = SSC(Images, Masks, Header, isCMG, ImageLabel, MaskLabel, 0, -1, -1, SSCweight)
            writeCBX(flabels, fvalues, outpath, imagefilename, r'cb4')


    end_ = time.time()
    print('Segmentation time: ' + str(end_ - start_))
    print('---------------------------------')