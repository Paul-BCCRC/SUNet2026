#Writen by Paul Gallagher (BC Cancer Research Centre)

if __name__ == '__main__': #REQUIRED TO PREVENT INFINITE FORKING

    import os
    import time
    #from UNET_Segmentation_Function import UNET_Segmentation_Function
    #from UNET_Segment_CMG_Function import UNET_Segment_CMG_Function
    from UNET_Segmentation_Sequence_Function import UNET_Segmentation_Function
    #from libtiff import TIFF



    imagefilepath=r'G:\DATA\SpilloverCompensationTEST\Out\2026_08_25TEST'
    outpath = r'G:\DATA\SpilloverCompensationTEST\Out\2026_08_25TEST'

    #unmixDir = r'TileUnmixing(2-Spectra_(2024-2-1)(12;45.7.042).xlsx)_Scalefactor=1.7 1.2'

    start_ = time.time()
    files = [fname for fname in os.listdir(imagefilepath) if os.path.isfile(os.path.join(imagefilepath, fname))]
    #files = [fname for fname in os.listdir(imagefilepath) if os.path.isdir(os.path.join(imagefilepath, fname))]
    for fname in files:
        seg_start_ = time.time()

        #-----Segment Folders------------------------------------------------------
        path = imagefilepath + r"\\"
        #filename = r"Thionin_ConcentrationMap(^).ome.TIFF"
        UNET_Segmentation_Function(path, fname)
        #--------------------------------------------------------------------------

        #-----Segment TIFFs--------------------------------------------------------
        # isTiff = False
        # UP = fname.upper()

        # if(UP.find(".TIF") != -1):
        #     isTiff = True;

        # if(isTiff):
        #     fname = '\\' + fname
        #     UNET_Segmentation_Function(imagefilepath, fname)
        #--------------------------------------------------------------------------


        #-----Segment CMGs---------------------------------------------------------
        # isFile = 0
        # UP = fname.upper()
        # if(len(UP) >= 4):
        #     if(UP[-5:].find(".CMG") != -1):
        #         isFile = 1;

        # if(isFile):
        #     UNET_Segment_CMG_Function(imagefilepath, fname[:-4], outpath)
        #--------------------------------------------------------------------------

        seg_end_ = time.time()
        print('Segmentation time: ' + str(seg_end_ - seg_start_))
        print('---------------------------------')



    end_ = time.time()
    print('Finished directory in: ' + str(end_ - start_))