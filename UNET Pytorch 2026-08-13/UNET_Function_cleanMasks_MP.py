#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np
import copy
import multiprocessing as mp
import scipy.ndimage as ndimage
from skimage import measure

def cleanMasks_MP(isRemove, X, Y, mask_Thresh, border_mirror, Raws, Masks, CenterCrops, xp , yp, numProc):

    N = len(Masks)
    if(numProc==0):
        numProc=12
        if(N<2000):
            numProc=8
        if(N>1000):
            numProc=4
        if(N<100):
            numProc=1

    out_Raws = mp.Queue()
    out_Masks = mp.Queue()
    out_CenterCrops = mp.Queue()
    out_xp = mp.Queue()
    out_yp = mp.Queue()

    procs = []
    B = np.int64( np.ceil(N/numProc) );
    cCell = 0;
    for i in range(0, numProc):
        if(i+1 != numProc):
            Rw = Raws[cCell:(cCell+B), :, :]
            Msk = Masks[cCell:(cCell+B), :, :]
            CCrp = CenterCrops[cCell:(cCell+B), :, :]
            xxp = xp[cCell:(cCell+B)]
            yyp = yp[cCell:(cCell+B)]
            cCell = cCell+B;

            p = mp.Process(
                target=cleanMask_Process,
                args=( isRemove, X, Y, mask_Thresh, border_mirror, Rw, Msk, CCrp, xxp , yyp, out_Raws, out_Masks, out_CenterCrops, out_xp , out_yp, i ) )


            procs.append(p)
            p.start()

        else:
            Rw = Raws[cCell:, :, :]
            Msk = Masks[cCell:, :, :]
            CCrp = CenterCrops[cCell:, :, :]
            xxp = xp[cCell:]
            yyp = yp[cCell:]


            p = mp.Process(
                target=cleanMask_Process,
                args=( isRemove, X, Y, mask_Thresh, border_mirror, Rw, Msk, CCrp, xxp , yyp, out_Raws, out_Masks, out_CenterCrops, out_xp , out_yp, i ) )

            procs.append(p)
            p.start()
            i=i

            resultdict_Raws = {}
            resultdict_Masks = {}
            resultdict_CenterCrops = {}
            resultdict_xp = {}
            resultdict_yp = {}
            for i in range(0, numProc):
                resultdict_Raws.update(out_Raws.get())
                resultdict_Masks.update(out_Masks.get())
                resultdict_CenterCrops.update(out_CenterCrops.get())
                resultdict_xp.update(out_xp.get())
                resultdict_yp.update(out_yp.get())

            for p in procs:
                p.join()


            for i in range(0, numProc):
                if(0==i):
                    rRaws = resultdict_Raws[i]
                    rMasks = resultdict_Masks[i]
                    rCenterCrops = resultdict_CenterCrops[i]
                    rxp = resultdict_xp[i]
                    ryp = resultdict_yp[i]
                else:
                    rRaws = np.append(rRaws, resultdict_Raws[i], axis=0)
                    rMasks = np.append(rMasks, resultdict_Masks[i], axis=0)
                    rCenterCrops = np.append(rCenterCrops, resultdict_CenterCrops[i], axis=0)
                    rxp = np.append(rxp, resultdict_xp[i], axis=0)
                    ryp = np.append(ryp, resultdict_yp[i], axis=0)

    return [rRaws, rMasks, rCenterCrops, rxp, ryp]

def cleanMask_Process(isRemove, X, Y, mask_Thresh, border_mirror, Raws, Masks, CenterCrops, xp , yp, out_Raws, out_Masks, out_CenterCrops, out_xp , out_yp, procID):

    #Suppress border predictions 3 pixels deep
    Masks[:, 0:3, :, :] = 0
    Masks[:, :, 0:3, :] = 0
    Masks[:, -3:, :, :] = 0
    Masks[:, :, -3:, :] = 0

    numImage = len(Masks)
    bit8Masks = np.uint8(Masks)

    Raws_dict = {}
    Masks_dict = {}
    CenterCrops_dict = {}
    xp_dict = {}
    yp_dict = {}

    removecount = 0;
    for n in range(0, numImage):
        nn = numImage - n -1
        mask = Masks[nn,:,:,0]
        H = mask.shape[0]
        W = mask.shape[1]

        #Set values > mask_Thresh to 1
        mask = mask >= mask_Thresh

        #Fill holes in masks
        mask = ndimage.binary_fill_holes(mask)

        if(np.sum(mask) >= H*W*0.9):
            mask = np.zeros((H,W), dtype='uint8')

        #Remove mask parts not connected to the center mask piece. If no mask at center use largest mask.
        [mask_label, numlabels] = measure.label(mask, return_num=True, connectivity=1 )
        center_label = mask_label[int((H/2)-1), int((W/2)-1)];
        if( center_label > 0):
            mask = mask_label == center_label
        else:
            mask = np.zeros((H,W), dtype='uint8')

        bit8Masks[nn,:,:,0] = mask

    Masks = copy.deepcopy(bit8Masks)
    del bit8Masks


    if(isRemove):

        for n in range(0, numImage):
            removeCell = 0
            nn = numImage - n -1
            mask = Masks[nn,:,:,0]

            #Delete nuclei with 10 or less pixels
            if(np.sum(mask) <= 10):
                removeCell = 1

            #Delete nuclei with masks 2 pixels from border
            if( (np.sum(Masks[nn,0:2,:,0]) > 0) or (np.sum(Masks[nn,:,0:2,0]) > 0) or (np.sum(Masks[nn,-2:,:,0]) > 0) or (np.sum(Masks[nn,:,-2:,0]) > 0) ):
                removeCell = 1

            #Delete nuclei Outside image border in mirrored area
            if(removeCell==0):
                maskH = np.sum(mask, 1);
                maskW = np.sum(mask, 0);

                maskH = maskH > 0
                maskW = maskW > 0
                bit = True
                for x in range(0, len(maskW)):

                    if(maskW[x]==True and bit==True):
                        left = x
                        bit = False
                    if(maskW[x]==False and bit==False):
                        right = x-1
                        bit = True
                if(bit == False):
                    right = x;
                bit = True
                for y in range(0, len(maskH)):

                    if(maskH[y]==True and bit==True):
                        top = y
                        bit = False
                    if(maskH[y]==False and bit==False):
                        bottom = y-1
                        bit = True
                if(bit == False):
                    bottom = y;

                #if( border_mirror >= (xp[nn]-(W/2))+left or border_mirror+X <= (xp[nn]-(W/2))+right or border_mirror >= (yp[nn]-(H/2))+top or border_mirror+Y <= (yp[nn]-(H/2))+bottom ):
                #if( border_mirror >= xp[nn]+left-W or border_mirror+X <= xp[nn]+right+W or border_mirror >= yp[nn]+top-H or border_mirror+Y <= yp[nn]+bottom+H ):
                if( border_mirror >= (xp[nn]-(W/2))+left-3 or border_mirror+X <= (xp[nn]-(W/2))+right+3 or border_mirror >= (yp[nn]-(H/2))+top-3 or border_mirror+Y <= (yp[nn]-(H/2))+bottom+3 ):
                    removeCell = 1


            #Remove the cell
            if(removeCell):
                Raws = np.delete(Raws, nn, 0)
                Masks = np.delete(Masks, nn, 0)
                CenterCrops = np.delete(CenterCrops, nn, 0)
                xp = np.delete(xp, nn, 0)
                yp = np.delete(yp, nn, 0)
                removecount=removecount+1;

        numImage = len(Masks)

    Raws_dict[procID] = Raws
    Masks_dict[procID] = Masks
    CenterCrops_dict[procID] = CenterCrops
    xp_dict[procID] = xp
    yp_dict[procID] = yp

    out_Raws.put(Raws_dict)
    out_Masks.put(Masks_dict)
    out_CenterCrops.put(CenterCrops_dict)
    out_xp.put(xp_dict)
    out_yp.put(yp_dict)

    return [Raws, Masks, CenterCrops, xp, yp]