#Writen by Paul Gallagher (BC Cancer Research Centre)

import numpy as np
import multiprocessing as mp

def proximityFilter_MP(XpredictB, YpredictB, CenterCrops, xp, yp, B_X, B_Y, prox_distance, iThresh, numProc):

    N, H, W, C = XpredictB.shape
    cellsum = np.sum(YpredictB, (1, 2, 3))
    plist = [[]] * N
    #removelist = []

    if(numProc==0):
        numProc=12
        if(N<2000):
            numProc=8
        if(N>1000):
            numProc=4
        if(N<100):
            numProc=1


    #Find cells within prox
    for n in range(0, N):
        xpos = xp[n]
        ypos = yp[n]

        xdiff = xp[(n+1):] - xpos
        ydiff = yp[(n+1):] - ypos

        prox = np.sqrt( (xdiff*xdiff) + (ydiff*ydiff) )

        I = np.array( np.where( prox < prox_distance ) )
        II = I+(n+1)

        if(II.shape[1]>=1):
            tlist = np.ones((II.shape[1], 2)) * n
            tlist[:, 1] = II
            plist[n] = tlist
            del tlist




    #-----MP-------------------------------------------------------------------
    out_removelist = mp.Queue()

    procs = []
    L = np.int64( np.ceil(N/numProc) );
    O = 0;
    for i in range(0, numProc):
        if(i+1 != numProc):

            p = mp.Process(
                target=proximityFilter_Process,
                args=( O, L, YpredictB, plist, cellsum, xp, yp, B_X, B_Y, iThresh, out_removelist, i ) )


            procs.append(p)
            p.start()
            O = O+L;
        else:

            LL = N - O

            p = mp.Process(
                target=proximityFilter_Process,
                args=( O, LL, YpredictB, plist, cellsum, xp, yp, B_X, B_Y, iThresh, out_removelist, i ) )

            procs.append(p)
            p.start()

            resultdict_Remove = {}

            for i in range(0, numProc):
                resultdict_Remove.update(out_removelist.get())

            for p in procs:
                p.join()


            for i in range(0, numProc):
                if(0==i):
                    removelist = resultdict_Remove[i]
                else:
                    removelist = np.append(removelist, resultdict_Remove[i], axis=0)

    #--------------------------------------------------------------------------



    sorted_removelist = list(dict.fromkeys(removelist))
    M = len(sorted_removelist)
    nXpredictB = np.zeros((N-M, H, W, C), dtype='uint8')
    nYpredictB = np.zeros((N-M, H, W, 1), dtype='uint8')
    nCenterCrops = np.zeros((N-M, H, W, 1), dtype='uint8')
    nxp = np.zeros(N-M, dtype='int64')
    nyp = np.zeros(N-M, dtype='int64')

    m=0
    for n in range(0, N):
        if(n in sorted_removelist):
            n=n#do nothing
        else:
            nXpredictB[m,:,:,:] = XpredictB[n,:,:,:]
            nYpredictB[m,:,:,0] = YpredictB[n,:,:,0]
            nCenterCrops[m,:,:,0] = CenterCrops[n,:,:,0]
            nxp[m] = xp[n]
            nyp[m] = yp[n]
            m=m+1


    return [nXpredictB, nYpredictB, nCenterCrops, nxp, nyp]

def proximityFilter_Process(O, L, YpredictB,  plist, cellsum, xp, yp, B_X, B_Y, iThresh, out_Remove, procID):
    removelist = []
    Remove_dict = {}

    #Determine which cells to remove
    for n in range(0, L):
        num = len(plist[O+n])
        if(num >= 1):
            nNum = np.int32(plist[O+n][0, 0])
            nScreenX = np.uint64( xp[nNum]-(B_X/2) )
            nScreenY = np.uint64( yp[nNum]-(B_Y/2) )
            nRangeX = [int(nScreenX), int(nScreenX+B_X)]
            nRangeY = [int(nScreenY), int(nScreenY+B_Y)]

            for m in range(0, num):
                mNum = np.int32(plist[O+n][m, 1])
                mScreenX = np.uint64( xp[mNum]-(B_X/2) )
                mScreenY = np.uint64( yp[mNum]-(B_Y/2) )
                mRangeX = [int(mScreenX), int(mScreenX+B_X)]
                mRangeY = [int(mScreenY), int(mScreenY+B_Y)]

                a_frame = np.zeros( (max( nRangeY[1], mRangeY[1]), max(nRangeX[1], mRangeX[1]) ), 'uint8');
                a_frame[ nRangeY[0]:nRangeY[1], nRangeX[0]:nRangeX[1] ] += YpredictB[nNum,:,:,0];
                a_frame[ mRangeY[0]:mRangeY[1], mRangeX[0]:mRangeX[1] ] += YpredictB[mNum,:,:,0];

                ind = np.where( a_frame > 1 )
                iSum = len(ind[0])
                nSum = cellsum[nNum]
                mSum = cellsum[mNum]
                nfrac = iSum/nSum
                mfrac = iSum/mSum

                if( (nfrac > iThresh) or (mfrac > iThresh) ):
                    if(nfrac >= mfrac):
                        removelist.append(plist[O+n][m, 0])
                    elif(nfrac < mfrac):
                        removelist.append(plist[O+n][m, 1])

            del a_frame

    Remove_dict[procID] = removelist
    out_Remove.put(Remove_dict)

    return removelist










def proximityFilter(XpredictB, YpredictB, CenterCrops, xp, yp, B_X, B_Y, prox_distance, iThresh):

    N, H, W, C = XpredictB.shape
    plist = [[]] * (N-1)
    removelist = []

    #Find cells within prox
    for n in range(0, (N-1)):
        xpos = xp[n]
        ypos = yp[n]

        xdiff = xp[(n+1):] - xpos
        ydiff = yp[(n+1):] - ypos

        prox = np.sqrt( (xdiff*xdiff) + (ydiff*ydiff) )

        I = np.array( np.where( prox < prox_distance ) )
        II = I+(n+1)

        if(II.shape[1]>=1):
            tlist = np.ones((II.shape[1], 2)) * n
            tlist[:, 1] = II
            plist[n] = tlist
            del tlist


    #Determine which cells to remove
    for n in range(0, (N-1)):
        num = len(plist[n])
        if(num >= 1):
            nNum = np.int32(plist[n][0, 0])
            nScreenX = np.uint64( xp[nNum]-(B_X/2) )
            nScreenY = np.uint64( yp[nNum]-(B_Y/2) )
            nRangeX = [int(nScreenX), int(nScreenX+B_X)]
            nRangeY = [int(nScreenY), int(nScreenY+B_Y)]

            for m in range(0, num):
                mNum = np.int32(plist[n][m, 1])
                mScreenX = np.uint64( xp[mNum]-(B_X/2) )
                mScreenY = np.uint64( yp[mNum]-(B_Y/2) )
                mRangeX = [int(mScreenX), int(mScreenX+B_X)]
                mRangeY = [int(mScreenY), int(mScreenY+B_Y)]

                nIOU_frame = np.zeros( (max( nRangeY[1], mRangeY[1]), max(nRangeX[1], mRangeX[1]) ), 'uint8');
                mIOU_frame = np.zeros( (max( nRangeY[1], mRangeY[1]), max(nRangeX[1], mRangeX[1]) ), 'uint8');

                nIOU_frame[ nRangeY[0]:nRangeY[1], nRangeX[0]:nRangeX[1] ] = YpredictB[nNum,:,:,0];
                mIOU_frame[ mRangeY[0]:mRangeY[1], mRangeX[0]:mRangeX[1] ] = YpredictB[mNum,:,:,0];


                a_frame = nIOU_frame + mIOU_frame;
                #iSum = sum(sum(a_frame-1));
                ind = np.where( a_frame >  1 )
                iSum = len(ind[0])
                nSum = sum(sum( np.uint32( YpredictB[nNum,:,:,0] ) ));
                mSum = sum(sum( np.uint32( YpredictB[mNum,:,:,0] ) ));

                if( (iSum/nSum > iThresh) or (iSum/mSum > iThresh) ):
                    if((iSum/nSum) >= (iSum/mSum)):
                        #removelist = cat(1, removelist, plist{n}(m, 1) );
                        removelist.append(plist[n][m, 0])
                    elif((iSum/nSum) < (iSum/mSum)):
                        #removelist = cat(1, removelist, plist{n}(m, 2) );
                        removelist.append(plist[n][m, 1])

            del a_frame, nIOU_frame, mIOU_frame

    sorted_removelist = list(dict.fromkeys(removelist))
    M = len(sorted_removelist)
    nXpredictB = np.zeros((N-M, H, W, C), dtype='uint8')
    nYpredictB = np.zeros((N-M, H, W, 1), dtype='uint8')
    nCenterCrops = np.zeros((N-M, H, W, 1), dtype='uint8')
    nxp = np.zeros(N-M, dtype='int64')
    nyp = np.zeros(N-M, dtype='int64')

    m=0
    for n in range(0, N):
        if(n in sorted_removelist):
            n=n#do nothing
        else:
            nXpredictB[m,:,:,:] = XpredictB[n,:,:,:]
            nYpredictB[m,:,:,0] = YpredictB[n,:,:,0]
            nCenterCrops[m,:,:,0] = CenterCrops[n,:,:,0]
            nxp[m] = xp[n]
            nyp[m] = yp[n]
            m=m+1


    #return -1
    #return [XpredictB, YpredictB, CenterCrops, xp, yp]
    return [nXpredictB, nYpredictB, nCenterCrops, nxp, nyp]