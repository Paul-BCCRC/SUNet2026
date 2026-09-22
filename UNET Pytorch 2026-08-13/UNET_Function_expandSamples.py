import numpy as np

def expandSamples(Data, T=5):
    
    dataSize = Data.shape
    Aug_Data = np.zeros((1, dataSize[1], dataSize[2]), dtype='uint8')
#--------------------------Add flips and rotations-----------------------------


    NewFrame = np.zeros(dataSize, dtype='uint8')
    for angle in range(0, 2):
        for flip in range(0, 4):
            for run in range(0, dataSize[0]):

                    
                frame = np.zeros((dataSize[1], dataSize[2]), dtype='uint8')
                frame = Data[run,:,:]
                    # if(run==0):
                    #     displayframe = np.clip(frame, 0, 255).astype('uint8')
                    #     img = Image.fromarray(displayframe)
                    #     img.show()                    
                    
                if angle==1:
                    frame = np.rot90(frame, k=-1)
                else:
                    frame = frame
                if flip==1:
                    frame = np.flip(frame, 0)
                elif flip==2:
                    frame = np.flip(frame, 1)
                elif flip==3:
                    frame = np.flip(frame, 0)
                    frame = np.flip(frame, 1)
                else:
                    frame = frame
                NewFrame[run, :, :] = frame
    #                displayframe = np.clip(frame, 0, 255).astype('uint8')
    #                img = Image.fromarray(displayframe)
    #                img.show()
            Aug_Data  = np.concatenate((Aug_Data , NewFrame), axis=0)
    
#------------------------------------------------------------------------------
    Aug_Data = np.delete(Aug_Data ,0 ,0)
#---------------------------Add Translations-----------------------------------
    # dataSize = Aug_Data.shape
    # X = dataSize[2];
    # Y = dataSize[1];
    # NewFrame = np.zeros(dataSize, dtype='uint8')
    # for n in range(1, 9):
    #     for run in range(0, dataSize[0]):
    #         frame = np.zeros((dataSize[1], dataSize[2]), dtype='uint8')
    #         if n==1:#Left
    #             Tdata = Aug_Data[run, :, T:X]
    #             Mdata = np.flip(Aug_Data[run, :, (X-T):X], 1)
    #             frame[:, 0:(X-T)] = Tdata;
    #             frame[:, (X-T):X] = Mdata;   
    #         elif n==5:#Right
    #             Tdata = Aug_Data[run, :, 0:(X-T)]
    #             Mdata = np.flip(Aug_Data[run, :, 0:T], 1)
    #             frame[:, T:X] = Tdata;  
    #             frame[:, 0:T] = Mdata;
    #         elif n==3:#Up                   
    #             Tdata = Aug_Data[run, T:Y, :]
    #             Mdata = np.flip(Aug_Data[run, (Y-T):Y, :], 0)
    #             frame[0:(Y-T), :] = Tdata;
    #             frame[(Y-T):Y, :] = Mdata;
    #         elif n==7:#Down
    #             Tdata = Aug_Data[run, 0:(Y-T), :]
    #             Mdata = np.flip(Aug_Data[run, 0:T, :], 0)
    #             frame[T:Y, :] = Tdata;
    #             frame[0:T, :] = Mdata;
    #         elif n==2:#UpLeft     
    #             Tdata = Aug_Data[run, T:Y, :]
    #             Mdata = np.flip(Aug_Data[run, (Y-T):Y, :], 0)
    #             frame[0:(Y-T), :] = Tdata;
    #             frame[(Y-T):Y, :] = Mdata;
                    
    #             Tdata = frame[:, T:X]  
    #             Mdata = np.flip(frame[:, (X-T):X], 1) 
    #             frame[:, 0:(X-T)] = Tdata;
    #             frame[:, (X-T):X] = Mdata;  
    #         elif n==4:#UpRight
    #             Tdata = Aug_Data[run, T:Y, :]
    #             Mdata = np.flip(Aug_Data[run, (Y-T):Y, :], 0)
    #             frame[0:(Y-T), :] = Tdata;
    #             frame[(Y-T):Y, :] = Mdata;
                    
    #             Tdata = frame[:, 0:(X-T)]
    #             Mdata = np.flip(frame[:, 0:T], 1)
    #             frame[:, T:X] = Tdata;
    #             frame[:, 0:T] = Mdata;  
    #         elif n==6:#DownRight  
    #             Tdata = Aug_Data[run, 0:(Y-T), :]
    #             Mdata = np.flip(Aug_Data[run, 0:T, :], 0)
    #             frame[T:Y, :] = Tdata;
    #             frame[0:T, :] = Mdata;
                    
    #             Tdata = frame[:, 0:(X-T)]
    #             Mdata = np.flip(frame[:, 0:T], 1)
    #             frame[:, T:X] = Tdata;
    #             frame[:, 0:T] = Mdata;  
    #         elif n==8:#DownLeft
    #             Tdata = Aug_Data[run, 0:(Y-T), :]
    #             Mdata = np.flip(Aug_Data[run, 0:T, :], 0)
    #             frame[T:Y, :] = Tdata;
    #             frame[0:T, :] = Mdata;
                    
    #             Tdata = frame[:, T:X]   
    #             Mdata = np.flip(frame[:, (X-T):X], 1) 
    #             frame[:, 0:(X-T)] = Tdata;
    #             frame[:, (X-T):X] = Mdata; 
    #         else:
    #             frame = frame
                    
    #         NewFrame[run, :, :] = frame

    #     Aug_Data  = np.concatenate((Aug_Data , NewFrame), axis=0)
            
#------------------------------------------------------------------------------

    return Aug_Data  