#Writen by Paul Gallagher (BC Cancer Research Centre)

#--------------------------------Initialize------------------------------------
from UNET_Function_readWriteTrainingSets import read_maskpairs_tiffs, write_tiffs, write_valtrain_tiffs, saveNames, loadNames
from UNET_Function_expandSamples import expandSamples
from UNET_Function_save_training_times import save_training_times
from UNET_Structures import UNet_B, accuracy

#import os
import time
import numpy as np
#from PIL import Image
import torch
#from torch.utils.data import DataLoader
#from torchvision import datasets, transforms
#from torch.nn.functional import leaky_relu, sigmoid
from torch import nn
from torch import optim
from tqdm import tqdm

import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#For deterministic behavoir lock all seeds
SEED = 42

np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)

# Force PyTorch to use deterministic algorithms where available
torch.use_deterministic_algorithms(True)

# Configure cuDNN to be deterministic
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False


total_start = time.time()

#--------------------------------Hyperparameters-------------------------------
print('Initializing Hyperparameters...')
X = 128
Y = 128

#in_data = r'A:\UNET_NucleiOverlap\UNET_Dataspace\Data\CenterPoints\64x64diamond'
in_data = r'G:\DATA\UNET_B_Trainingsets\Test'
outpath = r'G:\DATA\UNET_B_Trainingsets\Test'

initial_weights = r'G:\DATA\UNET_B_Trainingsets\UNET_B.pth'
save_weights = outpath + r'\UNET_B.pth'

#Training
train_epochs = 5
learning_rate = 3e-4
train_batch_size = 50
predict_batch_size = 50

#Tloss = 'mse'
Tloss = 'binary_crossentropy'
#Tloss = jaccard_distance_loss

Nfeatures = 32

isTrain = 1
loadTiffs = 0
useInitWeights = 1

isBackgroundDark = 0
SaveLossAccPlot = 1
isSGD = 0
isExpand = 0

#Check for GPU
device = (
    "cuda"
    if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu"
)
print(f"Using {device} device")

#Clear GPU memory from possible previous runs
torch.cuda.empty_cache()
#------------------------------------------------------------------------------

#-------------------------Load/Augment Training Samples------------------------
#[Xtrain, Ytrain, XVal, YVal, Xpredict] = read_tiffs(X, Y, r'C:\UNET_NucleiOverlap\UNET_Dataspace')
print('Loading Input Images...')
start_read = time.time()
if(not(loadTiffs)):
    Xtrain = np.load( in_data + r'\Xtrain.npy' )
    Ytrain = np.load( in_data + r'\Ytrain.npy' )
    XVal = np.load( in_data + r'\XVal.npy' )
    YVal = np.load( in_data + r'\YVal.npy' )
    Xpredict = np.load( in_data + r'\Xpredict.npy' )
    trainnames = loadNames( in_data + r'\trainnames.txt' )
    valnames = loadNames( in_data + r'\valnames.txt' )
    testnames = loadNames( in_data + r'\testnames.txt' )
else:
    [Xtrain, Ytrain, XVal, YVal, Xpredict, testnames, trainnames, valnames] = read_maskpairs_tiffs(X, Y, in_data, isBackgroundDark)
    np.save( in_data + r'\Xtrain.npy', Xtrain )
    np.save( in_data + r'\Ytrain.npy', Ytrain )
    np.save( in_data + r'\XVal.npy', XVal )
    np.save( in_data + r'\YVal.npy', YVal )
    np.save( in_data + r'\Xpredict.npy', Xpredict )
    saveNames( in_data + r'\trainnames.txt', trainnames )
    saveNames( in_data + r'\valnames.txt', valnames )
    saveNames( in_data + r'\testnames.txt', testnames )


#[Xtrain, Ytrain, XVal, YVal, Xpredict, testnames, trainnames, valnames] = read_maskpairs_tiffs(X, Y, in_data, isBackgroundDark)
end_read = time.time()
readtime = end_read - start_read
print('Load Time: ' + str(readtime))

Ytrain = np.uint8(Ytrain/255)
YVal = np.uint8(YVal/255)

print('Expanding Input Images...')
# dataSize = Xtrain.shape
if (isExpand & isTrain):
    eXtrain = expandSamples(Xtrain)
    eYtrain = expandSamples(Ytrain)
    eXVal = expandSamples(XVal)
    eYVal = expandSamples(YVal)
else:
    eXtrain = Xtrain
    eYtrain = Ytrain
    eXVal = XVal
    eYVal = YVal

NTrain, H, W = Xtrain.shape
NVal, H, W = XVal.shape
NPredict, H, W = Xpredict.shape
NeTrain, H, W = eXtrain.shape
NeVal, H, W = eXVal.shape

Xtrain = Xtrain.reshape((NTrain, 1, H, W))
Ytrain = Ytrain.reshape((NTrain, 1, H, W))
XVal = XVal.reshape((NVal, 1, H, W))
YVal = YVal.reshape((NVal, 1, H, W))
Xpredict = Xpredict.reshape((NPredict, 1, H, W))

eXtrain = eXtrain.reshape((NeTrain, 1, H, W))
eYtrain = eYtrain.reshape((NeTrain, 1, H, W))
eXVal = eXVal.reshape((NeVal, 1, H, W))
eYVal = eYVal.reshape((NeVal, 1, H, W))


NBatchTrain = np.int64(np.ceil(NeTrain/train_batch_size))
NBatchVal = np.int64(np.ceil(NeVal/train_batch_size))
NBatchPredict = np.int64(np.ceil(NPredict/predict_batch_size))
NPBatchTrain = np.int64(np.ceil(NTrain/predict_batch_size))
NPBatchVal = np.int64(np.ceil(NVal/predict_batch_size))
#----------------------------------------------------------------------------
model = UNet_B(Nfeatures).to(device)

if(isSGD):
    optimizer = optim.SGD(model.parameters(), lr=learning_rate)
else:
    optimizer = optim.AdamW(model.parameters(), lr=learning_rate)
#criterion = nn.BCEWithLogitsLoss()
criterion = nn.BCELoss()
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

start_train = time.time()

if(isTrain):
    if (useInitWeights):
        model.load_state_dict(torch.load(initial_weights, weights_only=True, map_location=torch.device(device)))

    train_losses = []
    train_dcs = []
    val_losses = []
    val_dcs = []

    for epoch in tqdm(range(0, train_epochs)):
        model.train()
        train_running_loss = 0
        train_running_dc = 0

        Tshuffler = np.random.permutation(len(eXtrain))
        eXtrain = eXtrain[Tshuffler]
        eYtrain = eYtrain[Tshuffler]

        Vshuffler = np.random.permutation(len(eXVal))
        eXVal = eXVal[Vshuffler]
        eYVal = eYVal[Vshuffler]

        tbi = 0
        vbi = 0

        for i in tqdm(range(0, NBatchTrain), position=0, leave=True):
            if(i == (NBatchTrain-1)):
                img = eXtrain[tbi:,[0],:,:]
                mask = eYtrain[tbi:,[0],:,:]
            else:
                img = eXtrain[tbi:(tbi+train_batch_size),[0],:,:]
                mask = eYtrain[tbi:(tbi+train_batch_size),[0],:,:]
                tbi=tbi+train_batch_size

            Timg = torch.from_numpy(img).float().to(device)
            Tmask = torch.from_numpy(mask).float().to(device)

            y_pred = model(Timg)
            optimizer.zero_grad()

            #dc = dice_coefficient(y_pred, Tmask)
            dc = accuracy(y_pred, Tmask)
            loss = criterion(y_pred, Tmask)

            train_running_loss += loss.item()
            train_running_dc += dc.item()

            loss.backward()
            optimizer.step()

        train_loss = train_running_loss / (i + 1)
        train_dc = train_running_dc / (i + 1)

        train_losses.append(train_loss)
        train_dcs.append(train_dc)

        model.eval()
        val_running_loss = 0
        val_running_dc = 0

        with torch.no_grad():
            for i in tqdm(range(0, NBatchVal), position=0, leave=True):
                if(i == (NBatchVal-1)):
                    img = eXVal[vbi:,[0],:,:]
                    mask = eYVal[vbi:,[0],:,:]
                else:
                    img = eXVal[vbi:(vbi+train_batch_size),[0],:,:]
                    mask = eYVal[vbi:(vbi+train_batch_size),[0],:,:]
                    vbi=vbi+train_batch_size

                Timg = torch.from_numpy(img).float().to(device)
                Tmask = torch.from_numpy(mask).float().to(device)

                y_pred = model(Timg)
                loss = criterion(y_pred, Tmask)
                #dc = dice_coefficient(y_pred, Tmask)
                dc = accuracy(y_pred, Tmask)

                val_running_loss += loss.item()
                val_running_dc += dc.item()

            val_loss = val_running_loss / (i + 1)
            val_dc = val_running_dc / (i + 1)

        val_losses.append(val_loss)
        val_dcs.append(val_dc)

        print("-" * 30)
        print(f"Training Acc EPOCH {epoch + 1}: {train_dc:.4f}")
        print(f"Validation Acc EPOCH {epoch + 1}: {val_dc:.4f}")
        #print("\n")
        print(f"Training Loss EPOCH {epoch + 1}: {train_loss:.4f}")
        print(f"Validation Loss EPOCH {epoch + 1}: {val_loss:.4f}")
        print("-" * 30)
        torch.save(model.state_dict(), save_weights)
    del Timg, Tmask, y_pred
    torch.cuda.empty_cache()
else:
    print('Loading weights...')
    model.load_state_dict(torch.load(initial_weights, map_location=torch.device(device)))
# Saving the model
#torch.save(model.state_dict(), wpath)
if(isTrain):
    torch.save(model.state_dict(), save_weights)
    print("Saved model weights to disk")

end_train = time.time()
traintime = end_train - start_train
print('Train Time: ' + str(traintime))

#trained_model = UNet_A(Nfeatures).to(device)
#trained_model.load_state_dict(torch.load(out_wpath, map_location=torch.device(device)))



start_predict = time.time()

Ypredict = np.zeros((NPredict, 1, X, Y),dtype='float')
Ypredict_train = np.zeros((NTrain, 1, X, Y),dtype='float')
Ypredict_val = np.zeros((NVal, 1, X, Y),dtype='float')
pbi=0
tbi=0
vbi=0
with torch.no_grad():
    for i in tqdm(range(0, NBatchPredict), position=0, leave=True):
        if(i == (NBatchPredict-1)):
            img = Xpredict[pbi:,[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict[pbi:,[0],:,:] = y_pred.cpu().detach().numpy()
        else:
            img = Xpredict[pbi:(pbi+predict_batch_size),[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict[pbi:(pbi+predict_batch_size),[0],:,:] = y_pred.cpu().detach().numpy()
            pbi=pbi+predict_batch_size
        #Ypredict = np.concatenate((Ypredict, y_pred.cpu().detach().numpy()), axis=0)

    for i in tqdm(range(0, NPBatchTrain), position=0, leave=True):
        if(i == (NPBatchTrain-1)):
            img = Xtrain[tbi:,[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict_train[tbi:,[0],:,:] = y_pred.cpu().detach().numpy()
        else:
            img = Xtrain[tbi:(tbi+predict_batch_size),[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict_train[tbi:(tbi+predict_batch_size),[0],:,:] = y_pred.cpu().detach().numpy()
            tbi=tbi+predict_batch_size
        #Ypredict_train = np.concatenate((Ypredict_train, y_pred.cpu().detach().numpy()), axis=0)

    for i in tqdm(range(0, NPBatchVal), position=0, leave=True):
        if(i == (NPBatchVal-1)):
            img = XVal[vbi:,[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict_val[vbi:,[0],:,:] = y_pred.cpu().detach().numpy()
        else:
            img = XVal[vbi:(vbi+predict_batch_size),[0],:,:]
            Timg = torch.from_numpy(img).float().to(device)
            y_pred = model.forward(Timg)
            Ypredict_val[vbi:(vbi+predict_batch_size),[0],:,:] = y_pred.cpu().detach().numpy()
            vbi=vbi+predict_batch_size
        #Ypredict_val = np.concatenate((Ypredict_val, y_pred.cpu().detach().numpy()), axis=0)


    del Timg, y_pred
    torch.cuda.empty_cache()

#Ypredict = np.delete(Ypredict ,0 ,0)
#Ypredict_train = np.delete(Ypredict_train ,0 ,0)
#Ypredict_val = np.delete(Ypredict_val ,0 ,0)

end_predict = time.time()
predicttime = end_predict - start_predict
print('Prediction Time: ' + str(predicttime))

print('Saving predicted...')
start_write = time.time()
write_tiffs(X, Y, Xpredict, Ypredict, testnames, in_data, isBackgroundDark)
write_valtrain_tiffs(X, Y, Ypredict_train, Xtrain , Ytrain, trainnames, True, in_data, isBackgroundDark)
write_valtrain_tiffs(X, Y, Ypredict_val, XVal, YVal, valnames, False, in_data, isBackgroundDark)

end_write = time.time()
writetime = end_write - start_write



if(SaveLossAccPlot & isTrain):
    epochs_list = list(range(1, train_epochs + 1))
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(epochs_list, train_losses, label='Training Loss')
    plt.plot(epochs_list, val_losses, label='Validation Loss')
    plt.xticks(ticks=list(range(1, train_epochs + 1, 1)))
    plt.title('Loss over epochs')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.grid()
    plt.tight_layout()

    plt.legend()


    plt.subplot(1, 2, 2)
    plt.plot(epochs_list, train_dcs, label='Training DICE')
    plt.plot(epochs_list, val_dcs, label='Validation DICE')
    plt.xticks(ticks=list(range(1, train_epochs + 1, 1)))
    plt.title('DICE Coefficient over epochs')
    plt.xlabel('Epochs')
    plt.ylabel('DICE')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    savename = in_data + r'\\TrainingPlots.tif'
    plt.savefig(savename, dpi=600, transparent=1)
    plt.show()


total_end = time.time()
totaltime = total_end - total_start

save_training_times( (outpath + r'\\Centers_A_Runtimes.xlsx'), totaltime,readtime, traintime, predicttime, writetime)

print("Done")

#------------------------------------------------------------------------------