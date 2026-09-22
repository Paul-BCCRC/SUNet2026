#Writen by Paul Gallagher (BC Cancer Research Centre)

import os
import numpy as np
import skimage.io as skio
from PIL import Image

def get_dir_png_names(batch_size, folder=r'A:\UNET_NucleiOverlap\UNET_Dataspace'):
    names = ['o'];
    imgpath = folder + r'\image\img';
    for fnameXtrain in os.listdir(imgpath):
        #imXtrain = TIFF.open(os.path.join(folder, fnameXtrain))
        names = names + [fnameXtrain[0:(len(fnameXtrain))]]

    names = names[1:]
    batches = np.uint32( np.ceil( len(names) / batch_size ) )
    #batchnames = ['o']

    if(len(names) > batch_size):
        for n in range(0, batches):
            if(n==0):
                batchnames = [names[0:batch_size]]
            else:
                if(len(names) > batch_size):
                    batchnames.append( names[0:batch_size] )
                else:
                    batchnames.append( names )

            names = names[batch_size+1:]
    else:
        batchnames = [names]

    return [batchnames, batches]


def read_tiffs(X, Y, folder=r'A:\UNET_NucleiOverlap\UNET_Dataspace'):
    dirXtrain = folder + r'\Train\X'
    dirYtrain = folder + r'\Train\Y'
    dirXval = folder + r'\Validation\X'
    dirYval = folder + r'\Validation\Y'
    dirXtest = folder + r'\Test\X'

    finalXtrain = np.zeros((1, X, Y))
    finalYtrain = np.zeros((1, X, Y))
    finalXval = np.zeros((1, X, Y))
    finalYval = np.zeros((1, X, Y))
    finalXtest = np.zeros((1, X, Y))

    for fnameXtrain in os.listdir(dirXtrain):
        #imXtrain = TIFF.open(os.path.join(dirXtrain, fnameXtrain))
        #imXarray = imXtrain.read_image()
        imXarray = skio.imread(os.path.join(dirXtrain, fnameXtrain))
#        imXarray = np.invert(np.array(imXtrain))
        imXvec = imXarray.reshape(1, X, Y)
        finalXtrain = np.concatenate((finalXtrain, imXvec), axis=0)
    for fnameYtrain in os.listdir(dirYtrain):
        #imYtrain = TIFF.open(os.path.join(dirYtrain, fnameYtrain))
        imYarray = skio.imread(os.path.join(dirYtrain, fnameYtrain))
#        imYarray = np.invert(np.array(imYtrain))
        imYvec = imYarray.reshape(1, X, Y)
        finalYtrain = np.concatenate((finalYtrain, imYvec), axis=0)
    for fnameXval in os.listdir(dirXval):
        #imXval = TIFF.open(os.path.join(dirXval, fnameXval))
        imXvalarray = skio.imread(os.path.join(dirXval, fnameXval))
#        imXvalarray = np.invert(np.array(imXval))
        imXvec = imXvalarray.reshape(1, X, Y)
        finalXval = np.concatenate((finalXval, imXvec), axis=0)
    for fnameYval in os.listdir(dirYval):
        #imYval = TIFF.open(os.path.join(dirYval, fnameYval))
        imYvalarray = skio.imread(os.path.join(dirYval, fnameYval))
#        imYvalarray = np.invert(np.array(imYval))
        imYvec = imYvalarray.reshape(1, X, Y)
        finalYval = np.concatenate((finalYval, imYvec), axis=0)
    for fnameXtest in os.listdir(dirXtest):
        #imXtest = TIFF.open(os.path.join(dirXtest, fnameXtest))
        imXtestarray = skio.imread(os.path.join(dirXtest, fnameXtest))
#        imXtestarray = np.invert(np.array(imXtest))
        imXvec = imXtestarray.reshape(1, X, Y)
        finalXtest = np.concatenate((finalXtest, imXvec), axis=0)

    finalXtrain = np.delete(finalXtrain ,0 ,0)
    finalYtrain = np.delete(finalYtrain ,0 ,0)
    finalXval = np.delete(finalXval ,0 ,0)
    finalYval = np.delete(finalYval ,0 ,0)
    finalXtest = np.delete(finalXtest ,0 ,0)

#    img = Image.fromarray(finalXtrain[0,:,:,0])
#    img.show()
    return [finalXtrain, finalYtrain, finalXval, finalYval, finalXtest]

def read_maskpairs_tiffs(X, Y, folder, isBackgroundDark):

    dirtrain = folder + r'\Train'
    dirtest = folder + r'\Test\X'
    dirval = folder + r'\Validation'


    finalXtrain = np.zeros((1, X, Y),dtype='uint8')
    finalYtrain = np.zeros((1, X, Y),dtype='uint8')
    finalXval = np.zeros((1, X, Y),dtype='uint8')
    finalYval = np.zeros((1, X, Y),dtype='uint8')
    finalXtest = np.zeros((1, X, Y),dtype='uint8')
    testnames = ['o'];
    trainnames = ['o'];
    valnames = ['o'];

    for fnametrain in os.listdir(dirtrain):
        #imtrain = TIFF.open(os.path.join(dirtrain, fnametrain))
        #imarray = imtrain.read_image()
        trainnames = trainnames + [fnametrain[1:(len(fnametrain))]]
        imarray = skio.imread(os.path.join(dirtrain, fnametrain))
        mask = np.array(imarray[:,:,0])
#        image = np.array(imarray[:,:,1])
        if(isBackgroundDark):
            image = np.array(imarray[:,:,1])
        else:
            image = np.invert(np.array(imarray[:,:,1]))
        imvecmask = mask.reshape(1, X, Y)
        imvec = image.reshape(1, X, Y)
        finalXtrain = np.concatenate((finalXtrain, imvec), axis=0)
        finalYtrain = np.concatenate((finalYtrain, imvecmask), axis=0)

    for fnameval in os.listdir(dirval):
        #imval = TIFF.open(os.path.join(dirval, fnameval))
        #imarray = imval.read_image()
        valnames = valnames + [fnameval[1:(len(fnameval))]]
        imarray = skio.imread(os.path.join(dirval, fnameval))
        mask = np.array(imarray[:,:,0])
#        image = np.array(imarray[:,:,1])
        if(isBackgroundDark):
            image = np.array(imarray[:,:,1])
        else:
            image = np.invert(np.array(imarray[:,:,1]))
        imvec = image.reshape(1, X, Y)
        imvecmask = mask.reshape(1, X, Y)
        finalXval = np.concatenate((finalXval, imvec), axis=0)
        finalYval = np.concatenate((finalYval, imvecmask), axis=0)

    for fnametest in os.listdir(dirtest):
        #imtest = TIFF.open(os.path.join(dirtest, fnametest))
        #imtestarray = np.invert(imtest.read_image())
        testnames = testnames + [fnametest[1:(len(fnametest))]]
#        imtestarray = np.array(imtest)
        if(isBackgroundDark):
            imtestarray = skio.imread(os.path.join(dirtest, fnametest))
        else:
            imtestarray = np.invert(skio.imread(os.path.join(dirtest, fnametest)))
        imvec = imtestarray.reshape(1, X, Y)
        finalXtest = np.concatenate((finalXtest, imvec), axis=0)

    finalXtrain = np.delete(finalXtrain ,0 ,0)
    finalYtrain = np.delete(finalYtrain ,0 ,0)
    finalXval = np.delete(finalXval ,0 ,0)
    finalYval = np.delete(finalYval ,0 ,0)
    finalXtest = np.delete(finalXtest ,0 ,0)

#    img = Image.fromarray(finalXtrain[0,:,:,0])
#    img.show()
    return [finalXtrain, finalYtrain, finalXval, finalYval, finalXtest, testnames, trainnames, valnames]




def write_tiffs(X, Y, X_tiff_data, Y_tiff_data, names, folder, isBackgroundDark):

    dirYtest = folder + r'\Test\Y'
    dirZtest = folder + r'\Test\Z'
    dataSize = Y_tiff_data.shape

    for i in range(0, dataSize[0]):
        rgb = np.zeros((X, Y, 3), dtype='uint8');

        rgb[:, :, 0] = np.uint8(Y_tiff_data[i, 0, :, :]*255)
        if(isBackgroundDark):
            rgb[:, :, 1] = np.uint8(X_tiff_data[i, 0, :, :])
            rgb[:, :, 2] = np.uint8(X_tiff_data[i, 0, :, :])
        else:
            rgb[:, :, 1] = np.invert(np.uint8(X_tiff_data[i, 0, :, :]))
            rgb[:, :, 2] = np.invert(np.uint8(X_tiff_data[i, 0, :, :]))
        im = Image.fromarray(rgb[:,:,:],'RGB')
        im2 = Image.fromarray(np.uint8(Y_tiff_data[i, 0, :, :]*255))
#        if(~usenames):
#            writepath = dirZtest + '\Z_'+ str(i+1) +'.TIFF'
#            writepath2 = dirYtest + '\Y_'+ str(i+1) +'.TIFF'
#        else:
        writepath = dirZtest + r'\Z' + str(names[i+1])
        writepath2 = dirYtest + r'\Y' + str(names[i+1])
        im.save(writepath)
        im2.save(writepath2)

def write_valtrain_tiffs(X, Y, RED, GREEN, BLUE, names, isTrain, folder, isBackgroundDark):

    if(isTrain):
        dir_path = folder + r'\Test\TrainPredict'
    else:
        dir_path = folder + r'\Test\ValPredict'

    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    dataSize = RED.shape

    for i in range(0, dataSize[0]):
        rgb = np.zeros((X, Y, 3), dtype='uint8')
        #rgb = np.concatenate((X_tiff_data[i, :, :, 0], Y_tiff_data[i, :, :, 0], X_tiff_data[i, :, :, 0]),axis=2)

        if(isBackgroundDark):
            nGREEN = np.float32(np.uint8(GREEN[i, 0, :, :]))
        else:
            nGREEN = np.float32(np.invert(np.uint8(GREEN[i, 0, :, :])))
        #nGREEN = np.float32(np.invert(np.uint8(GREEN[i, 0, :, :])))
        nRED = RED[i, 0, :, :]*255 #+nGREEN
        nBLUE = np.float32(BLUE[i, 0, :, :]*255) #+nGREEN
        rgb[:, :, 0] = np.uint8(np.clip(nRED+nGREEN,0,255))
        rgb[:, :, 1] = np.uint8(nGREEN)
        rgb[:, :, 2] = np.uint8(np.clip(nBLUE+nGREEN,0,255))
        im = []
        im.append(Image.fromarray(rgb[:,:,:],'RGB'))
        im.append(Image.fromarray(np.uint8(nGREEN)))
        im.append(Image.fromarray(np.uint8(nRED)))
        im.append(Image.fromarray(np.uint8(nBLUE)))
        #im = Image.fromarray(rgb[:,:,:],'RGB')

        file, ext = os.path.splitext(names[i])
        f = file.split(".")
        imageName = f[0] + r".TIFF"

        if(isTrain):
            writepath = dir_path + r'\Train' + str(imageName)
        else:
            writepath = dir_path + r'\Val' + str(imageName)

        im[0].save(writepath, save_all=True, append_images=im[1:])
        #im.save(writepath)

def write_tiffs2(X, Y, X_tiff_data, Y_tiff_data, folder):

    if not os.path.exists(folder + r'\Images'):
        os.mkdir(folder + r'\Images')

    dirXtest = folder + r'\Images\X'
    if not os.path.exists(dirXtest):
        os.mkdir(dirXtest)

    dirYtest = folder + r'\Images\Y'
    if not os.path.exists(dirYtest):
        os.mkdir(dirYtest)

    dirZtest = folder + r'\Images\Z'
    if not os.path.exists(dirZtest):
        os.mkdir(dirZtest)

    dataSize = Y_tiff_data.shape

    for i in range(0, dataSize[0]):
        rgb = np.zeros((X, Y, 3), dtype='float32');
        #rgb = np.concatenate((X_tiff_data[i, :, :, 0], Y_tiff_data[i, :, :, 0], X_tiff_data[i, :, :, 0]),axis=2)
        tframe = np.float32(np.invert(X_tiff_data[i, :, :, 0]))
        rgb[:, :, 0] = tframe + np.float32(Y_tiff_data[i, :, :, 0]*255)
        rgb[:, :, 1] = tframe
        rgb[:, :, 2] = tframe
        rgb = np.clip(rgb, 0, 255).astype('uint8')

        im = Image.fromarray(rgb[:,:,:],'RGB')
        im2 = Image.fromarray(np.uint8(Y_tiff_data[i, :, :, 0]*255))
        im3 = Image.fromarray(np.invert(np.uint8(X_tiff_data[i, :, :, 0])))

        writepath = dirZtest + r'\Z_'+ str(i+1) +'.TIFF'
        writepath2 = dirYtest + r'\Y_'+ str(i+1) +'.TIFF'
        writepath3 = dirXtest + r'\X_'+ str(i+1) +'.TIFF'

        im.save(writepath)
        im2.save(writepath2)
        im3.save(writepath3)

def loadNames(namepath):
    names=[]
    with open(namepath, 'r') as fp:
        for line in fp:
            x = line[:-1]
            names.append(x)
    return names

def saveNames(namepath, names):
    with open(namepath, 'w') as fp:
        for item in names:
            fp.write("%s\n" % item)