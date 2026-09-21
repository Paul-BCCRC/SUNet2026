#Writen by Paul Gallagher (BC Cancer Research Centre)

#import numpy as np
import torch
from torch import nn
#from torch.utils.data import DataLoader
#from torchvision import datasets, transforms
from torch.nn.functional import  leaky_relu, sigmoid #,relu
import torch.nn.init as init


A_useAttentionNet = 1
B_useAttentionNet = 1

ReLualpha = 0.1                           #Slope for negative portion of "Leaky" ReLu. A zero value is the same as normal ReLu.
isBias = 0                                #Use bias offsets as trainable weights

class AttentionBlock(nn.Module):
    def __init__(self, f_g, f_l, f_int):
        super().__init__()
        self.psi = nn.Sequential(
            nn.Conv2d(f_int, 1, kernel_size=1, stride=1, padding=0,  bias=True),
            #nn.BatchNorm2d(1),
            nn.Sigmoid(),
        )
        self.relu = nn.ReLU(inplace=True)
    def forward(self, g, x):
        psi = self.relu(g+x)
        psi = self.psi(psi)

        return psi*x

#-----UNET A-------------------------------------------------------------------
class UNet_A(nn.Module):
        def __init__(self, Nfeatures):
            super().__init__()

            self.e11 = nn.Conv2d(1, Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.e12 = nn.Conv2d(Nfeatures, Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.pool1 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e21 = nn.Conv2d(Nfeatures, 2*Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.e22 = nn.Conv2d(2*Nfeatures, 2*Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.pool2 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e31 = nn.Conv2d(2*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.e32 = nn.Conv2d(4*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.pool3 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e41 = nn.Conv2d(4*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.e42 = nn.Conv2d(8*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.pool4 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e51 = nn.Conv2d(8*Nfeatures, 16*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.e52 = nn.Conv2d(16*Nfeatures, 16*Nfeatures, kernel_size=3, padding=1, bias=isBias)


            self.upconv1 = nn.ConvTranspose2d(16*Nfeatures, 8*Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(A_useAttentionNet):
                self.att5 = AttentionBlock(f_g=8*Nfeatures, f_l=8*Nfeatures, f_int=8*Nfeatures)
            self.d11 = nn.Conv2d(16*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.d12 = nn.Conv2d(8*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)

            self.upconv2 = nn.ConvTranspose2d(8*Nfeatures, 4*Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(A_useAttentionNet):
                self.att4 = AttentionBlock(f_g=4*Nfeatures, f_l=4*Nfeatures, f_int=4*Nfeatures)
            self.d21 = nn.Conv2d(8*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.d22 = nn.Conv2d(4*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)

            self.upconv3 = nn.ConvTranspose2d(4*Nfeatures, 2*Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(A_useAttentionNet):
                self.att3 = AttentionBlock(f_g=2*Nfeatures, f_l=2*Nfeatures, f_int=2*Nfeatures)
            self.d31 = nn.Conv2d(4*Nfeatures, 2*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.d32 = nn.Conv2d(2*Nfeatures, 2*Nfeatures, kernel_size=3, padding=1, bias=isBias)

            self.upconv4 = nn.ConvTranspose2d(2*Nfeatures, Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(A_useAttentionNet):
                self.att2 = AttentionBlock(f_g=Nfeatures, f_l=Nfeatures, f_int=Nfeatures)
            self.d41 = nn.Conv2d(2*Nfeatures, Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.d42 = nn.Conv2d(Nfeatures, Nfeatures, kernel_size=3, padding=1, bias=isBias)

            # Output layer
            self.outconv = nn.Conv2d(Nfeatures, 1, kernel_size=1)

            self.apply(self._init_weights)


        def _init_weights(self, module):
            """Recursively checks the layer type and applies custom init routines."""

            if isinstance(module, nn.Conv2d):
                init.kaiming_normal_(module.weight, mode='fan_out', nonlinearity='leaky_relu')
                #init.ones_(module.weight)
                #init.constant_(module.weight, 1.0)
                if module.bias is not None:
                    init.constant_(module.bias, 0.0)

            # elif isinstance(module, nn.Linear):
            #     init.xavier_uniform_(module.weight)
            #     if module.bias is not None:
            #         init.constant_(module.bias, 0.0)

            # elif isinstance(module, nn.BatchNorm2d):
            #     init.constant_(module.weight, 1.0)
            #     init.constant_(module.bias, 0.0)

            #Configure network structure
        def forward(self, x):
            # Encoder
            xe11 = leaky_relu(self.e11(x), ReLualpha)
            xe12 = leaky_relu(self.e12(xe11), ReLualpha)
            xp1 = self.pool1(xe12)

            xe21 = leaky_relu(self.e21(xp1), ReLualpha)
            xe22 = leaky_relu(self.e22(xe21), ReLualpha)
            xp2 = self.pool2(xe22)

            xe31 = leaky_relu(self.e31(xp2), ReLualpha)
            xe32 = leaky_relu(self.e32(xe31), ReLualpha)
            xp3 = self.pool3(xe32)

            xe41 = leaky_relu(self.e41(xp3), ReLualpha)
            xe42 = leaky_relu(self.e42(xe41), ReLualpha)
            xp4 = self.pool4(xe42)

            xe51 = leaky_relu(self.e51(xp4), ReLualpha)
            xe52 = leaky_relu(self.e52(xe51), ReLualpha)

            # Decoder
            xu1 = self.upconv1(xe52)
            if(A_useAttentionNet):
                xe42 = self.att5(g=xu1, x=xe42)
            xu11 = torch.cat([xu1, xe42], dim=1)
            xd11 = leaky_relu(self.d11(xu11), ReLualpha)
            xd12 = leaky_relu(self.d12(xd11), ReLualpha)

            xu2 = self.upconv2(xd12)
            if(A_useAttentionNet):
                xe32 = self.att4(g=xu2, x=xe32)
            xu22 = torch.cat([xu2, xe32], dim=1)
            xd21 = leaky_relu(self.d21(xu22), ReLualpha)
            xd22 = leaky_relu(self.d22(xd21), ReLualpha)

            xu3 = self.upconv3(xd22)
            if(A_useAttentionNet):
                xe22 = self.att3(g=xu3, x=xe22)
            xu33 = torch.cat([xu3, xe22], dim=1)
            xd31 = leaky_relu(self.d31(xu33), ReLualpha)
            xd32 = leaky_relu(self.d32(xd31), ReLualpha)

            xu4 = self.upconv4(xd32)
            if(A_useAttentionNet):
                xe12 = self.att2(g=xu4, x=xe12)
            xu44 = torch.cat([xu4, xe12], dim=1)
            xd41 = leaky_relu(self.d41(xu44), ReLualpha)
            xd42 = leaky_relu(self.d42(xd41), ReLualpha)

            # Output layer
            out = sigmoid(self.outconv(xd42))

            return out
#-----UNET B-------------------------------------------------------------------
class UNet_B(nn.Module):
    def __init__(self, Nfeatures):
            super().__init__()

            self.e11 = nn.Conv2d(1, Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.e12 = nn.Conv2d(Nfeatures, Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.pool1 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e21 = nn.Conv2d(Nfeatures, 2*Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.e22 = nn.Conv2d(2*Nfeatures, 2*Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.pool2 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e31 = nn.Conv2d(2*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.e32 = nn.Conv2d(4*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.pool3 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e41 = nn.Conv2d(4*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.e42 = nn.Conv2d(8*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.pool4 = nn.MaxPool2d(kernel_size=2, stride=2)

            self.e51 = nn.Conv2d(8*Nfeatures, 16*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.e52 = nn.Conv2d(16*Nfeatures, 16*Nfeatures, kernel_size=3, padding=1, bias=isBias)


            self.upconv1 = nn.ConvTranspose2d(16*Nfeatures, 8*Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(B_useAttentionNet):
                self.att5 = AttentionBlock(f_g=8*Nfeatures, f_l=8*Nfeatures, f_int=8*Nfeatures)
            self.d11 = nn.Conv2d(16*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.d12 = nn.Conv2d(8*Nfeatures, 8*Nfeatures, kernel_size=3, padding=1, bias=isBias)

            self.upconv2 = nn.ConvTranspose2d(8*Nfeatures, 4*Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(B_useAttentionNet):
                self.att4 = AttentionBlock(f_g=4*Nfeatures, f_l=4*Nfeatures, f_int=4*Nfeatures)
            self.d21 = nn.Conv2d(8*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)
            self.d22 = nn.Conv2d(4*Nfeatures, 4*Nfeatures, kernel_size=3, padding=1, bias=isBias)

            self.upconv3 = nn.ConvTranspose2d(4*Nfeatures, 2*Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(B_useAttentionNet):
                self.att3 = AttentionBlock(f_g=2*Nfeatures, f_l=2*Nfeatures, f_int=2*Nfeatures)
            self.d31 = nn.Conv2d(4*Nfeatures, 2*Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.d32 = nn.Conv2d(2*Nfeatures, 2*Nfeatures, kernel_size=5, padding=2, bias=isBias)

            self.upconv4 = nn.ConvTranspose2d(2*Nfeatures, Nfeatures, kernel_size=2, stride=2, bias=isBias)
            if(B_useAttentionNet):
                self.att2 = AttentionBlock(f_g=Nfeatures, f_l=Nfeatures, f_int=Nfeatures)
            self.d41 = nn.Conv2d(2*Nfeatures, Nfeatures, kernel_size=5, padding=2, bias=isBias)
            self.d42 = nn.Conv2d(Nfeatures, Nfeatures, kernel_size=5, padding=2, bias=isBias)

            # Output layer
            self.outconv = nn.Conv2d(Nfeatures, 1, kernel_size=1)

        #Configure network structure
    def forward(self, x):
                # Encoder
                xe11 = leaky_relu(self.e11(x), ReLualpha)
                xe12 = leaky_relu(self.e12(xe11), ReLualpha)
                xp1 = self.pool1(xe12)

                xe21 = leaky_relu(self.e21(xp1), ReLualpha)
                xe22 = leaky_relu(self.e22(xe21), ReLualpha)
                xp2 = self.pool2(xe22)

                xe31 = leaky_relu(self.e31(xp2), ReLualpha)
                xe32 = leaky_relu(self.e32(xe31), ReLualpha)
                xp3 = self.pool3(xe32)

                xe41 = leaky_relu(self.e41(xp3), ReLualpha)
                xe42 = leaky_relu(self.e42(xe41), ReLualpha)
                xp4 = self.pool4(xe42)

                xe51 = leaky_relu(self.e51(xp4), ReLualpha)
                xe52 = leaky_relu(self.e52(xe51), ReLualpha)

                # Decoder
                xu1 = self.upconv1(xe52)
                if(B_useAttentionNet):
                    xe42 = self.att5(g=xu1, x=xe42)
                xu11 = torch.cat([xu1, xe42], dim=1)
                xd11 = leaky_relu(self.d11(xu11), ReLualpha)
                xd12 = leaky_relu(self.d12(xd11), ReLualpha)

                xu2 = self.upconv2(xd12)
                if(B_useAttentionNet):
                    xe32 = self.att4(g=xu2, x=xe32)
                xu22 = torch.cat([xu2, xe32], dim=1)
                xd21 = leaky_relu(self.d21(xu22), ReLualpha)
                xd22 = leaky_relu(self.d22(xd21), ReLualpha)

                xu3 = self.upconv3(xd22)
                if(B_useAttentionNet):
                    xe22 = self.att3(g=xu3, x=xe22)
                xu33 = torch.cat([xu3, xe22], dim=1)
                xd31 = leaky_relu(self.d31(xu33), ReLualpha)
                xd32 = leaky_relu(self.d32(xd31), ReLualpha)

                xu4 = self.upconv4(xd32)
                if(B_useAttentionNet):
                    xe12 = self.att2(g=xu4, x=xe12)
                xu44 = torch.cat([xu4, xe12], dim=1)
                xd41 = leaky_relu(self.d41(xu44), ReLualpha)
                xd42 = leaky_relu(self.d42(xd41), ReLualpha)

                # Output layer
                out = sigmoid(self.outconv(xd42))

                return out
#------------------------------------------------------------------------------

def dice_coefficient(prediction, target, epsilon=1e-07):
    prediction_copy = prediction.clone()

    prediction_copy[prediction_copy < 0] = 0
    prediction_copy[prediction_copy > 0] = 1

    intersection = abs(torch.sum(prediction_copy * target))
    union = abs(torch.sum(prediction_copy) + torch.sum(target))
    dice = (2. * intersection + epsilon) / (union + epsilon)

    return dice

def accuracy(prediction, target):
    N, C, H, W = target.shape
    #prediction_copy = prediction.clone()

    #prediction_copy[prediction_copy < 0.5] = 0
    #prediction_copy[prediction_copy >= 0.5] = 1

    #train_acc = torch.sum(prediction == target)
    train_acc = torch.sum(abs(prediction - target))
    acc = 1 - (train_acc / (H*W*N))
    return acc