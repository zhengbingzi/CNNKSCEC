# -*- coding: utf8 -*-
#model
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sys
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.metrics import precision_score, recall_score, accuracy_score, f1_score, precision_recall_curve, auc
import warnings
import random
from tqdm import tqdm

class GroupBatchnorm2d(nn.Module):
    def __init__(self, c_num:int, 
                 group_num:int = 16, 
                 eps:float = 1e-10
                 ):
        super(GroupBatchnorm2d,self).__init__()
        assert c_num    >= group_num
        self.group_num  = group_num
        self.weight     = nn.Parameter( torch.randn(c_num, 1, 1)    )
        self.bias       = nn.Parameter( torch.zeros(c_num, 1, 1)    )
        self.eps        = eps
    def forward(self, x):
        N, C, H, W  = x.size()
        x           = x.view(   N, self.group_num, -1   )
        mean        = x.mean(   dim = 2, keepdim = True )
        std         = x.std (   dim = 2, keepdim = True )
        x           = (x - mean) / (std+self.eps)
        x           = x.view(N, C, H, W)
        return x * self.weight + self.bias
    
class ECAAttention(nn.Module):
    def __init__(self, kernel_size=3):
        super(ECAAttention, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool2d(1)
        self.conv = nn.Conv1d(1, 1, kernel_size, padding=(kernel_size - 1) // 2, bias=False)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        y = self.avg_pool(x)  # [B, C, 1, 1]
        y = self.conv(y.squeeze(-1).transpose(-1, -2))  # [B, 1, C]
        y = self.sigmoid(y.transpose(-1, -2).unsqueeze(-1))  # [B, C, 1, 1]
        return x * y.expand_as(x)
    

    
# Convolutional Block Attention Module
class ChannelAttention(nn.Module):
    def __init__(self, in_planes, ratio=8):
        super(ChannelAttention, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool2d(1)
        self.max_pool = nn.AdaptiveMaxPool2d(1)

        self.fc = nn.Sequential(
            nn.Conv2d(in_planes, in_planes // ratio, 1, bias=False),
            nn.ReLU(),
            nn.Conv2d(in_planes // ratio, in_planes, 1, bias=False)
        )

        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        avg_out = self.fc(self.avg_pool(x))
        max_out = self.fc(self.max_pool(x))
        out = avg_out + max_out
        return x * self.sigmoid(out)

class SpatialAttention(nn.Module):
    def __init__(self, kernel_size=7):
        super(SpatialAttention, self).__init__()
        padding = (kernel_size - 1) // 2
        self.conv = nn.Conv2d(2, 1, kernel_size, padding=padding, bias=False)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        avg_out = torch.mean(x, dim=1, keepdim=True)
        max_out, _ = torch.max(x, dim=1, keepdim=True)
        x_cat = torch.cat([avg_out, max_out], dim=1)
        out = self.conv(x_cat)
        return x * self.sigmoid(out)

class CBAM(nn.Module):
    def __init__(self, in_planes, ratio=8, kernel_size=7):
        super(CBAM, self).__init__()
        self.ca = ChannelAttention(in_planes, ratio)
        self.sa = SpatialAttention(kernel_size)

    def forward(self, x):
        out = self.ca(x)
        out = self.sa(out)
        return out
    
class HybridAttention(nn.Module):
    def __init__(self, channel, eca_kernel_size=3, cbam_ratio=8):
        super(HybridAttention, self).__init__()
        self.eca = ECAAttention(kernel_size=eca_kernel_size)
        self.cbam = CBAM(channel, ratio=cbam_ratio)

    def forward(self, x):
        x = self.eca(x)
        x = self.cbam(x)
        return x


class SRU(nn.Module):
    def __init__(self,
                 oup_channels:int, 
                 group_num:int = 16,
                 gate_treshold:float = 0.5,
                 torch_gn:bool = True
                 ):
        super().__init__()
        
        self.gn             = nn.GroupNorm( num_channels = oup_channels, num_groups = group_num ) if torch_gn else GroupBatchnorm2d(c_num = oup_channels, group_num = group_num)
        self.gate_treshold  = gate_treshold
        self.sigomid        = nn.Sigmoid()

    def forward(self,x):
        gn_x        = self.gn(x)
        w_gamma     = self.gn.weight/sum(self.gn.weight)
        w_gamma     = w_gamma.view(1,-1,1,1)
        reweigts    = self.sigomid( gn_x * w_gamma )
        # Gate
        w1          = torch.where(reweigts > self.gate_treshold, torch.ones_like(reweigts), reweigts) # 大于门限值的设为1，否则保留原值
        w2          = torch.where(reweigts > self.gate_treshold, torch.zeros_like(reweigts), reweigts) # 大于门限值的设为0，否则保留原值
        x_1         = w1 * x
        x_2         = w2 * x 
        y           = self.reconstruct(x_1,x_2)
        return y
    
    def reconstruct(self,x_1,x_2):
        x_11,x_12 = torch.split(x_1, x_1.size(1)//2, dim=1)
        x_21,x_22 = torch.split(x_2, x_2.size(1)//2, dim=1)
        return torch.cat([ x_11+x_22, x_12+x_21 ],dim=1)


class CRU(nn.Module):
    '''
    alpha: 0<alpha<1
    '''
    def __init__(self, 
                 op_channel:int,
                 alpha:float = 1/2,
                 squeeze_radio:int = 2 ,
                 group_size:int = 2,
                 group_kernel_size:int = 3,
                 ):
        super().__init__()
        self.up_channel     = up_channel   =   int(alpha*op_channel)
        self.low_channel    = low_channel  =   op_channel-up_channel
        self.squeeze1       = nn.Conv2d(up_channel,up_channel//squeeze_radio,kernel_size=1,bias=False)
        self.squeeze2       = nn.Conv2d(low_channel,low_channel//squeeze_radio,kernel_size=1,bias=False)
        #up
        self.GWC            = nn.Conv2d(up_channel//squeeze_radio, op_channel,kernel_size=group_kernel_size, stride=1,padding=group_kernel_size//2, groups = group_size)
        self.PWC1           = nn.Conv2d(up_channel//squeeze_radio, op_channel,kernel_size=1, bias=False)
        #low
        self.PWC2           = nn.Conv2d(low_channel//squeeze_radio, op_channel-low_channel//squeeze_radio,kernel_size=1, bias=False)
        self.advavg         = nn.AdaptiveAvgPool2d(1)

    def forward(self,x):
        # Split
        up,low  = torch.split(x,[self.up_channel,self.low_channel],dim=1)
        up,low  = self.squeeze1(up),self.squeeze2(low)
        # Transform
        Y1      = self.GWC(up) + self.PWC1(up)
        Y2      = torch.cat( [self.PWC2(low), low], dim= 1 )
        # Fuse
        out     = torch.cat( [Y1,Y2], dim= 1 )
        out     = F.softmax( self.advavg(out), dim=1 ) * out
        out1,out2 = torch.split(out,out.size(1)//2,dim=1)
        return out1+out2


class ScConv(nn.Module):
    def __init__(self,
                op_channel:int,
                group_num:int = 4,
                gate_treshold:float = 0.5,
                alpha:float = 1/2,
                squeeze_radio:int = 2 ,
                group_size:int = 2,
                group_kernel_size:int = 3,
                 ):
        super().__init__()
        self.SRU = SRU( op_channel, 
                       group_num            = group_num,  
                       gate_treshold        = gate_treshold
                        #torch_gn = False 选择的是GroupBatchnorm2d
                         )
        self.CRU = CRU( op_channel, 
                       alpha                = alpha, 
                       squeeze_radio        = squeeze_radio ,
                       group_size           = group_size ,
                       group_kernel_size    = group_kernel_size )
    
    def forward(self,x):
        x = self.SRU(x)
        x = self.CRU(x)
        return x

class DeepCNN(nn.Module):
    """FPN for semantic segmentation"""
    def __init__(self, motiflen=7):
        super(DeepCNN, self).__init__()
        # encode process
        self.conv1 = nn.Conv2d(in_channels=2, out_channels=16, kernel_size=(5, 5), padding=2, bias=True)
        self.pool1 = nn.AvgPool2d(kernel_size=2, stride=2)
        self.conv2 = nn.Conv2d(in_channels=16, out_channels=32, kernel_size=(3, 3), padding=1, bias=True)
        self.pool2 = nn.AvgPool2d(kernel_size=2, stride=2)
        self.conv3 = nn.Conv2d(in_channels=32, out_channels=32, kernel_size=(3, 3), padding=1, bias=True)
        self.pool3 = nn.AvgPool2d(kernel_size=2, stride=2)
        
        
        self.conv1_ = nn.Conv2d(in_channels=2, out_channels=16, kernel_size=(3, 3), stride=2, bias=True)
        self.conv2_ = nn.Conv2d(in_channels=16, out_channels=32, kernel_size=(3, 3), stride=2,padding=1, bias=True)
        self.conv3_ = nn.Conv2d(in_channels=32, out_channels=32, kernel_size=(3, 3), stride=2, bias=True)
        
        self.scconv1 = ScConv(16)
        self.scconv2 = ScConv(32)
        self.scconv3 = ScConv(32)

        self.eca1 = HybridAttention(channel=16) 
        self.eca2 = HybridAttention(channel=32) 
        self.eca3 = HybridAttention(channel=32) 

        # classifier head
        c_in = 128 # 384
        self.bn = nn.BatchNorm2d(32)
        self.linear1 = nn.Linear(c_in, 64)
        self.drop = nn.Dropout(p=0.2)
        self.linear2 = nn.Linear(64, 1)
        # general functions
        self.sigmoid = nn.Sigmoid()
        self.relu = nn.ReLU(inplace=True)
        self.dropout = nn.Dropout(p=0.2)
        self._init_weights()

    def _init_weights(self):
        """Initialize the new built layers"""
        for layer in self.modules():
            if isinstance(layer, (nn.Conv2d, nn.Linear)):
                # nn.init.kaiming_uniform_(layer.weight, mode='fan_in', nonlinearity='relu')
                nn.init.xavier_uniform_(layer.weight)
                if layer.bias is not None:
                    nn.init.constant_(layer.bias, 0)
            elif isinstance(layer, nn.BatchNorm2d):
                nn.init.constant_(layer.weight, 1)
                nn.init.constant_(layer.bias, 0)


    def forward(self, data):
        """Construct a new computation graph at each froward"""
        b, _, _, _ = data.size()
        # encode process
        #print(data.shape) torch.Size([50, 2, 21, 21])
        out1 = self.conv1(data)
        #print(out1.shape)torch.Size([50, 16, 21, 21])
        out1 = self.relu(out1)
        #print(out1.shape)torch.Size([50, 16, 21, 21])
        out1 = self.pool1(out1) 
        #print(out1.shape)torch.Size([50, 16, 10, 10])
        out1_ = self.conv1_(data)
        #print(out1_.shape)torch.Size([50, 16, 10, 10])
        out1_ = self.relu(out1_)
        #print(out1_.shape)torch.Size([50, 16, 10, 10])
        out_e1 = (out1_+out1)/2 
        #print(out_e1.shape)torch.Size([50, 16, 10, 10])
        out_e1 = self.scconv1(out_e1) + out_e1
        #print(out_e1.shape)torch.Size([50, 16, 10, 10])
        out_e1 = self.eca1(out_e1) 

        out_e1 = self.dropout(out_e1)
        #print(out_e1.shape)torch.Size([50, 16, 10, 10])
        skip1 = out_e1
 
        out1 = self.conv2(out_e1)
        #print(out1.shape)torch.Size([50, 32, 10, 10])
        out1 = self.relu(out1)
        #print(out1.shape)torch.Size([50, 32, 10, 10])
        out1 = self.pool2(out1)
        #print(out1.shape)torch.Size([50, 32, 5, 5])
        out1_ = self.conv2_(out_e1)
        #print(out1_.shape)torch.Size([50, 32, 5, 5])
        out1_ = self.relu(out1_)
        #print(out1_.shape)torch.Size([50, 32, 5, 5])
        
        out_e2 = (out1_+out1)/2 
        #print(out_e2.shape)torch.Size([50, 32, 5, 5])
        out_e2 = self.scconv2(out_e2)+ out_e2
        #print(out_e2.shape)torch.Size([50, 32, 5, 5])

        out_e2 = self.eca2(out_e2) 

        out_e2 = self.dropout(out_e2)
        #print(out_e2.shape)torch.Size([50, 32, 5, 5])
        skip2 = out_e2
        
        
        out1 = self.conv3(out_e2)
        #print(out1.shape)torch.Size([50, 32, 5, 5])
        out1 = self.relu(out1)
        #print(out1.shape)torch.Size([50, 32, 5, 5])
        out1 = self.pool3(out1)
        #print(out1.shape)torch.Size([50, 32, 2, 2]) 
        kd_out1 = out1.view(b, -1)
        #print(kd_out1.shape)torch.Size([50, 128])

        out1_ = self.conv3_(out_e2)
        #print(out1_.shape)torch.Size([50, 32, 2, 2])
        out1_ = self.relu(out1_)
        #print(out1_.shape)torch.Size([50, 32, 2, 2])
        kd_out1_ = out1_.view(b, -1)
        #print(kd_out1_.shape)torch.Size([50, 128])
        
        out_e3 = (out1_+out1)/2 
        #print(out_e3.shape)torch.Size([50, 32, 2, 2])
        out_e3 = self.scconv3(out_e3)+ out_e3
        #print(out_e3.shape)torch.Size([50, 32, 2, 2])
        out_e3 = self.eca3(out_e3) 
        out_e3 = self.dropout(out_e3)
        
        kd_out_e3 = out_e3.view(b, -1)
        #print(kd_out_e3.shape)torch.Size([50, 128])
        
        
        #print(out_e3.shape)torch.Size([50, 32, 2, 2])
        skip3 = out_e3
        
        out1 = self.bn(out_e3)
        #print(out1.shape)torch.Size([50, 32, 2, 2])
        skip4 = out1

        # classifier
        out2 = skip4.view(b, -1)
        #print(out2.shape)torch.Size([50, 128])
        out2 = self.linear1(out2)
        #print(out2.shape)torch.Size([50, 64])
        out2 = self.relu(out2)
        #print(out2.shape)torch.Size([50, 64])
        out2 = self.drop(out2)
        #print(out2.shape)torch.Size([50, 64])
        out2 = self.linear2(out2)
        #print(out2.shape)torch.Size([50, 1])
        out_class = self.sigmoid(out2)
        #print(out_class.shape)torch.Size([50, 1])
        #print(kd_out_e3.shape, kd_out1_.shape, kd_out1.shape)torch.Size([50, 128]) torch.Size([50, 128]) torch.Size([50, 128]) 

        return out_class, [kd_out_e3, kd_out1_, kd_out1]#teacher,student1,student2
    

def KL(outputs, targets):
    log_softmax_outputs = F.log_softmax(outputs/4.0, dim=1)
    softmax_targets = F.softmax(targets/4.0, dim=1)
    return -(log_softmax_outputs * softmax_targets).sum(dim=1).mean()