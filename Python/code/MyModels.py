import torch.nn.functional as F
import torch
import torch.nn as nn
import scipy.io as scio

from torch.autograd import Variable
from torch.nn.parameter import Parameter
import math
import numpy as np
import torch
from collections import OrderedDict

gpu = 1

class RNN(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, tsize, classnum=2):
        super(RNN,self).__init__()

        self.hidden_size = hidden_size
        self.input_size = input_size
        self.num_layers = num_layers
        self.classnum = classnum

        self.rnn = torch.nn.LSTM(
            input_size,
            hidden_size,
            num_layers,
            batch_first=True
        )

        self.fct = nn.Linear(tsize*hidden_size,hidden_size)
        self.fc1 = nn.Linear(hidden_size,classnum)

        self.bn1 = nn.BatchNorm1d(tsize*hidden_size)
        self.bn2 = nn.BatchNorm1d(hidden_size)
        self.softmax = torch.nn.Softmax(dim=1)

    def forward(self,x):
        x = x.cuda()
        bsize=x.size(0)
        fsize=x.size(1)
        h0 = torch.zeros(self.num_layers,bsize,self.hidden_size).cuda()
        c0 = torch.zeros(self.num_layers,bsize,self.hidden_size).cuda()

        out_list, (h_n, c_n) = self.rnn(x, (h0,c0))

        r_out = torch.zeros([bsize,fsize,self.hidden_size]).cuda()

        for i in range(len(out_list)):
            r_out[i,:,:] = out_list[i]

        r_out = r_out.reshape([bsize,fsize*self.hidden_size])
        r_out = self.bn1(r_out)
        out = F.relu(self.fct(r_out))
        out = self.bn2(out)
        out = F.relu(self.fc1(out))
        out = self.softmax(out)

        return out
