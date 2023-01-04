#!/usr/bin/env python3

import xpore_readnames
import mbdmbd
import pandas as pd
import numpy as np

######## args for everything #########
xpath="/work/liam/nfruns/ADARknockdowncomparison/onB2/xpipor-pipe-pore/output/comparisons/mouse-antiadar-scramble/majority_direction_kmer_diffmod.table"
adir="/work/liam/nfruns/ADARknockdowncomparison/onB2/xpipor-pipe-pore/output/mouse_antiadar"
bdir="/work/liam/nfruns/ADARknockdowncomparison/onB2/xpipor-pipe-pore/output/mouse_scramble"
# compile args
arg_xporern = ["xpore_readnames.py","--xpore_table",xpath,"--sampleA_dir",adir,"--sampleB_dir",bdir]

### unmodified sample
bamdirA="/storage/common/rebasecall/p48_anti-adar_guppy6.4"
fdirA="/data/nanopore/data/p48-51_20210210_Adar-Abeta-Aza_R1_SQK_RNA002_HT22/P48_Anti-adar/20210210_2130_1-A1-D1_PAE86248_66c1fc9e/"
## modified sample
bamdirB="/storage/common/rebasecall/p50_Scramble_Si_guppy6.4"
fdirB="/data/nanopore/data/p48-51_20210210_Adar-Abeta-Aza_R1_SQK_RNA002_HT22/P50_Scramble-Si/20210210_2130_1-E9-H9_PAE86613_180257c4/"
# compile args
arg_kmertab_unmod = ['mbdmbd.py','--bamdir',bamdirA,'--fdir',fdirA]
arg_kmertab_mod = ['mbdmbd.py','--bamdir',bamdirB,'--fdir',fdirB]

########### run #############
xporereadnames = xpore_readnames.get_xpore_readnames(arg_xporern)
kmer_table_mod = mbdmbd.get_kmer_table_coord(xporereadnames, arg_kmertab_mod)
kmer_table_mod['call'] = 1
kmer_table_unmod = mbdmbd.get_kmer_table_coord(xporereadnames, arg_kmertab_unmod)
kmer_table_unmod['call'] = 0

# now that we have kmer, trying to do some torch
signal_table = pd.concat(
    [
    kmer_table_mod[['signal','call','kmer']],
    kmer_table_unmod[['signal','call','kmer']]
    ]
).reset_index()
del signal_table['index']

            #
          ## #
         ###
        ####
        #####
         #####
          ##
#######################
## THE TORCH ALIGHTS ##
#######################
# https://pytorch.org/tutorials/beginner/basics/quickstart_tutorial.html
import torch
from torch import nn
from torch.utils.data import DataLoader

# https://stackoverflow.com/questions/44429199/how-to-load-a-list-of-numpy-arrays-to-pytorch-dataset-loader
import torch
import numpy as np
from torch.utils.data import TensorDataset, DataLoader
######### data in
KMER = 'TTAAA'
data = signal_table[signal_table['kmer']==KMER]

######### split into train and test
p = 0.9
indices = np.random.choice(a=[True,False], size=len(data), p=[p,1-p])
train_data = data[indices]
test_data = data[[not i for i in indices]]
## torch dataloader shit
train_x = np.array(train_data['signal'].to_list()) # a list of numpy arrays
train_y = train_data['call'].to_list() # another list of numpy arrays (targets)
train_tensor_x = torch.Tensor(train_x) # transform to torch tensor
train_tensor_y = torch.Tensor(train_y)
train_dataset = TensorDataset(train_tensor_x,train_tensor_y) # create your datset
train_dataloader = DataLoader(train_dataset) # create your dataloader
test_x = np.array(test_data['signal'].to_list()) # a list of numpy arrays
test_y = test_data['call'].to_list() # another list of numpy arrays (targets)
test_tensor_x = torch.Tensor(test_x) # transform to torch tensor
test_tensor_y = torch.Tensor(test_y)
test_dataset = TensorDataset(test_tensor_x,test_tensor_y) # create your datset
test_dataloader = DataLoader(test_dataset) # create your dataloader

# Get cpu or gpu device for training.
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using {device} device")

# Define model
class NeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(50, 25),
            nn.ReLU(),
            nn.Linear(25, 12),
            nn.ReLU(),
            nn.Linear(12, 2)
        )
    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

model = NeuralNetwork().to(device)
print(model)

# Define loss function
loss_fn = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)

# Define training loop
def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (X, y) in enumerate(dataloader):
        X, y = X.to(device), y.to(device)

        # Compute prediction error
        pred = model(X)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch % 5 == 0:
            loss, current = loss.item(), batch * len(X)
            print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

def test(dataloader, model, loss_fn):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    model.eval()
    test_loss, correct = 0, 0
    with torch.no_grad():
        for X, y in dataloader:
            X, y = X.to(device), y.to(device)
            pred = model(X)
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) == y).type(torch.float).sum().item()
    test_loss /= num_batches
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")


epochs = 5
for t in range(epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(train_dataloader, model, loss_fn, optimizer)
    test(test_dataloader, model, loss_fn)
print("Done!")

