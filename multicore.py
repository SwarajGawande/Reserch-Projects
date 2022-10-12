# Import package
import numpy as np
import random
from numpy.core.numeric import ones_like
from torch_geometric.data import Data, DataLoader
import torch
from tqdm import tqdm
from tqdm import trange 
from multiprocessing import Pool
#from  utils.debruijn import DeBruijnGraph

DATALIST=[]
# Define tf and datapaths
TF = 'ctcf'
processes=5
PATH_ = "/final_training_short.csv"
# Using readlines()
file1 = open('training_inter.csv', 'r')
Lines = file1.readlines()
t=[]
for i in range(processes):
    if i==processes-1:
        t.append(tqdm(Lines[i*len(Lines)//processes:]))
    else:
        t.append(tqdm(Lines[i*len(Lines)//processes:(i+1)*len(Lines)//processes]))
def run(arg):
    kmer, DATALIST = 4, []
    for l in range(arg['from'],arg['to']):
        t[arg['tqdm']].update()
        #print(line)
        line=Lines[l]    
        onehot_x = []
        # get DNA sequence from the dummy file 
        s = line.strip().split(',')[0]
        d=DeBruijnGraph(s,kmer)    
        for node in d.x.flatten():
            one_hot_ = d.one_hot_encode(node).flatten()
            onehot_x.append(one_hot_.tolist())

        # Arrays to pytorch tensors
        onehot_x_tensor = torch.tensor(np.array(onehot_x), dtype=torch.float)
        onehot_edge_index_tensor = torch.tensor(d.edge_index, dtype=torch.long)
    
        # get the label from the data
        if line.strip().split(',')[1]=='U':
            y_tensor = torch.tensor(0)
        else:
            y_tensor = torch.tensor(1)

        # Add tensors to torch_geometric data object
        data = Data(x=onehot_x_tensor, edge_index=onehot_edge_index_tensor, y=y_tensor)

        DATALIST.append(data)
        l=l+1
    return DATALIST

#def merge(d):
#    for i in range(1,len(d)):
#        for j in range(len(d[i])):
#            d[0].append(d[i][j])
#        d[i].clear()
#    return d[0]

if __name__=='__main__':
    arg=[]
    for i in range(processes):
        if i==processes-1:
            arg.append({'from':i*len(Lines)//processes,'to':len(Lines),'tqdm':i})
        else:
            arg.append({'from':i*len(Lines)//processes,'to':(i+1)*len(Lines)//processes,'tqdm':i})
    pool = Pool(processes) 
    d=list(pool.map(run, arg))
    pool.close()
    pool.join()
    #pool1=Pool(2)
    #d=list(pool1.map(merge,[d[:len(d)//2],d[len(d)//2:]]))
    for i in trange(1,len(d),desc='1st loop'):
        for j in trange(len(d[i]),desc='2nd loop'):
            d[0].append(d[i][j])
        d[i].clear()
    DATALIST=d[0]
