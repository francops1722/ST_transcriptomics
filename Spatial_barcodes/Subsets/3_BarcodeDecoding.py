#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:01:45 2022

@author: fpomasot
"""
# DECODING Step 1. Move the code (assumed loaded above) and reads (from file) to the GPU.

from numpy import *
import pickle
import torch
import argparse
import re
import yaml
from RandomBar_Functions import encode, seqtomer, mertobin, mypopcount

parser = argparse.ArgumentParser()
parser.add_argument('--config',required=True, help='path to yaml config file')
yaml_file = parser.parse_args()

with open(yaml_file.config, "r") as fp:
        args = yaml.safe_load(fp)

class ApproximateLevenshtein :
    def __init__(s, M, N, Q, zsub, zins, zdel, zskew):
        torch.set_grad_enabled(False) # just in case not done elsewhere!
        s.M = M # length of seq1
        s.N = N # length of each seq2
        s.Q = Q # number of seq2s
        (s.zsub, s.zins, s.zdel, s.zskew) = (zsub, zins, zdel, zskew)
        s.tab = torch.zeros(N+1,Q, device=device)

    def __call__(s,seq1,seq2) :
        assert (len(seq1) == s.M) and (seq2.shape[1] == s.N) and (seq2.shape[0] == s.Q)
        s.tab[:,:] = (s.zskew * torch.arange(s.N+1., device=device)).unsqueeze(1) # force broadcast
        for i in range(1,s.M+1) :
            diag = s.tab[:-1,:] + torch.where(seq1[i-1] == seq2.t(), 0., s.zsub) # diagonal move
            s.tab[0,:] += s.zskew
            s.tab[1:,:] += s.zdel # down move
            s.tab[1:,:] = torch.minimum(s.tab[1:,:], diag) # or diag if better
            s.tab[1:,:] = torch.minimum(s.tab[1:,:], s.tab[:-1,:] + s.zins) # right move
            s.tab[1:,:] = torch.minimum(s.tab[1:,:], s.tab[:-1,:] + s.zins) # repeat (>= 0 times) as you can afford
           # N.B.: M >= N gives better approx than N > M, so change arg order accordingly
        return s.tab[s.N,:]

def allcoses(mer, tcosvecs) : # correlate a mer against all the cosine templates
    mmer = torch.LongTensor(mer).to(device)
    ncos = tcosvecs.shape[0]
    cosvec = torch.zeros(ncos, 64, dtype=torch.float, device=device)
    for k in range(ncos) :
        source = tcoses[k, torch.arange(len(mmer), dtype=torch.long, device=device)]
        cosvec[k,:].index_add_(0,mmer,source)
    return torch.sum(torch.unsqueeze(cosvec,dim=1)*tcosvecs,dim=2)

def prank(arr, descending=False) : # returns rank of each element in torch array
    argsrt = torch.argsort(arr, descending=descending)
    rank = torch.zeros(arr.shape, dtype=torch.float, device=device)
    rank[argsrt] = torch.arange(len(argsrt),dtype=torch.float,device=device)
    return rank


picklefilename = args["picklefilename"]#the barcodes pickle file used to deconvolute the reads
with open(picklefilename,'rb') as IN :
    pickledict = pickle.load(IN)
(N, M, allseqs, alltrimers, allbitmaps, coses, cosvecs) = \
    [pickledict[x] for x in ('N', 'M', 'allseqs', 'alltrimers', 'allbitmaps', 'coses', 'cosvecs')]
print(f"loaded code with {N} codewords of length {M} from {picklefilename}.")


readsfilename = args["readsfilename"] # user to set

if torch.cuda.is_available() :
    device = torch.device("cuda")
    cudaname = torch.cuda.get_device_name()
else :
    raise RuntimeError("Required GPU not found! Exiting.")
print(f"Using {device} device {cudaname}.")

with open(readsfilename, 'r') as IN:
    stringreads=IN.readlines()
reads = []
for read in stringreads :
    reads.append(encode(read[:-1])) # lose the \n
reads = array(reads)

torch.set_grad_enabled(False)
tallseqs = torch.tensor(array(allseqs), device=device)
talltrimers = torch.tensor(array(alltrimers), device=device)
tallbitmaps = torch.tensor(allbitmaps.astype(int64), dtype=torch.int64, device=device) # 
tcoses = torch.tensor(coses, dtype=torch.float, device=device)
tcosvecs = cosvecs.to(device)
print(f"found {reads.shape[0]} reads of length {reads.shape[1]}")

# DECODING Step 2 (main step). Primary and secondary triage, followed by Levenshtein

Qq = 1000000 # reads.shape[0] or smaller number, how many reads to do
Ntriage = 10000 # user set to number passed from triage to Levenshtein
Nthresh = args["Nthresh"] # user set to Levenshtein score greater than which is called an erasure

torch.set_grad_enabled(False)

mydist = ApproximateLevenshtein(M,M,Ntriage, 1.,1.,1.,1.)
#mydist = ParallelLevenshtein(M,M,Ntriage, 1.,1.,1.,1.)

Ncos = tcosvecs.shape[0]
dists = torch.zeros(Ncos+1, N, dtype=torch.float, device=device) # will hold distances for each read
allrank = torch.zeros(Ncos+1 ,N, dtype=torch.float, device=device)
best = torch.zeros(Qq, dtype=torch.long, device=device)

for j,seq in enumerate(reads[:Qq]) :
    # primary and secondary triage
    mer = seqtomer(seq)
    foo = int64(uint64(mertobin(mer))) # need to cast 64 bits to a type known to torch
    seqbin = torch.tensor(foo,dtype=torch.int64,device=device)
    xored = torch.bitwise_xor(seqbin,tallbitmaps)
    dists[0,:] = 64. - mypopcount(xored) # all Hamming distances
    cosvec = torch.zeros(Ncos, 64, dtype=torch.float, device=device)
    for k in range(Ncos) :
        cosvec[k,mer] =  tcoses[k, torch.arange(len(mer), dtype=torch.long, device=device)]
    dists[1:,:] = torch.sum(torch.unsqueeze(cosvec,dim=1)*tcosvecs,dim=2) # all cosine distances
    for k in range(Ncos+1) :
        allrank[k,:] = prank(dists[k,:], descending=True) # rank them all
    offset = 1.
    fm = torch.prod(offset+allrank,dim=0)
    fmargsort = torch.argsort(fm)
    # Levenshtein distance
    tseq1 = torch.tensor(seq,device=device)
    tseq2 = tallseqs[fmargsort[:Ntriage],:]
    ans = mydist(tseq1,tseq2)
    ia = torch.argmin(ans) # index into fmargsort of best
    ibest = fmargsort[ia] # index of best codeword in codes
    best[j] = (ibest if ans[ia] <= Nthresh else -1) # erasures returned as -1

best = best.cpu().numpy()

best_path = args["best_path"]

savetxt(best_path, best)

print(f"done with decoding with threshold: {Nthresh}")
