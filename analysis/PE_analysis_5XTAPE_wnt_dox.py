#!/usr/bin/env python 
# Author: Will Chen
import os,sys,csv
import re
from Bio import pairwise2
import pickle as pkl
import numpy as np
import pandas as pd
from umi_tools._dedup_umi import edit_distance

def min_edit(array,seq):
    m = 20
    ref=''
    for s in array.keys():
        dis = edit_distance(s.encode('UTF-8'),seq.encode('UTF-8'))
        if dis <= m:
            m = dis
            ref = s
    return [ref,m]

def merge_edit_1(counts):
    keep,collapse = [],counts[0:]
    while len(collapse)>0:
        ref = list(collapse[0]) # need to convert to list here otherwise will keep updating the value to new_cluster
        temp = []
        for s in collapse[1:]:
            s1 = ref[0]
            s2 = s[0]
            if abs(len(s1)-len(s2))<=1 and edit_distance(s1.encode('utf-8'),s2.encode('utf-8'))<2:
                ref[-1] += s[-1]
            else:
                temp.append(s)
        keep.append(ref)
        collapse = temp[0:]
    return np.array(keep,dtype=object)


sample  = sys.argv[1] 
workdir = sys.argv[2]
os.chdir(workdir)
TAPE = 'GGATGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCGCTTTAAGGCCGGTCCTAGCAA'
table = {}


regex = r".ACACC([ATCG]{8})(GGAT.*)"
pattern = re.compile(regex)


with open(workdir + sample + '.txt') as f:
    for line in f:
        rname,_,read = line.rstrip().split('\t')
        if len(read)>55:
            s = pattern.search(read)
            if s:
                table[rname] = [None]*10
                umi,read = s.groups()
                table[rname][0] = read
                table[rname][1] = TAPE
                table[rname][2] = umi
f.close()


temp_matrix = np.array([s for s in list(table.values())])
TAPE_6x = 'GGATGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCGCTTTAAGGCCGGTCCTAGCAA'
TAPE_4x = 'GGATGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCACGTGATGGTGAGCGCTTTAAGGCCGGTCCTAGCAA'
for i in range(len(temp_matrix)):
    s1,s2 = temp_matrix[i,0],temp_matrix[i,1]
    ali = pairwise2.align.globalms(s1, s2,5,-4,-10,-0.5,penalize_end_gaps=False)
    seq,ref,score = ali[-1][0][0:-20],ali[-1][1][0:-20],ali[-1][-3] # do not change start sequence. see below
    if '-----' in seq:
        s1,s2 = temp_matrix[i,0], TAPE_4x
        ali = pairwise2.align.globalms(s1, s2,5,-4,-10,-0.5,penalize_end_gaps=False)
        seq,ref,score = ali[-1][0][0:-20],ali[-1][1][0:-20],ali[-1][-3]
    elif '----------' in ref:
        s1,s2 = temp_matrix[i,0], TAPE_6x
        ali = pairwise2.align.globalms(s1, s2,5,-4,-10,-0.5,penalize_end_gaps=False)
        seq,ref,score = ali[-1][0][0:-20],ali[-1][1][0:-20],ali[-1][-3]
    temp_matrix[i,0:2] = seq,ref
    temp_matrix[i,-1]  = score
## Note: sometimes we will observe a tape missing from PCR: see 4 tapes instead of 5. alignment would be messy 

regex = r'CG------TGA'
p1 = re.compile(regex)
wl = {'GTT':'Dox','ACA':'Wnt'}
counts={i:{'GTT':0,'ACA':0} for i in range(1,6)}
combi = {a+'+'+b+' ' + str(i)+'to'+str(i+1) :0 for a in ['GTT','ACA'] for b in ['GTT','ACA'] for i in range(1,5)}
combi['GTT+ACA' + ' ' + str(1)+'to'+str(3)]=0 # add 1+3 position
combi['ACA+GTT' + ' ' + str(1)+'to'+str(3)]=0 # add 1+3 position 
## Some other error mode where the editing is not sequential, such as only first and last unit was edit. ignore it
insert_loc = {17:1,\
              31:2,37:2,\
              45:3,51:3,57:3,\
              59:4,65:4,71:4,77:4,\
              73:5,79:5,85:5,91:5,97:5} # this numbers are highly related to the start of reference sequence. it should start with GGATGA right after UMI. do not change
for s in temp_matrix:
    index = [(m.start(0)+2, m.end(0)-3) for m in p1.finditer(s[1])]
    insert = [min_edit(wl,s[0][k[0]:k[0]+3])[0] for k in index]     
    index_loc = []
    for ind in index:
        try:
            index_loc.append(insert_loc[ind[0]])
        except KeyError:
            pass # only 0.03% doesn't align well.
    if index and s[-1]>325 and index_loc:
        for i in range(len(index_loc)):
            counts[index_loc[i]][insert[i]] += 1
            try:
                combi[insert[i]+'+'+insert[i+1]+' '+str(index_loc[i])+'to'+str(index_loc[i+1])] += 1
                combi[insert[i]+'+'+insert[i+2]+' '+str(index_loc[i])+'to'+str(index_loc[i+2])] += 1
            except IndexError:
                pass
            except KeyError:
                pass

final_count=pd.concat([pd.DataFrame([(j, i, counts[i][j]) for i in counts.keys() for j in counts[i].keys() ], columns=['barcode','Position','count']),\
                       pd.DataFrame([a.split(' ')+[int(b)] for a,b in combi.items()],columns=['barcode','Position','count'])],axis=0)
final_count['ratio'] = final_count['count']/len(temp_matrix)*100


fname = workdir + sample
final_count.to_csv(fname+'_bc_count.csv',sep='\t',index=False)