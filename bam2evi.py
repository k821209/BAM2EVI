import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
## ... ... ...
import matplotlib
matplotlib.style.use('ggplot') 
import math
import sys
sys.path.append('/ref/analysis/pipelines/')
import kang

file_forward_npy = '/ref/analysis/Cre/braker/braker.try5_mario/intron3000.merge.sorted.bam.paired_end_continuity_forward.np.npy'
file_reverse_npy = '/ref/analysis/Cre/braker/braker.try5_mario/intron3000.merge.sorted.bam.paired_end_continuity_reverse.np.npy'
file_pk = '/ref/analysis/pipelines/pandas_df/Creinhardtii_281_v5.5.gene.gff3.intron.gff3.pandas.df.pk'
file_cm = '/ref/analysis/Cre/braker/braker.try5_mario/chromosome.map.txt'

succ_list = [x.strip() for x in open('mayra_succ_list.txt').readlines()]
nsucc_list = [x.strip() for x in open('mayra_no_succ_list.txt').readlines()]
polyA_list = [x.strip() for x in open('polyA.one.txt').readlines()]

df_cm = pd.read_csv(file_cm,sep='\t',header=None)
df_cm_ix_chr2n = df_cm.set_index(0)
df_cm_ix_n2chr = df_cm.set_index(1)
array_contiguity_p = np.load(file_forward_npy)
array_contiguity_m = np.load(file_reverse_npy)
df_gff_cre = pd.read_pickle(file_pk)

def get_block(array,depth_cut=0):
    lim_len_block = 100
    #depth_cut     = 0 # ... 10. ... .. .. ..
    block_list = []
    #print(len(np.shape(array)))
    if len(np.shape(array)) == 1:
        rows = 1
        block = []
        for n,j in enumerate(array):
            if j > depth_cut:
                block.append(n)
            else:
                if len(block) > lim_len_block:
                    block_list.append([block[0],block[-1]])
                    block = []
                else:
                    block = []
        if block != []:
            block_list.append([block[0],block[-1]])
    else: 
        rows, columns = np.shape(array)       
        for i in range(rows):
            earray = array[i]
            block = []
            for n,j in enumerate(earray):
                if j > depth_cut:
                    block.append(n)
                else:
                    if len(block) > lim_len_block:
                        block_list.append([i,block[0],block[-1]])
                        block = []
                    else:
                        block = []
    return block_list



dic = {'mRNA'       : [],
       'length'     : [],
       'contiguity.plus' : [],
       'ratio.zero.plus' : [],
       'total.depth.plus': [],
       'ratio.depth.plus': [],
       'blocks.plus'     : [],
       'coverage.plus'   : [],
       'start5.plus'     : [],
       'end5.plus'       : [],
       '#blocks.plus'    : [],
       'contiguity.minus' : [],
       'ratio.zero.minus' : [],
       'total.depth.minus': [],
       'ratio.depth.minus': [],
       'blocks.minus'     : [],
       'coverage.minus'   : [],
       'start5.minus'     : [],
       'end5.minus'       : [],
       '#blocks.minus'    : [],
       'strand'           : [],
       'PCR'              : [],
       'CDS.depth.plus' : [],
       'INTRON.depth.plus' : [],
       'CDS.depth.minus'    : [],
       'INTRON.depth.minus' : [],
       'CDS.cov.plus' : [],
       'INTRON.cov.plus' : [],
       'CDS.cov.minus'    : [],
       'INTRON.cov.minus' : [],
       'CDS.num'          : [],
       'INRTON.num'       : [],
       'polyAone'         : [],
      }


genelist = set([x for x,y in df_gff_cre.index])
for genename in genelist:
    try:
        if math.isnan(float(genename)):
            continue
    except ValueError:
        pass
    
    df      = df_gff_cre.loc[genename]
    
    # CDS support calc. 
    mask    = (df[2]=='CDS')
    df_mRNA = df[mask].loc['1']
    try:
        chromosome = df_mRNA[0].values[0]
    except AttributeError:
        chromosome = df_mRNA[0]
    echr          = df_cm_ix_chr2n.loc[chromosome][1]
    array         = df[mask][[3,4]].values
    totdepth_p    = 0
    coveredbase_p = 0
    totdepth_m    = 0
    coveredbase_m = 0
    length        = 0
    dic['CDS.num'].append( array.shape[0]) 
    for left, right in array:
        length        += right - left + 1 - 1
        contiguity_p   = array_contiguity_p[echr][left-1:right-1] # continuity value require minus 1 from right pos
        contiguity_m   = array_contiguity_m[echr][left-1:right-1]
        totdepth_p    += np.sum(contiguity_p)
        totdepth_m    += np.sum(contiguity_m)
        coveredbase_p += len((contiguity_p > 0).nonzero()[0])
        coveredbase_m += len((contiguity_m > 0).nonzero()[0])
    dic['CDS.depth.plus'].append(float(totdepth_p)/float(length))
    dic['CDS.depth.minus'].append(float(totdepth_m)/float(length))
    dic['CDS.cov.plus'].append(float(coveredbase_p)/float(length))
    dic['CDS.cov.minus'].append(float(coveredbase_m)/float(length))
    
    # intron support calc. 
    mask    = (df[2]=='intron')
    bintron = 1
    try:
        df_mRNA = df[mask].loc['1']
    except KeyError:
        dic['INRTON.num'].append(0)
        dic['INTRON.depth.plus'].append(0)
        dic['INTRON.depth.minus'].append(0)
        dic['INTRON.cov.plus'].append(0)
        dic['INTRON.cov.minus'].append(0)
        bintron = 0
        pass 
    if bintron == 1:
        try:
            chromosome = df_mRNA[0].values[0]
        except AttributeError:
            chromosome = df_mRNA[0]
        echr          = df_cm_ix_chr2n.loc[chromosome][1]
        array         = df[mask][[3,4]].values
        totdepth_p    = 0
        coveredbase_p = 0
        totdepth_m    = 0
        coveredbase_m = 0
        length        = 0 
        dic['INRTON.num'].append(array.shape[0])
        if array.shape[0] > 0:
            for left, right in array:
                length        += right - left + 1 - 1
                contiguity_p   = array_contiguity_p[echr][left-1:right-1] # continuity value require minus 1 from right pos
                contiguity_m   = array_contiguity_m[echr][left-1:right-1]
                totdepth_p    += np.sum(contiguity_p)
                totdepth_m    += np.sum(contiguity_m)
                coveredbase_p += len((contiguity_p > 0).nonzero()[0])
                coveredbase_m += len((contiguity_m > 0).nonzero()[0])
                print genename,left,right,totdepth_p,totdepth_m,coveredbase_p,coveredbase_m
            dic['INTRON.depth.plus'].append(float(totdepth_p)/float(length))
            dic['INTRON.depth.minus'].append(float(totdepth_m)/float(length))
            dic['INTRON.cov.plus'].append(float(coveredbase_p)/float(length))
            dic['INTRON.cov.minus'].append(float(coveredbase_m)/float(length))
        else:
            dic['INTRON.depth.plus'].append(None)
            dic['INTRON.depth.minus'].append(None)
            dic['INTRON.cov.plus'].append(None)
            dic['INTRON.cov.minus'].append(None)
    
    
    # Gene regions support calc.  
    mask    = (df[2]=='CDS')
    df_mRNA = df[mask].loc['1']
    try:
        chromosome = df_mRNA[0].values[0]
    except AttributeError:
        chromosome = df_mRNA[0]
    left       = np.min(df_mRNA[[3,4]].values) #int(df_mRNA[3])
    right      = np.max(df_mRNA[[3,4]].values) #int(df_mRNA[4])
    try:
        strand     = df_mRNA[6].values[0]
    except :
        strand     = df_mRNA[6]
    #print chromosome,left,right
    length     = right - left + 1 
    echr       = df_cm_ix_chr2n.loc[chromosome][1]
    #print echr,chromosome
    #print left,right
    contiguity_p = array_contiguity_p[echr][left-1:right-1] # continuity value require minus 1 from right pos
    contiguity_m = array_contiguity_m[echr][left-1:right-1]
    #if np.percentile(contiguity,1) - min(contiguity) > 100:
    #    contiguity = contiguity - np.percentile(contiguity,1)
    zeropart_p   = list(contiguity_p).count(0)
    zeropart_m   = list(contiguity_m).count(0)
    #if min(contiguity) == 0:
    if 1:
        blocks_p = get_block(contiguity_p)
        blocks_m = get_block(contiguity_m)
        if strand == '+':
            dic['strand'].append(1)
        else:
            dic['strand'].append(0)
        dic['mRNA'].append(genename)
        dic['contiguity.plus'].append(contiguity_p)
        dic['contiguity.minus'].append(contiguity_m)
        dic['ratio.zero.plus'].append(float(zeropart_p)/float(length))
        dic['ratio.zero.minus'].append(float(zeropart_m)/float(length))
        dic['length'].append(length)
        dic['total.depth.plus'].append(sum(contiguity_p))
        dic['total.depth.minus'].append(sum(contiguity_m))
        dic['ratio.depth.plus'].append(float(sum(contiguity_p))/float(length))
        dic['ratio.depth.minus'].append(float(sum(contiguity_m))/float(length))
        dic['blocks.plus'].append(blocks_p)
        dic['blocks.minus'].append(blocks_m)
        dic['#blocks.plus'].append(len(blocks_p))
        dic['#blocks.minus'].append(len(blocks_m))
        dic['coverage.plus'].append(1-float(zeropart_p)/float(length))
        dic['coverage.minus'].append(1-float(zeropart_m)/float(length))
        if sum(contiguity_p[0:5]) == 0:
            dic['start5.plus'].append(0)
        else:
            dic['start5.plus'].append(1)
        if sum(contiguity_m[0:5]) == 0:
            dic['start5.minus'].append(0)
        else:
            dic['start5.minus'].append(1)
        if sum(contiguity_p[-5:]) == 0:
            dic['end5.plus'].append(0)
        else:
            dic['end5.plus'].append(1)
        if sum(contiguity_m[-5:]) == 0:
            dic['end5.minus'].append(0)
        else:
            dic['end5.minus'].append(1)
        
        if genename.replace('.v5.5','') in set(succ_list):
            dic['PCR'].append(1)
        elif genename.replace('.v5.5','') in set(nsucc_list):
            dic['PCR'].append(0)
        else:
            dic['PCR'].append(None)
   
        if genename in set(polyA_list):
            dic['polyAone'].append(1)
        else:
            dic['polyAone'].append(0)
for key in dic:
    print(len(dic[key]))            
df_cont = pd.DataFrame(dic)
df_cont_ix = df_cont.set_index('mRNA')
df_cont.to_pickle('continuity_strand.df')
df_cont.to_csv('continuity_strand.txt',sep='\t')
