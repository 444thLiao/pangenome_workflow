"""
merged sequencing stats together into a table

"""
import pandas as pd
from glob import glob
import os
from os.path import *
import sys
sys.path.insert(0,__file__.rsplit('/',3)[0])
from toolkit.parse_kraken2 import merge_kraken2
# python ~/software/pangenome_workflow/toolkit/extra_toolkit_luolab/get_summary.py ./input.tab ./
if __name__=='__main__':
    infile, odir = sys.argv[1:]
# odir = './'
# infile = './input.tab'
    indf = pd.read_csv(infile,sep='\t',index_col=0)

    for postQC in glob(f'{odir}/seqtk_result/*_reads.summary'):
        if os.path.getsize(postQC)==0:
            continue    
        d = pd.read_csv(postQC,index_col=0,sep='\t')
        name = postQC.split('/')[-1].split('_reads')[0]
        indf.loc[name,'post QC/Mbp'] = round(d.loc[name,'metricsReads_Yield']/10**6,2)
        
    for postQC in glob(f'{odir}/seqtk_result/*_assembly.summary'):
        if os.path.getsize(postQC)==0:
            continue
        d = pd.read_csv(postQC,index_col=0,sep='\t')
        name = postQC.split('/')[-1].split('_assembly')[0]
        indf.loc[name,'assembly N50/Mbp'] = round(d.loc[name,'metricsContigs_N50']/10**6,2)
        indf.loc[name,'assembly #contigs'] = d.loc[name,'metricsContigs_no']
        indf.loc[name,'assembly length/Mbp'] = round(d.loc[name,'metricsContigs_bp']/10**6,2)
        
    indf = indf.iloc[:,4:]
    indf = indf.reindex(sorted(indf.index,key=lambda x: (x[:2],int(x[2:]))))
    taxdf = merge_kraken2(glob(f'{odir}/kraken2_report/*'))
    indf = indf.join(taxdf)
    indf = indf.fillna('Failed')
    indf.to_csv(f'{odir}/sequencing_stats.tsv',sep='\t')

