#################################################################################
#### For analysis some strain info, one should collect all reported data of this genus/species. This script is a guideline of this action.
####
####
#################################################################################
from Bio import SeqIO

all_pth = "/home/liaoth/data_bank/ncbi_Acinetobacter/raw_fna"

# 0. download required data

from subprocess import check_call

import pandas as pd
import os
from tqdm import tqdm
from glob import glob
import gzip
# df could download from https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/
# if download a whole genus, you could just downlaod the genus.
# Such as Acinetobacter. Although this catalogue just contains some unclassified species but belong to this genus(300+),but it will fetch all genome data from this genus.
df = pd.read_csv('/home/liaoth/temp_/Acinetobacter all genomes.txt', sep='\t', index_col=0)
odir = '/home/liaoth/data_bank/ncbi_Acinetobacter'

for idx, val in tqdm(df.iterrows(), total=df.shape[0]):
    strain = val['Strain']
    refseq_link = val["RefSeq FTP"]
    genbank_link = val["GenBank FTP"]
    cmd = "wget --quiet -N -w 1 --tries=10 -P  {odir} {path}"
    if refseq_link != '-':
        assembly_name = os.path.basename(refseq_link)
        path = os.path.join(refseq_link, "%s_genomic.fna.gz" % assembly_name)
    elif genbank_link != '-':
        assembly_name = os.path.basename(genbank_link)
        path = os.path.join(genbank_link, "%s_genomic.fna.gz" % assembly_name)
    else:
        continue
    new_cmd = cmd.format(odir=odir, path=path)
    check_call(new_cmd, shell=True, executable='/usr/bin/zsh')
# unzipped it
# 1. cut with 500bp and contigs number
df = pd.read_csv('/home/liaoth/temp_/Acinetobacter all genomes.txt', sep='\t', )
summary_info = pd.DataFrame()
for fna in tqdm(glob(os.path.join(all_pth, "*.fna"))):
    gid = '_'.join(os.path.basename(fna).split('_', maxsplit=2)[:2])
    genome = list(SeqIO.parse(gzip.open(fna,'rt'), format='fasta'))
    num_contig = len(genome)
    lens = [len(_) for _ in genome]
    contig_lt500 = len([_ for _ in lens if _ > 500])
    contig_lt1000 = len([_ for _ in lens if _ > 1000])

    summary_info = summary_info.append(pd.DataFrame([[num_contig, max(lens), sum(lens),
                                                      contig_lt500, contig_lt1000]],
                                                    index=[gid]))
    remained_contigs = [_ for _ in genome if len(_) > 500]
    with open(os.path.join(os.path.dirname(fna).replace('raw_fna',"preprocessed_fna"),
                           "%s_genomic.fna" % gid),'w') as f1:
        SeqIO.write(remained_contigs,f1,format='fasta')

basic_cols = ["num of contigs",
                        "length of largest contig",
                        "total length",
                        "number of contigs larger than 500",
                        "number of contigs larger than 1000", ]
summary_info.columns = basic_cols
fetch_cols = ['#Organism/Name', 'Strain', 'CladeID',
                    'BioSample', 'BioProject', 'GC%',
                    'Replicons', 'WGS', "Level"]
summary_info = summary_info.reindex(columns = basic_cols+fetch_cols)
for gid in tqdm(summary_info.index):
    # cannot use assembly id to match, maybe space after id.
    cache = df.loc[df.loc[:, "RefSeq FTP"].str.contains(gid),
                   fetch_cols]
    if cache.shape[0] == 0:
        cache = df.loc[df.loc[:, "GenBank FTP"].str.contains(gid),
                       fetch_cols]
    summary_info.loc[gid,fetch_cols] = [str(_).strip() if not pd.isna(_) else _ for _ in cache.values.ravel()]

summary_info.to_csv(os.path.join(os.path.dirname(all_pth),'summary_info.csv'),index=True)

# 2. annotated with prokka
check_call("python3 /home/liaoth/tools/genome_pipelines/toolkit/batch_prokka.py -i /home/liaoth/data_bank/ncbi_Acinetobacter/preprocessed_fna -o /home/liaoth/data_bank/ncbi_Acinetobacter/prokka_o")
