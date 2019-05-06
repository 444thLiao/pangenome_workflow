import os

import gffutils
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

tmp_dir = '/tmp'
def clean_db(tab,prokka_dir):
    data = pd.read_csv(tab, sep=',', index_col=0, low_memory=False)
    for name in data.loc[0, :][13:].index:
        gff_fn = os.path.join(prokka_dir, name, '%s.gff' % name)
        dbfn = os.path.join(tmp_dir, os.path.basename(gff_fn).split('.')[0])
        os.remove(dbfn)

def get_gff(gff_fn):
    dbfn = os.path.join(tmp_dir, os.path.basename(gff_fn).rsplit('.',1)[0])
    if os.path.isfile(dbfn):
        fn = gffutils.FeatureDB(dbfn, keep_order=True)
    else:
        fn = gffutils.create_db(gff_fn, dbfn=dbfn,merge_strategy='merge')
    return fn

def get_gene_with_regin(gff_fn,regions):
    gff_f = get_gff(gff_fn)
    all_genes = []
    for region in regions:
        for cds in gff_f.region(region=(region), completely_within=True):
            if "ID" in cds.attributes.keys():
                all_genes.append(cds["ID"][0])
            else:
                all_genes.append(cds["note"][0])
    return all_genes

if __name__ == '__main__':
    # todo: make it a api with output which describe the gene
    pass