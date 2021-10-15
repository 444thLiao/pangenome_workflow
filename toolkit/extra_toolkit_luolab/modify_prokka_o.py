"""
This is a script for modifying the output of prokka into a designed proteins sequence.
This should use `extract 16S.py` first and perform `16s annotation` first.
extract 16S.py is in the same directory
16s annotation is in the `script/EZbiocloudmatch.py`
"""
from Bio import SeqIO
from glob import glob
from os.path import *
from tqdm import tqdm
from collections import defaultdict
import pandas as pd
import os

# target_col1 = "16S_species"
# target_col2 = "gtdb genus"
# df = pd.read_excel('summary_sequencing_20200303.xlsx',index_col=0)
# g2name = df[target_col1].to_dict()
# for g,name in list(g2name.items()):
#     if name == 'NA' or pd.isna(name):
#         alt_name = df.loc[g,target_col2 ]
#         if pd.isna(alt_name):
#             alt_name = 'Unclassified '
#         else:
#             alt_name = alt_name+' sp. '
#         g2name.update({g:alt_name })


# def get_taxon_name(annotation_file):
#     a = pd.read_csv(annotation_file, sep='\t')
#     derep_a = a.drop_duplicates('sequence ID')

#     g2name = defaultdict(list)
#     genome2withMultiple16S = []
#     for _, row in derep_a.iterrows():
#         id = row['sequence ID']
#         gname = id.split('_')[0]
#         taxon_name = row['taxon_name']
#         strain_name = row['strain_name']
#         if gname in g2name:
#             print(f' {gname} might be wrongly assigned because of multiple 16S detected. please manually remove one of them. it will not change the output of prokka.')
#             genome2withMultiple16S.append(gname)
#             # g2name.pop(gname)
#         else:
#             g2name[gname] = (taxon_name.replace(' ', '_'),
#                              str(strain_name).replace(' ', '_'))
#     for _ in genome2withMultiple16S:
#         g2name.pop(_)
#     return g2name

def get_taxon_name(annotation_file):
    g2name = dict([row.split('\t') for row in  open(annotation_file).read().split('\n') if row])
    # e.g :  'Youngimonas_vesicularis_CC-AMW-E_SZCCSK-T6F11 ': ' Youngimonas_vesicularis_SZCCSK-T6F11',
    #g2name = {k.rpartition('_')[-1].strip():v.rpartition('_')[0].strip() for k,v in g2name.items()}
    g2name = {k.strip():v.strip().replace(' ','_') for k,v in g2name.items()}
    return g2name
    
def main(indir,odir,annotation_file,):
    g2name = get_taxon_name(annotation_file)

    if not exists(odir):
        os.makedirs(odir)
    for g in tqdm(glob(join(indir, '*', '*.gbk'))):
        records = SeqIO.parse(g, format='genbank')
        gname = basename(dirname(g))
        new_proteins = []
        annoate_name = g2name.get(gname)
        
        ofaa = join(odir, f"{gname}.faa")
        if exists(ofaa):
            tqdm.write(f"Bypass {gname} cause it exists")
            continue
        if  annoate_name is None:
            tqdm.write(f"{gname} doesn't have corresponding taxon assignment")
            continue
        for r in records:
            CDS_records = [_ for _ in r.features if _.type == 'CDS']
            for cds in CDS_records:
                _r = cds.translate(r)
                _r.id = _r.name = _r.description = ''
                q = cds.qualifiers
                locus_tag = q['locus_tag'][0]
                product = ' '.join(q['product'])
                start = cds.location.start.real
                end = cds.location.start.real
                pos = f"[{r.id}:c({start}..{end})]"
                _r.id = f"{annoate_name}_{gname}|{locus_tag} {pos} [{product}]"
                new_proteins.append(_r)
        with open(ofaa, 'w') as f1:
            SeqIO.write(new_proteins, f1, format='fasta-2line')

if __name__ == '__main__':
    # indir = './pipelines_o/prokka_o'
    # mode = 'prokka'
    # annotation_file = '/mnt/storage2/jjtao/jjtao_20200119/new_proteins/rename_20200119'
    # odir = './new_proteins'
    import sys
    assert len(sys.argv) == 4
    # modify_prokka_o.py ./pipelines_o/prokka_o ./new_proteins /mnt/storage2/jjtao/jjtao_20200119/new_proteins/rename_20200119
    def process_path(path):
        if not '/' in path:
            path = './' + path
        return abspath(path)
    indir = process_path(sys.argv[1])
    odir = process_path(sys.argv[2])
    annotate_16s = process_path(sys.argv[3])
    main(indir,odir,annotate_16s)