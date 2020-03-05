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

def main(indir,odir,annotate_16s,):
    a = pd.read_csv(annotate_16s, sep='\t')
    derep_a = a.drop_duplicates('sequence ID')

    g2name = defaultdict(list)
    genome2withMultiple16S = []
    for _, row in derep_a.iterrows():
        id = row['sequence ID']
        gname = id.split('_')[0]
        taxon_name = row['taxon_name']
        strain_name = row['strain_name']
        if gname in g2name:
            print(f'Wrong {gname} it detect multiple 16S. please manually remove one of them. it will not change the output of prokka.')
            genome2withMultiple16S.append(gname)
            # g2name.pop(gname)
        else:
            g2name[gname] = (taxon_name.replace(' ', '_'),
                             str(strain_name).replace(' ', '_'))
    for _ in genome2withMultiple16S:
        g2name.pop(_)

    if not exists(odir):
        os.makedirs(odir)
    for g in tqdm(glob(join(indir, '*', '*.gbk'))):
        records = SeqIO.parse(g, format='genbank')
        gname = basename(dirname(g))
        new_proteins = []
        annoate_name = '_'.join(g2name.get(gname, []))
        if not annoate_name:
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
        with open(join(odir, f"{gname}.faa"), 'w') as f1:
            SeqIO.write(new_proteins, f1, format='fasta-2line')

if __name__ == '__main__':
    # indir = './pipelines_o/prokka_o'
    # mode = 'prokka'
    # annotate_16s = './matchResults.csv'
    # odir = './new_proteins'
    import sys
    assert len(sys.argv) == 4
    # modify_prokka_o.py ./pipelines_o/prokka_o ./new_proteins ./matchResults.csv
    def process_path(path):
        if not '/' in path:
            path = './' + path
        return abspath(path)
    indir = process_path(sys.argv[1])
    odir = process_path(sys.argv[2])
    annotate_16s = process_path(sys.argv[3])
    main(indir,odir,annotate_16s)