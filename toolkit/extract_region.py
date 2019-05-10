import os

import gffutils
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

def get_gff(gff_fn):
    tmp_dir = '/tmp'
    dbfn = os.path.join(tmp_dir, os.path.basename(gff_fn).split('.')[0])
    if os.path.isfile(dbfn):
        fn = gffutils.FeatureDB(dbfn, keep_order=True)
    else:
        fn = gffutils.create_db(gff_fn, dbfn=dbfn)

    # cds.start
    return fn


def seq_between_gene1_gene2(tab, prokka_dir, outdir, g1, g2):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    data = pd.read_csv(tab, sep=',', index_col=0, low_memory=False)
    if g1 not in data.index or g2 not in data.index:
        raise Exception
    row1 = data.loc[g1, :][13:]
    row2 = data.loc[g2, :][13:]

    for name in tqdm(row1.index):
        gff_file = os.path.join(prokka_dir, name, '%s.gff' % name)
        fasta_file = os.path.join(prokka_dir, name, '%s.fna' % name)
        gff_obj = get_gff(gff_file)
        cds1 = gff_obj[row1[name]]
        cds2 = gff_obj[row2[name]]
        starts = [cds1.start, cds2.start, cds1.end, cds2.end]
        ends = [cds1.start, cds2.start, cds1.end, cds2.end]
        # if cds1.strand != cds2.strand:
        #     import pdb;pdb.set_trace()
        if cds1.chrom != cds2.chrom:
            print('%s: unmatch chrome' % name)
            continue
        start, end = min(starts)-10, max(ends) + 10
        fa = SeqIO.parse(open(fasta_file, 'r'), format='fasta')
        fa = [_ for _ in fa if _.description == cds1.chrom][0]
        subseq = fa[start:end]
        subseq.id = subseq.name = subseq.description
        subseq.id = "%s: gene region (%s-%s) (%s:%s-%s)" % (name, g1, g2, cds1.chrom, start, end)
        with open(os.path.join(outdir, name + ".fa"), 'w') as f1:
            SeqIO.write(subseq, f1, format='fasta')

if __name__ == '__main__':
    g1 = 'fkpA_1'
    g2 = 'lldP'
    tab = "/home/liaoth/project/shenzhen_actineto/roary_o/gene_presence_absence.csv"
    prokka_dir = "/home/liaoth/project/shenzhen_actineto/prokka_o"
    outdir = "/home/liaoth/project/shenzhen_actineto/KL_extracted_reads"
    base_dir = "%s-%s_region" % (g1, g2)

    # main(g1,g2,tab,prokka_dir,odir)
    seq_between_gene1_gene2(tab, prokka_dir, os.path.join(outdir, base_dir), g1, g2)

    g1 = 'fkpA_2'
    base_dir = "%s-%s_region" % (g1, g2)
    seq_between_gene1_gene2(tab, prokka_dir, os.path.join(outdir, base_dir), g1, g2)

    # main(args)