import os
import sys

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from BCBio import GFF
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from toolkit.get_gene_info import get_gff
from toolkit.utils import valid_path


def seq_between_gene1_gene2(roary_dir, prokka_dir, outdir, g1, g2):
    valid_path(outdir, check_odir=True)
    valid_path([roary_dir, prokka_dir], check_dir=True)

    input_tab = os.path.join(roary_dir, "gene_presence_absence.csv")
    data = pd.read_csv(input_tab, sep=',', index_col=0, low_memory=False)
    if (g1 not in data.index) or (g2 not in data.index):
        raise Exception("request gene %s or %s not found at %s" % (g1,
                                                                   g2,
                                                                   input_tab))

    row1 = data.loc[g1, :][13:]
    row2 = data.loc[g2, :][13:]

    for name in tqdm(row1.index):
        # prokka dir must list like this
        gff_file = os.path.join(prokka_dir,
                                name,
                                '%s.gff' % name)
        fasta_file = os.path.join(prokka_dir,
                                  name,
                                  '%s.fna' % name)

        gff_obj = get_gff(gff_file,mode="db")
        # get sqlite format gff obj
        # because it it faster
        cds1 = gff_obj[row1[name]]
        cds2 = gff_obj[row2[name]]
        starts = [cds1.start, cds2.start, cds1.end, cds2.end]
        ends = [cds1.start, cds2.start, cds1.end, cds2.end]
        record_dict = get_gff(gff_file,mode='bcbio')

        if cds1.chrom != cds2.chrom:
            print('%s: unmatch chrome' % name)
            print("It means that, sample %s got different pos of these genes" % name)
            continue

        start, end = min(starts) - 10, max(ends) + 10
        fa = SeqIO.parse(open(fasta_file, 'r'), format='fasta')
        fa = [_ for _ in fa if _.description == cds1.chrom][0]
        subseq = fa[start:end]
        contig_id = subseq.id
        gff_record = record_dict[contig_id]
        subseq.id = subseq.name = subseq.description = ''
        subseq.id = "%s:%s-%s" % (contig_id,
                                  start,
                                  end)
        sub_gff_record = gff_record[start:end]
        sub_gff_record.id = subseq.id
        with open(os.path.join(outdir, name + ".fa"), 'w') as f1:
            SeqIO.write(subseq, f1, format='fasta')
        with open(os.path.join(outdir, name + ".gff"), 'w') as f1:
            GFF.write([sub_gff_record], f1)

def main(args):
    gene1, gene2 = args.gene1, args.gene2
    roary_dir, prokka_dir, output_dir = (args.roary_dir,
                                         args.prokka_dir,
                                         args.output_dir)
    seq_between_gene1_gene2(roary_dir, prokka_dir, output_dir,
                            gene1, gene2)


if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("gene1",  type=str)
    parse.add_argument("gene2",  type=str)
    parse.add_argument("-r", "--roary_dir", help='')
    parse.add_argument("-p", "--prokka_dir", help='')
    parse.add_argument("-o", "--output_dir", help='')

    args = parse.parse_args()

    main(args)

    #python toolkit/extract_region.py fkpA_1 lldP -r /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi/all_roary_o -p /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi/prokka_o -o /home/liaoth/project/genome_pipelines/pipelines/test/test_toolkit/test_extract_region
    # g1 = 'fkpA_1'
    # g2 = 'lldP'
    # tab = "/home/liaoth/project/shenzhen_actineto/roary_o/gene_presence_absence.csv"
    # prokka_dir = "/home/liaoth/project/shenzhen_actineto/prokka_o"
    # outdir = "/home/liaoth/project/shenzhen_actineto/KL_extracted_reads"
    # base_dir = "%s-%s_region" % (g1, g2)
