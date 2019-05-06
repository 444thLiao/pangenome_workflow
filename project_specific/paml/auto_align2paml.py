import os
from glob import glob
from subprocess import check_call
from multiprocessing import Manager, Process, cpu_count, Pool
import pandas as pd
from Bio import AlignIO, SeqIO
from tqdm import tqdm

roary_dir = '/home/liaoth/project/shenzhen_actineto/21_outgroup_roary_o_test'


def get_locus2sample(roary_dir):
    ab_tab = os.path.join(roary_dir, "gene_presence_absence.csv")
    ab_df = pd.read_csv(ab_tab, sep=',', index_col=0)
    ab_df = ab_df.iloc[:,13:]
    query_dict = ab_df.to_dict('index')
    for group,adict in query_dict.items():
        new_dict = {}
        for k,v in adict.items():
            if not pd.isna(v):
                new_dict[v] = k
        query_dict[group] = new_dict
    return query_dict

def extract_core(roary_dir):
    ab_tab = os.path.join(roary_dir, "gene_presence_absence.Rtab")
    ab_df = pd.read_csv(ab_tab, sep='\t', index_col=0)
    core_genes = ab_df.index[ab_df.sum(1) == ab_df.shape[1]]
    return list(core_genes)


def _run_align(args):
    fa, outdir = args
    ofile = os.path.join(outdir, os.path.basename(fa) + '.aln')
    check_call("mafft --thread -1 --quiet %s > %s" % (fa, ofile), shell=True, executable='/usr/bin/zsh')

def run_alignment(indir, outdir):
    os.makedirs(outdir, exist_ok=True)

    differ_args = zip(glob(os.path.join(indir, '*.fa')),
                      [outdir]*len(glob(os.path.join(indir, '*.fa'))))
    with Pool(processes=20) as pool:
        results = list(tqdm(pool.imap(_run_align, differ_args),
                            total=len(glob(os.path.join(indir, '*.fa'))) ))

def remove_stop_codon(read):
    read.seq = read.seq.upper().ungap('-')
    for sc in ['TAA','TAG','TGA']:
        if read.seq.endswith(sc):
            read.seq = read.seq.rsplit(sc, 1)[0]
    return read


def get_aln_file(roary_dir, output_dir):
    core_genes = extract_core(roary_dir)
    os.makedirs(output_dir, exist_ok=True)
    for gene in core_genes:
        reads = []
        aln_file = os.path.join(roary_dir, 'pan_genome_sequences', '%s.fa.aln' % gene)
        if not os.path.isfile(aln_file):
            continue
        aln_obj = AlignIO.read(aln_file, format='fasta')
        for read in aln_obj:
            remove_stop_codon_read = remove_stop_codon(read)
            reads.append(remove_stop_codon_read)
            locus_id = remove_stop_codon_read.id
            remove_stop_codon_read.id = remove_stop_codon_read.description = remove_stop_codon_read.name = ''
            remove_stop_codon_read.id = query_dict[gene][locus_id]
        with open(os.path.join(output_dir, '%s.fa' % gene), 'w') as f1:
            SeqIO.write(reads, f1, format='fasta')

query_dict = get_locus2sample(roary_dir)
get_aln_file(roary_dir, "/home/liaoth/project/shenzhen_actineto/test_dnds/extract_read")

run_alignment("/home/liaoth/project/shenzhen_actineto/test_dnds/extract_read",
              "/home/liaoth/project/shenzhen_actineto/test_dnds/new_aln")
########################################################################################################################
