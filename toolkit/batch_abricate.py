import os
import re
import sys
from collections import defaultdict
from glob import glob
from multiprocessing import cpu_count

import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from tqdm import tqdm
from pandas.errors import EmptyDataError
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from toolkit.utils import get_locus2group, run_cmd, valid_path
from toolkit.get_gene_info import get_gff, add_fea2gff


def run_abricate(prokka_dir,
                 odir,
                 exe_path="abricate",
                 mincov=80,
                 thread=0,
                 dry_run=False,
                 log_file=None):
    """run abricate, and output to it dir with name like abricate_{db}.tab"""
    if thread == 0 or thread == -1:
        thread = cpu_count()
    samples2locus = defaultdict(list)
    for ffn_pth in tqdm(glob(os.path.join(prokka_dir,
                                          "*",
                                          '*.ffn'))):
        # ffn has been separated each locus into different record.
        sample_name = os.path.basename(os.path.dirname(ffn_pth))
        for db in ['card', 'argannot', 'ncbi', 'vfdb', 'vfdb_full', 'resfinder', 'plasmidfinder', 'victors']:
            cmd = "{exe_path} {infile} --threads {thread} --db {db} --mincov {mincov} --quiet > {ofile}"
            ofile = os.path.join(odir,
                                 sample_name,
                                 'abricate_{db}.tab'.format(db=db))
            valid_path(ofile, check_ofile=True)
            new_cmd = cmd.format(infile=ffn_pth,
                                 db=db,
                                 exe_path=exe_path,
                                 thread=thread,
                                 mincov=mincov,
                                 ofile=ofile)
            if not os.path.isfile(ofile):
                run_cmd(new_cmd,
                        dry_run=dry_run,
                        log_file=log_file)
            else:
                run_cmd("# %s has been analysised " % ofile,
                        dry_run=dry_run,
                        log_file=log_file)
            ############################################################
        samples2locus[sample_name] = [_.id for _ in SeqIO.parse(ffn_pth, format='fasta')]
    return samples2locus


def gene_process(gene):
    gene = gene.rsplit('.', 1)[0]
    if 'full' in gene:
        gene = gene.split('full_')[-1]
    if gene.endswith('_1') or gene.endswith('_2'):
        gene = re.split('_1', gene)[0]
        gene = re.split('_2', gene)[0]
    if gene.startswith('('):
        gene = gene.split(')')[1]
    if gene.startswith('bla'):
        gene = gene.split('bla')[1]
    return gene


def summary_abricate(odir):
    """summary output tab, and output locus2annotate and annotate2db"""
    locus2annotate = {}
    annotate2db = {}
    for tab in tqdm(glob(os.path.join(odir,
                                      '*',
                                      "abricate_*.tab"))):
        try:
            df = pd.read_csv(tab, sep='\t')
        except EmptyDataError:
            continue
        seq_id = list(df.loc[:, "SEQUENCE"])
        genes = [gene_process(_) for _ in df.loc[:, 'GENE']]
        db = os.path.basename(tab).split('_')[1].split('.')[0]
        annotate2db.update({g: db for g in genes})
        # only compare among one sample.
        prepare_update_dict = dict(zip(seq_id, genes))
        for k in seq_id:
            if k in locus2annotate.keys():
                ori = locus2annotate[k]
                aft = prepare_update_dict[k]
                if len(ori) < len(aft):
                    prepare_update_dict.pop(k)
        locus2annotate.update(prepare_update_dict)
    return locus2annotate, annotate2db


def refine_abricate_output(locus2annotate, roary_dir):
    """ given locus2annotate and roary_dir, it could use roary_cluster info to summary or unify the annotated genes."""
    locus2group = get_locus2group(roary_dir)
    group2annotate = defaultdict(list)
    for locus, annotate in locus2annotate.items():
        formatted_annotate = gene_process(annotate)
        if locus not in locus2group:
            # it mean the roary dir is not full
            print("unknown group for ", locus)
            pass
        else:
            group2annotate[locus2group[locus]].append(formatted_annotate)

    # for unify similary annotate genes
    rename_genes = defaultdict(list)
    # check for, exam one group match too many annotated genes
    for g, a in group2annotate.items():
        if len(set(a)) > 1:
            unify_name = min(set(a), key=lambda x: len(x))
            for _ in set(a):
                rename_genes[_] = unify_name
    return rename_genes


def output_abricate_result(samples2locus,
                           locus2annotate,
                           annotate2db,
                           rename_genes
                           ):
    ############################################################
    samples2genes = defaultdict(lambda: defaultdict(lambda: 0))
    locus2samples = {}
    for sample in samples2locus.keys():
        for locus in samples2locus[sample]:
            if locus in locus2annotate:
                annotate = locus2annotate[locus]
                formatted_annotate = rename_genes.get(annotate, annotate)
                samples2genes[sample][formatted_annotate] += 1
                locus2samples[locus] = sample
    for locus, annotate in locus2annotate.items():
        locus2annotate[locus] = rename_genes.get(annotate, annotate)
    ############################################################
    locus2annotate_df = pd.DataFrame.from_dict({'gene': locus2annotate})
    locus2annotate_df.loc[:, 'db'] = [annotate2db[_] for _ in locus2annotate_df.loc[:, 'gene']]
    locus2annotate_df.loc[:, 'sample'] = [locus2samples[_] for _ in locus2annotate_df.index]
    abricate_result = pd.DataFrame.from_dict(samples2genes,
                                             orient='index')
    abricate_result = abricate_result.fillna(0)
    abricate_result.loc['db', :] = [annotate2db[_]
                                    for _ in abricate_result.columns]

    return locus2annotate_df, abricate_result


if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='mainlty the prokka output, mean like indir/sample1/sample1.ffn; ffn is a fasta file(nucleotide) contain each locus.')
    parse.add_argument("-r", "--roary_dir", help='If you run roary yet, it could use cluster of roary to unify the result of abricate.')
    parse.add_argument("-o", "--outdir", help='output dir for abricate')
    parse.add_argument("--log", help="log file to write")
    parse.add_argument("-t", "--threads", help="how many threads you want to assign to abricate.", type=int)
    parse.add_argument("-db", "--database", help="the path of abricate database.")
    parse.add_argument("-mc", "--mincov", help="how many threads you want to assign to abricate.", default=80, type=float)
    parse.add_argument("-exe", "--abricate_path", help="the path of your abricate. Default is 'abricate' in you $PATH. ", default="abricate")
    parse.add_argument("--dry_run", help="Given this parameter, it will only log out to terminal or log file instead of running it.", action="store_true")

    args = parse.parse_args()
    indir = os.path.abspath(args.indir)
    roary_dir = os.path.abspath(args.roary_dir)
    odir = os.path.abspath(args.outdir)

    valid_path([indir, roary_dir], check_dir=True)
    valid_path([odir], check_odir=True)
    samples2locus = run_abricate(indir,
                                 odir=odir,
                                 exe_path=args.abricate_path,
                                 mincov=args.mincov,
                                 thread=args.threads,
                                 log_file=args.log,
                                 dry_run=args.dry_run
                                 )
    locus2annotate, annotate2db = summary_abricate(odir)
    if args.roary_dir:
        rename_genes = refine_abricate_output(locus2annotate, roary_dir)
    else:
        rename_genes = {}
    locus2annotate_df, abricate_result = output_abricate_result(samples2locus,
                                                                locus2annotate,
                                                                annotate2db,
                                                                rename_genes)
    l2a_pth = os.path.join(odir, 'locus2annotate.csv')
    s2g_pth = os.path.join(odir, "samples2annotate.csv")
    locus2annotate_df.to_csv(l2a_pth, index=1)
    abricate_result.to_csv(s2g_pth, index=1)
    ############################################################
    # reannotated gff
    prepare_dict = {}
    for sn in map(str, locus2annotate_df['sample'].unique()):
        gff_pth = os.path.join(indir, "{sn}", '{sn}.gff').format(sn=sn)
        gff_db = get_gff(gff_pth, mode='db')
        gff_obj = get_gff(gff_pth, mode='bcbio')
        # remove all existing features
        for record in gff_obj.values():
            record.features.clear()
            record.annotations['sequence-region'].clear()
            record.annotations['sequence-region'] = "%s %s %s" % (record.id, 1, len(record))
        prepare_dict[sn] = (gff_db, gff_obj)
    res_db = ["card", "ncbi", "resfinder", "argannot"]
    vf_db = ["vfdb_full", "vfdb", "victors"]
    for locus, vals in locus2annotate_df.iterrows():
        sn = str(vals["sample"])
        gff_db, gff_obj = prepare_dict[sn]
        locus_info = gff_db[locus]
        contig = locus_info.chrom
        db = vals["db"]
        type = "RES_gene" if db in res_db else "VF_gene"
        add_fea2gff(gff_obj[contig],
                    locus_info.start,
                    locus_info.end,
                    ID=vals["gene"],
                    strand=-1 if locus_info.strand == '-' else 1,
                    type=type,
                    source="abricate",
                    db=db)
    for sn, vals in prepare_dict.items():
        ogff = os.path.join(odir, sn, "%s.gff" % sn)
        with open(ogff, 'w') as f1:
            GFF.write(list(vals[1].values()), f1)
    # prokka_dir = "/home/liaoth/project/shenzhen_actineto/KL_extracted_reads/fkpA_1-lldP_region/prokka_o"
    # samples2locus = run_abricate(prokka_dir)
    # locus2annotate,annotate2db = summary_abricate(prokka_dir)
