from glob import glob
import os, re
from subprocess import check_output
from collections import defaultdict
from tqdm import tqdm
import pandas as pd

if __name__ == "__main__" and __package__ is None:
    import sys
    sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))

from toolkit.utils import run_cmd, get_locus2group, get_length_fasta
from toolkit.get_gene_info import get_gene_with_regin


def get_plasmids(indir):
    """ indir may end with assembly_o"""
    result_dict = {}
    for contig in tqdm(glob(os.path.join(indir, 'regular', '*', 'contigs.fa'))):
        run_cmd("bwa index %s" % contig, dry_run=False)
        plasmid_contig = contig.replace('/regular/', '/plasmidsSpades/') + 'sta'
        # regular is 'contigs.fa', plasmids must is 'contigs.fasta'
        if os.path.getsize(plasmid_contig) > 0:
            result = check_output("bwa mem -x intractg {regular_one} {plasmid_one}".format(regular_one=contig,
                                                                                           plasmid_one=plasmid_contig),
                                  shell=True,
                                  executable='/usr/bin/zsh')
            result = result.decode('utf-8')
            result = [_ for _ in result.split('\n') if not _.startswith('@') and _]
            match_plasmid_row = [row.split('\t')[0] for row in result]
            match_contig_row = [row.split('\t')[2] for row in result]
            match_leftmost_pos = [row.split('\t')[2] for row in result]  # 1-coordinate

            plasmid_count_dict = defaultdict(list)
            for p, c in zip(match_plasmid_row,
                            match_contig_row):
                num_p = re.findall("component_[0-9]+_pilon$", p)[0]  # get plasmids num/ID
                plasmid_count_dict[num_p].append(c)
            sample_name = os.path.basename(os.path.dirname(contig))
            result_dict[sample_name] = plasmid_count_dict
    return result_dict


def get_gene_in_plasmids(plasmids_dict, locus2group, prokka_dir):
    plasmids_genes = defaultdict(dict)
    for sn, p2contig in plasmids_dict.items():
        all_plasmid_r = [region for _v in p2contig.values() for region in _v]
        gff_p = os.path.join(prokka_dir, "{sn}/{sn}.gff")
        if not os.path.isfile(gff_p.format(sn=sn)):
            gff_p = os.path.join(prokka_dir, "{sn}.gff")
        if not os.path.isfile(gff_p.format(sn=sn)):
            raise Exception("weird prokka input")
        all_genes = get_gene_with_regin(gff_p.format(sn=sn),
                                        all_plasmid_r)
        for g in all_genes:
            g = locus2group.get(g, 'removed')  # todo: process these lost genes
            if g != 'removed':
                plasmids_genes[sn][g] = 1
    return plasmids_genes


def main(indir, roary_dir, prokka_dir, odir):
    plasmids_dict = get_plasmids(indir)
    locus2group = get_locus2group(roary_dir)
    plasmids_genes = get_gene_in_plasmids(plasmids_dict,
                                          locus2group,
                                          prokka_dir)
    summary_df = pd.DataFrame(
        columns=['total contigs',
                 'total CDS',
                 'total length',
                 'contigs belong to plasmid',
                 'CDS belong to plasmid',
                 'total length of plasmids',
                 'ratio of plasmid'])
    for sample_name in plasmids_genes.keys():
        contig_pth = os.path.join(indir, 'regular', sample_name, 'contigs.fa')
        num_contigs = int(check_output("grep -c '^>' %s " % contig_pth, shell=True))
        length_contigs = sum(get_length_fasta(contig_pth).values())
        gff_p = os.path.join(prokka_dir, "{sn}/{sn}.gff")
        if not os.path.isfile(gff_p.format(sn=sample_name)):
            gff_p = os.path.join(prokka_dir, "{sn}.gff")
        if not os.path.isfile(gff_p.format(sn=sample_name)):
            raise Exception("weird prokka input")
        gff_pth = gff_p.format(sn=sample_name)
        num_CDS = int(check_output("grep -c 'CDS' %s " % gff_pth, shell=True))
        plasmidcontig_pth = os.path.join(indir, 'plasmidsSpades', sample_name, 'contigs.fa')
        length_plasmidscontigs = sum(get_length_fasta(plasmidcontig_pth).values())
        num_contigs4plasmid = sum([len(_) for _ in plasmids_dict[sample_name].values()])
        _sub_df = pd.DataFrame(columns=summary_df.columns,
                               index=[sample_name],
                               data=[[num_contigs,
                                      num_CDS,
                                      length_contigs,
                                      num_contigs4plasmid,
                                      len(plasmids_genes[sample_name]),
                                      length_plasmidscontigs,
                                      length_plasmidscontigs / length_contigs
                                      ]])
        summary_df = summary_df.append(_sub_df)
    # todo: process output; except a statistic info(summary df), (plasmids_genes) also important.
    summary_pth = os.path.join(odir,"plasmid_summary.csv")
    summary_df.to_csv(summary_pth,sep=',')

if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-r", "--roary_dir", help='')
    parse.add_argument("-p", "--prokka_dir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    indir = os.path.abspath(args.indir)
    odir = os.path.abspath(args.outdir)
    roary_dir = os.path.abspath(args.roary_dir)
    prokka_dir = os.path.abspath(args.prokka_dir)
    # import pdb;pdb.set_trace()
    # import code;code.interact(local=vars())
    main(indir, roary_dir, prokka_dir, odir)
