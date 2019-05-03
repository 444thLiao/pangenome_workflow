from glob import glob
import os, re
from subprocess import check_call, check_output
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
from utils import get_length_fasta

indir = "/home/liaoth/project/shenzhen_actineto/shovill_output"

result_dict = {}
for contig in tqdm(glob(os.path.join(indir, 'regular', '*', 'contigs.fa'))):
    check_call("bwa index %s" % contig, shell=True, executable='/usr/bin/zsh')
    plasmid_contig = contig.replace('/regular/', '/plasmidsSpades/') + 'sta'
    if os.path.getsize(plasmid_contig) > 0:
        result = check_output("bwa mem -x intractg {regular_one} {plasmid_one}".format(regular_one=contig,
                                                                                       plasmid_one=plasmid_contig),
                              shell=True, executable='/usr/bin/zsh')
        result = result.decode('utf-8')
        result = [_ for _ in result.split('\n') if not _.startswith('@') and _]
        match_plasmid_row = [row.split('\t')[0] for row in result]
        match_contig_row = [row.split('\t')[2] for row in result]
        plasmid_count_dict = defaultdict(list)
        for p, c in zip(match_plasmid_row, match_contig_row):
            num_p = re.findall("component_[0-9]+_pilon$", p)[0]

            plasmid_count_dict[num_p].append(c)
        sample_name = os.path.basename(os.path.dirname(contig))
        result_dict[sample_name] = plasmid_count_dict

from toolkit.get_gene_info import get_gene_with_regin
from utils import get_locus2group

roary_dir = '/home/liaoth/project/shenzhen_actineto/roary_o'
locus2group = get_locus2group(roary_dir)

plasmids_genes = defaultdict(dict)
for k, vdict in result_dict.items():
    all_plasmid_r = [region for _v in vdict.values() for region in _v]
    gff_p = "/home/liaoth/project/shenzhen_actineto/prokka_o/{sn}/{sn}.gff"
    all_genes = get_gene_with_regin(gff_p.format(sn=k),
                                    all_plasmid_r)
    for g in all_genes:
        g = locus2group.get(g, 'removed')
        if g != 'removed':
            plasmids_genes[k][g] = 1
############################################################
# summary all require info
summary_df = pd.DataFrame(
    columns=['total contigs', 'total CDS', 'total length', 'contigs belong to plasmid', 'CDS belong to plasmid', 'total length of plasmids', 'ratio of plasmid'])
for sample_name in plasmids_genes.keys():
    contig_pth = os.path.join(indir, 'regular', sample_name, 'contigs.fa')
    num_contigs = int(check_output("grep -c '^>' %s " % contig_pth, shell=True))
    length_contigs = sum(get_length_fasta(contig_pth).values())
    gff_pth = "/home/liaoth/project/shenzhen_actineto/prokka_o/{sn}/{sn}.gff".format(sn=sample_name)
    num_CDS = int(check_output("grep -c 'CDS' %s " % gff_pth, shell=True))

    plasmidcontig_pth = os.path.join(indir, 'plasmidsSpades', sample_name, 'contigs.fa')
    num_contigs4plasmid = sum([len(_) for _ in result_dict[sample_name].values()])
    length_plasmidscontigs = sum(get_length_fasta(plasmidcontig_pth).values())
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
