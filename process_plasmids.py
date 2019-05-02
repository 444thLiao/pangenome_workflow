from glob import glob
import os, re
from subprocess import check_call, check_output
from collections import defaultdict
from tqdm import tqdm

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
roary_dir = '/home/liaoth/project/shenzhen_actineto/roary_o'
protein_cluster_info = os.path.join(roary_dir,'clustered_proteins')
locus2group = dict()
for row in open(protein_cluster_info).readlines():
    row = row.strip('\n')
    cluster_group = row.split(': ')[0]
    remained_locus = row.split(": ")[1]
    for locus in remained_locus.split('\t'):
        locus2group[locus] = cluster_group

plasmids_genes = defaultdict(dict)
for k, vdict in result_dict.items():
    all_plasmid_r = [region for _v in vdict.values() for region in _v]
    gff_p = "/home/liaoth/project/shenzhen_actineto/prokka_o/{sn}/{sn}.gff"
    all_genes = get_gene_with_regin(gff_p.format(sn=k),
                                    all_plasmid_r)
    for g in all_genes:
        g = locus2group.get(g,'removed')
        if g != 'removed':
            plasmids_genes[k][g] = 1


