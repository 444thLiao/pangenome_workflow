from glob import glob
import os, re
from subprocess import check_call, check_output
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
from toolkit.utils import get_length_fasta,get_locus2group,get_fasta_by_ID
from toolkit.get_gene_info import get_gene_with_regin

indir = "/home/liaoth/project/shenzhen_Acinetobacter/shovill_output"

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


roary_dir = '/home/liaoth/project/shenzhen_Acinetobacter/21_roary_o'
locus2group = get_locus2group(roary_dir)

plasmids_genes = defaultdict(dict)
for k, vdict in result_dict.items():
    all_plasmid_r = [region for _v in vdict.values() for region in _v]
    gff_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/{sn}.gff"
    all_genes = get_gene_with_regin(gff_p.format(sn=k),
                                    all_plasmid_r)
    for g in all_genes:
        g = locus2group.get(g, 'removed')
        if g != 'removed':
            plasmids_genes[k][g] = 1
############################################################

# summary all require info
samples21 = ['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
             '39502', '39528', '40067', '40116', '40615', '41194', '41296', '41833',
             '42121', '42268', '42349', '43458', '40283']
summary_df = pd.DataFrame(
    columns=['total contigs', 'total CDS', 'total length', 'contigs belong to plasmid', 'CDS belong to plasmid', 'total length of plasmids', 'ratio of plasmid'])
for sample_name in samples21:
    contig_pth = os.path.join(indir, 'regular', sample_name, 'contigs.fa')
    num_contigs = int(check_output("grep -c '^>' %s " % contig_pth, shell=True))
    length_contigs = sum(get_length_fasta(contig_pth).values())
    gff_pth = "/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/{sn}.gff".format(sn=sample_name)
    num_CDS = int(check_output("grep -c 'CDS' %s " % gff_pth, shell=True))

    plasmidcontig_pth = os.path.join(indir, 'plasmidsSpades', sample_name, 'contigs.fa')
    num_contigs4plasmid = sum([len(_) if sample_name in result_dict else 0 for _ in result_dict.get(sample_name,{}).values()])
    length_plasmidscontigs = sum(get_length_fasta(plasmidcontig_pth).values()) if sample_name in result_dict else 0
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

summary_df.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/plasmids_summary.csv",index=1)

plasmids_genes_presence_or_absecnces_df = pd.DataFrame.from_dict(plasmids_genes,orient='index')
from vis.heatmap import *

merged_df = pd.concat([plasmids_genes_presence_or_absecnces_df, total_df.loc[:, params.vf_cols]], axis=1)

tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/21_outgroup_roary_o/core_gene.newick"
rooted = 'AB030'
fig, mmatrix, amatrix = main(merged_df,
                             plasmids_genes_presence_or_absecnces_df.columns,
                             filter_sample=samples21,
                             tree_p=tree_p,
                             rooted=rooted,
                             accessory_cols=params.vf_cols,
                             width=2000, height=1000)
mmatrix.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/21_heatmap_with_metadata.csv", index=True)
plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/21_heatmap_with_metadata.html",
                    auto_open=False)
#
# all_results = []
# for col in tqdm(plasmids_genes_presence_or_absecnces_df.columns):
#     results = get_fasta_by_ID('/home/liaoth/data2/project/shenzhen_Acinetobacter/21_roary_o', col)
#     if len(results) != 1:
#         import pdb;pdb.set_trace()
#     else:
#         all_results += results
# from Bio import SeqIO
#
# with open('/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/IS_ref.fa', 'w') as f1:
#     SeqIO.write(all_results, f1, format='fasta')