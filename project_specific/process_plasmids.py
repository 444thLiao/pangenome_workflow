import os
import re
from collections import defaultdict
from glob import glob
from subprocess import check_call, check_output

import pandas as pd
import plotly.io as pio
from tqdm import tqdm

from toolkit.get_gene_info import get_gene_with_regin
from toolkit.utils import get_length_fasta, get_locus2group

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

roary_dir = '/home/liaoth/project/shenzhen_Acinetobacter/20_outgroup_roary_o'
locus2group = get_locus2group(roary_dir)

locus2annotate_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/all_65_annotated_locus2annotate.tab", sep='\t', index_col=0)
locus2annotate = dict(zip(locus2annotate_df.index, locus2annotate_df.loc[:, 'gene']))
samples2genes_inplasmid = defaultdict(dict)
for k, vdict in result_dict.items():
    all_plasmid_r = [region for _v in vdict.values() for region in _v]
    gff_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/{sn}.gff"
    all_genes = get_gene_with_regin(gff_p.format(sn=k),
                                    all_plasmid_r)
    for g in all_genes:
        if g in locus2annotate.keys():
            g = locus2annotate[g]
        else:
            g = locus2group.get(g, 'removed')
        if g != 'removed':
            samples2genes_inplasmid[k][g] = 1
############################################################

plasmids_genes_presence_or_absecnces_df = pd.DataFrame.from_dict(samples2genes_inplasmid, orient='index')
plasmids_genes_presence_or_absecnces_df = plasmids_genes_presence_or_absecnces_df.fillna(0)
from vis.heatmap import *

query_dict = params.res2fun.copy()
query_dict.update(params.vf2fun)
merged_df = pd.concat([plasmids_genes_presence_or_absecnces_df, total_df.loc[:, params.vf_cols]], axis=1)

tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/20_outgroup_roary_o/core_gene.newick"
rooted = 'AB030'
fig, mmatrix_plasmid, amatrix_plasmid = main(merged_df,
                                             plasmids_genes_presence_or_absecnces_df.columns,
                                             filter_sample=samples21,
                                             tree_p=tree_p,
                                             rooted=rooted,
                                             up_cluster=True,
                                             map_type=lambda x: query_dict.get(x, "unknown"),
                                             accessory_cols=params.vf_cols,
                                             width=2000, height=1000)
mmatrix_plasmid.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/20_heatmap_with_metadata.csv", index=True)
plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/20_heatmap_with_metadata.html",
                    auto_open=False)

# summary all require info
samples21 = ['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
             '39502', '39528', '40067', '40116', '40615', '41194', '41296', '41833',
             '42121', '42268', '42349', '43458', '40283']
summary_df = pd.DataFrame(
    columns=['total contigs',
             'total CDS',
             'total length',
             "number of plasmid",
             'contigs belong to plasmid',
             'CDS belong to plasmid',
             'total length of plasmids',
             'ratio of plasmid'])
for sample_name in samples21:
    contig_pth = os.path.join(indir, 'regular', sample_name, 'contigs.fa')
    num_contigs = int(check_output("grep -c '^>' %s " % contig_pth, shell=True))
    length_contigs = sum(get_length_fasta(contig_pth).values())
    gff_pth = "/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/{sn}.gff".format(sn=sample_name)
    num_CDS = int(check_output("grep -c 'CDS' %s " % gff_pth, shell=True))

    num_plasmid = len(set([_.split('_')[1] for _ in result_dict[sample_name].keys()])) if sample_name in result_dict.keys() else 0
    plasmidcontig_pth = os.path.join(indir, 'plasmidsSpades', sample_name, 'contigs.fa')
    num_contigs4plasmid = sum([len(_) if sample_name in result_dict else 0 for _ in result_dict.get(sample_name, {}).values()])
    length_plasmidscontigs = sum(get_length_fasta(plasmidcontig_pth).values()) if sample_name in result_dict else 0

    _sub_df = pd.DataFrame(columns=summary_df.columns,
                           index=[sample_name],
                           data=[[num_contigs,
                                  num_CDS,
                                  length_contigs,
                                  num_plasmid,
                                  num_contigs4plasmid,
                                  len(samples2genes_inplasmid[sample_name]),
                                  length_plasmidscontigs,
                                  length_plasmidscontigs / length_contigs
                                  ]])
    summary_df = summary_df.append(_sub_df)

_t = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/21_heatmap_with_metadata.csv", index_col=0)

summary_df.reindex(_t.index).to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/plasmids_summary.csv", index=1)

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

############################################################
# annotated RES html
locus2annotate_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/all_65_annotated_locus2annotate.tab", sep='\t', index_col=0)
locus2annotate = dict(zip(locus2annotate_df.index, locus2annotate_df.loc[:, 'gene']))

samples = samples20
mode = 'RES'
f_all = True
tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/20_outgroup_roary_o/core_gene.newick"
rooted = 'AB030'
if mode == 'VF':
    genes = vf_genes
    accessory_cols = params.vf_cols
    map_type = lambda x: params.vf2fun.get(x, "unknown")
elif mode == 'RES':
    genes = res_genes
    accessory_cols = params.res_cols
    map_type = lambda x: params.res2fun.get(x, "unknown")
fig, mmatrix, amatrix = main(total_df,
                             genes,
                             filter_sample=samples,
                             filter_fea_all=f_all,
                             tree_p=tree_p,
                             accessory_cols=accessory_cols,
                             map_type=map_type,
                             rooted=rooted,
                             up_cluster=True,
                             width=2000,
                             height=1000)
import plotly.graph_objs as go

xs = []
ys = []
for g in mmatrix_plasmid.columns:
    if g in mmatrix:
        samples_in_plasmid = mmatrix_plasmid.index[(mmatrix_plasmid.loc[:, g] != "nan<Br>NaN") & (~pd.isna(mmatrix_plasmid.loc[:, g]))]
        # y coordinate
        x_coord = g
        xs += [g] * len(samples_in_plasmid)
        ys += list(samples_in_plasmid)

convert_y = dict(zip(fig.layout.yaxis2.ticktext, fig.layout.yaxis2.tickvals))
convert_x = dict(zip(fig.layout.xaxis3.ticktext, fig.layout.xaxis3.tickvals))
trace1 = go.Scatter(
    x=list(map(lambda x: convert_x[x], xs)),
    y=list(map(lambda x: convert_y[x], ys)),
    mode='text',
    text='P',
    showlegend=False,
    textposition='middle center',
    textfont=dict(size=15)
)
fig.append_trace(trace1, 2, 3)
plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/20_RES_with_plasmid.html")
pio.write_image(fig, file="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/20_RES_with_plasmid.png")
############################################################

samples = samples20
mode = 'VF'
f_all = True
tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/20_outgroup_roary_o/core_gene.newick"
rooted = 'AB030'
if mode == 'VF':
    genes = vf_genes
    accessory_cols = params.vf_cols
    map_type = lambda x: params.vf2fun.get(x, "unknown")
elif mode == 'RES':
    genes = res_genes
    accessory_cols = params.res_cols
    map_type = lambda x: params.res2fun.get(x, "unknown")
fig, mmatrix, amatrix = main(total_df,
                             genes,
                             filter_sample=samples,
                             filter_fea_all=f_all,
                             tree_p=tree_p,
                             accessory_cols=accessory_cols,
                             map_type=map_type,
                             rooted=rooted,
                             up_cluster=True,
                             width=2000,
                             height=1000)
import plotly.graph_objs as go

xs = []
ys = []
for g in mmatrix_plasmid.columns:
    if g in mmatrix:
        samples_in_plasmid = mmatrix_plasmid.index[(mmatrix_plasmid.loc[:, g] != "nan<Br>NaN") & (~pd.isna(mmatrix_plasmid.loc[:, g]))]
        # y coordinate
        x_coord = g
        xs += [g] * len(samples_in_plasmid)
        ys += list(samples_in_plasmid)

convert_y = dict(zip(fig.layout.yaxis2.ticktext, fig.layout.yaxis2.tickvals))
convert_x = dict(zip(fig.layout.xaxis3.ticktext, fig.layout.xaxis3.tickvals))
trace1 = go.Scatter(
    x=list(map(lambda x: convert_x[x], xs)),
    y=list(map(lambda x: convert_y[x], ys)),
    mode='text',
    text='P',
    showlegend=False,
    textposition='middle center',
    textfont=dict(size=15)
)
fig.append_trace(trace1, 2, 3)
plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/20_VF_with_plasmid.html")
pio.write_image(fig, file="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/plasmids_genes/20_VF_with_plasmid.png")
