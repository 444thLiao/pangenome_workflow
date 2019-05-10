import itertools
import os
import re
from glob import glob
from subprocess import check_call

import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.io as pio
from Bio import AlignIO
from Bio.Phylo.PAML import yn00
from tqdm import tqdm

from toolkit.utils import get_locus2group


def ttt():
    indir = "/home/liaoth/project/shenzhen_actineto/test_dnds/new_aln"

    odir = '/home/liaoth/project/shenzhen_actineto/test_dnds/yn00_r'
    os.makedirs(odir, exist_ok=True)
    nuc_file = "/home/liaoth/temp_/test.nuc"

    count_error = 0
    for aln_file in tqdm(glob(os.path.join(indir, '*.fa.aln'))):
        gene_name = os.path.basename(aln_file).split('.fa')[0]
        aln_obj = AlignIO.read(aln_file, format='fasta')
        length_ = len(aln_obj[0])
        add_str = ''
        if length_ % 3 != 0:
            add_str = '-' * (3 - (length_ % 3))
            result = '%s %s\n' % (len(aln_obj), length_ + len(add_str))
        else:
            result = '%s %s\n' % (len(aln_obj), length_)

        for record in aln_obj:
            result += '%s\n' % (record.id.strip('-'))
            result += '%s\n' % (str(record.seq) + add_str)
        with open(nuc_file, 'w') as f1:
            f1.write(result)

        ctl = open('/home/liaoth/tools/paml4.9i/test.ctl').read()
        new_ctl = re.sub('outfile = (.*)', 'outfile = %s' % (os.path.join(odir, 'test_%s' % gene_name)), ctl)
        with open('/home/liaoth/tools/paml4.9i/test.ctl', 'w') as f1:
            f1.write(new_ctl)
        try:
            check_call('/home/liaoth/tools/paml4.9i/bin/yn00 /home/liaoth/tools/paml4.9i/test.ctl', shell=True, executable='/usr/bin/zsh')
        except:
            count_error += 1
            pass


def parse_yn00_r(pth):
    try:
        result = yn00.read(pth)
        return result
    except:
        return None


def get_all_result(indir):
    gene_dict = {}
    for result_pth in tqdm(glob(os.path.join(indir, 'test_*'))):
        gene_name = os.path.basename(result_pth).split('_', 1)[-1]
        result = parse_yn00_r(result_pth)
        if result is not None:
            gene_dict[gene_name] = result
    return gene_dict


locus2group = get_locus2group('/home/liaoth/project/shenzhen_actineto/21_outgroup_roary_o_test')
group2locus = {g: l for l, g in locus2group.items()}  # it will remove redundancy, it does not matter
locus2annotate_df = pd.read_csv("/home/liaoth/project/shenzhen_actineto/all_65_annotated_locus2annotate.tab", sep='\t', index_col=0)
locus2annotate = dict(zip(locus2annotate_df.index, locus2annotate_df.loc[:, '0']))

group1 = ['40615', '41833', '41194', '39058', '40116', '42121', '42268']
group2 = ['39292', '39502', '41296', '34960', '35082', '39054', '40067',
          '42349', '43458', '37706', '37502', '35989', '39528']
gene_dict = get_all_result('/home/liaoth/project/shenzhen_actineto/test_dnds/yn00_r')


def get_dnds_group(group1, result, group2=None):
    dnds_collect = []
    num_syn_collect = []
    num_nonsyn_collect = []
    if group2 is None:
        iter_obj = itertools.combinations(group1, 2)
    else:
        iter_obj = itertools.product(group1, group2)
    for g1, g2 in iter_obj:
        current_compare = result[g1][g2]
        if current_compare['YN00']['dS'] != 0:
            dnds_collect.append(current_compare['YN00']['omega'])
        num_syn_collect.append(current_compare['YN00']['dS'] * current_compare['YN00']['S'])
        num_nonsyn_collect.append(current_compare['YN00']['dN'] * current_compare['YN00']['N'])
    return dnds_collect, num_syn_collect, num_nonsyn_collect


def get_gene_info4group(group1, result_dict, genes=None, group2=None):
    gene_dict = {}
    if genes is None:
        iter_obj = tqdm(result_dict.keys())
    else:
        iter_obj = tqdm({gene: v for gene, v in result_dict.items()})
    for gene_name in iter_obj:
        result = result_dict[gene_name]
        if group2 is not None:
            dnds_c, num_syn, num_nonsyn = get_dnds_group(group1, result, group2=group2)
        else:
            dnds_c, num_syn, num_nonsyn = get_dnds_group(group1, result, group2=group2)
        gene_dict[gene_name] = (dnds_c,
                                num_syn,
                                num_nonsyn)
    return gene_dict


total_result = get_gene_info4group(group1 + group2, gene_dict)

mean_syn = np.mean([np.mean(_[1]) for _ in total_result.values()])
std_syn = np.std([np.mean(_[1]) for _ in total_result.values()])
mean_nonsyn = np.mean([np.mean(_[2]) for _ in total_result.values()])
std_nonsyn = np.std([np.mean(_[2]) for _ in total_result.values()])
mean_dnds = np.nanmean([np.nanmean(_[0]) for _ in total_result.values()])
std_dnds = np.nanstd([np.nanmean(_[0]) for _ in total_result.values()])

result1 = get_gene_info4group(group1, gene_dict)
result2 = get_gene_info4group(group2, gene_dict)
inter_result = get_gene_info4group(group1, gene_dict, group2=group2)

from scipy.stats import ranksums

collect_num_syn_diff = []
collect_dnds_diff = []
collect_a = []
for g in result1.keys():
    g1_dnds, g1_num_syn, g1_num_nonsyn = result1[g]
    g2_dnds, g2_num_syn, g2_num_nonsyn = result2[g]
    g3_dnds, g3_num_syn, g3_num_nonsyn = inter_result[g]

    MK_test_a = 1 - ((np.nanmean(g3_num_syn )* (np.nanmean(g2_num_nonsyn)+np.nanmean(g1_num_nonsyn))/2 )/ (np.nanmean(g3_num_nonsyn )* (np.nanmean(g2_num_syn)+np.nanmean(g1_num_syn))/2 ))
    collect_a.append((g,MK_test_a))
    if ranksums(g3_num_syn, g2_num_syn).pvalue < 0.05 / len(result1) and \
            ranksums(g3_num_syn, g1_num_syn).pvalue < 0.05 / len(result1):
        collect_num_syn_diff.append(g)
    if ranksums(g3_dnds, g2_dnds).pvalue < 0.05 / len(result1) and \
            ranksums(g3_dnds, g1_dnds).pvalue < 0.05 / len(result1):
        collect_dnds_diff.append(g)

print(len(collect_dnds_diff), len(collect_num_syn_diff))
idx = 0
print(np.mean(result1[collect_num_syn_diff[idx]][1]),
      np.mean(result2[collect_num_syn_diff[idx]][1]),
      np.mean(inter_result[collect_num_syn_diff[idx]][1]))

print(np.mean(result1[collect_dnds_diff[0]][0]),
      np.mean(result2[collect_dnds_diff[0]][0]),
      np.mean(inter_result[collect_dnds_diff[0]][0]))
############################################################
odir = "/home/liaoth/project/shenzhen_actineto/test_dnds/graph"
def plot_attr(name,title,idx):
    fig = plotly.tools.make_subplots(3, 1, shared_xaxes=True, subplot_titles=["Intra group1",
                                                                              "Inter groups",
                                                                              "Intra group2"])
    fig.append_trace(go.Scatter(x=[g for g, v in result1.items()],
                                y=[np.nanmean(v[idx]) for g, v in result1.items()],
                                mode='markers',
                                legendgroup='intragroup1',
                                name='intra group1'), 1, 1)
    fig.append_trace(go.Scatter(x=[list(result1.keys())[0],list(result1.keys())[-1]],
                                y=[np.nanmean([np.nanmean(v[idx]) for g, v in result1.items()])]*2,
                                mode='lines',
                                showlegend=False,
                                legendgroup='intragroup1',
                                name='intra group1'), 1, 1)
    fig.append_trace(go.Scatter(x=[g for g, v in result2.items()],
                                y=[np.nanmean(v[idx]) for g, v in result2.items()],
                                mode='markers',
                                legendgroup='intragroup2',
                                name='intra group2'), 3, 1)
    fig.append_trace(go.Scatter(x=[list(result2.keys())[0],list(result2.keys())[-1]],
                                y=[np.nanmean([np.nanmean(v[idx]) for g, v in result2.items()])]*2,
                                mode='lines',
                                showlegend=False,
                                legendgroup='intragroup2',
                                name='intra group2'), 3, 1)
    fig.append_trace(go.Scatter(x=[g for g, v in inter_result.items()],
                                y=[np.nanmean(v[idx]) for g, v in inter_result.items()],
                                mode='markers',
                                legendgroup='intergroup',
                                name='inter groups'), 2, 1)
    fig.append_trace(go.Scatter(x=[list(inter_result.keys())[0],list(inter_result.keys())[-1]],
                                y=[np.nanmean([np.nanmean(v[idx]) for g, v in inter_result.items()])]*2,
                                mode='lines',
                                showlegend=False,
                                legendgroup='intergroup',
                                name='inter groups'), 2, 1)
    fig.layout.title = title
    fig.layout.font.size = 15
    fig.layout.width = 1000
    fig.layout.height = 1000
    pio.write_image(fig, file=os.path.join(odir, '%s.png' % name))
    plotly.offline.plot(fig, filename=os.path.join(odir, '%s.html' % name), auto_open=False)
plot_attr('DnDs',"Compare of DnDs",0)
plot_attr('num_S','compare of number of Synonymous substitutions',1)
############################################################
fig = go.Figure()
collect_a = [_ for _ in collect_a if not np.isinf(_[1]) and not pd.isna(_[1])]
fig.add_bar(x=[_[0] for _ in sorted(collect_a,key=lambda x:x[1])],
            y=[_[1] for _ in sorted(collect_a,key=lambda x:x[1])],
            )
fig.layout.width = 2000
fig.layout.height = 1000
pio.write_image(fig, file=os.path.join(odir, 'WK_test.png'))
plotly.offline.plot(fig, filename=os.path.join(odir, 'WK_test.html'), auto_open=False)

############################################################
fig = plotly.tools.make_subplots(3, 1, shared_xaxes=True, subplot_titles=["Intra group1",
                                                                          "Inter groups",
                                                                          "Inter group2"])
fig.append_trace(go.Scatter(x=[g for g, v in result1.items()],
                            y=[np.nanmean(v[1]) for g, v in result1.items()],
                            mode='markers',
                            legendgroup='intragroup1',
                            name='intra group1'), 1, 1)
fig.append_trace(go.Scatter(x=[g for g, v in result2.items()],
                            y=[np.nanmean(v[1]) for g, v in result2.items()],
                            mode='markers',
                            legendgroup='intragroup2',
                            name='intra group2'), 3, 1)
fig.append_trace(go.Scatter(x=[g for g, v in inter_result.items()],
                            y=[np.nanmean(v[1]) for g, v in inter_result.items()],
                            mode='markers',
                            legendgroup='intergroup',
                            name='inter groups'), 2, 1)
fig.layout.title =
fig.layout.font.size = 15
fig.layout.width = 800
fig.layout.height = 800
pio.write_image(fig, file=os.path.join(odir, 'num_S.png'))
plotly.offline.plot(fig, filename=os.path.join(odir, 'num_S.html'), auto_open=False)

############################################################
for i in collect_num_syn_diff:
    locus = group2locus[i]
    if locus in locus2annotate:
        print(locus, locus2annotate[locus])

############################################################
# covert to df
fig = go.Figure()
for g in ds_gene:
    fig.add_bar(name=g,
                y=[np.mean(ds_gene[g])],
                showlegend=False)
plotly.offline.plot(fig, filename='/home/liaoth/temp_/test.html', auto_open=False)

############################################################
# plot result
samples21 = ['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
             '39502', '39528', '40067', '40116', '40615', '41194', '41296', '41833',
             '42121', '42268', '42349', '43458', '40283']
KL_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/fkpA_1-lldP_region/roary_o/gene_presence_absence.csv",
                    sep=',', index_col=0)
KL_df = KL_df.iloc[:, 13:].T
merged_df = pd.concat([KL_df, total_df.loc[:, params.vf_cols]], axis=1)
fig, mmatrix, amatrix = main(merged_df,
                             KL_df.columns,
                             filter_sample=samples21,
                             # tree_p="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/fkpA_1-lldP_region/roary_o/core_gene.newick",
                             accessory_cols=params.vf_cols,
                             width=2000, height=1000)
mmatrix.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/21_heatmap_with_metadata_KLsnp.csv", index=True)
plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/21_heatmap_with_metadata_KLsnp.html",
                    auto_open=False)
