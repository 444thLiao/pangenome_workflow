# -*- coding: UTF-8 -*-
import itertools
import re
from collections import defaultdict
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import plotly
import plotly.graph_objs as go
from ete3 import Tree

############################################################
med_data = pd.read_csv("/home/liaoth/data2/project/shenzhen_actineto/RS_judge/RS.csv", index_col=0)
med_data.index = med_data.index.astype(int).astype(str)
med_data = med_data.loc[:, [_ for _ in med_data.columns if _.startswith('RES_')]]
############################################################
total_csv = '/home/liaoth/data2/project/shenzhen_actineto/汇总结果/For phandango.csv'
total_df = pd.read_csv(total_csv, index_col=0)
res_name_convert = pd.read_csv('/home/liaoth/data2/project/shenzhen_actineto/gene_name_convert.csv', index_col=0)
res_name_convert = res_name_convert.to_dict()['0']
cols = [_ for _ in total_df.columns if _.startswith('abricate')]

med_data = med_data.reindex(total_df.index)
med_data = med_data.fillna('no')
med_data = med_data.loc[:, [col for col, val in med_data.iteritems() if not (val == 'no').all()]]

all_gene = {}
for col in cols:
    for db in ['argannot', 'card', 'ncbi', 'resfinder', 'vfdb', 'vfdb_full', 'victors']:
        if db in col:
            tmp = col.split(db + '_')[1]
            tmp = tmp.split('.')[0]
            if 'full' in tmp:
                tmp = tmp.split('full_')[-1]
            if 'resfinder' == db:
                tmp = re.split('_[0-9]', tmp)[0]
            if tmp.startswith('('):
                tmp = tmp.split(')')[1]
            if tmp.startswith('bla'):
                tmp = tmp.split('bla')[1]
            if col.endswith(tmp):
                all_gene[tmp.lower()] = col

drop_duplicated = ['abricate_ncbi_mph(E)',
                   'abricate_ncbi_msr(E)',
                   'abricate_ncbi_aac(3)-I',
                   'abricate_argannot_(Sul)SulI',
                   'abricate_argannot_(Sul)SulII',
                   'abricate_argannot_(Tet)Tet-39',
                   'abricate_argannot_(Tet)TetB',
                   'abricate_argannot_(AGly)Aac3-I',
                   "abricate_argannot_(AGly)Aph3''Ia"]

single_gene_cols = [_ for _ in all_gene.values() if _ not in drop_duplicated]

res_db = "_card|ncbi|resinder|argannot_"
vf_db = "_vfdb_full|vfdb|victors_"
for i in [_ for _ in single_gene_cols if re.search(vf_db, _)]:
    tmp = '_'.join(i.split('_')[2:])
    if tmp.startswith('full_'):
        tmp = tmp.split('full_')[1]
    res_name_convert[i] = tmp

res_genes = [res_name_convert.get(_, _) for _ in single_gene_cols if re.search(res_db, _)]
vf_genes = [res_name_convert.get(_, _) for _ in single_gene_cols if re.search(vf_db, _)]
total_df.columns = [res_name_convert.get(_, _) for _ in total_df.columns]
res_genes = [_ for _ in res_genes if pd.isnull(total_df.loc[:, _]).any() and not pd.isnull(total_df.loc[:, _]).all()]
vf_genes = [_ for _ in vf_genes if pd.isnull(total_df.loc[:, _]).any() and not pd.isnull(total_df.loc[:, _]).all()]

############################################################


# absolute difference
# group_info = {'collect': lambda x: bool(re.findall('^[0-9]*$', x)),
#               'ST457': lambda x: len(re.findall('^[0-9]*$', x)) == 0,
#               }
# sub_df = total_df.loc[:, single_gene_cols]
# for g, f in group_info.items():
#     match_rows = sub_df.loc[map(f, sub_df.index), :]
#     mismatch_rows = sub_df.loc[sub_df.index.difference(match_rows.index), :]
#     for col in sub_df.columns:
#         match_all_miss = pd.isnull(match_rows.loc[:, col]).all()
#         match_all_own = (~pd.isnull(match_rows.loc[:, col])).all()
#         mismatch_all_miss = pd.isnull(mismatch_rows.loc[:, col]).all()
#         mismatch_all_own = (~pd.isnull(mismatch_rows.loc[:, col])).all()
#         if match_all_miss and mismatch_all_own:
#             print(col)
#         elif match_all_own and mismatch_all_miss:
#             print(col)

# hamming distance and clustering
from scipy.spatial.distance import pdist, squareform
import plotly.figure_factory as ff


# target_gens = vf_genes
def draw_heatmap(target_gens, filter_sample=True, use_core=False, dist=None):
    sub_df = total_df.loc[:, target_gens]
    sub_df = sub_df.applymap(lambda x: abs(int(pd.isna(x)) - 1))
    if filter_sample:
        sub_df = sub_df.loc[map(lambda x: bool(re.findall('^[0-9]*$', x)), sub_df.index), :]
    dist = pd.DataFrame(squareform(pdist(sub_df)), index=sub_df.index, columns=sub_df.index)
    sub_df = sub_df.loc[:, sub_df.sum(0) != 0]
    col_dist = pd.DataFrame(squareform(pdist(sub_df.T)), index=sub_df.columns, columns=sub_df.columns)

    if use_core:
        t = Tree(open('/home/liaoth/data2/project/shenzhen_actineto/roary_o/core_gene.newick').read())
        n_dict = {}
        for n in t.traverse():
            if n.is_leaf():
                n_dict[n.name] = n

        pairwise_matrix = defaultdict(dict)
        for t1, t2 in itertools.product(total_df.index, total_df.index):
            pairwise_matrix[t1][t2] = n_dict[t1].get_distance(n_dict[t2])
        dist = pd.DataFrame(pairwise_matrix)
        dist = dist.loc[sub_df.index, sub_df.index]
        dist = pd.DataFrame(dist, index=dist.index, columns=dist.index)
    # prepare dist and data
    ############################################################
    # create fig
    fig = plotly.tools.make_subplots(2, 3,
                                     specs=[[None, None, {}],
                                            [{}, {}, {}]],
                                     shared_yaxes=True,
                                     shared_xaxes=True)
    # This is the format of your plot grid:
    #     (empty)          (empty)      [ (1,3) x3,y1 ]
    # [ (2,1) x1,y2 ]  [ (2,2) x2,y2 ]  [ (2,3) x3,y2 ]

    # create up side dendrogram
    up_dendro = ff.create_dendrogram(col_dist.values,
                                     orientation='bottom',
                                     labels=col_dist.index)
    sub_df = sub_df.reindex(columns=up_dendro.layout.xaxis.ticktext)
    for i in up_dendro.data:
        i.showlegend = False
        fig.append_trace(i, 1, 3)

    # Create left Side Dendrogram
    side_dendro = ff.create_dendrogram(dist.values,
                                       orientation='right',
                                       labels=dist.index)
    for i in side_dendro.data:
        i.showlegend = False
        fig.append_trace(i, 2, 1)
    sub_df = sub_df.reindex(side_dendro.layout.yaxis.ticktext)

    ############################################################
    # draw two heatmap
    heatmap = go.Heatmap(
        x=sub_df.index,
        y=sub_df.columns,
        z=sub_df.values,
        colorscale='Earth',
        reversescale=True,
        showscale=False
    )
    heatmap['x'] = up_dendro.layout.xaxis.tickvals
    heatmap['y'] = side_dendro.layout.yaxis.tickvals
    fig.append_trace(heatmap, 2, 3)

    ############################################################
    _sub_df = total_df.reindex(side_dendro.layout.yaxis.ticktext).loc[:, ["Sp_krkn_FinalCall",
                                                                          "Oxf_ST.1",
                                                                          "Oxf_ST.2",
                                                                          "Pas_ST.1",
                                                                          "year",
                                                                          "source",
                                                                          "gender",
                                                                          "科室",
                                                                          "RIS_IPM",
                                                                          "RIS_MEM"]]
    new_df = _sub_df.copy()
    for idx, col in _sub_df.iteritems():
        le = LabelEncoder()
        col = col.fillna('0')
        col = col.astype(str)
        if (col=='0').any():
            new_df.loc[:, idx] = le.fit_transform(col)
        else:
            new_df.loc[:, idx] = le.fit_transform(col) + 1
    _sub_df = med_data.reindex(side_dendro.layout.yaxis.ticktext).loc[:, ['RES_AMK', 'RES_AMP', 'RES_CAZ', 'RES_CIP', 'RES_CRO', 'RES_CSL',
                                                                           'RES_FEP', 'RES_GEN', 'RES_IPM', 'RES_LVX', 'RES_MEM', 'RES_MNO',
                                                                           'RES_PIP', 'RES_SAM', 'RES_SXT', 'RES_TGC', 'RES_TZP']]
    new_df = _sub_df.applymap(lambda x: {'I': 2, 'R': 3, 'S': 1, 'no': 0}.get(x))

    heatmap = go.Heatmap(
        x=list(_sub_df.columns),
        y=side_dendro['layout']['yaxis']['tickvals'],
        z=new_df.values,
        text=_sub_df.values,
        hoverinfo='all',
        colorscale='Earth',
        reversescale=True,
        showscale=False
    )
    fig.append_trace(heatmap, 2, 2)
    ############################################################
    fig.layout.xaxis1.domain = [0, 0.05]
    fig.layout.xaxis2.domain = [0.05, 0.1]
    fig.layout.xaxis3.domain = [0.1, 1]
    fig.layout.yaxis1.domain = [0.9, 1]
    fig.layout.yaxis2.domain = [0.0, 0.9]
    fig.layout.yaxis2.ticktext = side_dendro.layout.yaxis.ticktext
    fig.layout.yaxis2.tickvals = side_dendro.layout.yaxis.tickvals
    fig.layout.xaxis3.ticktext = up_dendro.layout.xaxis.ticktext
    fig.layout.xaxis3.tickvals = up_dendro.layout.xaxis.tickvals
    # fig.layout.xaxis.anchor = 'x2'
    fig.layout.margin.b = 250
    fig.layout.width = len(target_gens) * 23
    fig.layout.height = 1000
    fig.layout.xaxis3.tickangle = 30
    fig.layout.xaxis1.showgrid = False
    fig.layout.yaxis2.showgrid = False
    fig.layout.xaxis3.showgrid = False
    fig.layout.yaxis1.showgrid = False
    fig.layout.xaxis1.showticklabels = False
    fig.layout.yaxis1.showticklabels = False
    return fig


vf_fig = draw_heatmap(vf_genes)
res_fig = draw_heatmap(res_genes)

vf_coreT_fig = draw_heatmap(vf_genes, use_core=True)
res_coreT_fig = draw_heatmap(res_genes, use_core=True)

vf_withall_fig = draw_heatmap(vf_genes, filter_sample=False)
res_withall_fig = draw_heatmap(res_genes, filter_sample=False)

vf_coreT_withall_fig = draw_heatmap(vf_genes, use_core=True, filter_sample=False)
res_coreT_withall_fig = draw_heatmap(res_genes, use_core=True, filter_sample=False)
import os
import plotly.io as pio
base_dir = "/home/liaoth/data2/project/shenzhen_actineto/汇总结果/heatmap2"
os.makedirs(base_dir, exist_ok=True)
for t in ["vf_{tag}fig", "res_{tag}fig"]:
    for tag in ['', 'coreT_', 'withall_', 'coreT_withall_']:
        exec("plotly.offline.plot(%s,filename=os.path.join(base_dir,'%s.html'),auto_open=False)" % (t.format(tag=tag),
                                                                                                                           t.format(tag=tag).replace('_fig', '')))
        exec("#pio.write_image(%s,file=os.path.join(base_dir,'%s.pdf'),format='pdf')" % (t.format(tag=tag),
                                                                                        t.format(tag=tag).replace('_fig', '')))
