# -*- coding: UTF-8 -*-
import itertools
import os
import re
import sys
from collections import defaultdict

import pandas as pd
from ete3 import Tree
from sklearn.preprocessing import LabelEncoder
############################################################
import plotly
from project_specific import params
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

############################################################
pandoo_o_dir = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pandoo_result"
merged_df_path = os.path.join(pandoo_o_dir, 'For phandango.csv')
total_df = pd.read_csv(merged_df_path, index_col=0)

res_name_convert = pd.read_csv('/home/liaoth/data2/project/shenzhen_Acinetobacter/gene_name_convert.csv', index_col=0)
res_name_convert = res_name_convert.to_dict()['0']

cols = [_ for _ in total_df.columns if _.startswith('abricate')]

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
print("original remove name redundant, left:", len(single_gene_cols))
res_db = "_card|ncbi|resinder|argannot_"
vf_db = "_vfdb_full|vfdb|victors_"
for i in [_ for _ in single_gene_cols if re.search(vf_db, _)]:
    tmp = '_'.join(i.split('_')[2:])
    if tmp.startswith('full_'):
        tmp = tmp.split('full_')[1]
    res_name_convert[i] = tmp

res_genes = [res_name_convert.get(_, _) for _ in single_gene_cols if re.search(res_db, _)]
vf_genes = [res_name_convert.get(_, _) for _ in single_gene_cols if re.search(vf_db, _)]
print("original remove name redundant, left:", len(res_genes), len(vf_genes))
total_df.columns = [res_name_convert.get(_, _) for _ in total_df.columns]
# remove gene which all samples contain this gene & gene which all samples doesn't have this gene.
res_genes = [_ for _ in res_genes if pd.isnull(total_df.loc[:, _]).any() and not pd.isnull(total_df.loc[:, _]).all()]
vf_genes = [_ for _ in vf_genes if pd.isnull(total_df.loc[:, _]).any() and not pd.isnull(total_df.loc[:, _]).all()]
print("original remove all have/none, left:", len(res_genes), len(vf_genes))
############################################################
# target_gens = vf_genes
# eu distance and clustering
from scipy.spatial.distance import pdist, squareform
from vis.heatmap_api import create_heatmap


def process_df(total_df, target_gens, filter_sample=True, map_type='normal'):
    sub_df = total_df.loc[:, target_gens]

    if filter_sample:
        sub_df = sub_df.loc[map(lambda x: bool(re.findall('^[0-9]*$', x)), sub_df.index), :]
        sub_df = sub_df.loc[:, (pd.isna(sub_df).any(0)) & (~pd.isna(sub_df).all(0))]

    if map_type == "normal":
        sub_df = sub_df.applymap(lambda x: abs(int(pd.isna(x)) - 1))
        sub_df_text = sub_df.copy()
    elif "__call__" in dir(map_type):
        sub_df = sub_df.T
        sub_df = sub_df.apply(lambda row:[sub_df.index[idx] if not pd.isna(v) else "NaN" for idx,v in enumerate(row)] ,axis=0)
        sub_df = sub_df.T
        sub_df = sub_df.applymap(lambda x: map_type(x))
        sub_df_text = sub_df.copy()
        empty_list = sorted(set(sub_df.values.ravel()))
        le_dict = dict(zip(empty_list, range(1, len(empty_list) + 1)))
        le_dict.update({'nan':0})
        sub_df = sub_df.applymap(lambda x: le_dict[x])

    sub_df = sub_df.loc[:, sub_df.sum(0) != 0]
    return sub_df, sub_df_text


def main(total_df, target_gens, filter_sample=True, tree_p=None, accessory_cols=None, map_type="normal", up_cluster=None, **kwargs):
    main_matrix, main_matrix_text = process_df(total_df, target_gens, filter_sample, map_type=map_type)


    s_df = main_matrix.T
    if up_cluster is not None:
        up_dist = pd.DataFrame(squareform(pdist(s_df.values)), index=s_df.index, columns=s_df.index)
    if tree_p is not None:
        t = Tree(open(tree_p).read())
        n_dict = {}
        for n in t.traverse():
            if n.is_leaf():
                n_dict[n.name] = n

        pairwise_matrix = defaultdict(dict)
        for t1, t2 in itertools.product(main_matrix.index, main_matrix.index):
            pairwise_matrix[t1][t2] = n_dict[t1].get_distance(n_dict[t2])
        dist = pd.DataFrame(pairwise_matrix)
        dist = dist.loc[main_matrix.index, main_matrix.index]
        left_dist = pd.DataFrame(dist, index=dist.index, columns=dist.index)
    else:
        left_dist = pd.DataFrame(squareform(pdist(main_matrix.values)), index=main_matrix.index, columns=main_matrix.index)
    if accessory_cols is not None:
        accessory_matrix = total_df.copy()
        accessory_matrix = accessory_matrix.loc[:, accessory_cols]
        accessory_matrix_text = accessory_matrix.loc[:, accessory_cols]
        accessory_matrix = accessory_matrix.fillna('unknown').astype(str)
        le = LabelEncoder()
        empty_list = []
        empty_list += [i for idx,_ in accessory_matrix.iteritems() for i in set(_)]
        le.fit(np.array(empty_list).reshape(-1, 1))
        accessory_matrix = accessory_matrix.applymap(lambda x: le.transform([[x]])[0])
        fig = create_heatmap(left_dist=left_dist,
                             main_matrix=main_matrix,
                             main_matrix_text=main_matrix_text,
                             accessory_matrix=accessory_matrix,
                             accessory_matrix_text=accessory_matrix_text,
                             up_dist=up_dist if up_cluster is not None else None,
                             **kwargs)
    else:
        fig = create_heatmap(left_dist=left_dist,
                             main_matrix=main_matrix,
                             up_dist=up_dist if up_cluster is not None else None,
                             **kwargs)
    return fig
if __name__ == '__main__':



    vf_fig = main(total_df,
                  vf_genes,
                  filter_sample=True,
                  accessory_cols=params.vf_cols,
                  map_type=lambda x: params.vf2fun.get(x, "unknown"),
                  width=2000, height=1000)
    plotly.offline.plot(vf_fig,
                        filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pandoo_result/VF_with_genes.html",
                        auto_open=False)

    res_fig = main(total_df,
                   res_genes,
                   filter_sample=True,
                   accessory_cols=params.res_cols,
                   map_type=lambda x: params.res2fun.get(x, "unknown"),
                   width=2000, height=1000)
    plotly.offline.plot(res_fig,
                        filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pandoo_result/Res_with_genes.html",
                        auto_open=False)
    ############################################################
    full_fig = main(total_df,
                  total_df.columns,
                  filter_sample=False,
                  accessory_cols=params.res_cols,
                  map_type='normal',
                  width=2500, height=1000)
    plotly.offline.plot(full_fig,
                        filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/whole_region_roary/all_with_metadata.html",
                        auto_open=False)


    ############################################################
    # KL specific
    KL_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/fkpA_1-lldP_region/roary_o/gene_presence_absence.csv",
                        sep=',', index_col=0)
    KL_df = KL_df.iloc[:, 13:].T
    merged_df = pd.concat([KL_df, total_df.loc[:, params.vf_cols]], axis=1)

    fig = main(merged_df,
               KL_df.columns,
               filter_sample=True,
               tree_p='/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/fkpA_1-lldP_region/roary_o/accessory_binary_genes.fa.newick',
               accessory_cols=params.vf_cols,
               width=2000, height=1000)
    plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/heatmap_with_metadata.html",
                        auto_open=False)
