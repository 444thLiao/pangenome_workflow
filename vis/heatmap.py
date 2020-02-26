# -*- coding: UTF-8 -*-

import numpy as np
import pandas as pd
############################################################
import plotly
from sklearn.preprocessing import LabelEncoder

from project_specific import params
from vis.draw_dendrogram import main as c_dendrogram

pandoo_df = pd.read_csv('/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pandoo_result/For phandango.csv',
                        sep=',', index_col=0)
############################################################
abricate_data_path = "/home/liaoth/data2/project/shenzhen_Acinetobacter/21_outgroup_roary_o/abricate_result.tab"
total_df = pd.read_csv(abricate_data_path, sep='\t', index_col=0)
total_df = total_df.fillna(0)

res_db = ["card", "ncbi", "resfinder", "argannot"]
vf_db = ["vfdb", "victors"]
res_df = total_df.loc[total_df.index[:-1], total_df.loc['db', :].isin(res_db)].astype(float).astype(int)
vf_df = total_df.loc[total_df.index[:-1], total_df.loc['db', :].isin(vf_db)].astype(float).astype(int)

res_genes = list(res_df.columns)
vf_genes = list(vf_df.columns)
total_df = pd.concat([res_df, vf_df, pandoo_df], axis=1)
# drop_duplicated = ['abricate_ncbi_mph(E)',
#                    'abricate_ncbi_msr(E)',
#                    'abricate_ncbi_aac(3)-I',
#                    'abricate_argannot_(Sul)SulI',
#                    'abricate_argannot_(Sul)SulII',
#                    'abricate_argannot_(Tet)Tet-39',
#                    'abricate_argannot_(Tet)TetB',
#                    'abricate_argannot_(AGly)Aac3-I',
#                    "abricate_argannot_(AGly)Aph3''Ia"]
#
# unique_genes = list(set([_ for v in all_gene.values() for _ in v]))
#
#
# print("original remove name redundant, left:", len(unique_genes))
#
# res_genes = list(set([_ for _ in unique_genes if  _ in res_genes]))
# vf_genes = list(set([_ for _ in unique_genes if _ in vf_genes]))
# print("original remove name redundant, left:", len(res_genes), len(vf_genes))
#
# # remove gene which all samples contain this gene & gene which all samples doesn't have this gene.
# res_genes = [_ for _ in res_genes if pd.isnull(interested_df.loc[:, _]).any() and not pd.isnull(interested_df.loc[:, _]).all()]
# vf_genes = [_ for _ in vf_genes if pd.isnull(interested_df.loc[:, _]).any() and not pd.isnull(interested_df.loc[:, _]).all()]
# print("original remove all have/none, left:", len(res_genes), len(vf_genes))
############################################################
# target_gens = vf_genes
# eu distance and clustering
from scipy.spatial.distance import pdist, squareform
from vis.heatmap_api import create_heatmap


def process_df_for_str(total_df, target_gens, filter_fea_all=True, filter_sample=True, map_type='normal'):
    """ For dataframe contains str/float/obj, like metadata."""
    sub_df = total_df.loc[:, target_gens]

    if type(filter_sample) == list:
        sub_df = sub_df.loc[filter_sample, :]
        if filter_fea_all:
            sub_df = sub_df.loc[:, (pd.isnull(sub_df).any(0)) & (~ pd.isnull(sub_df).all(0))]
        else:
            sub_df = sub_df.loc[:, (~ pd.isnull(sub_df).all(0))]

    if map_type == "normal":
        sub_df = sub_df.fillna('NaN')
        set_values = set(sub_df.values.ravel())
        set_num = range(len(set_values))
        set_dict = dict(zip(set_values, set_num))
        sub_df_text = sub_df.copy()
        sub_df = sub_df.applymap(lambda x: set_dict[x])

    elif "__call__" in dir(map_type):
        sub_df = sub_df.T
        sub_df = sub_df.apply(lambda row: [sub_df.index[idx] if v > 0 else "NaN" for idx, v in enumerate(row)], axis=0)
        sub_df = sub_df.T
        sub_df_text = sub_df.applymap(lambda x: map_type(x) + '<Br>' + x)
        sub_df = sub_df.applymap(lambda x: map_type(x))
        empty_list = sorted(set(sub_df.values.ravel()))
        le_dict = dict(zip(empty_list, range(1, len(empty_list) + 1)))
        le_dict.update({'nan': 0})
        sub_df = sub_df.applymap(lambda x: le_dict[x])

    sub_df = sub_df.loc[:, sub_df.sum(0) != 0]
    return sub_df, sub_df_text


def process_df(total_df, target_gens, filter_fea_all=True, filter_sample=True, map_type='normal'):
    """ For 1/0 dataframe create by roary. like gene_presence_absence.Rtab (!!!instead of csv)"""
    sub_df = total_df.loc[:, target_gens]

    if type(filter_sample) == list:
        sub_df = sub_df.loc[filter_sample, :]
        if filter_fea_all:
            sub_df = sub_df.loc[:, ((sub_df.sum(0) != sub_df.shape[0]) & (sub_df.sum(0) != 0))]
        else:
            sub_df = sub_df.loc[:, (sub_df.sum(0) != 0)]
    if map_type == "normal":
        sub_df = sub_df.applymap(lambda x: abs(int(pd.isna(x)) - 1))
        sub_df_text = sub_df.copy()

    elif "__call__" in dir(map_type):
        sub_df = sub_df.T
        sub_df = sub_df.apply(lambda row: [sub_df.index[idx] if v > 0 else "NaN" for idx, v in enumerate(row)], axis=0)
        sub_df = sub_df.T
        sub_df_text = sub_df.applymap(lambda x: map_type(x) + '<Br>' + x)
        sub_df = sub_df.applymap(lambda x: map_type(x))
        empty_list = sorted(set(sub_df.values.ravel()))
        le_dict = dict(zip(empty_list,
                           range(1, len(empty_list) + 1)))
        le_dict.update({'nan': 0})
        sub_df = sub_df.applymap(lambda x: le_dict[x])

    sub_df = sub_df.loc[:, sub_df.sum(0) != 0]
    return sub_df, sub_df_text


def main(total_df, target_gens, filter_fea_all=True, filter_sample=True, tree_p=None, accessory_cols=None, map_type="normal", up_cluster=None, rooted='midpoint', **kwargs):
    main_matrix, main_matrix_text = process_df(total_df,
                                               target_gens,
                                               filter_sample=filter_sample,
                                               map_type=map_type,
                                               filter_fea_all=filter_fea_all)

    s_df = main_matrix.T
    if up_cluster is not None:
        up_dist = pd.DataFrame(squareform(pdist(s_df.values)), index=s_df.index, columns=s_df.index)
    if tree_p is not None:
        left_dist = c_dendrogram(tree_p, rooted=rooted)
    else:
        left_dist = pd.DataFrame(squareform(pdist(main_matrix.values)), index=main_matrix.index, columns=main_matrix.index)
    if accessory_cols is not None:
        accessory_matrix = total_df.copy()
        accessory_matrix = accessory_matrix.loc[:, accessory_cols]
        accessory_matrix_text = accessory_matrix.loc[:, accessory_cols]
        accessory_matrix = accessory_matrix.fillna('unknown').astype(str)
        le = LabelEncoder()
        empty_list = []
        empty_list += [i for idx, _ in accessory_matrix.iteritems() for i in set(_)]
        le.fit(np.array(empty_list).reshape(-1, 1))
        accessory_matrix = accessory_matrix.applymap(lambda x: le.transform([[x]])[0])
        fig, main_matrix, accessory_matrix = create_heatmap(left_dist=left_dist,
                                                            main_matrix=main_matrix,
                                                            main_matrix_text=main_matrix_text,
                                                            accessory_matrix=accessory_matrix,
                                                            accessory_matrix_text=accessory_matrix_text,
                                                            up_dist=up_dist if up_cluster is not None else None,
                                                            return_matrix=True,
                                                            **kwargs)
    else:
        fig, main_matrix, accessory_matrix = create_heatmap(left_dist=left_dist,
                                                            main_matrix=main_matrix,
                                                            main_matrix_text=main_matrix_text,
                                                            up_dist=up_dist if up_cluster is not None else None,
                                                            return_matrix=True,
                                                            **kwargs)
    return fig, main_matrix, accessory_matrix


samples20 = ['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
             '39502', '39528', '40067', '40116', '40615', '41194', '41296', '41833',
             '42121', '42268', '42349', '43458']

samples21 = ['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
             '39502', '39528', '40067', '40116', '40615', '41194', '41296', '41833',
             '42121', '42268', '42349', '43458', '40283']

samples64 = ['34960', '35082', '35989', '37502', '37706', '39054', '39058', '39292',
             '39502', '39528', '40067', '40116', '40283', '40615', '41194', '41296',
             '41833', '42121', '42268', '42349', '43458', 'AB1', 'AB10', 'AB11',
             'AB12', 'AB13', 'AB14', 'AB15', 'AB16', 'AB17', 'AB18', 'AB19', 'AB2',
             'AB20', 'AB21', 'AB22', 'AB23', 'AB24', 'AB25', 'AB26', 'AB27', 'AB28',
             'AB29', 'AB3', 'AB30', 'AB31', 'AB32', 'AB33', 'AB34', 'AB35', 'AB36',
             'AB37', 'AB38', 'AB39', 'AB4', 'AB40', 'AB41', 'AB42', 'AB5', 'AB6',
             'AB7', 'AB8', 'AB9', 'XH860']
if __name__ == '__main__':
    import os

    odir = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pandoo_result"
    for samples in [samples20, samples21, samples64]:
        if len(samples) == 64:
            tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/%s_roary_o/core_gene.newick" % len(samples)
            rooted = 'midpoint'
        elif len(samples) == 21:
            tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/21_outgroup_roary_o/core_gene.newick"
            rooted = 'AB030'
        elif len(samples) == 20:
            tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/20_outgroup_roary_o/core_gene.newick"
            rooted = 'AB030'

        for mode in ['VF', 'RES']:
            if mode == 'VF':
                genes = vf_genes
                accessory_cols = params.vf_cols
                map_type = lambda x: params.vf2fun.get(x, "unknown")
            elif mode == 'RES':
                genes = res_genes
                accessory_cols = params.res_cols
                map_type = lambda x: params.res2fun.get(x, "unknown")

            for f_all in [True, False]:
                if f_all:
                    mmatrix_name = os.path.join(odir, "%s_%s_based_SNP.csv" % (len(samples), mode))
                    fig_name = os.path.join(odir, "%s_%s_based_SNP.html" % (len(samples), mode))
                else:
                    mmatrix_name = os.path.join(odir, "%s_%s_FULL_based_SNP.csv" % (len(samples), mode))
                    fig_name = os.path.join(odir, "%s_%s_FULL_based_SNP.html" % (len(samples), mode))

                fig, mmatrix, amatrix = main(total_df,
                                             genes,
                                             filter_sample=samples,
                                             filter_fea_all=f_all,
                                             tree_p=tree_p,
                                             accessory_cols=accessory_cols,
                                             map_type=map_type,
                                             rooted=rooted,
                                             up_cluster=True,
                                             width=2000, height=1000)
                mmatrix.to_csv(mmatrix_name, index=True)

                plotly.offline.plot(fig,
                                    filename=fig_name,
                                    auto_open=False)

    ############################################################
