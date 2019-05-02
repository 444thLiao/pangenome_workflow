import os

import pandas as pd
from sklearn.decomposition import PCA
import re

data = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/roary_o/gene_presence_absence.csv",
                   sep=',',
                   index_col=0, low_memory=False)

data = data.iloc[:, 13:]
# data = data.apply(lambda x: [data.index[idx] if not pd.isna(v) else 0 for idx,v in enumerate(x)],axis=0)
data = data.apply(lambda x: [1 if not pd.isna(v) else 0 for idx, v in enumerate(x)], axis=0)
data = data.T

############################################################
pandoo_o_dir = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pandoo_result"
merged_df_path = os.path.join(pandoo_o_dir, 'For phandango.csv')
metadata = pd.read_csv(merged_df_path, index_col=0, low_memory=False)
metadata = metadata.reindex(data.index)
############################################################
import plotly
import plotly.graph_objs as go
import plotly.io as pio


def draw_PCA(data, col=None):
    pca = PCA()
    pca_r = pca.fit_transform(data)
    fig = go.Figure()
    if col is not None:
        for val in set(col):
            fig.add_scatter(x=pca_r[col == val, 0],
                            y=pca_r[col == val, 1],
                            text=data.index[col == val],
                            marker=dict(size=20, opacity=0.6),
                            name=val,
                            mode="markers", )
    else:
        fig.add_scatter(x=pca_r[:, 0],
                        y=pca_r[:, 1],
                        text=data.index,
                        marker=dict(size=20, opacity=0.7),
                        mode="markers", )
    fig.layout.hovermode = "closest"
    fig.layout.xaxis.title = "PC1({:.2f}%)".format(pca.explained_variance_ratio_[0] * 100)
    fig.layout.yaxis.title = "PC2({:.2f}%)".format(pca.explained_variance_ratio_[1] * 100)
    fig.layout.font.size = 15
    return fig


def extract_samples(data, metadata):
    subdata = data.loc[map(lambda x: bool(re.findall('^[0-9]*$', x)), data.index), :]
    subdata = subdata.loc[:, data.sum(0) != 0]
    submetadata = metadata.reindex(subdata.index)
    return subdata, submetadata


def iter_output_PCA(data, metadata, odir, label):
    os.makedirs(os.path.join(odir, label), exist_ok=True)
    for fea in ["Sp_krkn_FinalCall", 'source', '科室', 'year', 'Oxf_ST.1'] + list([_ for _ in metadata.columns if _.startswith('RES_')]):
        fig = draw_PCA(data, metadata.loc[:, fea])
        pio.write_image(fig, file=os.path.join(odir, label, "%s.png" % fea))


odir = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/pca"

iter_output_PCA(data, metadata, odir, label='all_samples_based_all_genes')
############################################################
subdata, submetadata = extract_samples(data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='21_samples_based_all_genes')
############################################################
# input other data
from vis.heatmap import process_df, total_df, vf_genes, res_genes

vf_data, vf_data_text = process_df(total_df, vf_genes, filter_sample=True)
subdata, submetadata = extract_samples(vf_data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='21_samples_based_VF_genes')

res_data, res_data_text = process_df(total_df, res_genes, filter_sample=True)
subdata, submetadata = extract_samples(res_data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='21_samples_based_RES_genes')

vf_data, vf_data_text = process_df(total_df, vf_genes, filter_sample=False)
subdata, submetadata = extract_samples(vf_data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='all_samples_based_VF_genes')

res_data, res_data_text = process_df(total_df, res_genes, filter_sample=False)
subdata, submetadata = extract_samples(res_data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='all_samples_based_RES_genes')
############################################################
# For merge OXA-like
from project_specific.params import OXAdict
from collections import defaultdict
res_data, res_data_text = process_df(total_df, res_genes, filter_sample=False)
rev_OXAdict = defaultdict(list)
for k,v in OXAdict.items():
    rev_OXAdict[v].append(k)
for merged_OXA,before_oxa in rev_OXAdict.items():
    _sub_data = res_data.loc[:,before_oxa]
    _sub_data = _sub_data.any(1).astype(int)
    res_data.drop(before_oxa,axis=1,inplace=True)
    res_data.loc[:,merged_OXA] = _sub_data
subdata, submetadata = extract_samples(res_data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='all_samples_based_merged_OXA_RES_genes')
res_data, res_data_text = process_df(total_df, res_genes, filter_sample=True)
rev_OXAdict = defaultdict(list)
for k,v in OXAdict.items():
    rev_OXAdict[v].append(k)
for merged_OXA,before_oxa in rev_OXAdict.items():
    before_oxa = [_ for _ in before_oxa if _ in res_data.columns]
    _sub_data = res_data.loc[:,before_oxa]
    _sub_data = _sub_data.any(1).astype(int)
    res_data.drop(before_oxa,axis=1,inplace=True)
    res_data.loc[:,merged_OXA] = _sub_data
subdata, submetadata = extract_samples(res_data, metadata)
iter_output_PCA(subdata, submetadata, odir, label='21_samples_based_merged_OXA_RES_genes')
############################################################
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from scipy.spatial.distance import squareform,pdist
import numpy as np
importr("vegan")
dist_m = pd.DataFrame(squareform(pdist(data.values)),index=data.index,columns=data.index)
dist_m.fillna(0).to_csv('/tmp/test_dist.csv',index=True)
vf_data, vf_data_text = process_df(total_df, vf_genes, filter_sample=False)
vf_data.loc[:,vf_genes].fillna(0).to_csv('/tmp/test_metadata.csv',index=True)
def envfit_metadata( metadata_path, dist_path, n_iter=500, return_ord=False):
    rcode = """
    metadata <- read.csv('{path_metadata}',row.names = 1,check.names=FALSE)
    dist <- read.csv('{path_dist}',row.names = 1,check.names=FALSE)
    dist <- as.dist(dist)
    ord <- capscale(dist ~ -1)
    """.format(
               path_dist=dist_path,
               path_metadata=metadata_path)
    robjects.r(rcode)

    envfit_result = robjects.r(
        """
        fit <- envfit(ord,metadata,permutations = {n_iter})
        fit$vectors
        """.format(n_iter=n_iter))
    R_ord = robjects.r('summary(ord)$sites')
    R_pro_X_df = pd.DataFrame(data=np.array(R_ord), index=R_ord.rownames, columns=R_ord.colnames)
    fit_result = pd.DataFrame(columns=["r2", "pvals", "Source", "End"]
                              , index=envfit_result[envfit_result.names.index("arrows")].rownames)
    fit_result.loc[:, "r2"] = envfit_result[envfit_result.names.index("r")]
    fit_result.loc[:, "pvals"] = envfit_result[envfit_result.names.index("pvals")]
    fit_result.loc[:, ["Source", "End"]] = np.array(envfit_result[envfit_result.names.index("arrows")])
    if return_ord:
        return fit_result, R_pro_X_df
    else:
        return fit_result
fit_result, R_pro_X_df = envfit_metadata('/tmp/test_metadata.csv','/tmp/test_dist.csv',return_ord=True)
from vis.envfit_plot import draw_PCA_with_envfit_arrow
fig = draw_PCA_with_envfit_arrow(R_pro_X_df,fit_result)
plotly.offline.plot(fig,filename="vis/test.html",auto_open=False)



