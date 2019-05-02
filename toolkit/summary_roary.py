import itertools
import os
from collections import defaultdict

import pandas as pd
import plotly.io as pio
from ete3 import Tree
from scipy.spatial.distance import pdist, squareform

import sys
sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
from vis.heatmap_api import create_heatmap


def main(indir, odir, tree=None):
    # tree for samples
    core_file = os.path.join(indir, "gene_presence_absence.csv")
    data = pd.read_csv(core_file, index_col=0, sep=',',low_memory=False)
    sub_data = data.iloc[:, 13:]
    sparse_matrix_df = sub_data.applymap(lambda x: 0 if pd.isna(x) else 1).T
    sample_dist = pd.DataFrame(squareform(pdist(sparse_matrix_df)),
                               index=sparse_matrix_df.index,
                               columns=sparse_matrix_df.index)
    if tree:
        tree_obj = Tree(open(tree).read())
        n_dict = {}
        for n in tree_obj.traverse():
            if n.is_leaf():
                n_dict[n.name] = n

        pairwise_matrix = defaultdict(dict)
        for t1, t2 in itertools.product(sample_dist.index, sample_dist.index):
            pairwise_matrix[t1][t2] = n_dict[t1].get_distance(n_dict[t2])
        dist = pd.DataFrame(pairwise_matrix)
        dist = dist.loc[sample_dist.index, sample_dist.index]
        dist = pd.DataFrame(dist, index=dist.index, columns=dist.index)
        fig = create_heatmap(left_dist=dist,
                             # up_dist=gene_dist,
                             main_matrix=sparse_matrix_df,
                             width=1700,
                             height=1500)
    else:
        fig = create_heatmap(left_dist=sample_dist,
                             # up_dist=gene_dist,
                             main_matrix=sparse_matrix_df,
                             width=1700,
                             height=1500)

    new_data = sub_data.apply(lambda x: [sub_data.index[idx] if not pd.isna(_) else '-' for idx, _ in enumerate(x)])
    new_data = new_data.reindex(columns=fig.layout.yaxis2.ticktext)
    gene_info = data.loc[:, "Annotation"]
    new_data = pd.concat([new_data,gene_info],axis=1)
    new_data = new_data.T
    new_data.to_csv(os.path.join(odir,"heatmap_data.csv"))
    #todo: implement a more usefull output instead of a useless `heatmap_data.csv`
    heatmap_path = os.path.join(odir,"heatmap_with_tree.png")
    pio.write_image(fig, file=heatmap_path)



if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    indir = os.path.abspath(args.indir)
    odir = os.path.abspath(args.outdir)
    os.makedirs(odir,exist_ok=True)
    main(indir,odir )
