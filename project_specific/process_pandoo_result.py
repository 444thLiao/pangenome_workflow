# -*- coding: UTF-8 -*-
import os
from glob import glob

import pandas as pd
import sys
sys.path.insert(0,os.path.dirname(os.path.dirname(__file__)))
from process_mlst import main as process_mlst


def get_and_remove(data, match_part, keep_name):
    cols = [_ for _ in data.columns if match_part in _ and keep_name not in _]
    data = data.drop(cols, axis=1)
    return data


def filter_out(data):
    dropped_col = []
    for idx, col in data.iteritems():
        tmp = col
        unique_vals = set([_ for _ in col.values if not pd.isna(_)])
        for val in list(unique_vals):
            if val and str(val).startswith('COV'):
                coverage = float(str(val).split('_')[1])
                identity = float(str(val).split('_')[-1])
                if coverage <= 90:
                    unique_vals.remove(val)
                    tmp[tmp == val] = ''
                    data.loc[:, idx] = tmp
                if identity <= 60:
                    unique_vals.remove(val)
                    tmp[tmp == val] = ''
                    data.loc[:, idx] = tmp
        if len(unique_vals) == 0:
            dropped_col.append(idx)
    data = data.drop(dropped_col, axis=1)
    data = data.loc[:, ~data.isna().all(0)]
    return data


def merge_for_phandango(analysis_dir):
    final_df = []
    basic_info = os.path.join(analysis_dir, 'basic_info.csv')
    data = pd.read_csv(basic_info, index_col=0)
    share_index = list(data.index)
    final_df.append(data.loc[:, ['Sp_krkn_FinalCall', 'metricsContigs_N50']])
    ############################################################
    output_mlst = glob(os.path.join(analysis_dir, 'mlst', '*_mlst.csv'))

    for mlst_file in output_mlst:
        data = pd.read_csv(mlst_file, index_col=0)
        final_df.append(data.loc[:, [_ for _ in data.columns if '_ST' in _]])
    ############################################################
    basic_info = os.path.join(analysis_dir, 'metadata.csv')
    data = pd.read_csv(basic_info, index_col=0)
    data = data.loc[share_index, :]
    final_df.append(data)
    ############################################################
    abricate_info = os.path.join(analysis_dir, 'abricate.csv')
    data = pd.read_csv(abricate_info, index_col=0)
    data = data.loc[share_index, :]
    data = filter_out(data)
    final_df.append(data)
    ############################################################
    # ariba_info = os.path.join(analysis_dir, 'ariba.csv')
    # data = pd.read_csv(ariba_info, index_col=0)
    # data = data.loc[share_index,:]
    # data = filter_out(data)
    # final_df.append(data)
    ############################################################
    final_df = pd.concat(final_df, axis=1)
    return final_df


def main(indir, outdir):
    output_mlst = os.path.join(outdir, 'mlst', '{scheme}_mlst.csv')
    os.makedirs(os.path.dirname(output_mlst), exist_ok=True)
    ############################################################

    # mlst part
    mlst_files = glob(os.path.join(indir, '*/mlst.txt'))
    df_list = []
    for mlst_file in mlst_files:
        df_list.append(pd.read_csv(mlst_file, sep=',', index_col=0))
    mlst_df = pd.concat(df_list, axis=0)
    process_mlst(mlst_df, output_mlst)
    ############################################################
    summary_file = os.path.join(indir, 'isolates_metadataAll.csv')
    data = pd.read_csv(summary_file, index_col=0)
    # remove version column and mlst
    data = data.loc[:, [_ for _ in data.columns if 'version' not in _.lower()]]
    data = data.drop([_ for _ in data.columns if 'mlst' in _.lower()], axis=1)
    # summarize kraken
    data = get_and_remove(data, "Sp_krkn_", "krkn_FinalCall")
    # summarize abricate
    abricate_file = os.path.join(outdir, 'abricate.csv')
    abricate_part = [_ for _ in data.columns if 'abricate' in _]
    abricate_data = data.loc[:, abricate_part]
    data = data.drop(abricate_part, axis=1)
    abricate_data.to_csv(abricate_file, index=True)
    # summarize abricate
    ariba_file = os.path.join(outdir, 'ariba.csv')
    ariba_part = [_ for _ in data.columns if 'ariba' in _]
    ariba_data = data.loc[:, ariba_part]
    data = data.drop(ariba_part, axis=1)
    ariba_data.to_csv(ariba_file, index=True)
    # remove path
    basic_info = os.path.join(outdir, 'basic_info.csv')
    useless_info = [_ for _ in data.columns if 'software' in _ or 'path' in _]
    data = data.drop(useless_info, axis=1)
    data.to_csv(basic_info, index=True)
    ############################################################
    # merge for phandango
    merged_path = os.path.join(outdir, 'For phandango.csv')
    merged_df = merge_for_phandango(outdir)
    merged_df.to_csv(merged_path, index=1)
    ############################################################
    from pathlib import Path
    infile = os.path.join(indir, 'mashtree_temp_isolates.tre')
    outtre = os.path.join(indir, 'mashtree_isolates.tre')
    treestring = Path(infile).read_text()
    from ete3 import Tree
    tree = Tree(treestring, format=1)
    for i in tree.traverse():
        i.name = i.name.split('_')[0]
    tree.set_outgroup(tree.get_midpoint_outgroup())
    tree.write(outfile=outtre, format=1)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", help='')
    parser.add_argument("-o", "--odir", help='')

    args = parser.parse_args()
    indir = args.indir
    odir = args.odir

    main(str(indir), outdir=odir)
