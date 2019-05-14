import os
from collections import defaultdict
from glob import glob

import pandas as pd
from tqdm import tqdm

from toolkit.get_gene_info import get_gff
from toolkit.utils import get_locus2group


def get_IS_CDS(ori_gff, IS_gff, locus2annotate, locus2group):
    IS_gff_info = get_gff(IS_gff)

    gff_info = get_gff(ori_gff)
    IS2CDS = defaultdict(list)
    IS2INFO = defaultdict(dict)
    for IS in list(IS_gff_info.features_of_type("insertion_sequence")):
        start = IS.start
        end = IS.end
        seqid = IS.seqid
        strand = IS.strand
        IS_id = IS.attributes['ID'][0]
        get_region = gff_info.region(region="%s:%s-%s" % (seqid, start, end),
                                     strand=strand)
        notis_gene = []
        for cds in get_region:

            if 'ISfinder' not in cds.attributes['inference'][-1]:
                locus_id = cds.attributes['locus_tag'][0]
                if locus_id in locus2annotate.keys():
                    group = locus2annotate[locus_id]
                else:
                    group = locus2group.get(locus_id, 'unknown')
                if group != 'unknown':
                    notis_gene.append(group)
        if notis_gene:
            IS2CDS[IS_id] = notis_gene
            IS2INFO[IS_id].update(IS.attributes)
    return IS2CDS, IS2INFO


def batch_get(prokka_dir, IS_pattern, locus2annotate, locus2group):
    final_r = {}
    for ori_gff in tqdm(glob(os.path.join(prokka_dir, '*', '*.gff'))):
        sn = os.path.basename(ori_gff).split('.')[0]
        # ISgff = os.path.join(IS_dir,'%s.gff' % sn)
        ISgff = IS_pattern.format(sample_name=sn)
        if not os.path.isfile(ISgff):
            continue
        IS2CDS, IS2INFO = get_IS_CDS(ori_gff, ISgff, locus2annotate, locus2group)
        final_r[sn] = (IS2CDS, IS2INFO)
    # sample name: (IS2CDS,IS2INFO)
    return final_r


if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-r", "--roary_dir", help='')
    parse.add_argument("-p", "--prokka_dir", help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    indir = os.path.abspath(args.indir)
    odir = os.path.abspath(args.outdir)
    roary_dir = os.path.abspath(args.roary_dir)
    prokka_dir = os.path.abspath(args.prokka_dir)

    locus2group = get_locus2group(roary_dir)

    locus2annotate_df = pd.read_csv(abricate_result_file, sep='\t', index_col=0)
    locus2annotate = dict(zip(locus2annotate_df.index, locus2annotate_df.loc[:, 'gene']))

    final_r = batch_get(prokka_dir=prokka_dir,
                        IS_dir=IS_dir)
    ready2df = {k: {group: _v.count(group) for _v in v[0].values() for group in _v}
                for k, v in final_r.items()}
    #
    result_df = pd.DataFrame.from_dict(ready2df, orient='index')
    result_df = result_df.fillna(0)
    from vis.heatmap import *

    query_dict = params.res2fun.copy()
    query_dict.update(params.vf2fun)
    merged_df = pd.concat([result_df, total_df.loc[:, params.vf_cols]], axis=1)

    # todo: directly draw heatmap
    # tree_p = os.path.join(roary_dir,"core_gene.newick")
    # rooted = 'midpoint'
    # fig, mmatrix, amatrix = main(merged_df,
    #                              result_df.columns,
    #                              filter_fea_all=False,
    #                              tree_p=tree_p,
    #                              rooted=rooted,
    #                              map_type=lambda x: query_dict.get(x, "unknown"),
    #                              accessory_cols=params.vf_cols,
    #                              width=2000, height=1000)
    # mmatrix.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/IS_region/21_heatmap_with_metadata.csv", index=True)
    # plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/IS_region/21_heatmap_with_metadata.html",
    #                     auto_open=False)
