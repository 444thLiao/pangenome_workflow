#################################################################################
#### For post process the output of ISEscan with
#### annotated genes(RES,VF etc) and cluster genes(cd-hit of roary)
#################################################################################
import argparse
import os
import sys
from collections import defaultdict
from glob import glob

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipelines.tasks import valid_path
from toolkit.get_gene_info import get_gff
from toolkit.utils import get_locus2group


def get_IS_CDS(ori_gff, IS_gff, locus2annotate, locus2group):
    IS_gff_dict = get_gff(IS_gff,mode="bcbio")

    gff_dict = get_gff(ori_gff,mode="bcbio")
    IS2CDS = defaultdict(list)
    IS2INFO = defaultdict(dict)

    for contig, IS in IS_gff_dict.items():
        for each_IS in IS.features:
            if each_IS.type != "insertion_sequence":
                continue
            ori_chrom = gff_dict[contig]
            start = each_IS.location._start
            end = each_IS.location._end

            IS_id = each_IS.qualifiers['ID'][0]

            get_ori_region = ori_chrom[start:end]

            notis_gene = []
            for cds in get_ori_region.features:
                inference = cds.qualifiers["inference"][-1]
                if 'ISfinder' not in inference:
                    locus_id = cds.id
                    if locus_id in locus2annotate.keys():
                        new_name = locus2annotate[locus_id]
                        # try to find annotate
                    else:
                        new_name = locus2group.get(locus_id, 'unknown')
                        # if failed, it will find in group
                    if new_name != 'unknown':
                        notis_gene.append(new_name)
            if notis_gene:
                IS2CDS[IS_id] = notis_gene
                IS2INFO[IS_id].update(IS.attributes)
    return IS2CDS, IS2INFO


def batch_get(prokka_dir,
              IS_pattern,
              locus2annotate,
              locus2group):
    final_r = {}
    for ori_gff in tqdm(glob(os.path.join(prokka_dir, '*', '*.gff'))):
        sn = os.path.basename(ori_gff).split('.')[0]
        # ISgff = os.path.join(IS_dir,'%s.gff' % sn)
        ISgff = IS_pattern.format(sample_name=sn)
        if '*' in ISgff:
            ISgff = glob(ISgff)[0]
        valid_path(ISgff, check_size=1)
        IS2CDS, IS2INFO = get_IS_CDS(ori_gff,
                                     ISgff,
                                     locus2annotate,
                                     locus2group)
        final_r[sn] = (IS2CDS, IS2INFO)
    # sample name: (IS2CDS,IS2INFO)
    return final_r


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--ISdir", help='')
    parse.add_argument("--l2a", help='')
    parse.add_argument("-r", "--roary_dir", help='')
    parse.add_argument("-p", "--prokka_dir", help='')
    parse.add_argument("-o", "--outfile", help='')

    args = parse.parse_args()
    ISdir = os.path.abspath(args.ISdir)
    ofile = os.path.abspath(args.outfile)
    roary_dir = os.path.abspath(args.roary_dir)
    prokka_dir = os.path.abspath(args.prokka_dir)
    abricate_result_file = os.path.abspath(args.l2a)

    valid_path([ISdir, roary_dir, prokka_dir], check_dir=1)
    valid_path(abricate_result_file, check_size=1)
    valid_path(ofile, check_ofile=1)

    locus2group = get_locus2group(roary_dir)

    locus2annotate_df = pd.read_csv(abricate_result_file, sep=',', index_col=0)
    locus2annotate = locus2annotate_df.loc[:, 'gene'].to_dict()

    final_r = batch_get(prokka_dir=prokka_dir,
                        IS_pattern=os.path.join(ISdir,
                                                "{sample_name}",
                                                "*.gff"),
                        locus2annotate=locus2annotate,
                        locus2group=locus2group)

    ready2df = {sn: {ISgroup: ISgroups.count(ISgroup)
                     for ISgroups in info[0].values()
                     for ISgroup in ISgroups}
                for sn, info in final_r.items()}
    result_df = pd.DataFrame.from_dict(ready2df, orient='index')
    result_df = result_df.fillna(0)
    result_df.to_csv(ofile, sep=',', index=1)

    # python /home/liaoth/project/genome_pipelines/toolkit/process_IS.py -i /home/liaoth/data2/project/genome_pipelines/pipelines/test/test_luigi/ISscan_result --l2a /home/liaoth/data2/project/genome_pipelines/pipelines/test/test_luigi/abricate_result/locus2annotate.csv -r /home/liaoth/data2/project/genome_pipelines/pipelines/test/test_luigi/all_roary_o -p /home/liaoth/data2/project/genome_pipelines/pipelines/test/test_luigi/prokka_o -o /home/liaoth/data2/project/genome_pipelines/pipelines/test/test_toolkit/test_IS_process/result.csv

    #
    # from vis.heatmap import *
    #
    # query_dict = params.res2fun.copy()
    # query_dict.update(params.vf2fun)
    # merged_df = pd.concat([result_df, total_df.loc[:, params.vf_cols]], axis=1)

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
