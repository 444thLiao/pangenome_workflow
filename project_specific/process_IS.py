from collections import defaultdict
from toolkit.get_gene_info import get_gff
from glob import glob
import os
from tqdm import tqdm
import pandas as pd
from toolkit.utils import get_locus2group,get_fasta_by_ID

locus2group = get_locus2group('/home/liaoth/data2/project/shenzhen_Acinetobacter/64_roary_o')

def get_IS_CDS(ori_gff,IS_gff):
    IS_gff_info = get_gff(IS_gff)

    gff_info = get_gff(ori_gff)
    IS2CDS = defaultdict(list)
    IS2INFO = defaultdict(dict)
    for IS in list(IS_gff_info.features_of_type("insertion_sequence")):
        start = IS.start
        end = IS.end
        seqid = IS.seqid
        strand = IS.strand
        id = IS.attributes['ID'][0]
        get_region = gff_info.region(region="%s:%s-%s" % (seqid, start, end),
                                     strand=strand)
        notis_gene = []
        for cds in get_region:
            try:
                if 'ISfinder' not in cds.attributes['inference'][-1]:

                    group = locus2group.get(cds.attributes['locus_tag'][0],'unknown')
                    if group != 'unknown':
                        notis_gene.append(group)
            except:
                import pdb;pdb.set_trace()
        if notis_gene:
            IS2CDS[id] = notis_gene
            IS2INFO[id].update(IS.attributes)
    return IS2CDS,IS2INFO

def batch_get():
    final_r = {}
    for ori_gff in tqdm(glob("/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/*.gff")):
        sn = os.path.basename(ori_gff).split('.')[0]
        ISgff = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ready_prokka/ISscan_result/ready_prokka/%s.fasta.gff" % sn
        if not os.path.isfile(ISgff):
            continue
        IS2CDS, IS2INFO = get_IS_CDS(ori_gff, ISgff)
        final_r[sn] = (IS2CDS,IS2INFO)
    return final_r

if __name__ == '__main__':
    ISgff = "/home/liaoth/data2/project/shenzhen_Acinetobacter/ready_prokka/ISscan_result/ready_prokka/34960.fasta.gff"
    ori_gff = "/home/liaoth/data2/project/shenzhen_Acinetobacter/prokka_html/prokka_gff/34960.gff"
    IS2CDS,IS2INFO = get_IS_CDS(ori_gff,ISgff)

    final_r = batch_get()
    ready2df = {k:{group:_v.count(group) for _v in v[0].values() for group in _v} for k,v in final_r.items()}
    result_df = pd.DataFrame.from_dict(ready2df,orient='index')

    from vis.heatmap import *
    merged_df = pd.concat([result_df, total_df.loc[:, params.vf_cols]], axis=1)

    tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/21_outgroup_roary_o/core_gene.newick"
    rooted = 'AB030'
    fig, mmatrix, amatrix = main(merged_df,
                                 result_df.columns,
                                 filter_sample=samples21,
                                 tree_p=tree_p,
                                 rooted=rooted,
                                 accessory_cols=params.vf_cols,
                                 width=2000, height=1000)
    mmatrix.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/IS_region/21_heatmap_with_metadata.csv", index=True)
    plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/IS_region/21_heatmap_with_metadata.html",
                        auto_open=False)

    all_results = []
    for col in tqdm(result_df.columns):
        results = get_fasta_by_ID('/home/liaoth/data2/project/shenzhen_Acinetobacter/64_roary_o',col)
        if len(results) != 1:
            import pdb;pdb.set_trace()
        else:
            all_results += results
    from Bio import SeqIO
    with open('/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/IS_region/IS_ref.fa','w') as f1:
        SeqIO.write(all_results,f1,format='fasta')