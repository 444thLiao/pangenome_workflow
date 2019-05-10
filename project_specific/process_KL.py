from vis.heatmap import *
############################################################
# KL specific

locus2annotate_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/all_65_annotated_locus2annotate.tab",sep='\t',index_col=0)
locus2annotate = dict(zip(locus2annotate_df.index,locus2annotate_df.loc[:,'gene']))

KL_df = pd.read_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/fkpA_1-lldP_region/roary_o/gene_presence_absence.csv",
                    sep=',', index_col=0)
KL_df = KL_df.iloc[:, 13:].T

new_cols = {}
dropped_col = []
for col,vals in KL_df.iteritems():
    new_vals = []
    for v in list(vals):
        if v in locus2annotate.keys():
            g = locus2annotate[v]
            dropped_col.append(col)
            col = g
        else:
            g = 1 if not pd.isna(v) else 0
        new_vals.append(g)
    new_cols[col] = new_vals

KL_df = KL_df.apply(lambda x: [1 if not pd.isna(v) else 0 for idx,v in enumerate(x)],axis=0,broadcast='list')

merged_df = pd.concat([KL_df, total_df.loc[:, params.vf_cols]], axis=1)
for samples in [samples20, samples21]:
    if len(samples) == 20:
        tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/20_outgroup_roary_o/core_gene.newick"
        rooted = 'AB030'
    else:
        tree_p = "/home/liaoth/data2/project/shenzhen_Acinetobacter/21_outgroup_roary_o/core_gene.newick"
        rooted = 'AB030'
    fig, mmatrix, amatrix = main(merged_df,
                                 KL_df.columns,
                                 filter_sample=samples,
                                 tree_p=tree_p,
                                 rooted=rooted,
                                 accessory_cols=params.vf_cols,
                                 width=2000, height=1000)
    mmatrix.to_csv("/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/%s_heatmap_with_metadata.csv" % len(samples), index=True)
    plotly.offline.plot(fig, filename="/home/liaoth/data2/project/shenzhen_Acinetobacter/ad_analysis/KL_extracted_reads/%s_heatmap_with_metadata.html" % len(samples),
                        auto_open=False)
