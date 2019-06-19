import pandas as pd

# 1. Percentage of fragments covered by the clade rooted at this taxon
# 2. Number of fragments covered by the clade rooted at this taxon
# 3. Number of fragments assigned directly to this taxon
# 4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
#    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
#    Taxa that are not at any of these 10 ranks have a rank code that is
#    formed by using the rank code of the closest ancestor rank with
#    a number indicating the distance from that rank.  E.g., "G2" is a
#    rank code indicating a taxon is between genus and species and the
#    grandparent taxon is at the genus rank.
# 5. NCBI taxonomic ID number
# 6. Indented scientific name

kraken2_header = ["percentage_frag",
                  "num frag",
                  "num assigned frag",
                  "rank code",
                  "NCBI taxid",
                  "scientific name"]


def parse_kraken2(infile,output_df=False):
    # todo: check some abnormal situations.
    df = pd.read_csv(infile, sep='\t', header=None)
    df.columns = kraken2_header
    df.loc[:, "scientific name"] = [_.strip() for _ in df.loc[:, "scientific name"]]
    sorted_df = df.sort_values("num assigned frag", ascending=False)

    if output_df:
        return sorted_df

    if sorted_df.iloc[0, 0] >= 90:
        return sorted_df.iloc[0, -1]
    else:
        return "unclassified, closely is %s" % sorted_df.iloc[0, -1]


def merge_kraken2(infiles):
    assembly_files = [pd.read_csv(_.path, index_col=0, sep='\t')
                      for _ in infiles
                      if _.endswith("_assembly.k2report")]
    reads_files = [pd.read_csv(_.path, index_col=0, sep='\t')
                   for _ in infiles
                   if _.endswith("_reads.k2report")]
    df_list = [parse_kraken2(_,
                             output_df=True)
               for _ in infiles]
