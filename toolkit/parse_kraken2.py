import os
from collections import defaultdict

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


def parse_kraken2(infile, output_df=False):
    # todo: check some abnormal situations.
    df = pd.read_csv(infile, sep='\t', header=None)
    df.columns = kraken2_header
    df.loc[:, "scientific name"] = [_.strip() for _ in df.loc[:, "scientific name"]]
    sorted_df = df.sort_values("num assigned frag", ascending=False)

    if output_df:
        return sorted_df

    get_lt90_df = sorted_df.loc[sorted_df.iloc[:, 0] >= 90, :]
    if get_lt90_df.shape[0] != 0:

        return (get_lt90_df.iloc[0, -1], get_lt90_df.iloc[0, 0])

    else:
        return ("closely like %s" % sorted_df.iloc[0, -1],
                sorted_df.iloc[0, 0])


def merge_kraken2(infiles):
    merged_df = pd.DataFrame(columns=["read_kraken2_classify",
                                      "read_kraken2 (%)",
                                      "assembly_kraken2_classify",
                                      "assembly_kraken2 (%)",
                                      ])
    s2paths = defaultdict(list)
    for _ in infiles:
        sample = os.path.basename(_).split('_assembly')[0]
        sample = sample.split("_reads")[0]
        s2paths[sample].append(_)

    for s, paths in s2paths.items():
        data = []
        for _ in paths:
            if _.endswith("_reads.k2report"):
                result = parse_kraken2(_)
                data += result
        if not data:
            data += ["",""]
        for _ in paths:
            if _.endswith("_assembly.k2report"):
                result = parse_kraken2(_)
                data += result
        if not data:
            data += ["",""]
        try:
            assert len(data) == 4
        except:
            import pdb;pdb.set_trace()

        merged_df = merged_df.append(pd.DataFrame([data],
                                                  columns=merged_df.columns,
                                                  index=[s]))

    return merged_df
