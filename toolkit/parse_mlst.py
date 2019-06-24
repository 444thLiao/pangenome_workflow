import pandas as pd

def parse_mlst(infiles):

    collect_df = []
    for infile in infiles:
        indf = pd.read_csv(infile,sep='\t',index_col=None,header=None)
        header = ["raw_path","Scheme","ST"]
        header += ["Locus%s" % (_+1)
                   for _ in range(indf.shape[1]-len(header))]
        indf.columns = header
        collect_df.append(indf)
    final_df = pd.concat(collect_df,axis=0)
    return final_df
