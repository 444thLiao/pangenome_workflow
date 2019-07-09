import io
import os
import sys
from glob import glob
from os.path import join, dirname, basename

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipelines import run_cmd
from pipelines.constant_str import assembly_summary_header


def batch_mash(prokka_dir, db='', thread=40):
    db = os.path.abspath(db)
    header = ["reference-ID", "query-ID", "distance", "p-value", "shared-hashes"]
    aln_dict = {}

    for fna in tqdm(glob(join(prokka_dir, '*', '*.fna'))):
        sid = basename(dirname(fna))
        output = run_cmd(f"mash dist {db} {fna} -p {thread} -d 0.03 -v 0.05",
                         get_output=True,
                         dry_run=False)
        try:
            a = pd.read_csv(io.StringIO(output), sep='\t', header=None, )
            # 97% for species level identity.
            a.columns = header
            a = a.sort_values("distance")
            # similar_id = ['_'.join(_.split('_')[:2]) for _ in similar_fna]
            result = ("97% cutoff most likely", a)
        except pd.errors.EmptyDataError:
            output = run_cmd(f"mash dist {db} {fna} -p {thread} -v 0.05",
                             get_output=True,
                             dry_run=False)
            a = pd.read_csv(io.StringIO(output), sep='\t', header=None, )
            # read in a total dataframe instead of one name
            a.columns = header
            result = ("no cutoff", a)
        aln_dict[sid] = result
    return aln_dict


def merge_mash(aln_dict, refseq_summary=''):
    refseq_df = pd.read_csv(refseq_summary, sep='\t', index_col=None,
                            header=None, low_memory=False, comment="#")
    refseq_df.columns = assembly_summary_header
    aid2name = dict(zip(list(refseq_df.loc[:, 'assembly_accession']),
                        list(refseq_df.loc[:, 'organism_name'])))
    aln2org = {}
    aln2org_dis = {}
    for sid, (annotated, df) in aln_dict.items():

        sorted_df = df.sort_values('distance')
        rids = ['_'.join(_.split('_')[:2])
                for _ in list(sorted_df.iloc[:, 0])]
        orgs = [aid2name.get(_, '')
                for _ in rids]
        org2rids = dict(zip(orgs, rids))
        if annotated == "no cutoff":
            if len(orgs) >= 5:
                aln2org[sid] = max(orgs[:5], key=lambda x: orgs[:5].count(x))
            else:
                aln2org[sid] = max(orgs, key=lambda x: orgs.count(x))
            # get the mean distance of annotated organism
        else:
            aln2org[sid] = max(orgs, key=lambda x: orgs.count(x))

        aln2org_dis[sid] = sorted_df.loc[sorted_df.iloc[:, 0].str.contains(org2rids[aln2org[sid]]),
                                         "distance"].mean()
    result_df = pd.DataFrame(columns=["mash distance", "annotated organism"])
    for sid in aln2org:
        result_df.loc[sid, :] = [aln2org_dis[sid],
                                 aln2org[sid]]

    return result_df


def parse_batch_result(infiles, refseq_summary):
    # For pipelines required
    refseq_df = pd.read_csv(refseq_summary, sep='\t', index_col=None,
                            header=None, low_memory=False, comment="#")
    refseq_df.columns = assembly_summary_header
    aid2name = dict(zip(list(refseq_df.loc[:, 'assembly_accession']),
                        list(refseq_df.loc[:, 'organism_name'])))
    ###################
    header = ["reference-ID", "query-ID", "distance", "p-value", "shared-hashes"]
    aln2org = {}
    aln2org_dis = {}
    for path in infiles:
        sid = os.path.basename(path).split('_')[0]
        a = pd.read_csv(path, sep='\t', header=None, )
        a.columns = header

        filter_df = a.loc[a.distance <= 0.03, :]
        # get result within species level
        # if within 0.03, take all result for annotating species
        # else: take top 5 or less.
        if filter_df.shape[0] == 0:
            filter_df = a.sort_values('distance')
            if len(filter_df.shape[0]) >= 5:
                filter_df = filter_df.iloc[:5, :]

        rids = ['_'.join(_.split('_')[:2])
                for _ in list(filter_df.iloc[:, 0])]
        orgs = [aid2name.get(_, '')
                for _ in rids]
        org2rids = dict(zip(orgs, rids))
        aln2org[sid] = max(orgs, key=lambda x: orgs.count(x))
        aln2org_dis[sid] = filter_df.loc[filter_df.iloc[:, 0].str.contains(org2rids[aln2org[sid]]),
                                         "distance"].mean()
        # mean all(not just top5) distance with annotated result.

    result_df = pd.DataFrame(columns=["mash distance", "annotated organism"])
    for sid in aln2org:
        result_df.loc[sid, :] = [aln2org_dis[sid],
                                 aln2org[sid]]
    return result_df
