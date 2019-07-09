# pre roary based on pariwised mash distance
# For pangenome workflow pipelines

import io
import os
import sys
from os.path import basename, join

import click
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from pipelines import run_cmd, valid_path
from pipelines.soft_db_path import mash_path


@click.command()
@click.argument('infiles', nargs=-1)
@click.option("-o", "--odir", help='output directory')
@click.option("--db", help='ref database for pairwise mash, could be none. It will auto named')
@click.option("--thread", help='number of threads to perform mash analysis.')
@click.option("--cutoff", help="cutoff values which taken as bound of each cluster.")
def pairwise_mash(infiles,
                  odir,
                  thread,
                  org_annotated,
                  db=None,
                  cutoff=0.05,
                  force_cmd=False):
    if db is None:
        db = join(odir, 'pairwise_ref.msh')
    if not db.endswith('.msh'):
        db += '.msh'
    valid_path(db, check_ofile=1)
    valid_path(odir, check_odir=1)

    infiles_str = ' '.join(infiles)
    if (not os.path.exists(db)) or force_cmd:
        run_cmd(f"{mash_path} sketch -o {db} {infiles_str} -p {thread}",
                dry_run=False)
    # read in organism annotated dataframe
    aln2org = pd.read_csv(org_annotated, index_col=0, dtype=str)
    assert len(infiles) <= aln2org.shape[0]
    cache_result = run_cmd(f"{mash_path} dist {db} {infiles_str} -p {thread}",
                           dry_run=False, get_output=1)
    parsed_df = pd.read_csv(io.StringIO(cache_result), sep='\t', header=None)
    parsed_df.columns = ["reference-ID", "query-ID", "distance", "p-value", "shared-hashes"]

    # assume infiles formatted like blablabla/sid.fna; sid mean sample ID
    parsed_df.loc[:, 'reference-ID'] = [basename(_).split('.')[0]
                                        for _ in parsed_df.loc[:, 'reference-ID']]
    parsed_df.loc[:, "query-ID"] = [basename(_).split('.')[0]
                                    for _ in parsed_df.loc[:, "query-ID"]]
    parsed_df = parsed_df.iloc[:, [0, 1, 2]]
    pairwise_df = parsed_df.pivot('reference-ID', "query-ID")
    pairwise_df.columns = pairwise_df.index
    pairwise_df.to_csv(join(odir,
                            "pairwise_mash.dist"),
                       index=1, index_label=pairwise_df.index.name)
    model = AgglomerativeClustering(distance_threshold=cutoff,
                                    n_clusters=None,
                                    affinity='precomputed',
                                    linkage='complete')
    clustering_r = model.fit_predict(pairwise_df)
    aln2org.loc[:, "clustering"] = 0
    for g in set(clustering_r):
        sids = pairwise_df.index[clustering_r == g]
        aln2org.loc[sids, "clustering"] = g
    ofile = join(odir, 'clustered_sid.csv')
    aln2org.to_csv(ofile, index=True, index_label=aln2org.index.name)

    return aln2org

if __name__ == '__main__':
    pairwise_mash()
