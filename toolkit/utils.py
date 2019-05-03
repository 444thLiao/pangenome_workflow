from subprocess import check_call
from Bio import SeqIO
import os
from ete3 import Tree
from collections import defaultdict
import itertools
import pandas as pd


def run_cmd(cmd, dryrun=False, **kwargs):
    if dryrun:
        print(cmd)
    else:
        try:
            check_call(cmd, shell=True, executable="/usr/bin/zsh", **kwargs)
        except:
            print('##' * 50, '\n', cmd, '\n', "##" * 50)


def get_length_fasta(fasta):
    fh = SeqIO.parse(fasta, format='fasta')
    return {_.id: len(_.seq) for _ in fh}


def get_locus2group(roary_dir):
    protein_cluster_info = os.path.join(roary_dir, 'clustered_proteins')
    locus2group = dict()
    for row in open(protein_cluster_info).readlines():
        row = row.strip('\n')
        cluster_group, remained_locus = row.split(': ')
        for locus in remained_locus.split('\t'):
            locus2group[locus] = cluster_group
    return locus2group

def get_distance(tree_pth):
    t = Tree(open(tree_pth).read())

    pairwise_matrix = defaultdict(dict)
    all_leaves = t.get_leaf_names()
    for t1, t2 in itertools.combinations(all_leaves, 2):
        n1,n2 = t.get_leaves_by_name(t1)[0],t.get_leaves_by_name(t2)[0]
        pairwise_matrix[t1][t2] =pairwise_matrix[t2][t1] = n1.get_distance(n2)
    dist = pd.DataFrame.from_dict(pairwise_matrix)
    dist = pd.DataFrame(dist, index=dist.index, columns=dist.index)
    dist = dist.fillna(0)
    return dist