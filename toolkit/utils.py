from subprocess import check_call
from Bio import SeqIO,Phylo
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

def get_tree(tree_pth,rooted=False):
    t = Phylo.read(tree_pth, 'newick')
    if rooted is False:
        return t
    elif rooted == 'midpoint':
        t.root_at_midpoint()
        return t
    elif rooted in [_.name for _ in t.get_terminals()]:
        rooted_node = [_ for _ in t.get_terminals() if _.name == rooted]
        t.root_with_outgroup(rooted_node[0])
        return t
    return t