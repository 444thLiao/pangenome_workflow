"""
parse the shovil log and extract estimated coverage
"""
import click
from glob import glob
from os.path import *
from tqdm import tqdm

def get_depth(f):
    depth = 0
    for row in open(f):
        row = row.strip('\n')
        if 'Estimated sequencing depth' in row:
            depth = row.split(' ')[-2]
    return int(depth)

def get_coverage(f):
    coverage = 0
    for row in open(f):
        row = row.strip('\n')
        if 'Mean total coverage:' in row:
            coverage = row.split(' ')[-1]
    return int(coverage)

def main(indir,):
    all_logs = glob(join(indir,'*','*.log'))
    id2coverage = {}
    id2depth = {}
    for l in tqdm(all_logs):
        name = basename(dirname(l))
        coverage = get_coverage(l)
        id2coverage[name] = coverage
        id2depth[name] = get_depth(l)
    return id2depth,id2coverage

id2depth,id2coverage = main("/mnt/storage2/jjtao/jjtao_20200119/pipelines_o/assembly_o/regular/")
ofile = '/home-user/thliao/tmp/jjtao_20200119.depth.txt'
with open(ofile,'w') as f1:
    for id,c in id2depth.items():
        f1.write(f"{id}\t{c}\n")


def cli():
    pass