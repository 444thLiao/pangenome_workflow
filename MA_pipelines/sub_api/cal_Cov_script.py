import sys
from os.path import dirname

sys.path.insert(0, dirname(dirname(__file__)))
import argparse
import os

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

from parse_file_name import fileparser
import multiprocessing as mp

"""
final edit at 20190617: add multiprocessing exec
Function to cal each base info from given bed file and bam file.
with same bed file and same fasta file

:type Cal_Function
"""


def bam2info(bam_path, output_cov, bed_file, REF_file):
    print('Start loading required file......')
    bed_file = pd.read_csv(bed_file,
                           index_col=False,
                           sep='\t',
                           header=None)
    try:
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
    except:
        raise IOError("bamefile need to cal doesn't exist")
    fastafile = pysam.FastaFile(filename=REF_file)
    f1 = open(output_cov, 'w')
    # define the header
    header = ['Gene',
              'Chr',
              'Position',
              'Reference',
              'base',
              'A',
              'C',
              'G',
              'T',
              'A_Rate',
              'C_Rate',
              'G_Rate',
              'T_Rate',
              '1',
              '2',
              '3',
              '4']
    f1.write('\t'.join(header) + '\n')
    f1.flush()
    # define columns.
    # print 'Loading required file. using %d' % (time.time()-t1)
    print('Iterating regions within bed file')
    for idx, row in tqdm(bed_file.iterrows()):
        Chr, start, end, Gene_name = row[:4]
        start = min(int(start), int(end))
        end = max(int(start), int(end))
        # fetch basic info.
        coverage_ACGT = bamfile.count_coverage(Chr,
                                               start,
                                               end,
                                               read_callback='nofilter',
                                               quality_threshold=0)
        ###fixed coordinate
        coverage_ACGT = np.array(coverage_ACGT)
        ref_base_str = fastafile.fetch(Chr, start, end)
        # return str
        if coverage_ACGT.sum().sum() != 0:
            # total base for all ACGT doesn't equal to 0
            for relative_pos, real_pos in enumerate(range(start, end)):
                depth_ = coverage_ACGT[:, relative_pos].sum()
                # total base num at this pos
                num_A, num_C, num_G, num_T = coverage_ACGT[:, relative_pos]

                if depth_ == 0:
                    # if position here didn't have any base, then pass this pos.
                    continue
                row = [str(Gene_name) if str(Gene_name) == 'nan' else '',
                       Chr,
                       str(real_pos),
                       ref_base_str[relative_pos].upper(),
                       str(depth_),
                       str(num_A),
                       str(num_C),
                       str(num_G),
                       str(num_T)]  # left an empty line
                for num_base in [num_A, num_C, num_G, num_T]:
                    row.append(str(
                        round(float(num_base) / depth_,
                              4)))
                    # wirte the rate of A/T/C/G RATE
                ranked_base = sorted(list("ACGT"),
                                     key=lambda x: coverage_ACGT["ACGT".index(x),
                                                                 relative_pos],
                                     reverse=True)
                row += list(ranked_base)
                # construct a list to sort A/T/C/G each base num.
                row = '\t'.join(row) + '\n'
                f1.write(row)
                f1.flush()


def exec_fun(args):
    func, kwargs = args
    func(**kwargs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input_tab', required=True, type=str)
    parser.add_argument('-o', dest='output', required=True, type=str)
    parser.add_argument('-B', '--bed', dest='bed_path', type=str, help='bed file path')
    parser.add_argument('-r', '--ref', dest='ref_fasta', type=str, help='reference fasta file path', default='/home/liaoth/data/hg19/ucsc.hg19.fasta')
    parser.add_argument('-t', '--thread', dest='thread', type=int, help='how many process you want to used as calculating coverage.', default=0)
    args = parser.parse_args()

    input_tab = os.path.abspath(args.input_tab)
    odir = os.path.abspath(args.output)
    bed_path = os.path.abspath(args.bed_path)
    ref_path = os.path.abspath(args.ref_fasta)
    thread = mp.cpu_count() if args.thread == 0 else args.thread
    df = fileparser(input_tab)
    sample_dict = df.get_full_info(odir)

    params = []
    for sid in sample_dict.keys():
        infodict = sample_dict[sid]
        for alignment, output_cov in zip([infodict["sorted_bam"],
                                          infodict["recal_bam"]],
                                         [infodict["sorted_cov"],
                                          infodict["recal_cov"]]):
            params.append((bam2info, dict(bam_path=alignment,
                                          output_cov=output_cov,
                                          bed_file=bed_path,
                                          REF_file=ref_path)
                           ))
    with mp.Pool(processes=thread) as tp:
        tp.map(exec_fun, params)
