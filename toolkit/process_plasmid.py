import os
import re
import sys
from collections import defaultdict
from glob import glob
from subprocess import check_output

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from toolkit.utils import run_cmd, valid_path


def align_plasmid(indir):
    # indir may the assembly_o which contains both regular and plasmidspades
    # parse the plasmid result into a unify format of data
    # unify format:
    # {sample_name : [(chrom:start-end, other info),
    #                 (chrom2:start2-end2, other info)]
    #                 }
    aligned_info = {}
    for contig in tqdm(glob(os.path.join(indir,
                                         'regular',
                                         '*',
                                         'contigs.fa'))):
        run_cmd("bwa index %s" % contig, dry_run=False, log_file='/tmp/tmp.log')
        plasmid_contig = contig.replace('/regular/',
                                        '/plasmidsSpades/') + 'sta'
        osam = plasmid_contig.replace('.fasta', '.sam')
        # regular is 'contigs.fa', plasmids must is 'contigs.fasta'
        if not os.path.isfile(plasmid_contig):
            # it may the SE data, so it should not process plasmid
            continue
        plasmid_count_dict = defaultdict(list)
        sample_name = os.path.basename(os.path.dirname(contig))
        if os.path.getsize(plasmid_contig) > 0:
            if not os.path.isfile(osam):
                sam_content = check_output("bwa mem -x intractg {regular_one} {plasmid_one} > {ofile}".format(regular_one=contig,
                                                                                                              plasmid_one=plasmid_contig,
                                                                                                              ofile=osam),
                                           shell=True,
                                           executable='/usr/bin/zsh')
                sam_content = sam_content.decode('utf-8')
            else:
                sam_content = open(osam, 'r').read()
            validated_row_of_sam = [_ for _ in sam_content.split('\n') if (not _.startswith('@')) and (_) and (not _.startswith('*'))]

            matched_plasmid_belong = [re.findall("component_([0-9]+)$",
                                                 row.split('\t')[0])[0]
                                      for row in validated_row_of_sam]
            # parse the contig name of plasmid to get belong plasmid id
            match_region = ["%s:%s-%s" % (row.split('\t')[2],
                                          int(row.split('\t')[3]) - 1,  # convert 0-coord
                                          int(row.split('\t')[3]) - 1 + len(row.split('\t')[9]))
                            for row in validated_row_of_sam]
            aligned_info[sample_name] = list(zip(match_region,
                                                 matched_plasmid_belong))

    return aligned_info


def main(indir, ofile):
    # indir may the assembly_o which contains both regular and plasmidspades
    aligned_plasmid_info = align_plasmid(indir)
    summary_df = pd.DataFrame(
        columns=["sample",
                 "region",
                 "other info"])

    for sample_name, infos in aligned_plasmid_info.items():
        count = 0
        for region, plasmid_id in infos:
            formatted_ID = 'plasmid%s_piece%s' % (str(plasmid_id),
                                                  str(count))
            summary_df = summary_df.append(pd.DataFrame([[sample_name,
                                                          region,
                                                          ''
                                                          ]],
                                                        columns=summary_df.columns,
                                                        index=[formatted_ID]))
            count += 1
    summary_df.to_csv(ofile, sep='\t', index=1)


if __name__ == '__main__':
    import argparse

    """
    This file main to aligned assembly contigs which from plasmidSpades and corresponding contigs which from regular Spades.
    For coventinal usage, it is a script embed at pangenome_pipeliens, so the file structure is fixed(Especially the plasmidSpades and regular Spades place.) 
    Finally, it will output a summary file. Of course it will output sam file into plasmidSpades dir. 
    """
    parse = argparse.ArgumentParser()
    parse.add_argument("-i", "--indir", help='')
    parse.add_argument("-o", "--ofile", help='')

    args = parse.parse_args()
    indir = os.path.abspath(args.indir)
    ofile = os.path.abspath(args.ofile)

    valid_path(ofile, check_ofile=1)
    valid_path(indir, check_dir=1)
    main(indir, ofile)
