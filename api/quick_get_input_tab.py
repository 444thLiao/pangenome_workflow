import os
from glob import glob
from os.path import *

import click,re

template_path = join(dirname(dirname(__file__)), 'toolkit', 'data_input.template')

header = open(template_path).read()

@click.command()
@click.option("-p", "pattern", help='*.R1.fastq.gz')
@click.option("-o", "ofile", help='', default=None)
@click.option('-s', "name_pattern", default='_L1_')
@click.option('-PE', "pair_format", default='_1.fq')
def cli(ofile, pattern, name_pattern,pair_format):
    rows = [header]
    if ofile is None:
        ofile = './data_input.tab'
    else:
        if not exists(dirname(abspath(ofile))):
            os.makedirs(dirname(abspath(ofile)))
    fs = glob(pattern)
    for r1 in fs:
        name = re.findall(name_pattern,basename(r1))[0]
        R2 = r1.replace(pair_format, pair_format.replace('1','2'))
        rows.append('\t'.join([name, r1, R2]))
    with open(ofile, 'w') as f1:
        f1.write('\n'.join(rows))


if __name__ == '__main__':
    cli()
