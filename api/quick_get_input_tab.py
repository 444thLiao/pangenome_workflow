import os
from glob import glob
from os.path import *

import click

template_path = join(dirname(dirname(__file__)), 'toolkit', 'data_input.template')

header = open(template_path).read()




@click.command()
@click.option("-p", "pattern", help='*.R1.fastq.gz')
@click.option("-o", "ofile", help='', default=None)
@click.option('-s', "simple_split", default='_L1_')
def cli(ofile, pattern, simple_split):
    rows = [header]
    if ofile is None:
        ofile = './data_input.tab'
    else:
        if not exists(dirname(abspath(ofile))):
            os.makedirs(dirname(abspath(ofile)))
    fs = glob(pattern)
    for r1 in fs:
        name = basename(r1).split(simple_split)[0]
        R2 = r1.replace('.R1.', '.R2.')
        rows.append('\t'.join([name, r1, R2]))
    with open(ofile, 'w') as f1:
        f1.write('\n'.join(rows))


if __name__ == '__main__':
    cli()
