from glob import glob
from tqdm import tqdm
from subprocess import check_call
import os

def main(idir,odir):
    for seq in tqdm(glob(os.path.join(idir,'*.f*'))):
        sample_name = os.path.basename(seq).split('.')[0]
        cmd = "prokka {infile} --outdir {odir} --prefix {sn} --force --quiet"

        check_call(cmd.format(infile=seq,
                              odir=odir,
                              sn=sample_name),executable="/usr/bin/zsh",shell=True)


if __name__ == '__main__':
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument("-i","--indir",help='')
    parse.add_argument("-o", "--outdir", help='')

    args = parse.parse_args()
    main(args.indir,args.outdir)