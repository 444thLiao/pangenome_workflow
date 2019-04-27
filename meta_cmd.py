from glob import glob
from tqdm import tqdm
import os
from subprocess import check_call

odir = os.path.abspath('/home/liaoth/project/shenzhen_myco/shovill_output/')
for fq1 in tqdm(glob("/home/liaoth/project/shenzhen_myco/data/raw/cleandata/*/*_1.clean.fq.gz")):
    sample_name = os.path.basename(os.path.dirname(fq1).strip('-1a'))
    real_odir = os.path.join(odir,sample_name)
    fq2 = fq1.replace('_1','_2')
    cmd = "/home/liaoth/anaconda3/envs/shovill/bin/shovill --outdir {odir} --ram 100 --R1 {R1} --R2 {R2}".format(odir=real_odir,
                                                                                     R1 = fq1,
                                                                                     R2=fq2)
    # print(cmd)
    check_call(cmd.split(' '))

for fa in tqdm(glob("/home/liaoth/project/shenzhen_myco/shovill_output/cleandata_o/*/contigs.fa")):

    dirname = os.path.basename(os.path.dirname(fa))
    nname = dirname.split('_')[0].strip('S') +'.fa'
    os.system("cp %s %s" % (fa,os.path.join(os.path.dirname(fa),nname)))


for fa in tqdm(glob("/home/liaoth/project/shenzhen_myco/ready_roary_input/*.f*")):
    cmd = "source /tools/.bashrc; /tools/prokka/bin//prokka %s --outdir %s --cpus 20" % (fa,os.path.basename(fa).split('.')[0])
    os.system(cmd)