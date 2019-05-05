import os

from glob import glob
from utils import run_cmd

quast_path = "/home/liaoth/.local/bin/quast.py"
quast_cmd = """{exe_path} {contig} -r {ref} -g {gff} --circos --gene-finding -1 {r1} -2 {r2} --threads {threads} --no-check -o {odir} """
ref = "/home/liaoth/project/shenzhen_actineto/reference/GCF_000746645.1_ASM74664v1_genomic.fna"
gff = "/home/liaoth/project/shenzhen_actineto/reference/GCF_000746645.1_ASM74664v1_genomic.gff"

def main(indir,odir):
    os.makedirs(odir,exist_ok=True)
    for sample_dir in glob(os.path.join(indir,'*')):
        sample_name = os.path.basename(sample_dir)
        if os.path.isdir(sample_dir):
            quast_dir = os.path.join(odir, sample_name)
            os.makedirs(quast_dir,exist_ok=True)
            contig_fn = os.path.join(sample_dir,'contigs.fasta')
            r1 = glob(os.path.join("/home/liaoth/project/shenzhen_actineto/data/after_qc/","S%s*_1.clean.fq.gz" % sample_name))[0]
            r2 = glob(os.path.join("/home/liaoth/project/shenzhen_actineto/data/after_qc/","S%s*_2.clean.fq.gz" % sample_name))[0]
            cmd = quast_cmd.format(exe_path=quast_path,
                               contig=contig_fn,
                               r1=r1,
                               r2=r2,
                               ref=ref,
                               gff=gff,
                               threads=30,
                               odir=quast_dir)
            run_cmd(cmd, dryrun=False)

if __name__ == '__main__':
    #todo: make it a api
    main("/home/liaoth/project/shenzhen_actineto/shovill_output/regular",
         "/home/liaoth/project/shenzhen_actineto/shovill_output/quast_all")

