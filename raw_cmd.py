import os
import sys
from glob import glob
from subprocess import check_call

import pandas as pd


def run_cmd(cmd, dryrun=False, **kwargs):
    if dryrun:
        print(cmd)
    else:
        try:
            check_call(cmd, shell=True, executable="/usr/bin/zsh", **kwargs)
        except:
            print('##' * 50, '\n', cmd, '\n', "##" * 50)


# todo: make it a pipelines instead a bunch of commands
#####################################################################
# qc
fastqc_path = "/tools/FastQC/fastqc"
fastqc_cmd = "{exe_path} {in_pattern} -t {threads} -o {odir} --quiet"
fastqc_cmd.format(exe_path=fastqc_path,
                  in_pattern='',
                  threads=40,
                  odir='')

#####################################################################
# multiqc
multiqc_path = "/home/liaoth/.local/bin/multiqc"
multiqc_cmd = "{exe_path} {indir} --outdir {odir} --filename {fn} --force -q"
multiqc_cmd.format(exe_path=multiqc_path,
                   indir='',
                   odir='',
                   fn='')
#####################################################################
# trimmomatic
trimmomatic_setting = "ILLUMINACLIP:/home/liaoth/tools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5"
jar_path = "/home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar"
trimmomatic_cmd = """java -jar {exe_path} PE -threads {threads} {R1} {R2} -trimlog {log} {clean_r1} {unpaired_r1} {clean_r2} {unpaired_r2} {params}"""


def get_trimmomatic_cmd(R1, R2, odir, threads, split_str='.f'):
    _R1 = os.path.basename(R1)
    _R2 = os.path.basename(R2)
    clean_r1 = os.path.join(odir, _R1.split(split_str)[0]) + '.clean.fq.gz'
    unpaired_r1 = os.path.join(odir, _R1.split(split_str)[0]) + '.unpaired.fq.gz'
    clean_r2 = os.path.join(odir, _R2.split(split_str)[0]) + '.clean.fq.gz'
    unpaired_r2 = os.path.join(odir, _R2.split(split_str)[0]) + '.unpaired.fq.gz'
    log = os.path.join(odir, _R1.split(split_str)[0]) + '.log'
    cmd = trimmomatic_cmd.format(exe_path=jar_path,
                                 threads=threads,
                                 R1=R1,
                                 R2=R2,
                                 log=log,
                                 clean_r1=clean_r1,
                                 unpaired_r1=unpaired_r1,
                                 clean_r2=clean_r2,
                                 unpaired_r2=unpaired_r2,
                                 params=trimmomatic_setting)
    return cmd


#####################################################################
# spades/shovill
shovill_path = "/tools/anaconda3/envs/shovill/bin/shovill"
shovill_cmd = """source /tools/.bashrc;source activate shovill;{exe_path} --outdir {odir} --ram 100 --R1 {r1} --R2 {r2} --depth {depth} --cpus {thread}"""
spades_extra_indicator = "--opts"
spades_extra_options = ""


def get_shovill_cmd(r1, r2, odir, depth, thread, spades_extra_options, extra_option=''):
    if spades_extra_options:
        extra_str = ' ' + spades_extra_indicator + ' ' + spades_extra_options
    else:
        extra_str = ''
    cmd = shovill_cmd.format(exe_path=shovill_path,
                             r1=r1,
                             r2=r2,
                             odir=odir,
                             depth=depth,
                             thread=thread)
    cmd = ' '.join([cmd, extra_option, extra_str, ]) + ';source deactivate'
    return cmd


#####################################################################
# prokka
prokka_path = "/tools/prokka/bin/prokka"
prokka_cmd = "{exe_path} {in_fasta} --outdir {odir} --cpus {threads} --prefix {sample_name} --force "
prokka_cmd.format(exe_path=prokka_path,
                  in_fasta='',
                  odir='',
                  threads=30,
                  sample_name='')
#####################################################################
# roary
roary_path = "/usr/local/bin/roary"
roary_cmd = "{exe_path} -v -e --mafft -p {threads} -f {odir} -r {gff_pattern}"
roary_cmd.format(exe_path=roary_path,
                 threads=35,
                 odir='',
                 gff_pattern='')

#########################################
# config
pandoo_path = "/tools/anaconda3/envs/pandoo/bin/pandoo"
ariba_db = "/home/db_public/ariba_db"
abricate_db = "/tools/abricate/db/"

ariba_db_list = ["CARD", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_full", "virulencefinder"]
ariba_str = str({k: "{db}/{name}".format(db=ariba_db, name=k.lower()) for k in ariba_db_list})
abricate_db_list = ['card', "argannot", "ncbi", "vfdb", "vfdb_full", "resfinder", "plasmidfinder", "victors"]
abricate_str = str({k: "{db}".format(db=abricate_db) for k in abricate_db_list})
#########################################
# preset cmd

pandoo_cmd = """source /tools/.bashrc;source activate pandoo;{exe_path} run -i {input} -o {odir} -t -r -b "{ariba_str}" -a "{abricate_str}" ;source deactivate """
# -t: --infer_tree_on; -r: --ariba_on; -b: --ariba_dbs;

pandoo_cmd.format(exe_path=pandoo_path,
                  ariba_str=ariba_str,
                  abricate_str=abricate_str,
                  input="isolates_pre.tab",
                  odir='')
# print(pandoo_cmd)


if __name__ == '__main__':
    dryrun = False
    isolates_df = pd.DataFrame()
    indir = os.path.abspath('/home/liaoth/project/shenzhen_actineto/data/raw/cleandata')
    base_odir = os.path.abspath("/home/liaoth/project/shenzhen_actineto/")
    global_err_log = os.path.join(base_odir, 'global.err')

    reads_pattern = os.path.join(indir, '*', '*.gz')
    r1_pattern = os.path.join(indir, '*', '*_1.clean.fq.gz')
    with open(global_err_log, 'w') as sys.stderr:
        #####################################################################
        # QC
        qc_before = os.path.join(base_odir, "fastqc", "fastqc_before")
        qc_after = os.path.join(base_odir, 'fastqc', "fastqc_after")
        os.makedirs(qc_before, exist_ok=True)
        os.makedirs(qc_after, exist_ok=True)
        cmd = fastqc_cmd.format(exe_path=fastqc_path,
                                in_pattern=reads_pattern,
                                threads=40,
                                odir=qc_before)
        run_cmd(cmd, dryrun=dryrun)
        #####################################################################
        # trimmomatic
        after_qc_seq_dir = os.path.join(base_odir, 'data', 'after_qc')
        os.makedirs(after_qc_seq_dir, exist_ok=True)
        for r1 in glob(r1_pattern):
            cmd = get_trimmomatic_cmd(R1=r1,
                                      R2=r1.replace('_1', '_2'),
                                      odir=after_qc_seq_dir,
                                      threads=30,
                                      split_str='.clean')
            run_cmd(cmd, dryrun=dryrun)
        #####################################################################
        # fastqc after file
        cmd = fastqc_cmd.format(exe_path=fastqc_path,
                                in_pattern=os.path.join(after_qc_seq_dir, '*clean.fq.gz'),
                                threads=40,
                                odir=qc_after)
        run_cmd(cmd, dryrun=dryrun)
        #####################################################################
        # merge all fastqc into two file
        outdir = os.path.join(base_odir, 'fastqc')
        cmd = multiqc_cmd.format(exe_path=multiqc_path,
                                 indir=qc_before,
                                 odir=outdir,
                                 fn='beforeQC')
        run_cmd(cmd, dryrun=dryrun)
        cmd = multiqc_cmd.format(exe_path=multiqc_path,
                                 indir=qc_after,
                                 odir=outdir,
                                 fn='afterQC')


        run_cmd(cmd, dryrun=dryrun)

        #####################################################################
        # accessment option #todo: split it into accessory function.
        def access_assembly(r1, r2, ref, gff, test_dir, dryrun=True):
            os.makedirs(test_dir, exist_ok=True)
            contigs_list = []
            for depth in [100, 200, 300]:
                odir = os.path.join(test_dir, 'shovill_%s' % depth)
                # os.makedirs(odir, exist_ok=True)
                if not os.path.isdir(odir):
                    contigs_list.append(os.path.join(odir, 'contigs.fasta'))
                    cmd = get_shovill_cmd(r1, r2, odir, depth, 30, "", "--minlen 500 --force")
                    run_cmd(cmd, dryrun=dryrun)
            odir = os.path.join(test_dir, 'full_spades_%s' % depth)
            # os.makedirs(odir, exist_ok=True)
            spades_path = "/tools/SPAdes-3.13.0-Linux/bin/spades.py"
            spades_cmd = """{spades} --pe1-1 {r1} --pe1-2 {r2} --threads {threads} --memory 150 -o {odir}"""
            cmd = spades_cmd.format(spades=spades_path,
                                    r1=r1,
                                    r2=r2,
                                    threads=35,
                                    odir=odir)
            if not os.path.isdir(odir):
                contigs_list.append(os.path.join(odir, 'scaffolds.fasta'))
                run_cmd(cmd, dryrun=dryrun)

            #####################################################################
            quast_dir = os.path.join(test_dir, 'quast_all')
            quast_path = "/home/liaoth/.local/bin/quast.py"
            quast_cmd = """{exe_path} {contig} -r {ref} -g {gff} --circos --gene-finding -1 {r1} -2 {r2} --threads {threads} --no-check -o {odir} """
            for contig_fn in contigs_list:
                dir_name = os.path.basename(os.path.dirname(contig_fn))
                cmd = quast_cmd.format(exe_path=quast_path,
                                       contig=contig_fn,
                                       r1=r1,
                                       r2=r1,
                                       ref=ref,
                                       gff=gff,
                                       threads=30,
                                       odir=os.path.join(quast_dir, dir_name))
                run_cmd(cmd, dryrun=dryrun)


        #####################################################################
        # regular shovill part and find plasmids
        shovill_output = os.path.join(base_odir, "shovill_output/")
        os.makedirs(shovill_output, exist_ok=True)

        for r1 in glob(os.path.join(after_qc_seq_dir, '*_1.clean.fq.gz')):
            r2 = r1.replace('_1.clean', "_2.clean")
            base_name = os.path.basename(r1)
            sample_name = base_name.split('_')[0].strip('S')
            plasmids_odir = os.path.join(shovill_output, 'plasmidsSpades', sample_name)
            regular_odir = os.path.join(shovill_output, 'regular', sample_name)
            # os.makedirs(regular_odir, exist_ok=True)
            # os.makedirs(plasmids_odir, exist_ok=True)
            cmd1 = get_shovill_cmd(r1, r2, plasmids_odir, 100, 30, " '--plasmid' ", "--nocorr")
            cmd2 = get_shovill_cmd(r1, r2, regular_odir, 100, 30, "", "--minlen 500")

            run_cmd(cmd1, dryrun=dryrun)
            run_cmd(cmd2, dryrun=dryrun)

            isolates_df = isolates_df.append(pd.DataFrame(data=[[sample_name,
                                                                 os.path.join(regular_odir, 'contigs.fasta'),
                                                                 r1,
                                                                 r2]], index=[0]))
            # if '41833' in r1:
            #     access_assembly(r1, r2,
            #                     ref='/home/liaoth/project/shenzhen_actineto/reference/GCF_000746645.1_ASM74664v1_genomic.fna.gz',
            #                     gff='/home/liaoth/project/shenzhen_actineto/reference/GCF_000746645.1_ASM74664v1_genomic.gff.gz',
            #                     test_dir='/home/liaoth/project/shenzhen_actineto/test_access/41833',
            #                     dryrun=True)
        #####################################################################
        # ln source
        all_contigs_dir = os.path.join(base_odir, 'ready_prokka')
        os.makedirs(all_contigs_dir, exist_ok=True)
        for contig in glob(os.path.join(shovill_output, 'regular', '*', 'contigs.fa')):
            # focus: the name of contig must short......using 'contigs.fasta' will raise error.
            cmd = 'ln -s {source} {target}'
            sample_name = os.path.basename(os.path.dirname(contig))
            target = os.path.join(all_contigs_dir, sample_name + '.fasta')
            if os.path.isfile(target):
                os.remove(target)

            cmd = cmd.format(source=contig,
                             target=target)
            run_cmd(cmd, dryrun=dryrun)
        #####################################################################
        # prokka
        prokka_dir = os.path.join(base_odir, "prokka_o")
        for contig in glob(os.path.join(all_contigs_dir, '*.fasta')):
            sample_name = os.path.basename(contig).split('.')[0]
            odir = os.path.join(prokka_dir, sample_name)
            os.makedirs(odir, exist_ok=True)
            cmd = prokka_cmd.format(exe_path=prokka_path,
                                    in_fasta=contig,
                                    odir=odir,
                                    threads=35,
                                    sample_name=sample_name)
            if sample_name not in list(isolates_df.iloc[:, 0]):
                isolates_df = isolates_df.append(pd.DataFrame(data=[[sample_name,
                                                                     contig,
                                                                     '',
                                                                     '']], index=[0]))
            if not os.path.isfile(os.path.join(odir,sample_name+'.gff')):
                run_cmd(cmd, dryrun=dryrun)
        #####################################################################
        # roary
        roary_dir = os.path.join(base_odir, "roary_o")
        cmd = roary_cmd.format(exe_path=roary_path,
                               threads=35,
                               odir=roary_dir,
                               gff_pattern=os.path.join(prokka_dir, '*', '*.gff'))
        run_cmd(cmd, dryrun=dryrun)
        #####################################################################
        # pandoo
        pandoo_input = os.path.join(base_odir, 'pandoo_o', 'isolates.tab')
        pandoo_output_dir = os.path.join(base_odir, 'pandoo_o', 'result')
        os.makedirs(pandoo_output_dir, exist_ok=True)
        isolates_df.to_csv(pandoo_input, index=False, header=False, sep='\t')
        cmd = pandoo_cmd.format(exe_path=pandoo_path,
                                ariba_str=ariba_str,
                                abricate_str=abricate_str,
                                input=pandoo_input,
                                odir=pandoo_output_dir)
        run_cmd(cmd, dryrun=dryrun)
        ####################################################################
        # fasttree
        fasttree_path = "/usr/bin/fasttreeMP"
        fasttree_cmd = "{exe_path} -nt -gtr {roary_dir}/core_gene_alignment.aln > {roary_dir}/core_gene.newick"
        cmd = fasttree_cmd.format(exe_path = fastqc_path,
                                  roary_dir=roary_dir)
        run_cmd(cmd,dryrun=dryrun)
        ############################################################
        #
