import os
import sys
from glob import glob
import pandas as pd


from toolkit.utils import run_cmd, valid_path
from .constant_str import *


def check_exe():
    def check_exists(exe_file):
        if os.path.isfile(exe_file):
            return exe_file
        else:
            return

    check_list = dict(fastqc=fastqc_path,
                      multiqc=multiqc_path,
                      trimmomatic=trimmomatic_jar_path,
                      shovill=shovill_path,
                      prokka=prokka_path,
                      roary=roary_path,
                      quast=quast_path,
                      pandoo=pandoo_path,
                      fasttree=fasttree_path,
                      isescan=ISEscan_path,
                      abricate=abricate_path,
                      abricate_py=abricate_py_path,
                      )
    # todo: uncompleted check need to check it useful or not
    for software in sorted(check_list):
        path = check_exists(check_list[software])

        if path is not None:
            print(software.ljust(11) + ':\tok\t' + path.ljust(11), file=sys.stderr)
        else:
            print('Dependency ' + software + ' is not existed.  Please ' + \
                  'check constrant_str',
                  file=sys.stderr)


def auto_update_exe():
    # todo: use which or others to auto implement exe
    check_list = dict(fastqc=fastqc_path,
                      multiqc=multiqc_path,
                      trimmomatic="java -jar " + trimmomatic_jar_path,
                      shovill=shovill_path,
                      prokka=prokka_path,
                      roary=roary_path,
                      quast=quast_path,
                      pandoo=pandoo_path,
                      fasttree=fasttree_path,
                      isescan=ISEscan_path,
                      abricate=abricate_path,
                      abricate_py=abricate_py_path,
                      )

    pass


def run_fastqc(in_files,
               odir,
               exe_path=fastqc_path,
               dry_run=False,
               log_file=None):
    if not dry_run:
        valid_path(in_files, check_size=True)
    valid_path(odir, check_odir=True)
    cmd = fastqc_cmd.format(exe_path=exe_path,
                            in_files=' '.join(in_files),
                            odir=odir)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_multiqc(in_dir,
                odir,
                fn,
                extra_str='',
                exe_path=multiqc_path,
                dry_run=False,
                log_file=None):
    if not dry_run:
        valid_path(in_dir, check_dir=True)
    valid_path(odir, check_odir=True)
    cmd = multiqc_cmd.format(exe_path=exe_path,
                             indir=in_dir,
                             odir=odir,
                             fn=fn,
                             extra_str=extra_str)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_trimmomatic(R1, R2, odir,
                    sample_name,
                    thread=0,
                    exe_path=trimmomatic_jar_path,
                    dry_run=False,
                    log_file=None):
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if not dry_run:
        valid_path([R1, R2], check_size=True)
    valid_path(odir, check_odir=True)
    sample_name = str(sample_name)
    clean_r1 = os.path.join(odir, sample_name + '_R1.clean.fq.gz')
    unpaired_r1 = os.path.join(odir, sample_name + '_R1.unpaired.fq.gz')
    clean_r2 = os.path.join(odir, sample_name + '_R2.clean.fq.gz')
    unpaired_r2 = os.path.join(odir, sample_name + '_R2.unpaired.fq.gz')
    log = os.path.join(odir, sample_name + '.log')
    cmd = trimmomatic_cmd.format(exe_path=exe_path,
                                 threads=thread,
                                 R1=R1,
                                 R2=R2,
                                 log=log,
                                 clean_r1=clean_r1,
                                 unpaired_r1=unpaired_r1,
                                 clean_r2=clean_r2,
                                 unpaired_r2=unpaired_r2,
                                 params=trimmomatic_setting)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_shovill(R1, R2, odir,
                spades_extra_options=None,
                thread=0,
                extra_option='',
                ram=100,  # unit G
                depth=100,
                exe_path=shovill_path,
                dry_run=False,
                log_file=None):
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if spades_extra_options is not None:
        extra_str = ' --opts ' + spades_extra_options
    else:
        extra_str = ''
    if not dry_run:
        valid_path([R1, R2], check_size=True)
    valid_path(odir, check_odir=True)
    cmd = shovill_cmd.format(exe_path=exe_path,
                             R1=R1,
                             R2=R2,
                             odir=odir,
                             depth=depth,
                             thread=thread,
                             ram=ram)
    cmd = ' '.join([cmd, extra_option, extra_str])
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_prokka(infile, odir,
               dry_run=False,
               log_file=None):
    if not dry_run:
        valid_path(infile, check_size=True)
    valid_path(odir, check_odir=True)
    sample_name = os.path.basename(odir)
    cmd = prokka_cmd.format(exe_path=prokka_path,
                            infile=infile,
                            odir=odir,
                            sn=sample_name,
                            )
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_roary(indir, odir,
              thread=0,
              dry_run=False,
              log_file=None):
    if glob(os.path.join(indir, '*.gff')):
        gff_pattern = os.path.join(indir, '*.gff')
    else:
        gff_pattern = os.path.join(indir, '*', '*.gff')
    if not dry_run:
        valid_path(gff_pattern, check_glob=True)
    valid_path(os.path.dirname(odir), check_odir=True)
    # roary would not overwrite existing directory... so we don't create it directly, we just check the parent directory.
    cmd = roary_cmd.format(exe_path=roary_path,
                           thread=thread,
                           odir=odir,
                           gff_pattern=gff_pattern)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_abricate(indir,
                 roary_dir,
                 odir,
                 thread=0,
                 mincov=80,
                 dry_run=False,
                 log_file=None):
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if not dry_run:
        valid_path(indir, check_dir=True)
    valid_path(odir, check_odir=True)
    inside_log_file = os.path.join(odir, 'abricate.log')
    extra_str = " --log %s" % inside_log_file
    if dry_run:
        extra_str += ' --dry_run'
    cmd = abricate_cmd.format(py_path=abricate_py_path,
                              exe_path=abricate_path,
                              indir=indir,
                              roary_dir=roary_dir,
                              odir=odir,
                              db=abricate_db,
                              mincov=mincov,
                              thread=thread,
                              extra_str=extra_str)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_pandoo(in_file,
               odir,
               thread=0,
               dry_run=False,
               log_file=None):
    """ run without abricate"""
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if not dry_run:
        valid_path(odir, check_odir=True)
    valid_path(in_file, check_size=True)
    cmd = pandoo_cmd.format(exe_path=pandoo_path,
                            ariba_str=ariba_str,
                            input=in_file,  # "isolates_pre.tab",
                            thread=thread,
                            odir=odir)
    cmd += """ -a "{}" """
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_fasttree(in_file,
                 ofile,
                 dry_run=False,
                 log_file=None
                 ):
    if not dry_run:
        valid_path(in_file, check_size=True)
    cmd = fasttree_cmd.format(exe_path=fasttree_path,
                              in_aln=in_file,
                              o_newick=ofile)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_quast(contig,
              R1,
              R2,
              odir,
              thread=0,
              ref=None,
              gff=None,
              dry_run=False,
              log_file=None):
    """ not loop inside a function in order to parallel"""
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if not dry_run:
        valid_path([contig, R1, R2, ref, gff], check_size=True)
    valid_path(odir, check_odir=True)
    extra_str = ''
    if ref is not None:
        extra_str = """ -r "%s" """ % ref
    if gff is not None:
        extra_str += """ -g "%s" """ % gff
    cmd = quast_cmd.format(exe_path=quast_path,
                           contig=contig,
                           R1=R1,
                           R2=R2,
                           extra_str=extra_str,
                           threads=thread,
                           odir=odir)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_ISEscan(infile,
                odir,
                sample_name,
                dry_run=False,
                log_file=None):
    if not dry_run:
        valid_path(infile, check_size=True)
    valid_path(odir, check_odir=True)
    cmd = isescan_cmd.format(exe_path=ISEscan_path,
                             infile=infile,
                             proteome_dir=os.path.join(odir, "proteome"),
                             hmm_dir=os.path.join(odir, "hmm"),
                             odir=odir,
                             sn=sample_name)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_plasmid_detect(indir,
                       roary_dir,
                       prokka_dir,
                       odir,
                       dry_run=False,
                       log_file=None):
    if not dry_run:
        valid_path(indir, check_dir=True)
    valid_path(odir, check_odir=True)

    cmd = plasmid_detect_cmd.format(py_path=plasmid_detect_path,
                                    indir=indir,
                                    roary_dir=roary_dir,
                                    prokka_dir=prokka_dir,
                                    odir=odir)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


# todo: add more post-analysis to it.
# checkM (completeness of assembly data) https://github.com/Ecogenomics/CheckM
# popins (new IS) https://github.com/bkehr/popins

# todo: add auto generated heatmap for this pipelines.

def run_phigaro(infile,
                ofile,
                thread=0,
                dry_run=False,
                log_file=None):
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if not dry_run:
        valid_path(infile, check_size=1)
    valid_path(ofile, check_ofile=1)
    cmd = phigaro_cmd.format(exe_path=phigaro_path,
                             infile=infile,
                             ofile=ofile,
                             thread=thread,)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


def run_gubbins(infile,
                oprefix,
                thread=0,
                dry_run=False,
                log_file=None):
    if thread == 0 or thread == -1:
        thread = cpu_count()
    if not dry_run:
        valid_path(infile, check_size=1)
    valid_path(oprefix, check_ofile=1)

    cmd = gubbins_cmd.format(exe_path=gubbins_path,
                             infile=infile,
                             oprefix=oprefix,
                             thread=thread)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)

############################################################
def post_analysis(workflow_task):
    # only ofr workflow post analysis
    odir = workflow_task.odir
    summary_odir = os.path.join(odir,"pipelines_summary")
    new_gff_odir = os.path.join(odir,"pipelines_summary","new_gff")
    valid_path(new_gff_odir,check_odir=1)
    ############################################################
    # mlst
    pandoo_ofile = workflow_task.output()["pandoo"].path
    pandoo_df = pd.read_csv(pandoo_ofile, index_col=0, )
    mlst_df = pandoo_df.loc[:, pandoo_df.columns.str.contains("MLST")]
    from toolkit.process_mlst import main as process_mlst
    output_mlst_df = process_mlst(mlst_df)
    for scheme,mlst_df in output_mlst_df:
        with open(os.path.join(summary_odir,"%s_mlst.csv" % scheme),'w') as f1:
            mlst_df.to_csv(f1, index=1)
    ############################################################
    # abricate
    abricate_ofile = workflow_task.output()["abricate"].path
    abricate_dir = os.path.dirname(abricate_ofile)
    new_abricate_dir = os.path.join(new_gff_odir,"abricate_annotated_gff")
    valid_path(new_abricate_dir,check_odir=1)

    os.system("cp %s %s" % (abricate_ofile,
                            summary_odir))
    os.system("cp %s %s" % (os.path.join(abricate_dir,"samples2annotate.csv"),
                            summary_odir))
    os.system("cp %s %s" % (os.path.join(abricate_dir,"*","*.gff"),
                            new_abricate_dir))  # todo: not yet
    ############################################################
    # fasttree
    core_gene_tree = workflow_task.output()["fasttree"].path
    os.system("cp %s %s" % (core_gene_tree,
                            summary_odir))
    ############################################################
    # plasmid
    plasmid_summary_file = workflow_task.output()["detect_plasmid"].path
    plasmid_dir = os.path.dirname(plasmid_summary_file)
    os.system("cp %s %s" % (os.path.join(plasmid_dir,"fullWithAnnotated"),
                            new_gff_odir))
    ############################################################
    # IS
    IS_summary_file = workflow_task.output()["ISEscan_summary"].path
    IS_dir = os.path.dirname(IS_summary_file)
    new_IS_dir = os.path.join(new_gff_odir, "IS_gff")
    valid_path(new_IS_dir, check_odir=1)
    os.system("cp %s %s" % (os.path.join(IS_dir,"*","*.gff"),
                            new_IS_dir))
    ############################################################
    #



############################################################
def access_assembly(r1, r2, ref, gff, test_dir, dryrun=True):
    os.makedirs(test_dir, exist_ok=True)
    contigs_list = []
    for depth in [100, 200, 300]:
        odir = os.path.join(test_dir, 'shovill_%s' % depth)
        # os.makedirs(odir, exist_ok=True)
        if not os.path.isdir(odir):
            contigs_list.append(os.path.join(odir, 'contigs.fasta'))
            # cmd = get_shovill_cmd(r1, r2, odir, depth, 30, "", "--minlen 500 --force")
            # run_cmd(cmd, dryrun=dryrun)
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
