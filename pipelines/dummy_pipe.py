#################################################################################
####  it should sequential concat all tasks.
####  it could be taken as test module.
####  obsolete file (190716)
#################################################################################
from pipelines.tasks import *


def main(indir, odir, dry_run=False):
    log_file = os.path.join(odir, "pipelines.log")
    log_file_stream = open(log_file, 'w')
    ############################################################
    fastqc_o_dir = os.path.join(odir,
                                "fastqc_before")
    fastqc_in_pattern = os.path.join(indir,
                                     '*',
                                     "*.fq.gz")
    run_fastqc(in_files=glob(fastqc_in_pattern),
               odir=fastqc_o_dir,
               dry_run=dry_run,
               log_file=log_file_stream)
    run_multiqc(in_dir=fastqc_o_dir,
                odir=fastqc_o_dir,
                fn=os.path.basename(fastqc_o_dir),
                dry_run=dry_run,
                log_file=log_file_stream)
    ############################################################
    R1s = [_ for _ in glob(os.path.join(indir,
                                        '*',
                                        "*.fq.gz")) if '_1' in _]
    R2s = [_.replace('_1.fq.gz', '_2.fq.gz') for _ in R1s]
    sample_names = [str(os.path.basename(_).split('_')[0]).strip('S') for _ in R1s]
    for R1, R2, sn in zip(R1s, R2s, sample_names):
        run_trimmomatic(R1=R1,
                        R2=R2,
                        odir=os.path.join(odir, "cleandata"),
                        sample_name=sn,
                        thread=7,
                        dry_run=dry_run,
                        log_file=log_file_stream
                        )
    ############################################################
    fastqc_o_dir = os.path.join(odir,
                                "fastqc_aft")
    fastqc_in_pattern = os.path.join(odir,
                                     "cleandata",
                                     "*.clean.fq.gz")
    run_fastqc(in_files=glob(fastqc_in_pattern),
               odir=fastqc_o_dir,
               dry_run=dry_run,
               log_file=log_file_stream)
    run_multiqc(in_dir=fastqc_o_dir,
                odir=fastqc_o_dir,
                fn=os.path.basename(fastqc_o_dir),
                dry_run=dry_run,
                log_file=log_file_stream)
    ############################################################
    # regular shovill
    sample_names = [str(os.path.basename(_).split('_')[0]).strip('S') for _ in R1s]
    for sn in sample_names:
        R1 = os.path.join(odir, "cleandata", sn + '_R1.clean.fq.gz')
        R2 = os.path.join(odir, "cleandata", sn + '_R2.clean.fq.gz')
        run_shovill(R1=R1,
                    R2=R2,
                    odir=os.path.join(odir, "assembly_o", "regular", sn),
                    thread=7,
                    ram=15,
                    dry_run=dry_run,
                    log_file=log_file_stream
                    )
    # plasmid shovill
    for sn in sample_names:
        R1 = os.path.join(odir, "cleandata", sn + '_R1.clean.fq.gz')
        R2 = os.path.join(odir, "cleandata", sn + '_R2.clean.fq.gz')
        run_shovill(R1=R1,
                    R2=R2,
                    odir=os.path.join(odir, "assembly_o", "plasmids", sn),
                    thread=7,
                    ram=15,
                    depth=10,
                    dry_run=dry_run,
                    log_file=log_file_stream
                    )
    ############################################################
    ref = '/home/liaoth/data2/project/shenzhen_Acinetobacter/reference/baumannii ATCC 19606.fna'
    gff = '/home/liaoth/data2/project/shenzhen_Acinetobacter/reference/baumannii ATCC 19606.gff'
    for sn in sample_names:
        R1 = os.path.join(odir, "cleandata", sn + '_R1.clean.fq.gz')
        R2 = os.path.join(odir, "cleandata", sn + '_R2.clean.fq.gz')
        prokka_in_file = os.path.join(odir, "assembly_o", "regular", sn, "contigs.fasta")
        run_quast(contig=prokka_in_file,
                  R1=R1,
                  R2=R2,
                  ref=ref,
                  gff=gff,
                  odir=os.path.join(odir, 'assembly_o', "regular_quast", sn),
                  thread=7,
                  dry_run=dry_run,
                  log_file=log_file_stream)
    ############################################################
    # prokka_pattern = os.path.join(odir, "assembly_o", "regular", "*", "contigs.fasta")
    for sn in sample_names:
        prokka_in_file = os.path.join(odir, "assembly_o", "regular", sn, "contigs.fasta")
        run_prokka(infile=prokka_in_file,
                   odir=os.path.join(odir, "prokka_o", sn),
                   dry_run=dry_run,
                   log_file=log_file_stream)
    ############################################################
    run_roary(os.path.join(odir, "prokka_o"),
              os.path.join(odir, "all_roary_o"),
              thread=7,
              dry_run=dry_run,
              log_file=log_file_stream)
    ############################################################
    run_fasttree(os.path.join(odir, "all_roary_o", "core_gene_alignment.aln"),
                 os.path.join(odir, "all_roary_o", "core_gene.newick"),
                 dry_run=dry_run,
                 log_file=log_file_stream)
    ############################################################
    run_abricate(os.path.join(odir, "prokka_o"),
                 roary_dir=os.path.join(odir, "all_roary_o"),
                 odir=os.path.join(odir, "abricate_result"),
                 thread=7,
                 mincov=80,
                 dry_run=dry_run,
                 log_file=log_file_stream)
    from toolkit.utils import construct_pandoo_table
    pandoo_infile = os.path.join(odir, "pandoo_input.tab")
    isolates_table = construct_pandoo_table(sample_names,
                                            odir)
    with open(pandoo_infile, 'w') as f1:
        isolates_table.to_csv(f1, sep='\t', index=False, header=None)
    run_pandoo(in_file=pandoo_infile,
               odir=os.path.join(odir, "pandoo_o"),
               thread=7,
               dry_run=dry_run,
               log_file=log_file_stream)


if __name__ == '__main__':
    indir = "/home/liaoth/data2/pangenome/shenzhen_Acinetobacter"
    odir = "/home/liaoth/project/genome_pipelines/pipelines/test/test_output"
    import sys
    if sys.argv[1] == 'check':
        check_exe()
    elif sys.argv[1] == 'test':
        main(indir, odir, dry_run=True)
    elif sys.argv[1] == 'run':
        main(indir, odir, dry_run=False)



