import copy
import os
import sys
from glob import glob

import pandas as pd
from BCBio import GFF

from toolkit.process_region_annotated import *
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
                       ofile,
                       dry_run=False,
                       log_file=None):
    if not dry_run:
        valid_path(indir, check_dir=1)
    valid_path(ofile, check_ofile=1)

    cmd = plasmid_detect_cmd.format(py_path=plasmid_detect_path,
                                    indir=indir,
                                    ofile=ofile)
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)


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
                             thread=thread, )
    run_cmd(cmd, dry_run=dry_run, log_file=log_file)

# todo: add more post-analysis to it.
# checkM (completeness of assembly data) https://github.com/Ecogenomics/CheckM
# popins (new IS) https://github.com/bkehr/popins

# todo: add auto generated heatmap for this pipelines.

def run_gubbins(infile,
                oprefix,
                thread=0,
                dry_run=False,
                log_file=None):
    # just like roary
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


def run_seqtk(infile, outfile, isolate):
    '''
    # fixme:
    Runs SeqTK on the contigs.fa (assemblies).  Gathers the fasta metrics.
    '''
    if len(infile) == 0:
        metrics = {}
    else:
        os.system('seqtk comp ' + ''.join(infile) + ' > ' + outfile)
        df1 = pd.read_csv(outfile, header=None, index_col=0, sep='\t',
                          names=['chr', 'length', '#A', '#C', '#G', '#T', '#2',
                                 '#3', '#4', '#CpG', '#tv', '#ts', '#CpG-ts'])
        contig_lengths = df1['length'].tolist()
        contig_lengths.sort(), contig_lengths.reverse()
        contig_lengths_prime = [[i] * i for i in contig_lengths]
        metrics = {}
        pfx = 'metricsContigs_'
        metrics[pfx + 'N50'] = int(np.median([i for j in contig_lengths_prime
                                              for i in j]))
        bps = sum(contig_lengths)
        nns = sum(df1['#4'].tolist())
        metrics[pfx + 'Ns'] = nns
        metrics[pfx + 'no'] = len(df1.index.values)
        metrics[pfx + 'bp'] = bps
        metrics[pfx + 'avg'] = bps / len(contig_lengths)
        metrics[pfx + 'Max'] = max(contig_lengths)
        metrics[pfx + 'Min'] = min(contig_lengths)
        metrics[pfx + 'ok'] = bps - nns
    metrics['sotwareSeqTKversion_comp'] = seqtk_version()
    seqtk_comp_df = create_pandas_df(metrics, isolate)
    # NB: this will overwrite the outfile that was read in at the first step.
    write_pandas_df(outfile, seqtk_comp_df)

def run_seqtk2(infile, outfile, isolate):
    '''
    fixme
    Runs SeqTK fqchk on the reads.  Gathers fastq metrics.
    '''
    pfx = 'metricsReads_'
    if len(infiles) == 0:
        metrics = {}
    else:
        os.system('seqtk fqchk ' + ' '.join(infiles) + ' > ' + outfile)
        df_all = pd.read_csv(outfile, index_col=0, sep='\t',
                             skiprows=1)
        with open(outfile, 'r') as outf:
            metrics = outf.readline()
            metrics = dict([j.split(':') for j in [i.replace(' ', '') for i in
                                                   metrics.split(';')[0:3]]])
            for key, value in list(metrics.items()):
                metrics[pfx + key] = metrics.pop(key)
        metrics[pfx + 'AvgQual'] = round(df_all.loc['ALL', 'avgQ'], 2)
        metrics[pfx + 'GeeCee'] = round(df_all.loc['ALL', '%C'] +
                                        df_all.loc['ALL', '%G'], 2)
        metrics[pfx + 'Yield'] = int(df_all.loc['ALL', '#bases'] * 2)

        # Count the number of reads in the infiles.
        cmd = 'zcat ' + ' '.join(infiles)
        cmd2 = 'wc -l'
        args = shlex.split(cmd)
        args2 = shlex.split(cmd2)
        proc = Popen(args, stdout=PIPE, stderr=PIPE)
        proc2 = Popen(args2, stdin=proc.stdout, stdout=PIPE)
        output = proc2.communicate()[0].decode('UTF-8')
        metrics[pfx + 'nReads'] = int(int(output) / 4)

        # Get the mode read length.
        n_bases = df_all['#bases'].values[1:]
        lengths = []
        pos = 0
        while pos < len(n_bases) - 1:
            lengths.append(n_bases[pos] - n_bases[pos + 1])
            pos += 1
        else:
            lengths.append(n_bases[pos])
        metrics[pfx + 'ModeLen'] = lengths.index(max(lengths)) + 1
    # Add the version, create the pandas object, write it to file.
    metrics['softwareSeqTKversion_fqchk'] = seqtk_version()
    metrics_df = create_pandas_df(metrics, isolate)
    # NB: this will overwrite the outfile that was read in at the first step.
    write_pandas_df(outfile, metrics_df)



def run_kraken(infile, outfile, fmt, isolate, dbase, threads):
    '''
    Run Kraken on the infile.
    '''

    def do_kraken(cmd_kraken):
        '''
        This function exists so do_kraken() can be bypassed if there are no
        infiles.
        '''
        #         cmd_krk_r = 'kraken-report --db '+dbase
        cmd_grep = "grep -P '\\tS\\t'"
        cmd_sed = "sed -e 's/%//g'"
        cmd_sort = 'sort -k 1 -g -r'
        cmd_head = 'head -3'
        cmd_full = f'{cmd_kraken} | {cmd_grep} | {cmd_sed} | {cmd_sort} | {cmd_head}'
        sys.stderr.write(cmd_full + '\n')

        output = check_output(cmd_full, shell=True)
        output2 = output.decode('UTF-8').split('\n')
        kraken = [line.strip().split('\t') for line
                  in [_f for _f in output2 if _f]]
        print(f"{isolate} {kraken}")
        return kraken

    if fmt == 'reads':
        if len(infile) == 2:
            infiles = ' '.join(infile)
            compression = ''
            for i in infile:
                # Compression test based on file extension using linux 'file'.
                args = shlex.split('file ' + i)
                proc = Popen(args, stdout=PIPE)
                f_fmt = proc.communicate()[0].decode('UTF-8').rstrip().split()
                if 'gzip' in f_fmt:
                    compression = '--gzip-compressed '
                    break
                if 'bzip2' in f_fmt:
                    compression = '--bzip2-compressed '
                    break
            cmd_kraken = f"kraken2 --threads {threads} --db {dbase} {compression} --paired --report {outfile} --output - --memory-mapping {infiles} && cat {outfile}"
            #             cmd_kraken = 'kraken --threads '+str(threads)+' --db '+dbase +\
            #                          ' --fastq-input '+compression +\
            #                          '--paired --check-names '+infiles
            kraken = do_kraken(cmd_kraken)
        else:
            # If no read pairs in list, kraken is an empty list.
            kraken = []

    if fmt == 'contigs':
        if len(infile) == 1:
            infile = ''.join(infile)
            cmd_kraken = f"kraken2 --threads {threads} --db {dbase} --report {outfile} --output - --memory-mapping {infile} && cat {outfile}"
            #
            #             cmd_kraken = 'kraken --threads '+str(threads)+' --db '+dbase +\
            #                          ' --fasta-input '+infile
            kraken = do_kraken(cmd_kraken)
        else:
            kraken = []

    def kraken_results_df_creator(kraken_hits, fmt):
        '''
        Take the 2D array from a kraken search and return result as a
        dictionary.
        '''
        dict_hits = {}
        k = 1
        for i in range(0, len(kraken_hits)):
            dict_hits['Sp_krkn_' + fmt + '_' + str(k)] = \
                kraken_hits[i][5].lstrip()
            dict_hits['Sp_krkn_' + fmt + '_' + str(k) + '_pc'] = kraken_hits[i][0]
            k += 1
        return dict_hits

    if len(kraken) > 0:
        krk_df = kraken_results_df_creator(kraken, fmt)
    else:
        krk_df = {}

    def get_krkn_version():
        '''
        Get the Kraken software version and Kraken Database path.
        '''
        args = shlex.split('kraken2 -v')
        proc = Popen(args, stdout=PIPE)
        version = proc.communicate()[0].decode('UTF-8').rstrip().split('\n')[0]
        return {'softwareKrakenVersion_' + fmt: version,
                'softwareKrakenDB_' + fmt: dbase}

    krk_vers_no = create_pandas_df(get_krkn_version(), isolate)
    krk_result = create_pandas_df(krk_df, isolate)
    kraken_df = pd.concat([krk_result, krk_vers_no], axis=1)
    write_pandas_df(outfile, kraken_df)



def run_mlst(assembly, outfile, isolate, species):
    '''
    Run Torsten's MLST program on the assembly.
    '''

    def parse_MLST_output(output):
        ncol = len(output)
        mlst_formatted_dict = {'MLST_Scheme': output[0],
                               'MLST_ST': output[1]}
        k = 1
        for i in range(2, ncol):
            mlst_formatted_dict['MLST_Locus' + str(k)] = output[i]
            k += 1
        return mlst_formatted_dict

    mlst_df_collect = []
    if len(assembly) == 0:
        mlst_formatted_dict = {}
        mlst_df = create_pandas_df(mlst_formatted_dict, isolate)
        mlst_df_collect.append(mlst_df)

    if len(assembly) == 1:
        assembly = assembly[0]
        sp_scheme = None
        if isinstance(species, str):
            if species in FORCE_MLST_SCHEME:
                sp_scheme = FORCE_MLST_SCHEME[species]
            elif species.split(' ')[0] in FORCE_MLST_SCHEME:
                sp_scheme = FORCE_MLST_SCHEME[species.split(' ')[0]]
            else:
                sp_scheme = FORCE_MLST_SCHEME["Acinetobacter baumannii"]
        if sp_scheme is not None:
            if type(sp_scheme) == str:
                sp_scheme = [sp_scheme]

            for _sp_scheme in sp_scheme:
                cmd = 'mlst --scheme ' + _sp_scheme + ' --quiet ' + \
                      assembly
                args_mlst = shlex.split(cmd)
                proc = Popen(args_mlst, stdout=PIPE)
                output = proc.communicate()[0].decode('UTF-8')
                out = output.rstrip().split('\t')[1:]
                mlst_formatted_dict = parse_MLST_output(out)
                mlst_df = create_pandas_df(mlst_formatted_dict, isolate)
                mlst_df_collect.append(mlst_df)
        else:
            cmd = 'mlst --quiet ' + assembly
            args_mlst = shlex.split(cmd)
            proc = Popen(args_mlst, stdout=PIPE)
            output = proc.communicate()[0].decode('UTF-8')
            out = output.rstrip().split('\t')[1:]
            mlst_formatted_dict = parse_MLST_output(out)
            mlst_df = create_pandas_df(mlst_formatted_dict, isolate)
            mlst_df_collect.append(mlst_df)

    def get_mlst_version():
        '''
        Get the MLST software version.
        '''
        args = shlex.split('mlst -v')
        proc = Popen(args, stdout=PIPE)
        version = proc.communicate()[0].decode('UTF-8').rstrip().split('\t')[0]
        return {'softwareMLSTversion': version}

    mlst_version = create_pandas_df(get_mlst_version(), isolate)
    mlst_df = pd.concat(mlst_df_collect + [mlst_version], axis=1)
    write_pandas_df(outfile, mlst_df)


############################################################

def post_analysis(workflow_task):
    # only ofr workflow post analysis
    roary_dir = os.path.dirname(workflow_task.input()["fasttree"].path)
    prokka_o = os.path.join(workflow_task.odir,
                            'prokka_o')
    abricate_file = workflow_task.input()["abricate"].path

    summary_odir = workflow_task.output().path
    valid_path(summary_odir, check_odir=1)
    ############################################################
    # For summary the region based annotated result
    # 1. write info to empty gff for each sample
    # 2. extract each region with annotated info for each sample
    # 3. summary into matrix
    # 4. summary statistic info
    # prepare the accessory obj
    locus2group, locus2annotate, sample2gff = get_accessory_obj(roary_dir,
                                                                abricate_file,
                                                                prokka_o)
    empty_sample2gff = {sn: vals[1] for sn, vals in sample2gff.items()}
    ori_sample2gff = {sn: vals[2] for sn, vals in sample2gff.items()}

    merged_locus2annotate = locus2group.copy()
    merged_locus2annotate.update(locus2annotate)
    # merged together
    ############
    summary_task_tags = ["detect_prophage",
                         "detect_plasmid",
                         "ISEscan_summary"]
    summary_task_source = ["phigaro",
                           "plasmidSpades+bwa",
                           "isescan"]
    names = ["Prophage", "Plasmid", "IS"]
    annotated_sample2gff = copy.deepcopy(ori_sample2gff)
    for record in [record
                   for contig2record in annotated_sample2gff.values()
                   for record in contig2record.values()]:
        for fea in record.features:
            if fea.type == 'CDS':
                locus_id = fea.id
                annotated = merged_locus2annotate.get(locus_id, locus_id)
                fea.qualifiers["ID"] = fea.qualifiers['locus_tag'] = [annotated]
                fea.id = annotated
    ###########
    for tag, source, name in zip(summary_task_tags,
                                 summary_task_source,
                                 names):
        pth = workflow_task.input()[tag].path

        annotated_sample2gff_records = write_new_gff(pth,
                                                     empty_sample2gff,
                                                     source=source)
        subset_sample2gff_records = cut_old_gff(pth,
                                                annotated_sample2gff,
                                                source=source)
        samples2annotated_df = summary_into_matrix(subset_sample2gff_records,
                                                   unique_by=None)
        summary_df = summary_statistic(ori_sample2gff,
                                       subset_sample2gff_records,
                                       name
                                       )
        # output
        full_gff_with_region_dir = os.path.join(summary_odir,
                                                "gff_with_%s" % name)
        subset_gff_with_annotated_dir = os.path.join(summary_odir,
                                                     "gff_of_%s" % name)
        annotated_gff_odir = os.path.join(summary_odir,
                                          "annotated_gff")
        valid_path([full_gff_with_region_dir,
                    subset_gff_with_annotated_dir,
                    annotated_gff_odir
                    ], check_odir=1)
        for sn, records in annotated_sample2gff_records.items():
            filename = "%s.gff" % sn
            with open(os.path.join(full_gff_with_region_dir,
                                   filename), 'w') as f1:
                GFF.write(records, f1)
            with open(os.path.join(subset_gff_with_annotated_dir,
                                   filename), 'w') as f1:
                GFF.write(subset_sample2gff_records[sn], f1)
            with open(os.path.join(annotated_gff_odir,
                                   filename), 'w') as f1:
                GFF.write(list(annotated_sample2gff[sn].values()), f1)
        with open(os.path.join(summary_odir, "%s_annotated_matrix.csv" % name), 'w') as f1:
            samples2annotated_df.to_csv(f1, sep=',', index=1)
        with open(os.path.join(summary_odir, "%s_statistic.csv" % name), 'w') as f1:
            summary_df.to_csv(f1, sep=',', index=1)

    ############################################################
    # mlst
    pandoo_ofile = workflow_task.input()["pandoo"].path
    pandoo_df = pd.read_csv(pandoo_ofile, index_col=0, )
    mlst_df = pandoo_df.loc[:, pandoo_df.columns.str.contains("MLST")]
    from toolkit.process_mlst import main as process_mlst
    output_mlst_df = process_mlst(mlst_df)
    for scheme, mlst_df in output_mlst_df.items():
        with open(os.path.join(summary_odir, "%s_mlst.csv" % scheme), 'w') as f1:
            mlst_df.to_csv(f1, index=1)
    ############################################################
    # abricate
    abricate_ofile = workflow_task.input()["abricate"].path
    abricate_dir = os.path.dirname(abricate_ofile)
    os.system("cp %s %s" % (abricate_ofile,
                            summary_odir))
    os.system("cp %s %s" % (os.path.join(abricate_dir,
                                         "samples2annotate.csv"),
                            summary_odir))
    abricate_gff_dir = os.path.join(summary_odir,
                                      "annotated_gff_simplified")
    valid_path([abricate_gff_dir], check_odir=1)
    os.system("cp %s %s" % (os.path.join(abricate_dir,
                                         "*",
                                         "*.gff"),
                            abricate_gff_dir))
    ############################################################
    # fasttree
    core_gene_tree = workflow_task.input()["fasttree"].path
    os.system("cp %s %s" % (core_gene_tree,
                            summary_odir))


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
