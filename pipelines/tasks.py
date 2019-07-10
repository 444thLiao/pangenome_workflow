from glob import glob

import numpy as np
from BCBio import GFF

from toolkit.get_version import *
from toolkit.process_region_annotated import *
from toolkit.utils import valid_path, write_pandas_df
from .constant_str import *
from .soft_db_path import *


def check_exe():
    # todo: use it to check the validation of application
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


def run_roary(indir,  # indir/ infiles
              odir,
              thread=0,
              dry_run=False,
              log_file=None):
    if type(indir) == str:
        if glob(os.path.join(indir, '*.gff')):
            gff_pattern = os.path.join(indir, '*.gff')
        else:
            gff_pattern = os.path.join(indir, '*', '*.gff')
        if not dry_run:
            valid_path(gff_pattern, check_glob=True)
    else:
        gff_pattern = ' '.join(indir)

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
        paths = [_
                 for _ in [contig, R1, R2, ref, gff]
                 if not pd.isna(_)]
        valid_path(paths, check_size=True)
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
                             thread=thread,
                             phigaro_config=phigaro_config)
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
    # todo:
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


def run_seqtk_contig(infile, outfile, isolate,
                     dry_run=False,
                     log_file=None):
    '''
    Runs SeqTK on the contigs.fa (assemblies).
    Gathers the fasta metrics.
    process only one infile each time.
    '''
    assert os.path.isfile(infile)
    seqtk = seqtk_path
    run_cmd(f'{seqtk} comp ' + infile + ' > ' + outfile,
            dry_run=dry_run,
            log_file=log_file)
    df1 = pd.read_csv(outfile,
                      header=None,
                      index_col=0,
                      sep='\t',
                      names=['chr', 'length', '#A', '#C', '#G', '#T', '#2',
                             '#3', '#4', '#CpG', '#tv', '#ts', '#CpG-ts'])
    contig_lengths = df1['length'].tolist()
    contig_lengths.sort(reverse=True)
    # 从大到小
    contig_lengths_prime = [[i] * i for i in contig_lengths]
    metrics = {}
    pfx = 'metricsContigs_'
    metrics[pfx + 'N50'] = int(np.median([i for j in contig_lengths_prime
                                          for i in j]))
    # confirmed...
    bps = sum(contig_lengths)
    nns = sum(df1['#4'].tolist())
    metrics[pfx + 'Ns'] = nns
    metrics[pfx + 'no'] = len(df1.index.values)
    metrics[pfx + 'bp'] = bps
    metrics[pfx + 'avg'] = bps / len(contig_lengths)
    metrics[pfx + 'Max'] = max(contig_lengths)
    metrics[pfx + 'Min'] = min(contig_lengths)
    metrics[pfx + 'lt1K'] = len([_ for _ in contig_lengths if _ > 1000])
    metrics[pfx + 'lt5K'] = len([_ for _ in contig_lengths if _ > 5000])
    metrics[pfx + 'ok'] = bps - nns

    metrics['sotwareSeqTKversion_comp'] = seqtk_version()
    seqtk_comp_df = pd.DataFrame.from_dict({isolate: metrics},
                                           orient='index')
    # become (1,8) df with isolate as only index
    write_pandas_df(outfile, seqtk_comp_df)


def run_seqtk_reads(infile, outfile, isolate,
                    dry_run=False,
                    log_file=None):
    '''
    Runs SeqTK fqchk on the reads.
    If pair end sequencing fastq, input one.
    Gathers fastq metrics.
    '''
    assert os.path.isfile(infile)
    outfile = os.path.realpath(outfile)
    pfx = 'metricsReads_'
    seqtk = seqtk_path
    run_cmd(f'{seqtk} fqchk ' + infile + ' > ' + outfile,
            dry_run=dry_run,
            log_file=log_file)
    df_all = pd.read_csv(outfile,
                         index_col=0,
                         sep='\t',
                         skiprows=1)

    with open(outfile, 'r') as outf:
        metrics = next(outf)
        # 获取首行
        metrics = dict([j.split(': ')
                        for j in [i.strip()
                                  for i in metrics.split(';')[0:3]]])
        for key, value in list(metrics.items()):
            metrics[pfx + key] = metrics.pop(key)
        # 更新字典
    metrics[pfx + 'AvgQual'] = round(df_all.loc['ALL', 'avgQ'], 2)
    metrics[pfx + 'GeeCee'] = round(df_all.loc['ALL', '%C'] +
                                    df_all.loc['ALL', '%G'], 2)
    metrics[pfx + 'Yield'] = int(df_all.loc['ALL', '#bases'] * 2)

    # Count the number of reads in the infiles.
    cmd = "zgrep -c '^+$' " + infile
    output = run_cmd(cmd,
                     get_output=True,
                     dry_run=False)
    metrics[pfx + 'nReads(single file)'] = int(output)

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
    metrics_df = pd.DataFrame.from_dict({isolate: metrics},
                                        orient='index')
    # NB: this will overwrite the outfile that was read in at the first step.
    write_pandas_df(outfile, metrics_df)


def run_mash(infile,
             outfile,
             thread,
             db=mash_db,
             dry_run=False,
             log_file=None
             ):
    if db.endswith('.msh'):
        pass
    run_cmd(f"{mash_path} dist {db} {infile} -p {thread} -v 0.05 > {outfile}",
            dry_run=dry_run, log_file=log_file)


def run_kraken(infile,
               outfile,
               fmt,
               dbase,
               threads,
               dry_run=False,
               log_file=None
               ):
    '''
    infile must be a list (even is a list with length==1)
    Run Kraken on the infile.
    '''
    assert type(infile) != str

    cmd_kraken = ''
    if fmt == 'reads':
        assert len(infile) == 2
        infiles = ' '.join(infile)
        compression = ''
        for i in infile:
            # todo: that is a good idea to get file format!!!
            # todo: write a unify function/module to get the format of any files.
            # Compression test based on file extension using linux 'file'.
            f_fmt = run_cmd('file ' + i,
                            dry_run=False,
                            get_output=True).rstrip().split()
            if 'gzip' in f_fmt:
                compression = '--gzip-compressed'
                break
            if 'bzip2' in f_fmt:
                compression = '--bzip2-compressed'
                break

        cmd_kraken = "{kraken2_path} --threads {threads} --db {dbase} {compression} --paired --report {outfile} --memory-mapping {infiles} ".format(
            kraken2_path=kraken2_path,
            threads=threads,
            dbase=dbase,
            compression=compression,
            outfile=outfile,
            infiles=infiles)

    if fmt == 'contigs':
        assert len(infile) == 1
        infile = infile[0]
        cmd_kraken = "{kraken2_path} --threads {threads} --db {dbase} --report {outfile} --memory-mapping {infile} ".format(
            kraken2_path=kraken2_path,
            threads=threads,
            dbase=dbase,
            outfile=outfile,
            infile=infile)

    run_cmd(cmd_kraken, dry_run=dry_run, log_file=log_file)


def run_mlst(assembly,
             outfile,
             species,
             return_opath=False,
             dry_run=False,
             log_file=None):
    '''
    Run Torsten's MLST program on the assembly.
    For multiple input...
    assembly must be a list even the length is 1
    species is also a list with equal length as assembly
    '''

    assert type(assembly) != str
    sp_scheme = defaultdict(list)
    for _seq, _species in zip(assembly, species):
        if _species in FORCE_MLST_SCHEME:
            _scheme = FORCE_MLST_SCHEME[_species]
            sp_scheme[_scheme].append(_seq)
        elif _species.split(' ')[0] in FORCE_MLST_SCHEME:
            _scheme = FORCE_MLST_SCHEME[_species.split(' ')[0]]
            sp_scheme[_scheme].append(_seq)
        else:
            # todo: exception? or just print?
            raise Exception("%s not have scheme" % _species)

    if return_opath:
        dry_run = True
        log_file = None
    ofiles = []

    for idx, (_scheme, all_files) in enumerate(sp_scheme.items()):
        if type(_scheme) == str:
            cmd = '{mlst_path} --scheme {scheme} --quiet {infile} > {ofile}.{idx}'.format(
                mlst_path=mlst_path,
                scheme=_scheme,
                infile=' '.join(all_files),
                ofile=outfile,
                idx=idx)
            run_cmd(cmd,
                    dry_run=dry_run,
                    log_file=log_file)
            ofiles.append(f"{outfile}.{idx}")
        else:
            for sub_scheme in _scheme:
                cmd = '{mlst_path} --scheme {scheme} --quiet {infile} > {ofile}.{idx}.{scheme}'.format(
                    mlst_path=mlst_path,
                    scheme=sub_scheme,
                    infile=' '.join(all_files),
                    ofile=outfile,
                    idx=idx)
                run_cmd(cmd,
                        dry_run=dry_run,
                        log_file=log_file)
                ofiles.append(f"{outfile}.{idx}.{sub_scheme}")
    return ofiles


############################################################

def post_analysis(workflow_task):
    # only ofr workflow post analysis
    if workflow_task.dry_run:
        print("Dry run complete without post-analysis function.")
        return
    list_roary_dir = glob(os.path.join(workflow_task.odir,'*_roary_o'))
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
    for roary_dir in list_roary_dir:
        set_name = os.path.basename(roary_dir).split('_')[0]
        locus2group, locus2annotate, sample2gff = get_accessory_obj(roary_dir,
                                                                    abricate_file,
                                                                    prokka_o)
        # get locus annotation(abricate)/group(roary) from different files.
        # sample2gff contains three objs:
        #   gff_db(for query), empty_gff_obj(gff obj but removed features), gff_obj(complete gff obj)
        empty_sample2gff = {sn: vals[1]
                            for sn, vals in sample2gff.items()}
        ori_sample2gff = {sn: vals[2]
                          for sn, vals in sample2gff.items()}

        merged_locus2annotate = locus2group.copy()
        merged_locus2annotate.update(locus2annotate)
        # use annotate to overlap group info.
        # merged them
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
        # main part for generate multiple GFF files
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
                                                    set_name,
                                                    "gff_with_%s" % name)
            subset_gff_with_annotated_dir = os.path.join(summary_odir,
                                                         set_name,
                                                         "gff_of_%s" % name)
            annotated_gff_odir = os.path.join(summary_odir,
                                              set_name,
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
    mlst_ofile = workflow_task.input()["mlst_summary"].path
    pandoo_df = pd.read_csv(mlst_ofile, index_col=0)
    type_ST = list(pandoo_df.columns)
    type_ST = set([_.split('.')[0] for _ in type_ST])

    for scheme in type_ST:
        with open(os.path.join(summary_odir,
                               "%s_mlst.csv" % scheme.replace('_ST', '')), 'w') as f1:
            pandoo_df.loc[:, pandoo_df.columns.str.startswith(scheme)].to_csv(f1, index=1)
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
    # core_gene_tree = workflow_task.input()["fasttree"].path
    # os.system("cp %s %s" % (core_gene_tree,
    #                         summary_odir))
    ############################################################
    # roary plot
    # cmdline = "cd {roarydir}; {roary_plot} {core_gene_tree} {ab_csv}".format(roarydir=roary_dir,
    #                                                                          roary_plot=roary_plot_path,
    #                                                                          core_gene_tree=core_gene_tree,
    #                                                                          ab_csv=os.path.join(roary_dir,
    #                                                                                              "gene_presence_absence.csv"))
    # run_cmd(cmdline,
    #         dry_run=workflow_task.dry_run,
    #         log_file=workflow_task.log_path)


############################################################
def accessment_assembly(r1, r2, ref, gff, test_dir, dryrun=True):
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
