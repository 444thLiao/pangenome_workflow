from os.path import join, dirname

import luigi

from pipelines import constant_str as constant
from pipelines.tasks import *


class base_luigi_task(luigi.Task):
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)
    thread = luigi.IntParameter(default=constant.total_thread)

    def get_log_path(self):
        base_log_path = self.log_path
        if base_log_path is not None:
            return base_log_path

    def get_kwargs(self):
        return dict(odir=self.odir,
                    dry_run=self.dry_run,
                    log_path=self.log_path,
                    thread=self.thread)


# give some default parameter
class fastqc(base_luigi_task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter(default=None)
    sample_name = luigi.Parameter(default=None)
    status = luigi.Parameter(default="before")

    def requires(self):
        kwargs = self.get_kwargs()
        if self.status == 'before':
            "before trimmomatic/other QC"
            "it doesn't need to require to any tasks"
            return
        elif self.status == 'after':
            return trimmomatic(R1=self.R1,
                               R2=self.R2,
                               sample_name=self.sample_name,
                               **kwargs)

    def output(self):
        odir = os.path.join(str(self.odir),
                            "fastqc_%s" % self.status)
        if self.status == 'before':
            ofiles = ["%s_fastqc.zip" % os.path.basename(str(self.R1)).rsplit('.', maxsplit=2)[0],
                      "%s_fastqc.zip" % os.path.basename(str(self.R2)).rsplit('.', maxsplit=2)[0]]
        elif self.status == 'after':
            ofiles = ["%s_fastqc.zip" % os.path.basename(self.input()[0].path).rsplit('.', maxsplit=2)[0],
                      "%s_fastqc.zip" % os.path.basename(self.input()[1].path).rsplit('.', maxsplit=2)[0], ]
        else:
            raise Exception
        # ofiles is a list of ouput of R1 & R2
        return [luigi.LocalTarget(os.path.join(odir, f)) for f in ofiles]

    def run(self):
        if self.status == 'before':
            R1 = self.R1
            R2 = self.R2
        elif self.status == 'after':
            R1 = self.input()[0].path
            R2 = self.input()[1].path

        run_fastqc(in_files=[R1, R2],
                   odir=os.path.dirname(self.output()[0].path),
                   dry_run=self.dry_run,
                   log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class multiqc(base_luigi_task):
    PE_data = luigi.TupleParameter()
    status = luigi.Parameter()
    other_info = luigi.DictParameter(default={})

    def requires(self):
        kwargs = self.get_kwargs()
        if self.status == 'before' or self.status == 'after':
            return [fastqc(R1=_R1,
                           R2=_R2,
                           status=self.status,
                           sample_name=sn,
                           **kwargs)
                    for sn, _R1, _R2 in self.PE_data]
        elif self.status == "quast":
            return [quast(R1=_R1,
                          R2=_R2,
                          sample_name=sn,
                          other_info=self.other_info[sn],
                          **kwargs)
                    for sn, _R1, _R2 in self.PE_data]
        else:
            raise Exception

    def output(self):
        if self.status == 'before' or self.status == 'after':
            indir = os.path.dirname(self.input()[0][0].path)  # any one is ok
        elif self.status == "quast":
            indir = os.path.dirname(os.path.dirname(self.input()[0].path))  # any one is ok
        else:
            raise Exception
        filename = os.path.basename(indir)
        target_file = os.path.join(indir, filename + '.html')
        return luigi.LocalTarget(target_file)

    def run(self):
        # if self.status == 'before' or self.status == 'after':
        #     indir = os.path.dirname(self.input()[0].path)  # any one is ok
        # elif self.status == "quast":
        if self.status == 'quast':
            extra_str = " -dd 1"
        else:
            extra_str = ''
        indir = os.path.dirname(self.output().path)  # any one is ok
        filename = os.path.basename(indir)
        run_multiqc(in_dir=indir,
                    odir=indir,
                    fn=filename,
                    extra_str=extra_str,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


class preprocess_SE(base_luigi_task):
    """
    For create symblioc link to cleandata directory.
    """
    R1 = luigi.Parameter()
    sample_name = luigi.Parameter()

    def output(self):
        formatted_file = os.path.join(str(self.odir),
                                      "cleandata",
                                      "%s.fasta" % self.sample_name)
        return luigi.LocalTarget(formatted_file)

    def run(self):
        formatted_file = os.path.join(str(self.odir),
                                      "cleandata",
                                      "%s.fasta" % self.sample_name)
        if not os.path.isfile(str(self.R1)):
            raise Exception
        valid_path(os.path.dirname(formatted_file), check_odir=True)

        run_cmd("ln -s '{ori}' {new}".format(ori=self.R1,
                                             new=self.output().path),
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


class trimmomatic(base_luigi_task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    sample_name = luigi.Parameter()

    def output(self):
        # two files which are clean R1,R2.
        # return [luigi.LocalTarget(f) for f in ofiles]
        ofile1 = os.path.join(str(self.odir),
                              "cleandata",
                              str(self.sample_name) + '_R1.clean.fq.gz')
        ofile2 = os.path.join(str(self.odir),
                              "cleandata",
                              str(self.sample_name) + '_R2.clean.fq.gz')
        return [luigi.LocalTarget(ofile1),
                luigi.LocalTarget(ofile2)]

    def run(self):
        run_trimmomatic(R1=self.R1,
                        R2=self.R2,
                        odir=os.path.join(str(self.odir),
                                          "cleandata"),
                        sample_name=self.sample_name,
                        thread=self.thread - 1,  # todo: determine the thread
                        dry_run=self.dry_run,
                        log_file=self.get_log_path()
                        )
        # remove the log file, too waste
        try:
            pth = self.output()[0].path.replace('_R1.clean.fq.gz', '.log')
            os.remove(pth)
        except:
            pass
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class quast(base_luigi_task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    other_info = luigi.Parameter(default=None)
    sample_name = luigi.Parameter()

    def requires(self):
        kwargs = self.get_kwargs()
        return [shovill(R1=self.R1,
                        R2=self.R2,
                        sample_name=self.sample_name,
                        status="regular",
                        **kwargs),

                trimmomatic(R1=self.R1,
                            R2=self.R2,
                            sample_name=self.sample_name,
                            **kwargs)]

    def output(self):
        # todo: formatted the output directory
        odir = os.path.join(str(self.odir),
                            'assembly_o',
                            "regular_quast",
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              "report.html"))

    def run(self):
        if self.other_info is not None:
            ref = self.other_info["ref"]
            # todo: formatted header of dataframe
            gff = self.other_info["gff"]
            if pd.isna(ref) or str(ref) == "nan":
                ref = None
            if pd.isna(gff) or str(gff) == "nan":
                gff = None
        else:
            ref = gff = None

        run_quast(contig=self.input()[0].path,
                  R1=self.input()[1][0].path,
                  R2=self.input()[1][1].path,
                  ref=ref,
                  gff=gff,
                  odir=os.path.dirname(self.output().path),
                  thread=self.thread - 1,  # todo: determine the thread
                  dry_run=self.dry_run,
                  log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class shovill(base_luigi_task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    sample_name = luigi.Parameter()
    status = luigi.Parameter()

    def requires(self):
        kwargs = self.get_kwargs()
        return trimmomatic(R1=self.R1,
                           R2=self.R2,
                           sample_name=self.sample_name,
                           **kwargs)

    def output(self):
        if self.status == 'plasmid':
            dirname = "plasmidsSpades"
            ofile = "contigs.fasta"
        elif self.status == 'regular':
            dirname = "regular"
            ofile = "contigs.fa"
        else:
            raise Exception
        odir = os.path.join(str(self.odir),
                            "assembly_o",
                            dirname,
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              ofile))

    # .fasta header like >NODE_1_length_292801_cov_27.424733_pilon
    # .fa header like >contig00001 len=292801 cov=27.4 corr=0 spades=NODE_1_length_292801_cov_27.424733_pilon
    # plasmid pipelines need to use .fasta(which include index of plasmid info),
    # regular pipelines need to use .fa(which use contig0000x as formatted name)

    def run(self):
        if self.status == 'plasmid':
            spades_extra_options = '--plasmid'
            extra_option = "--nocorr"
        elif self.status == 'regular':
            spades_extra_options = None
            extra_option = ''
        else:
            raise Exception
        run_shovill(R1=self.input()[0].path,
                    R2=self.input()[1].path,
                    odir=os.path.dirname(self.output().path),
                    thread=self.thread - 1,  # todo: determine the thread
                    ram=constant.ram_shovill,  # todo
                    spades_extra_options=spades_extra_options,
                    extra_option=extra_option,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path()
                    )
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class prokka(base_luigi_task):
    """
    representation of task which process assembly contigs or single raw reads.
    Could used for inheritance
    """
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    sample_name = luigi.Parameter()

    def requires(self):
        kwargs = self.get_kwargs()
        if not self.R2:
            return preprocess_SE(R1=self.R1,
                                 sample_name=self.sample_name,
                                 **kwargs)
        elif self.R2:
            return shovill(R1=self.R1,
                           R2=self.R2,
                           sample_name=self.sample_name,
                           status='regular',
                           **kwargs
                           )

    def output(self):
        odir = os.path.join(str(self.odir),
                            "prokka_o",
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              str(self.sample_name) + '.gff'))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        prokka_in_file = self.input().path

        run_prokka(infile=prokka_in_file,
                   odir=os.path.dirname(self.output().path),
                   dry_run=self.dry_run,
                   log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


## pangenome analysis part
class pre_roary(base_luigi_task):
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        tasks = {}
        tasks["prokka"] = [prokka(R1=_R1,
                                  R2=_R2,
                                  sample_name=sn,
                                  **kwargs)
                           for sn, _R1, _R2 in self.PE_data] + \
                          [prokka(R1=_R1,
                                  R2='',
                                  sample_name=sn,
                                  **kwargs)
                           for sn, _R1 in self.SE_data]
        tasks["annotated_species"] = species_annotated_summary(PE_data=self.PE_data,
                                                               SE_data=self.SE_data,
                                                               **kwargs)
        return tasks

    def output(self):
        necessary_file = join(self.odir,
                              summary_dir,
                              "all_set_roary.done")
        odir = os.path.join(str(self.odir),
                            summary_dir,
                            "pairwise_mash")
        ofile = os.path.join(odir, 'pairwise_mash.dist')
        valid_path(ofile, check_ofile=1)
        return [luigi.LocalTarget(ofile),
                luigi.LocalTarget(necessary_file)]

    def run(self):
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)
        else:
            from toolkit.pre_roary import pairwise_mash
            kwargs = self.get_kwargs()
            infiles = [_.path for _ in self.input()["prokka"]]
            odir = dirname(self.output()[0].path)
            db = os.path.join(odir, 'pairwise_ref.msh')
            thread = int(self.thread) - 1
            org_annotated = self.input()["annotated_species"].path

            cluster_df = pairwise_mash(infiles=infiles,
                                       odir=odir,
                                       db=db,
                                       thread=thread,
                                       org_annotated=org_annotated,
                                       log_file=self.get_log_path())
            sids = list(cluster_df.index)
            names = list(cluster_df.loc[:, "clustering"])
            set2sids = defaultdict(list)
            for n, sid in zip(names, sids):
                annotated_org = cluster_df.loc[cluster_df.loc[:, "clustering"] == n,
                                               "annotated organism"]
                most_org = annotated_org.value_counts().index[0]
                most_species = most_org.split()[-1]
                set2sids["set%s_%s" % (n, most_species)].append(sid)
            yield batch_roary(set2sids=set2sids,
                              **kwargs)


class batch_roary(base_luigi_task):
    set2sids = luigi.DictParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        all_sids = []
        tasks = []
        for name, sids in self.set2sids.items():
            all_sids += sids
            if len(sids) > 5:
                # pass sample ID need to perform pangenome analysis
                tasks.append(fasttree(sids=sids,
                                      name=name,
                                      **kwargs))
        all_sids = list(set(sids))
        tasks.append(fasttree(sids=all_sids,
                              name="all",
                              **kwargs))
        return tasks

    def output(self):
        necessary_file = join(str(self.odir),
                              summary_dir,
                              "all_set_roary.done")
        # touch a file to indicate the status of batch_roary
        return luigi.LocalTarget(necessary_file)

    def run(self):
        # do roary_plot?
        for _ in [self.output()]:
            run_cmd("touch %s" % _.path, dry_run=False)


class roary(base_luigi_task):
    # comparative
    sids = luigi.TupleParameter()
    name = luigi.Parameter()

    # requiresment
    # has been meet at `pre_roary` tasks
    def output(self):

        odir = os.path.join(str(self.odir),
                            "%s_roary_o" % self.name)
        return [luigi.LocalTarget(os.path.join(odir, "core_gene_alignment.aln")),
                luigi.LocalTarget(os.path.join(odir, "clustered_proteins"))]

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)
        prokka_dir = os.path.join(self.odir, 'prokka_o')
        gff_files = glob(join(prokka_dir, '*', '*.gff'))

        if self.name != 'all':
            processed_sid = self.sids
            gff_files = [_
                         for _ in gff_files
                         if os.path.basename(os.path.dirname(_)) in processed_sid]
            # filter out gff files which samples not in given processed_sid

        run_roary(gff_files,
                  os.path.dirname(self.output()[0].path),
                  thread=self.thread - 1,  # todo: determine the thread
                  dry_run=self.dry_run,
                  log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class fasttree(base_luigi_task):
    # comparative
    sids = luigi.TupleParameter()
    name = luigi.Parameter()

    def requires(self):
        kwargs = self.get_kwargs()
        return roary(sids=self.sids,
                     name=self.name,
                     **kwargs)

    def output(self):
        roary_dir = dirname(self.input()[0].path)
        aln_file = self.input()[0].path
        png_file = join(roary_dir, "pangenome_matrix.png")

        return [luigi.LocalTarget(aln_file.replace('.aln', '.newick')),
                luigi.LocalTarget(png_file)]

    def run(self):
        if not os.path.exists(self.output()[0].path):
            run_fasttree(self.input()[0].path,
                         self.output()[0].path,
                         dry_run=self.dry_run,
                         log_file=self.get_log_path())
        roary_dir = dirname(self.output()[0].path)
        cmdline = "cd {roarydir}; {roary_plot} {core_gene_tree} {ab_csv}".format(roarydir=roary_dir,
                                                                                 roary_plot=roary_plot_path,
                                                                                 core_gene_tree=self.output()[0].path,
                                                                                 ab_csv=os.path.join(roary_dir,
                                                                                                     "gene_presence_absence.csv"))
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


## quality accessment
class seqtk_tasks(prokka):
    def output(self):
        odir = os.path.join(str(self.odir), "seqtk_result")
        valid_path(odir, check_odir=1)
        ofiles = []

        ofiles.append(os.path.join(odir,
                                   "%s_assembly.summary" % self.sample_name))
        if self.R2 is not None:
            ofiles.append(os.path.join(odir,
                                       "%s_reads.summary" % self.sample_name))

        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        odir = os.path.join(str(self.odir), "seqtk_result")
        valid_path(odir, check_odir=1)
        if self.dry_run:
            for _ in self.output():
                run_cmd("touch %s" % _,
                        dry_run=False)
            return

        kwargs = dict(isolate=self.sample_name,
                      dry_run=self.dry_run,
                      log_file=self.log_path)

        run_seqtk_contig(infile=self.input().path,
                         outfile=os.path.join(odir,
                                              "%s_assembly.summary" % self.sample_name),
                         **kwargs
                         )

        if self.R2 is not None:
            # input is PE fa, need to access the quality before and after assembly
            run_seqtk_reads(infile=self.R1,
                            outfile=os.path.join(odir,
                                                 "%s_reads.summary" % self.sample_name),
                            **kwargs
                            )


class seqtk_summary(base_luigi_task):
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        required_tasks = []
        required_tasks += [seqtk_tasks(R1=_R1,
                                       R2=_R2,
                                       sample_name=sn,
                                       **kwargs)
                           for sn, _R1, _R2 in self.PE_data]
        required_tasks += [seqtk_tasks(R1=_R1,
                                       sample_name=sn,
                                       **kwargs)
                           for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        # merged reads.summary into single file
        # merged assembly.summary into single file
        ofiles = [os.path.join(str(self.odir),
                               constant.summary_dir,
                               "seqtk_assembly_accessment.csv"),
                  os.path.join(str(self.odir),
                               constant.summary_dir,
                               "seqtk_reads_accessment.csv")
                  ]
        valid_path(ofiles, check_ofile=1)
        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        if self.dry_run:
            for _ in self.ofiles:
                run_cmd("touch %s" % _,
                        dry_run=False)
            return
        input_files = [_
                       for each_list in self.input()
                       for _ in each_list]
        assembly_files = [pd.read_csv(_.path, index_col=0, sep='\t')
                          for _ in input_files
                          if _.path.endswith("_assembly.summary")]
        reads_files = [pd.read_csv(_.path, index_col=0, sep='\t')
                       for _ in input_files
                       if _.path.endswith("_reads.summary")]
        if assembly_files:
            total_assembly_f = pd.concat(assembly_files, axis=0)
            total_assembly_f.to_csv(self.output()[0].path, index=1)

        if reads_files:
            total_reads_f = pd.concat(reads_files, axis=0)
            total_reads_f.to_csv(self.output()[1].path, index=1)
        else:
            run_cmd("touch %s" % self.output()[1].path,
                    dry_run=False)


## gene annotated
class abricate(base_luigi_task):
    # todo: separate into independent task and summary task
    # single and joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        all_sids = [sn for sn, _R1, _R2 in self.PE_data] + \
                   [sn for sn, _R1 in self.SE_data]
        require_tasks = []
        require_tasks.append(roary(sids=all_sids,
                                   name='all',
                                   **kwargs))
        require_tasks += [prokka(R1=_R1,
                                 R2=_R2,
                                 sample_name=sn,
                                 **kwargs) for sn, _R1, _R2 in self.PE_data]
        require_tasks += [prokka(R1=_R1,
                                 R2='',
                                 sample_name=sn,
                                 **kwargs) for sn, _R1 in self.SE_data]
        return require_tasks

    def output(self):
        odir = os.path.join(str(self.odir),
                            constant.summary_dir)
        ofile = os.path.join(odir,
                             "locus2annotate.csv")
        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        prokka_o = os.path.dirname(os.path.dirname(self.input()[1].path))
        roary_dir = os.path.dirname(self.input()[0][0].path)
        # add(instead of annotated) annotated features into gff
        run_abricate(prokka_o,
                     roary_dir=roary_dir,
                     odir=os.path.dirname(self.output().path),
                     thread=int(self.thread) - 1,
                     mincov=constant.mincov_abricate,
                     dry_run=self.dry_run,
                     log_file=self.get_log_path())
        source = os.path.join(str(self.odir),
                              "abricate_result",
                              "locus2annotate.csv")
        target = os.path.join(str(self.odir),
                              constant.summary_dir)
        run_cmd(f"mv {source} {target}", dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


## mlst typing
class mlst_task(prokka):
    def output(self):
        odir = os.path.join(str(self.odir),
                            "mlst_o",
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              str(self.sample_name)))

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        mlst_in_file = self.input().path

        run_mlst(assembly=[mlst_in_file],
                 outfile=self.output().path,
                 # it just a prefix not a file name
                 species=[constant.specific_species],
                 dry_run=self.dry_run,
                 log_file=self.get_log_path())
        from toolkit.process_mlst import parse_mlst
        ofiles = glob(self.output().path + '*')
        merged_df = parse_mlst(ofiles, sample_name=self.sample_name)
        merged_df.to_csv(self.output().path, index=1)
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


class mlst_summary(base_luigi_task):
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        required_tasks = []
        required_tasks += [mlst_task(R1=_R1,
                                     R2=_R2,
                                     sample_name=sn,
                                     **kwargs)
                           for sn, _R1, _R2 in self.PE_data]
        required_tasks += [mlst_task(R1=_R1,
                                     sample_name=sn,
                                     **kwargs)
                           for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        ofile = os.path.join(str(self.odir),
                             constant.summary_dir,
                             "mlst_all.csv", )
        return luigi.LocalTarget(ofile)

    def run(self):
        infile_paths = [_.path for _ in self.input()]
        from toolkit.process_mlst import main as merge_mlst
        all_df = [pd.read_csv(_, index_col=0)
                  for _ in infile_paths]
        merged_df = pd.concat(all_df, axis=0)
        mlst_df, merged_df = merge_mlst(merged_df)
        for s, df in mlst_df.items():
            df.to_csv(os.path.join(str(self.odir),
                                   "mlst_o",
                                   "total_%s.csv" % s), index=1)
        merged_df.to_csv(self.output().path, index=1)
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)


## annotated species
class kraken2_tasks(prokka):

    def output(self):
        odir = os.path.join(str(self.odir), "kraken2_report")
        valid_path(odir, check_odir=1)
        ofiles = []

        ofiles.append(os.path.join(odir,
                                   "%s_assembly.k2report" % self.sample_name))
        if self.R2 is not None:
            ofiles.append(os.path.join(odir,
                                       "%s_reads.k2report" % self.sample_name))

        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        odir = os.path.join(str(self.odir), "kraken2_report")
        valid_path(odir, check_odir=1)
        if self.dry_run:
            for _ in self.output():
                run_cmd("touch %s" % _.path,
                        dry_run=False)
        kwargs = dict(dbase=kraken2_db,
                      threads=self.thread - 1,
                      dry_run=self.dry_run,
                      log_file=self.log_path)

        run_kraken(infile=[self.input().path],
                   fmt="contigs",
                   outfile=os.path.join(odir,
                                        "%s_assembly.k2report" % self.sample_name),
                   **kwargs
                   )

        if self.R2 is not None:
            # input is PE fa, need to access the quality before and after assembly
            run_kraken(infile=[self.R1, self.R2],
                       fmt="reads",
                       outfile=os.path.join(odir,
                                            "%s_reads.k2report" % self.sample_name),
                       **kwargs
                       )


class kraken2_summary(base_luigi_task):
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        required_tasks = []
        required_tasks += [kraken2_tasks(R1=_R1,
                                         R2=_R2,
                                         sample_name=sn,
                                         **kwargs)
                           for sn, _R1, _R2 in self.PE_data]
        required_tasks += [kraken2_tasks(R1=_R1,
                                         sample_name=sn,
                                         **kwargs)
                           for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        ofile = os.path.join(str(self.odir),
                             constant.summary_dir,
                             "kraken2_report.csv", )
        return luigi.LocalTarget(ofile)

    def run(self):
        from toolkit.parse_kraken2 import merge_kraken2
        infile_paths = [_.path
                        for _l in self.input()
                        for _ in _l]
        merged_kraken2_report = merge_kraken2(infile_paths)
        merged_kraken2_report.to_csv(self.output().path, index=1)


class mash_tasks(prokka):
    def output(self):
        odir = os.path.join(str(self.odir), "mash_report")
        valid_path(odir, check_odir=1)
        ofiles = []

        ofiles.append(os.path.join(odir,
                                   "%s_contig.mash_report" % self.sample_name))

        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        odir = os.path.join(str(self.odir), "mash_report")
        valid_path(odir, check_odir=1)
        if self.dry_run:
            for _ in self.output():
                run_cmd("touch %s" % _.path,
                        dry_run=False)

        run_mash(infile=self.input().path,
                 outfile=os.path.join(odir,
                                      "%s_contig.mash_report" % self.sample_name),
                 thread=self.thread - 1,
                 db=mash_db,
                 dry_run=self.dry_run,
                 log_file=self.get_log_path())


class species_annotated_summary(base_luigi_task):
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        required_tasks = {}
        required_tasks["k2_summary"] = kraken2_summary(PE_data=self.PE_data,
                                                       SE_data=self.SE_data,
                                                       **kwargs)
        required_tasks["mash"] = []
        required_tasks["mash"] += [mash_tasks(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              **kwargs)
                                   for sn, _R1, _R2 in self.PE_data]
        required_tasks["mash"] += [mash_tasks(R1=_R1,
                                              sample_name=sn,
                                              **kwargs)
                                   for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        ofile = os.path.join(str(self.odir),
                             constant.summary_dir,
                             "species_annotated.csv")
        return luigi.LocalTarget(ofile)

    def run(self):
        if self.dry_run:
            for _ in [self.output()]:
                run_cmd("touch %s" % _.path,
                        dry_run=False)
        else:
            from toolkit.process_mash import parse_batch_result
            kraken2_summary_path = self.input()['k2_summary'].path
            kraken2_df = pd.read_csv(kraken2_summary_path, index_col=0)

            mash_paths = [otarget.path
                          for req_task in self.input()["mash"]
                          # get tasks which is mash_tasks
                          for otarget in req_task
                          # get tasks.output which is mash_tasks.output
                          if otarget.path.endswith(".mash_report")]
            mash_df = parse_batch_result(mash_paths, mash_db_summary)
            annotated_df = pd.concat([kraken2_df, mash_df], axis=1)
            annotated_df.to_csv(self.output().path, index=1, index_label="sample ID")


## transposed element annotated
class ISEscan(base_luigi_task):
    # single
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    sample_name = luigi.Parameter()

    def requires(self):
        kwargs = dict(R1=self.R1,
                      odir=self.odir,
                      dry_run=self.dry_run,
                      sample_name=self.sample_name,
                      log_path=self.log_path,
                      thread=self.thread)
        if not self.R2:
            return preprocess_SE(**kwargs)
        elif self.R2:
            return shovill(R2=self.R2,
                           status='regular',
                           **kwargs)

    def output(self):
        final_name = "%s.gff" % self.sample_name
        ofile = os.path.join(str(self.odir),
                             "ISscan_result",
                             str(self.sample_name),
                             final_name)

        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        infile_pth = self.input().path

        run_ISEscan(infile=infile_pth,
                    odir=os.path.dirname(os.path.dirname(self.output().path)),
                    sample_name=self.sample_name,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path())
        # For rename
        odir = os.path.dirname(self.output().path)
        for ori_file in glob(os.path.join(odir, '*')):
            suffix = os.path.basename(ori_file).split('.')[-1]
            os.renames(ori_file,
                       os.path.join(odir,
                                    str(self.sample_name) + '.%s' % suffix))
        if not glob(os.path.join(odir, '*')):
            # it mean that no IS could be detected.
            run_cmd("touch %s" % self.output().path, dry_run=False)

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class ISEscan_summary(base_luigi_task):
    # joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        kwargs = self.get_kwargs()
        required_tasks = {"IS_scan": []}
        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              **kwargs
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2='',
                                              sample_name=sn,
                                              **kwargs) for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        ofile = os.path.join(str(self.odir),
                             constant.summary_dir,
                             'IS_summary.tab')

        return luigi.LocalTarget(ofile)

    def run(self):
        from toolkit.get_gene_info import get_gff
        total_samples = [sn for sn, _R1, _R2 in self.PE_data] + \
                        [sn for sn, _R1 in self.SE_data]
        summary_df = pd.DataFrame(
            columns=["sample",
                     "region",
                     "other info"])
        for IS_gff, sample_name in zip(self.input()["IS_scan"],
                                       total_samples):
            if os.path.getsize(os.path.abspath(IS_gff.path)) == 0:
                # it mean no IS could be detected.
                continue
            gff_dict = get_gff(IS_gff.path,
                               mode='bcbio')
            for contig, record in gff_dict.items():
                IS_id = region = 'None'
                other_info = {}
                for fea in record.features:
                    if fea.type == "insertion_sequence":
                        IS_id = fea.id
                        region = "%s:%s-%s" % (record.id,
                                               fea.location.start.real,
                                               fea.location.end.real)
                        other_info = {}
                        other_info["family"] = fea.qualifiers["family"][0]
                        other_info["cluster"] = fea.qualifiers["cluster"][0]
                    else:
                        if fea.id == IS_id + '_TIR':
                            other_info["length_TIR"] = len(fea)
                    if IS_id not in summary_df.index and IS_id != 'None':
                        # if new IS_id,it must have region
                        summary_df = summary_df.append(pd.DataFrame([[sample_name,
                                                                      region,
                                                                      json.dumps(other_info)]],
                                                                    columns=summary_df.columns,
                                                                    index=[IS_id + '_' + sample_name]))
        summary_df.to_csv(self.output().path, index=1, sep='\t')
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class detect_plasmid(base_luigi_task):
    # joint
    PE_data = luigi.TupleParameter()

    def requires(self):
        required_tasks = {"shovill": []}

        required_tasks["shovill"] += [shovill(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              status='plasmid',
                                              log_path=self.log_path,
                                              thread=self.thread
                                              ) for sn, _R1, _R2 in self.PE_data]

        return required_tasks

    def output(self):
        return luigi.LocalTarget(os.path.join(str(self.odir),
                                              constant.summary_dir,
                                              "plasmid_summary.tab"))

    def run(self):
        assembly_odir = str(self.input()["shovill"][0].path).rsplit('/', maxsplit=3)[0]

        run_plasmid_detect(indir=assembly_odir,  # todo: maybe not compatible to windows OS.
                           ofile=self.output().path,
                           dry_run=self.dry_run,
                           log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class phigaro(ISEscan):
    # single

    def output(self):
        final_name = "%s.out" % self.sample_name
        ofile = os.path.join(str(self.odir),
                             "phage_detect_result",
                             str(self.sample_name),
                             final_name)

        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path, check_ofile=1)

        infile_pth = self.input().path

        run_phigaro(infile_pth,
                    ofile=self.output().path,
                    thread=self.thread - 1,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path()
                    )
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class phigaro_summary(ISEscan_summary):
    # regions specific annotated
    def requires(self):
        kwargs = self.get_kwargs()
        required_tasks = {"phigaro": []}
        required_tasks["phigaro"] += [phigaro(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              **kwargs
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["phigaro"] += [phigaro(R1=_R1,
                                              R2='',
                                              sample_name=sn,
                                              **kwargs) for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        ofile = os.path.join(str(self.odir),
                             constant.summary_dir,
                             'phage_summary.tab')

        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        if not self.dry_run:
            total_samples = [sn for sn, _R1, _R2 in self.PE_data] + \
                            [sn for sn, _R1 in self.SE_data]
            summary_df = pd.DataFrame(
                columns=["sample",
                         "region",
                         "other info"])
            for phigaro_tab_pth, sample_name in zip(self.input()["phigaro"],
                                                    total_samples):
                # output to same path of phigaro_tab
                phigaro_tab_pth = phigaro_tab_pth.path
                phigaro_tab = pd.read_csv(phigaro_tab_pth, sep='\t')
                for idx, vals in phigaro_tab.iterrows():
                    region = "%s:%s-%s" % (vals["scaffold"],
                                           vals["begin"],
                                           vals["end"])
                    other_info = vals.loc[set(vals.index).difference({"scaffold",
                                                                      "begin",
                                                                      "end"})].to_dict()
                    summary_df = summary_df.append(pd.DataFrame([[sample_name,
                                                                  region,
                                                                  json.dumps(other_info)]],
                                                                index=[sample_name],
                                                                columns=summary_df.columns))
            summary_df.to_csv(self.output().path, sep='\t', index=1)

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


## transposed element annotated <end>
#
# class pandoo(base_luigi_task):
#     # single and joint
#     PE_data = luigi.TupleParameter()
#     SE_data = luigi.TupleParameter()
#
#     def requires(self):
#         required_tasks = {}
#
#         required_tasks["PE"] = [shovill(R1=_R1,
#                                         R2=_R2,
#                                         sample_name=sn,
#                                         odir=self.odir,
#                                         dry_run=self.dry_run,
#                                         status='regular',
#                                         log_path=self.log_path) for sn, _R1, _R2 in self.PE_data]
#         required_tasks["SE"] = [preprocess_SE(R1=_R1,
#                                               odir=self.odir,
#                                               dry_run=self.dry_run,
#                                               sample_name=sn,
#                                               log_path=self.log_path) for sn, _R1 in self.SE_data]
#         return required_tasks
#
#     def output(self):
#         odir = os.path.join(str(self.odir),
#                             "pandoo_o")
#         ofile = os.path.join(odir,
#                              "pandoo_input_metadataAll.csv")
#         return luigi.LocalTarget(ofile)
#
#     def run(self):
#         valid_path(self.output().path, check_ofile=1)
#         pandoo_tab = pd.DataFrame()
#         for idx in range(len(self.PE_data)):
#             sn, _R1, _R2 = self.PE_data[idx]
#             contig_pth = self.input()["PE"][idx].path
#             pandoo_tab = pandoo_tab.append(pd.DataFrame([[contig_pth,
#                                                           _R1,
#                                                           _R2,
#                                                           ]], index=[sn]))
#         for idx in range(len(self.SE_data)):
#             sn, _R1 = self.SE_data[idx]
#             formatted_pth = self.input()["SE"][idx].path
#             pandoo_tab = pandoo_tab.append(pd.DataFrame([[formatted_pth,
#                                                           '',
#                                                           '',
#                                                           ]], index=[sn]))
#         pandoo_file = os.path.join(str(self.odir), "pandoo_input.tab")
#         with open(pandoo_file, 'w') as f1:
#             pandoo_tab.to_csv(f1, sep='\t', header=None, index=1)
#         run_pandoo(in_file=pandoo_file,
#                    odir=os.path.dirname(self.output().path),
#                    thread=constant.p_pandoo,  # todo: determine the thread
#                    dry_run=self.dry_run,
#                    log_file=self.get_log_path())
#         if self.dry_run:
#             for _o in [self.output()]:
#                 run_cmd("touch %s" % _o.path, dry_run=False)
#

if __name__ == '__main__':
    luigi.run()

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # local cmd

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/tools/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/tools/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # server cmd
