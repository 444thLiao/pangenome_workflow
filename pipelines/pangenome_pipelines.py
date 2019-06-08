import time

import luigi

from pipelines import constant_str as constant
from pipelines.tasks import *

class base_luigi_task(luigi.Task):
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)

    def get_log_path(self):
        base_log_path = self.log_path
        if base_log_path is not None:
            return base_log_path


# give some default parameter
class fastqc(base_luigi_task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter(default=None)
    sample_name = luigi.Parameter(default=None)
    status = luigi.Parameter(default="before")

    def requires(self):
        if self.status == 'before':
            return
        elif self.status == 'after':
            return trimmomatic(R1=self.R1,
                               R2=self.R2,
                               odir=self.odir,
                               sample_name=self.sample_name,
                               dry_run=self.dry_run,
                               log_path=self.log_path)

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
        if self.status == 'before' or self.status == 'after':
            return [fastqc(R1=_R1,
                           R2=_R2,
                           status=self.status,
                           sample_name=sn,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           log_path=self.log_path) for sn, _R1, _R2 in self.PE_data]
        elif self.status == "quast":
            return [quast(R1=_R1,
                          R2=_R2,
                          odir=self.odir,
                          sample_name=sn,
                          other_info=self.other_info,
                          dry_run=self.dry_run,
                          log_path=self.log_path) for sn, _R1, _R2 in self.PE_data]
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
                        thread=constant.p_trimmomatic,  # todo: determine the thread
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
        return [shovill(R1=self.R1,
                        R2=self.R2,
                        odir=self.odir,
                        dry_run=self.dry_run,
                        sample_name=self.sample_name,
                        status="regular",
                        log_path=self.log_path),

                trimmomatic(R1=self.R1,
                            R2=self.R2,
                            odir=self.odir,
                            sample_name=self.sample_name,
                            dry_run=self.dry_run,
                            log_path=self.log_path)]

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
            other_info = json.loads(self.other_info)
            ref = other_info[str(self.sample_name)]["ref"]
            # todo: formatted header of dataframe
            gff = other_info[str(self.sample_name)]["gff"]
        else:
            ref = gff = None

        run_quast(contig=self.input()[0].path,
                  R1=self.input()[1][0].path,
                  R2=self.input()[1][1].path,
                  ref=ref,
                  gff=gff,
                  odir=os.path.dirname(self.output().path),
                  thread=constant.p_quast,  # todo: determine the thread
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
        return trimmomatic(R1=self.R1,
                           R2=self.R2,
                           odir=self.odir,
                           sample_name=self.sample_name,
                           dry_run=self.dry_run,
                           log_path=self.log_path)

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
                    thread=constant.p_shovill,  # todo: determine the thread
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
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    sample_name = luigi.Parameter()

    def requires(self):
        if not self.R2:
            return preprocess_SE(R1=self.R1,
                                 odir=self.odir,
                                 dry_run=self.dry_run,
                                 sample_name=self.sample_name,
                                 log_path=self.log_path)
        elif self.R2:
            return shovill(R1=self.R1,
                           R2=self.R2,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           sample_name=self.sample_name,
                           status='regular',
                           log_path=self.log_path)

    def output(self):
        odir = os.path.join(str(self.odir),
                            "prokka_o",
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              str(self.sample_name) + '.gff'))

    def run(self):
        valid_path(self.output().path,check_ofile=1)
        prokka_in_file = self.input().path

        run_prokka(infile=prokka_in_file,
                   odir=os.path.dirname(self.output().path),
                   dry_run=self.dry_run,
                   log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)

class roary(base_luigi_task):
    # comparative
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        return [prokka(R1=_R1,
                       R2=_R2,
                       sample_name=sn,
                       odir=self.odir,
                       dry_run=self.dry_run,
                       log_path=self.log_path)
                for sn, _R1, _R2 in self.PE_data] + \
               [prokka(R1=_R1,
                       R2='',
                       sample_name=sn,
                       odir=self.odir,
                       dry_run=self.dry_run,
                       log_path=self.log_path)
                for sn, _R1 in self.SE_data]

    def output(self):
        odir = os.path.join(str(self.odir), "all_roary_o")
        return [luigi.LocalTarget(os.path.join(odir, "core_gene_alignment.aln")),
                luigi.LocalTarget(os.path.join(odir, "clustered_proteins"))]

    def run(self):
        valid_path(self.output()[0].path,check_ofile=1)
        run_roary(os.path.dirname(os.path.dirname(self.input()[0].path)),
                  os.path.dirname(self.output()[0].path),
                  thread=constant.p_roary,  # todo: determine the thread
                  dry_run=self.dry_run,
                  log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)

class fasttree(base_luigi_task):
    # comparative
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        return roary(PE_data=self.PE_data,
                     SE_data=self.SE_data,
                     odir=self.odir,
                     dry_run=self.dry_run,
                     log_path=self.log_path)

    def output(self):
        aln_file = self.input()[0].path
        return luigi.LocalTarget(aln_file.replace('.aln', '.newick'))

    def run(self):
        run_fasttree(self.input()[0].path,
                     self.output().path,
                     dry_run=self.dry_run,
                     log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)

class pandoo(base_luigi_task):
    # todo: dissect pandoo into two.
    # single and joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        required_tasks = {}

        required_tasks["PE"] = [shovill(R1=_R1,
                                        R2=_R2,
                                        sample_name=sn,
                                        odir=self.odir,
                                        dry_run=self.dry_run,
                                        status='regular',
                                        log_path=self.log_path) for sn, _R1, _R2 in self.PE_data]
        required_tasks["SE"] = [preprocess_SE(R1=_R1,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              sample_name=sn,
                                              log_path=self.log_path) for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        odir = os.path.join(str(self.odir),
                            "pandoo_o")
        ofile = os.path.join(odir,
                             "pandoo_input_metadataAll.csv")
        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path, check_ofile=1)
        pandoo_tab = pd.DataFrame()
        for idx in range(len(self.PE_data)):
            sn, _R1, _R2 = self.PE_data[idx]
            contig_pth = self.input()["PE"][idx].path
            pandoo_tab = pandoo_tab.append(pd.DataFrame([[contig_pth,
                                                          _R1,
                                                          _R2,
                                                          ]], index=[sn]))
        for idx in range(len(self.SE_data)):
            sn, _R1 = self.SE_data[idx]
            formatted_pth = self.input()["SE"][idx].path
            pandoo_tab = pandoo_tab.append(pd.DataFrame([[formatted_pth,
                                                          '',
                                                          '',
                                                          ]], index=[sn]))
        pandoo_file = os.path.join(str(self.odir), "pandoo_input.tab")
        with open(pandoo_file, 'w') as f1:
            pandoo_tab.to_csv(f1, sep='\t', header=None, index=1)
        run_pandoo(in_file=pandoo_file,
                   odir=os.path.dirname(self.output().path),
                   thread=constant.p_pandoo,  # todo: determine the thread
                   dry_run=self.dry_run,
                   log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)

class abricate(base_luigi_task):
    # single and joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()
    def requires(self):
        require_tasks = []
        require_tasks.append(roary(PE_data=self.PE_data,
                                   SE_data=self.SE_data,
                                   odir=self.odir,
                                   dry_run=self.dry_run,
                                   log_path=self.log_path))
        require_tasks += [prokka(R1=_R1,
                                 R2=_R2,
                                 sample_name=sn,
                                 odir=self.odir,
                                 dry_run=self.dry_run,
                                 log_path=self.log_path) for sn, _R1, _R2 in self.PE_data]
        require_tasks += [prokka(R1=_R1,
                                 R2='',
                                 sample_name=sn,
                                 odir=self.odir,
                                 dry_run=self.dry_run,
                                 log_path=self.log_path) for sn, _R1 in self.SE_data]
        return require_tasks

    def output(self):
        odir = os.path.join(str(self.odir),
                            "abricate_result")
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
                     thread=constant.p_abricate,  # todo: determine the thread
                     mincov=constant.mincov_abricate,  # todo
                     dry_run=self.dry_run,
                     log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)

class ISEscan(base_luigi_task):
    # single
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    sample_name = luigi.Parameter()

    def requires(self):
        if not self.R2:
            return preprocess_SE(R1=self.R1,
                                 odir=self.odir,
                                 dry_run=self.dry_run,
                                 sample_name=self.sample_name,
                                 log_path=self.log_path)
        elif self.R2:
            return shovill(R1=self.R1,
                           R2=self.R2,
                           sample_name=self.sample_name,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           status='regular',
                           log_path=self.log_path)

    def output(self):
        final_name = "%s.gff" % self.sample_name
        ofile = os.path.join(str(self.odir),
                             "ISscan_result",
                             str(self.sample_name),
                             final_name)

        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path,check_ofile=1)
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
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)

class ISEscan_summary(base_luigi_task):
    # joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()

    def requires(self):
        required_tasks = {"IS_scan": []}
        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              log_path=self.log_path
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2='',
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              log_path=self.log_path) for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        odir = os.path.dirname(os.path.dirname(self.input()["IS_scan"][0].path))
        ofile = os.path.join(odir,
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
        for IS_gff, sample_name in zip(self.input()["IS_scan"], total_samples):
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
                                                                    index=[IS_id]))
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
                                              log_path=self.log_path
                                              ) for sn, _R1, _R2 in self.PE_data]

        return required_tasks

    def output(self):
        return luigi.LocalTarget(os.path.join(str(self.odir),
                                              "plasmid_summary",
                                              "plasmid_summary.tab"))

    def run(self):
        assembly_odir = str(self.input()["shovill"][0].path).rsplit('/', maxsplit=3)[0]

        run_plasmid_detect(indir=assembly_odir,  # todo: maybe not compatiable to windows OS.
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
        valid_path(self.output().path,check_ofile=1)

        infile_pth = self.input().path

        run_phigaro(infile_pth,
                    ofile=self.output().path,
                    thread=constant.p_phigaro,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path()
                    )
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class phigaro_summary(ISEscan_summary):
    # regions specific annotated
    def requires(self):
        required_tasks = {"phigaro": []}
        required_tasks["phigaro"] += [phigaro(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              log_path=self.log_path
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["phigaro"] += [phigaro(R1=_R1,
                                              R2='',
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              log_path=self.log_path) for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        odir = os.path.dirname(os.path.dirname(self.input()["phigaro"][0].path))
        ofile = os.path.join(odir,
                             'phage_summary.tab')

        return luigi.LocalTarget(ofile)

    def run(self):
        valid_path(self.output().path,check_ofile=1)
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



if __name__ == '__main__':
    luigi.run()

    # luigi.build([workflow(tab="/home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab",
    #                       odir="/home/liaoth/project/genome_pipelines/pipelines/test/test_luigi",
    #                       dry_run=False)],
    #             workers=5,
    #             local_scheduler=True
    #             )
    # _log_stream.close()
    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # local cmd

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/tools/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/tools/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # server cmd