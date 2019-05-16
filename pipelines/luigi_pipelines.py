import json

import luigi

from pipelines import constant_str as constant
from pipelines.tasks import *
from toolkit.utils import validate_table


# give some default parameter

class fastqc(luigi.Task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter(default=None)
    sample_name = luigi.Parameter(default=None)
    status = luigi.Parameter(default="before")

    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        if self.status == 'before':
            return
        elif self.status == 'after':
            return trimmomatic(R1=self.R1,
                               R2=self.R2,
                               odir=self.odir,
                               sample_name=self.sample_name,
                               dry_run=self.dry_run)

    def output(self):
        odir = os.path.join(self.odir, "fastqc_%s" % self.status)
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
                   log_file=log_file_stream)


class multiqc(luigi.Task):
    PE_data = luigi.TupleParameter()
    status = luigi.Parameter()
    odir = luigi.Parameter()
    other_info = luigi.DictParameter(default={})
    dry_run = luigi.BoolParameter()

    def requires(self):
        if self.status == 'before' or self.status == 'after':
            return [fastqc(R1=_R1,
                           R2=_R2,
                           status=self.status,
                           sample_name=sn,
                           odir=self.odir,
                           dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]
        elif self.status == "quast":
            return [quast(R1=_R1,
                          R2=_R2,
                          odir=self.odir,
                          sample_name=sn,
                          other_info=self.other_info,
                          dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]
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
                    log_file=log_file_stream)


class preprocess_SE(luigi.Task):
    R1 = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    sample_name = luigi.Parameter()

    def output(self):
        formatted_file = os.path.join(self.odir,
                                      "cleandata",
                                      "%s.fasta" % self.sample_name)
        return luigi.LocalTarget(formatted_file)

    def run(self):
        formatted_file = os.path.join(self.odir,
                                      "cleandata",
                                      "%s.fasta" % self.sample_name)
        if not os.path.isfile(self.R1):
            raise Exception
        valid_path(os.path.dirname(formatted_file), check_odir=True)

        run_cmd("ln -s '{ori}' {new}".format(ori=self.R1,
                                             new=self.output().path),
                dry_run=self.dry_run,
                log_file=log_file_stream)


class trimmomatic(luigi.Task):
    dry_run = luigi.BoolParameter()
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    odir = luigi.Parameter()
    sample_name = luigi.Parameter()

    def output(self):
        # two files which are clean R1,R2.
        # return [luigi.LocalTarget(f) for f in ofiles]
        ofile1 = os.path.join(self.odir,
                              "cleandata",
                              str(self.sample_name) + '_R1.clean.fq.gz')
        ofile2 = os.path.join(self.odir,
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
                        log_file=log_file_stream
                        )
        # remove the log file, too waste
        try:
            pth = self.output()[0].path.replace('_R1.clean.fq.gz', '.log')
            os.remove(pth)
        except:
            pass


class quast(luigi.Task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    other_info = luigi.Parameter(default=None)
    sample_name = luigi.Parameter()

    def requires(self):
        return [shovill(R1=self.R1,
                        R2=self.R2,
                        odir=self.odir,
                        dry_run=self.dry_run,
                        sample_name=self.sample_name,
                        status="regular"),

                trimmomatic(R1=self.R1,
                            R2=self.R2,
                            odir=self.odir,
                            sample_name=self.sample_name,
                            dry_run=self.dry_run)]

    def output(self):
        # todo:
        odir = os.path.join(self.odir,
                            'assembly_o',
                            "regular_quast",
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              "report.html"))

    def run(self):
        if self.other_info is not None:
            other_info = json.loads(self.other_info)
            ref = other_info[str(self.sample_name)]["ref"]  # todo: formatted header of dataframe
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
                  log_file=log_file_stream)


class shovill(luigi.Task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    sample_name = luigi.Parameter()
    status = luigi.Parameter()

    def requires(self):
        return trimmomatic(R1=self.R1,
                           R2=self.R2,
                           odir=self.odir,
                           sample_name=self.sample_name,
                           dry_run=self.dry_run)

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
                    log_file=log_file_stream
                    )


class prokka(luigi.Task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    sample_name = luigi.Parameter()

    def requires(self):
        if not self.R2:
            return preprocess_SE(R1=self.R1,
                                 odir=self.odir,
                                 dry_run=self.dry_run,
                                 sample_name=self.sample_name)
        elif self.R2:
            return shovill(R1=self.R1,
                           R2=self.R2,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           sample_name=self.sample_name,
                           status='regular')

    def output(self):
        odir = os.path.join(self.odir,
                            "prokka_o",
                            str(self.sample_name))
        return luigi.LocalTarget(os.path.join(odir,
                                              str(self.sample_name) + '.gff'))

    def run(self):
        prokka_in_file = self.input().path

        run_prokka(infile=prokka_in_file,
                   odir=os.path.dirname(self.output().path),
                   dry_run=self.dry_run,
                   log_file=log_file_stream)


class roary(luigi.Task):
    # comparative
    PE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        return [prokka(R1=_R1,
                       R2=_R2,
                       sample_name=sn,
                       odir=self.odir,
                       dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]

    def output(self):
        odir = os.path.join(self.odir, "all_roary_o")
        return [luigi.LocalTarget(os.path.join(odir, "core_gene_alignment.aln")),
                luigi.LocalTarget(os.path.join(odir, "clustered_proteins"))]

    def run(self):
        run_roary(os.path.dirname(os.path.dirname(self.input()[0].path)),
                  os.path.dirname(self.output()[0].path),
                  thread=constant.p_roary,  # todo: determine the thread
                  dry_run=self.dry_run,
                  log_file=log_file_stream)


class fasttree(luigi.Task):
    # comparative
    PE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        return roary(PE_data=self.PE_data,
                     odir=self.odir,
                     dry_run=self.dry_run)

    def output(self):
        aln_file = self.input()[0].path
        return luigi.LocalTarget(aln_file.replace('.aln', '.newick'))

    def run(self):
        run_fasttree(self.input()[0].path,
                     self.output().path,
                     dry_run=self.dry_run,
                     log_file=log_file_stream)


class pandoo(luigi.Task):
    # single and joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        required_tasks = {}

        required_tasks["PE"] = [shovill(R1=_R1,
                                        R2=_R2,
                                        sample_name=sn,
                                        odir=self.odir,
                                        dry_run=self.dry_run,
                                        status='regular') for sn, _R1, _R2 in self.PE_data]
        required_tasks["SE"] = [preprocess_SE(R1=_R1,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              sample_name=sn) for sn, _R1 in self.SE_data]
        return required_tasks

    def output(self):
        odir = os.path.join(str(self.odir),
                            "pandoo_o")
        ofile = os.path.join(odir,
                             "pandoo_input_metadataAll.csv")
        return luigi.LocalTarget(ofile)

    def run(self):
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
                   log_file=log_file_stream)


class abricate(luigi.Task):
    # single and joint
    PE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        require_tasks = []
        require_tasks.append(roary(PE_data=self.PE_data,
                                   odir=self.odir,
                                   dry_run=self.dry_run))
        require_tasks += [prokka(R1=_R1,
                                 R2=_R2,
                                 sample_name=sn,
                                 odir=self.odir,
                                 dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]
        return require_tasks

    def output(self):
        odir = os.path.join(str(self.odir),
                            "abricate_result")
        ofile = os.path.join(odir,
                             "locus2annotate.csv")
        return luigi.LocalTarget(ofile)

    def run(self):
        prokka_o = os.path.dirname(os.path.dirname(self.input()[1].path))
        roary_dir = os.path.dirname(self.input()[0][0].path)
        # add(instead of annotated) annotated features into gff
        run_abricate(prokka_o,
                     roary_dir=roary_dir,
                     odir=os.path.dirname(self.output().path),
                     thread=constant.p_abricate,  # todo: determine the thread
                     mincov=constant.mincov_abricate,  # todo
                     dry_run=self.dry_run,
                     log_file=log_file_stream)


class ISEscan(luigi.Task):
    # single
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    sample_name = luigi.Parameter()

    def requires(self):
        if not self.R2:
            return preprocess_SE(R1=self.R1,
                                 odir=self.odir,
                                 dry_run=self.dry_run,
                                 sample_name=self.sample_name)
        elif self.R2:
            return shovill(R1=self.R1,
                           R2=self.R2,
                           sample_name=self.sample_name,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           status='regular')

    def output(self):
        final_name = "%s.gff" % self.sample_name
        ofile = os.path.join(str(self.odir),
                             "ISscan_result",
                             str(self.sample_name),
                             final_name)

        return luigi.LocalTarget(ofile)

    def run(self):
        infile_pth = self.input().path

        run_ISEscan(infile=infile_pth,
                    odir=os.path.dirname(os.path.dirname(self.output().path)),
                    sample_name=self.sample_name,
                    dry_run=self.dry_run,
                    log_file=log_file_stream)
        # For rename
        odir = os.path.dirname(self.output().path)
        for ori_file in glob(os.path.join(odir, '*')):
            suffix = os.path.basename(ori_file).split('.')[-1]
            os.renames(ori_file,
                       os.path.join(odir,
                                    str(self.sample_name) + '.%s' % suffix))


class ISEscan_summary(luigi.Task):
    # joint
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        required_tasks = {"IS_scan": [],
                          "prokka": [],
                          "roary": '',
                          "abricate": ''}
        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2='',
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run) for sn, _R1 in self.SE_data]
        required_tasks["prokka"] += [prokka(R1=_R1,
                                            R2=_R2,
                                            sample_name=sn,
                                            odir=self.odir,
                                            dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]
        required_tasks["prokka"] += [prokka(R1=_R1,
                                            R2='',
                                            sample_name=sn,
                                            odir=self.odir,
                                            dry_run=self.dry_run) for sn, _R1 in self.SE_data]
        required_tasks["abricate"] = abricate(PE_data=self.PE_data,
                                              odir=self.odir,
                                              dry_run=self.dry_run)
        required_tasks["roary"] = roary(PE_data=self.PE_data,
                                        odir=self.odir,
                                        dry_run=self.dry_run)
        # fixed... "Self adjusted ISEscan" would work.....
        return required_tasks

    def output(self):
        odir = os.path.dirname(os.path.dirname(self.input()["IS_scan"][0].path))
        ofile = os.path.join(odir,
                             'IS_result.sum')

        return luigi.LocalTarget(ofile)

    def run(self):
        total_samples = [sn for sn, _R1, _R2 in self.PE_data] + \
                        [sn for sn, _R1 in self.SE_data]
        from toolkit.process_IS import get_IS_CDS, get_locus2group
        roary_dir = os.path.dirname(self.input()["roary"][0].path)
        locus2group = get_locus2group(roary_dir)

        abricate_file = self.input()["abricate"].path
        locus2annotate_df = pd.read_csv(abricate_file, sep=',', index_col=0)
        locus2annotate = locus2annotate_df.loc[:, 'gene'].to_dict()
        final_r = {}
        for IS_gff, ori_gff, sample_name in zip(self.input()["IS_scan"],
                                                self.input()["prokka"],
                                                total_samples
                                                ):
            ori_gff = ori_gff.path
            IS_gff = IS_gff.path
            IS2CDS, IS2INFO = get_IS_CDS(ori_gff, IS_gff, locus2annotate, locus2group)
            final_r[sample_name] = (IS2CDS, IS2INFO)
        # {sample: ({IS_id: [group1,group2]},
        #           {IS_id: IS_info_dict})
        ready2df = {sn: {ISgroup: ISgroups.count(ISgroup)
                         for ISgroups in info[0].values()
                         for ISgroup in ISgroups}
                    for sn, info in final_r.items()}
        result_df = pd.DataFrame.from_dict(ready2df, orient='index')
        result_df = result_df.fillna(0)
        result_df.to_csv(self.output().path, sep=',', index=1)


class detect_plasmid(luigi.Task):
    # joint
    PE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        required_tasks = {"IS_scan": [],
                          "prokka": [],
                          "roary": '',
                          "shovill": []}

        required_tasks["IS_scan"] += [ISEscan(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["prokka"] += [prokka(R1=_R1,
                                            R2=_R2,
                                            sample_name=sn,
                                            odir=self.odir,
                                            dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]
        required_tasks["roary"] = roary(PE_data=self.PE_data,
                                        odir=self.odir,
                                        dry_run=self.dry_run)
        required_tasks["shovill"] += [shovill(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              status='plasmid'
                                              ) for sn, _R1, _R2 in self.PE_data]

        return required_tasks

    def output(self):
        return luigi.LocalTarget(os.path.join(str(self.odir),
                                              "plasmid_summary",
                                              "plasmid_summary.csv"))

    def run(self):
        run_plasmid_detect(indir=str(self.input()["shovill"][0].path).rsplit('/', maxsplit=3)[0],  # todo: maybe not compatiable to windows OS.
                           roary_dir=os.path.dirname(self.input()["roary"][0].path),
                           prokka_dir=os.path.dirname(os.path.dirname(self.input()["prokka"][0].path)),
                           odir=os.path.dirname(self.output().path),
                           dry_run=self.dry_run,
                           log_file=log_file_stream)


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
        infile_pth = self.input().path

        run_phigaro(infile_pth,
                    ofile=self.output().path,
                    thread=constant.p_phigaro,
                    dry_run=self.dry_run,
                    log_file=log_file_stream
                    )


class phigaro_summary(ISEscan_summary):

    def requires(self):
        required_tasks = {"phigaro": [],
                          "prokka": [],
                          "roary": '',
                          "abricate": ''}
        required_tasks["phigaro"] += [phigaro(R1=_R1,
                                              R2=_R2,
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run,
                                              ) for sn, _R1, _R2 in self.PE_data]
        required_tasks["phigaro"] += [phigaro(R1=_R1,
                                              R2='',
                                              sample_name=sn,
                                              odir=self.odir,
                                              dry_run=self.dry_run) for sn, _R1 in self.SE_data]
        required_tasks["prokka"] += [prokka(R1=_R1,
                                            R2=_R2,
                                            sample_name=sn,
                                            odir=self.odir,
                                            dry_run=self.dry_run) for sn, _R1, _R2 in self.PE_data]
        required_tasks["prokka"] += [prokka(R1=_R1,
                                            R2='',
                                            sample_name=sn,
                                            odir=self.odir,
                                            dry_run=self.dry_run) for sn, _R1 in self.SE_data]
        required_tasks["abricate"] = abricate(PE_data=self.PE_data,
                                              odir=self.odir,
                                              dry_run=self.dry_run)
        required_tasks["roary"] = roary(PE_data=self.PE_data,
                                        odir=self.odir,
                                        dry_run=self.dry_run)
        return required_tasks

    def output(self):
        odir = os.path.dirname(os.path.dirname(self.input()["phigaro"][0].path))
        ofile = os.path.join(odir,
                             'phage_summary.tab')

        return luigi.LocalTarget(ofile)

    def run(self):
        import json
        from toolkit.process_phigaro import write_gff

        total_samples = [sn for sn, _R1, _R2 in self.PE_data] + \
                        [sn for sn, _R1 in self.SE_data]
        summary_df = pd.DataFrame(
            columns=['region',
                     "other info"])

        for phigaro_tab_pth, ori_gff, sample_name in zip(self.input()["phigaro"],
                                                         self.input()["prokka"],
                                                         total_samples):
            # write new phage annotated gff
            write_gff(phigaro_tab_pth,
                      ori_gff)
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
                summary_df.append(pd.DataFrame([[region,
                                                 json.dumps(other_info)]], index=[0]))
        summary_df.to_csv(self.output().path, sep='\t', index=0)


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    log_file = luigi.Parameter(default=None)

    def output(self):
        pass

    def requires(self):
        input_df = pd.read_csv(self.tab, sep='\t', index_col=0, dtype=str)
        # index is samples ID, following is [R1,R2,ref_fna,ref_gff]
        validate_table(input_df)

        PE_rows = input_df.loc[(~input_df.iloc[:, 0].isna()) & (~input_df.iloc[:, 1].isna()), :]
        # R1 and R2 are both not null
        Single_rows = input_df.loc[~input_df.index.isin(PE_rows.index), :]
        # except the PE_rows
        valid_path(self.odir, check_odir=True)
        require_tasks = {}
        pairreads = tuple(zip(PE_rows.index,
                              PE_rows.iloc[:, 0],
                              PE_rows.iloc[:, 1], ))
        singlereads = tuple(zip(Single_rows.index,
                                Single_rows.iloc[:, 0]))
        if PE_rows.shape[1] > 2:
            other_info = PE_rows.to_dict(orient='index')
        else:
            other_info = None

        global log_file_stream
        if self.log_file is None:
            log_file = os.path.join(str(self.odir), "pipelines.log")
        else:
            log_file = os.path.abspath(self.odir)

        log_file_stream = None
        if os.path.isfile(log_file):
            if os.path.getsize(log_file) == 0:
                log_file_stream = open(log_file, 'a')
        if log_file_stream is None:
            log_file_stream = open(log_file, 'w')

        require_tasks["fastqc_before"] = multiqc(PE_data=pairreads,
                                                 status='before',
                                                 odir=self.odir,
                                                 dry_run=self.dry_run)
        require_tasks["fastqc_after"] = multiqc(PE_data=pairreads,
                                                status='after',
                                                odir=self.odir,
                                                dry_run=self.dry_run,
                                                )
        require_tasks["fastqc_quast"] = multiqc(PE_data=pairreads,
                                                status='quast',
                                                other_info=json.dumps(other_info),
                                                odir=self.odir,
                                                dry_run=self.dry_run,
                                                )
        require_tasks["abricate"] = abricate(PE_data=pairreads,
                                             odir=self.odir,
                                             dry_run=self.dry_run)
        require_tasks["fasttree"] = fasttree(PE_data=pairreads,
                                             odir=self.odir,
                                             dry_run=self.dry_run)
        require_tasks["ISEscan_summary"] = ISEscan_summary(PE_data=pairreads,
                                                           SE_data=singlereads,
                                                           odir=self.odir,
                                                           dry_run=self.dry_run)
        require_tasks["pandoo"] = pandoo(PE_data=pairreads,
                                         SE_data=singlereads,
                                         odir=self.odir,
                                         dry_run=self.dry_run)
        require_tasks["detect_plasmid"] = detect_plasmid(PE_data=pairreads,
                                                         odir=self.odir,
                                                         dry_run=self.dry_run)
        require_tasks["detect_prophage"] = phigaro_summary(PE_data=pairreads,
                                                           SE_data=singlereads,
                                                           odir=self.odir,
                                                           dry_run=self.dry_run)
        return require_tasks

    def run(self):
        # post pipelines
        post_analysis(self)


if __name__ == '__main__':
    luigi.run()
    # luigi.build([workflow(tab="/home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab",
    #                       odir="/home/liaoth/project/genome_pipelines/pipelines/test/test_luigi",
    #                       dry_run=False)],
    #             workers=5,
    #             local_scheduler=True
    #             )
    # log_file_stream.close()
    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # local cmd

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/tools/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/tools/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # server cmd
