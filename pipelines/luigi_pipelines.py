import json

import luigi
import pandas as pd

from pipelines.tasks import *


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
            ofiles = ["%s_fastqc.zip" % os.path.basename(self.R1).rsplit('.', maxsplit=2)[0],
                      "%s_fastqc.zip" % os.path.basename(self.R2).rsplit('.', maxsplit=2)[0]]
        elif self.status == 'after':
            ofiles = ["%s_fastqc.zip" % os.path.basename(self.input()[0].path).rsplit('.', maxsplit=2)[0],
                      "%s_fastqc.zip" % os.path.basename(self.input()[1].path).rsplit('.', maxsplit=2)[0], ]
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

    def output(self):
        indir = os.path.dirname(self.input()[0][0].path)  # any one is ok
        filename = os.path.basename(indir)
        target_file = os.path.join(indir, filename + '.html')
        return luigi.LocalTarget(target_file)

    def run(self):
        # if self.status == 'before' or self.status == 'after':
        #     indir = os.path.dirname(self.input()[0].path)  # any one is ok
        # elif self.status == "quast":
        indir = os.path.dirname(self.input()[0][0].path)  # any one is ok
        filename = os.path.basename(indir)
        run_multiqc(in_dir=indir,
                    odir=indir,
                    fn=filename,
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
                        odir=os.path.join(self.odir, "cleandata"),
                        sample_name=self.sample_name,
                        thread=1,  # todo
                        dry_run=self.dry_run,
                        log_file=log_file_stream
                        )


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
        return [luigi.LocalTarget(os.path.join(odir, "report.html"))]

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
                  odir=os.path.dirname(self.output()[0].path),
                  thread=1,  # todo
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
        odir = os.path.join(self.odir,
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
                    thread=0,  # todo
                    ram=2,  # todo
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
            return
        elif not self.R2:
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
        if not self.R2:
            prokka_in_file = self.R1
        elif not self.R2:
            prokka_in_file = self.input().path
        run_prokka(infile=prokka_in_file,
                   odir=os.path.dirname(self.output().path),
                   dry_run=self.dry_run,
                   log_file=log_file_stream)


class roary(luigi.Task):
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
                  thread=7,  # todo
                  dry_run=self.dry_run,
                  log_file=log_file_stream)


class fasttree(luigi.Task):
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
    PE_data = luigi.TupleParameter()
    SE_data = luigi.TupleParameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()

    def requires(self):
        return [shovill(R1=_R1,
                        R2=_R2,
                        sample_name=sn,
                        odir=self.odir,
                        dry_run=self.dry_run,
                        status='regular') for sn, _R1, _R2 in self.PE_data]

    def output(self):
        odir = os.path.join(self.odir, "pandoo_o")
        ofile = os.path.join(odir, "isolates_metadataAll.csv")
        return luigi.LocalTarget(ofile)

    def run(self):
        pandoo_tab = pd.DataFrame()
        for idx in range(len(self.PE_data)):
            sn, _R1, _R2 = self.PE_data[idx]
            pandoo_tab = pandoo_tab.append(pd.DataFrame([[self.input()[idx].path,
                                                          _R1,
                                                          _R2,
                                                          ]], index=[sn]))
        for idx in range(len(self.SE_data)):
            sn, _R1 = self.SE_data[idx]
            formatted_name = os.path.join(self.odir,
                                          "cleandata",
                                          "%s.fasta" % sn)
            if not os.path.isfile(formatted_name):
                os.system("ln -s {ori} {new}".format(ori=_R1,
                                                     new=formatted_name))
            inpth = formatted_name
            pandoo_tab = pandoo_tab.append(pd.DataFrame([[inpth,
                                                          '',
                                                          '',
                                                          ]], index=[sn]))
        pandoo_file = os.path.join(self.odir, "pandoo_input.tab")
        with open(pandoo_file, 'w') as f1:
            pandoo_tab.to_csv(f1, sep='\t', header=None)
        run_pandoo(in_file=pandoo_file,
                   odir=os.path.dirname(self.output().path),
                   thread=7,  # todo
                   dry_run=self.dry_run,
                   log_file=log_file_stream)


class abricate(luigi.Task):
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
        run_abricate(prokka_o,
                     roary_dir=roary_dir,
                     odir=os.path.dirname(self.output().path),
                     thread=7,  # todo
                     mincov=80,  # todo
                     dry_run=self.dry_run,
                     log_file=log_file_stream)


class ISEscan(luigi.Task):
    R1 = luigi.Parameter()
    R2 = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    sample_name = luigi.Parameter()

    def requires(self):
        if not self.R2:
            return
        elif self.R2:
            return shovill(R1=self.R1,
                           R2=self.R2,
                           sample_name=self.sample_name,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           status='regular')

    def output(self):
        ofile = os.path.join(str(self.odir),
                         "ISscan_result",
                         str(self.sample_name),
                         'contigs.fa.gff')

        return luigi.LocalTarget(ofile)

    def run(self):
        if not self.R2:
            formatted_name = os.path.join(self.odir,
                                          "cleandata",
                                          "%s.fasta" % self.sample_name)
            if not os.path.isfile(formatted_name):
                os.system("ln -s '{ori}' {new}".format(ori=self.R1,
                                                     new=formatted_name))
            infile_pth = formatted_name
        elif self.R2:
            infile_pth = self.input().path
        else:
            raise Exception

        run_ISEscan(infile=infile_pth,
                    odir=os.path.dirname(os.path.dirname(self.output().path)),
                    sample_name=self.sample_name,
                    dry_run=self.dry_run,
                    log_file=log_file_stream)


class ISEscan_summary(luigi.Task):
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
        # todo: SE data may go wrong, because of the name of directory is unknown.
        return required_tasks

    def output(self):
        ofile = os.path.join(str(self.odir),
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
        locus2annotate = dict(zip(locus2annotate_df.index,
                                  locus2annotate_df.loc[:, 'gene']))
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
                                              status='regular'
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


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    log_file = luigi.Parameter(default=None)

    def requires(self):

        input_df = pd.read_csv(self.tab, sep='\t', index_col=0, dtype=str)
        # index is samples ID, following is [R1,R2,ref_fna,ref_gff]
        validate_table(input_df)

        PE_rows = input_df.loc[(~input_df.iloc[:, 0].isna()) & (~input_df.iloc[:, 1].isna()), :]
        # R1 and R2 are both not null
        Single_rows = input_df.loc[~input_df.index.isin(PE_rows.index), :]
        # except the PE_rows
        valid_path(self.odir, check_odir=True)
        require_tasks = []
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
        log_file_stream = open(log_file, 'w')

        require_tasks.append(multiqc(PE_data=pairreads,
                                     status='before',
                                     odir=self.odir,
                                     dry_run=self.dry_run,
                                     ))
        require_tasks.append(multiqc(PE_data=pairreads,
                                     status='after',
                                     odir=self.odir,
                                     dry_run=self.dry_run,
                                     ))
        require_tasks.append(multiqc(PE_data=pairreads,
                                     status='quast',
                                     other_info=json.dumps(other_info),
                                     odir=self.odir,
                                     dry_run=self.dry_run,
                                     ))
        require_tasks.append(abricate(PE_data=pairreads,
                                      odir=self.odir,
                                      dry_run=self.dry_run))
        require_tasks.append(fasttree(PE_data=pairreads,
                                      odir=self.odir,
                                      dry_run=self.dry_run))
        require_tasks.append(ISEscan_summary(PE_data=pairreads,
                                             SE_data=singlereads,
                                             odir=self.odir,
                                             dry_run=self.dry_run))
        require_tasks.append(pandoo(PE_data=pairreads,
                                    SE_data=singlereads,
                                    odir=self.odir,
                                    dry_run=self.dry_run))
        require_tasks.append(detect_plasmid(PE_data=pairreads,
                                            odir=self.odir,
                                            dry_run=self.dry_run))

        return require_tasks


if __name__ == '__main__':
    luigi.run()
    # luigi.build([workflow(tab="/home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab",
    #                       odir="/home/liaoth/project/genome_pipelines/pipelines/test/test_luigi",
    #                       dry_run=False)],
    #             workers=5,
    #             local_scheduler=True
    #             )
    log_file_stream.close()
    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # local cmd

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/tools/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/tools/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # server cmd
