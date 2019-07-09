import sys
import time
from os.path import dirname

import click

sys.path.insert(0, dirname(dirname(__file__)))
from pipelines import *


@click.group()
def cli():
    pass


@cli.command()
@click.argument('cmd', nargs=-1)
def run(cmd):
    luigi.run(cmdline_args=cmd)


@cli.command(help="Cleaning some output for performing sub-samples analysis. (just run it luigi pipelines again instead of run all steps again)")
@click.option("-i", "--input_dir", "indir",
              help="output directory which passed to `run` command")
def clean(indir):
    "Cleaning some output for performing sub-samples analysis. (just run it luigi pipelines again instead of run all steps again)"
    # todo:
    pass


@cli.command(help="Just like the `clean` command instead of deleting ")
@click.option("-i", "--input_dir", "indir",
              help="output directory which passed to `run` command")
@click.option("-o", "--output_dir", "odir",
              help="A big directory for archive useless files")
@click.option("-n", "--name", "name",
              help="Given a name for identify this archived files...")
def archive(indir, odir, name=None):
    # Just like the `clean` command, but it will archive required result into a named directory.
    if name is None:
        name = str(int(time.time()))
    output_directory = os.path.join(odir, "archived", name)
    valid_path(output_directory, check_odir=1)
    cmd_template = "mv {source} {target_dir};"
    cmd = ''
    # for workflow
    cmd += cmd_template.format(source=os.path.join(indir,
                                                   "pipelines_summary"),
                               target_dir=output_directory)
    # for roary
    cmd += cmd_template.format(source=os.path.join(indir,
                                                   "all_roary_o"),
                               target_dir=output_directory)
    # for phigaro, seqtk, mlst, plasmid, IS,
    cmd += cmd_template.format(source=os.path.join(indir,
                                                   constant.summary_dir),
                               target_dir=output_directory)
    # for abricate
    cmd += cmd_template.format(source=os.path.join(indir,
                                                   "abricate_result",
                                                   "locus2annotate.csv"),
                               target_dir=output_directory)
    # for multiqc
    cmd += cmd_template.format(source=os.path.join(indir, "fastqc_after", "fastqc_after*"),
                               target_dir=output_directory)
    cmd += cmd_template.format(source=os.path.join(indir, "fastqc_before", "fastqc_before*"),
                               target_dir=output_directory)
    cmd += cmd_template.format(source=os.path.join(indir,
                                                   "assembly_o",
                                                   "regular_quast",
                                                   "regular_quast*"),
                               target_dir=output_directory)
    run_cmd(cmd, dry_run=False)


@cli.command()
@click.option("-i", "--input_dir", "indir",
              help="output directory which passed to `run` command")
@click.option("-n", "--name", "name",
              help="Given a name for identify this archived files...")
def recovery(indir, name=None):
    # Just like the `recovery` command, but it will recovery archived files into it original directory.
    # it will overwrite required files
    def rev_mv(source, target, rm=True):
        cmd = ''
        if rm:
            cmd = "rm -r %s; " % target
            target = os.path.dirname(target)
        cmd += "mv %s %s ;" % (source, target)
        return cmd

    if name is None:
        name = [os.path.basename(_)
                for _ in glob(os.path.join(indir, "archived", "*"))
                ]
        name = [_ for _ in name
                if _.isnumeric()]
        name = str(max(name))

    in_directory = os.path.join(indir, "archived", name)
    cmd = ''
    # for workflow
    cmd += rev_mv(source=os.path.join(in_directory,
                                      "pipelines_summary"),
                  target=os.path.join(indir,
                                      "pipelines_summary"))
    # for roary
    cmd += rev_mv(source=os.path.join(in_directory,
                                      "all_roary_o"),
                  target=os.path.join(indir,
                                      "all_roary_o"))
    # for phigaro, seqtk, mlst, plasmid, IS,
    cmd += rev_mv(source=os.path.join(in_directory,
                                      constant.summary_dir),
                  target=os.path.join(indir,
                                      constant.summary_dir))
    # for abricate
    cmd += rev_mv(source=os.path.join(in_directory,
                                      "locus2annotate.csv"),
                  target=os.path.join(indir,
                                      "abricate_result"),
                  rm=False)
    # for multiqc
    cmd += rev_mv(source=os.path.join(in_directory,
                                      "fastqc_after*"),
                  target=os.path.join(indir,
                                      "fastqc_after"),
                  rm=False)
    cmd += rev_mv(source=os.path.join(in_directory,
                                      "fastqc_before*"),
                  target=os.path.join(indir,
                                      "fastqc_before"),
                  rm=False)
    cmd += rev_mv(source=os.path.join(in_directory,
                                      "regular_quast*", ),
                  target=os.path.join(indir,
                                      "assembly_o",
                                      "regular_quast"),
                  rm=False)
    print(cmd)
    input_cmd = input(
        f"Are you sure to recovery this file from {in_directory} \n It will first delete a lot of some directory at {indir}. such as `pipelines_summary`, `all_roary_o`, `summary_output`",
    )

    if str(input_cmd).lower() in ["y", "yes"]:
        run_cmd(cmd, dry_run=False)


@cli.command(help="analysis with test dataset, need to assign a output directory.")
@click.option("-o", "--odir", help="output directory for testing ...")
def testdata(odir):
    project_root_path = dirname(dirname(__file__))
    run_cmd(
        f"python3 {project_root_path}/pipelines/main.py run -- workflow --tab {project_root_path}/pipelines/test/test_input.tab --odir {odir} --workers 2 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command(help="test the dependency of required softwares. ")
@click.option("-u", "--update", help="updated the soft_db_path.py or not ")
def check(update):
    required_soft = {"fastqc", "multiqc", "shovill", "prokka", "roary", "quast", "pandoo", "fasttree", "ISEscan", "abricate", "phigaro",
                     "mlst", "kraken2", "seqtk"}
    for s in required_soft:
        path_exist = eval("os.path.exists(ori_path.%s_path)" % s)
        if path_exist:
            print("\033[1;32;40m {:<10}: software exists".format(s))
        else:
            print("\033[1;31;40m {:<10}: no requested files".format(s))


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)
    thread = luigi.IntParameter(default=constant.total_thread)

    def output(self):
        return luigi.LocalTarget(os.path.join(str(self.odir),
                                              "pipelines_summary"))

    def requires(self):
        inputdata = fileparser(self.tab)
        # it will validate the input df
        require_tasks = {}
        pairreads = inputdata.get_PE_info()
        singlereads = inputdata.get_SE_info()
        other_info = inputdata.get_full_PE_info()
        ############################################################
        valid_path(self.odir, check_odir=True)
        unify_kwargs = dict(odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path,
                            thread=self.thread,
                            PE_data=pairreads, )

        require_tasks["fastqc_before"] = multiqc(status='before',
                                                 **unify_kwargs
                                                 )
        require_tasks["fastqc_after"] = multiqc(status='after',
                                                **unify_kwargs
                                                )
        require_tasks["fastqc_quast"] = multiqc(status='quast',
                                                other_info=other_info,
                                                **unify_kwargs
                                                )
        require_tasks["abricate"] = abricate(SE_data=singlereads,
                                             **unify_kwargs)
        require_tasks["pre_roary"] = pre_roary(SE_data=singlereads,
                                             **unify_kwargs)
        # fixme
        require_tasks["seqtk"] = seqtk_summary(SE_data=singlereads,
                                               **unify_kwargs)
        require_tasks["species_annotated"] = species_annotated_summary(SE_data=singlereads,
                                                             **unify_kwargs)
        require_tasks["ISEscan_summary"] = ISEscan_summary(SE_data=singlereads,
                                                           **unify_kwargs)
        require_tasks["mlst_summary"] = mlst_summary(SE_data=singlereads,
                                                     **unify_kwargs)
        require_tasks["detect_plasmid"] = detect_plasmid(**unify_kwargs)
        require_tasks["detect_prophage"] = phigaro_summary(SE_data=singlereads,
                                                           **unify_kwargs)
        return require_tasks

    def run(self):
        # post pipelines
        post_analysis(self)


if __name__ == '__main__':
    cli()
    # luigi.run()

    # python  pipelines.luigi_pipelines.py workflow --tab /home/liaoth/project/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/project/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # local cmd

    # python -m luigi --module pipelines.luigi_pipelines workflow --tab /home/liaoth/tools/genome_pipelines/pipelines/test/test_input.tab --odir /home/liaoth/tools/genome_pipelines/pipelines/test/test_luigi  --parallel-scheduling --workers 12
    # server cmd
