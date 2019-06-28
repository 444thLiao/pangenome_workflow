import sys
from os.path import dirname
import time
import click

sys.path.insert(0, dirname(dirname(__file__)))
from pipelines import *


@click.group()
def cli():
    pass


@cli.command()
@click.argument('cmd',nargs=-1)
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
    "Just like the `clean` command, but it will archive these same result into a named directory."
    if name is None:
        name = str(int(time.time()))



@cli.command(help="analysis with test dataset, need to assign a output directory.")
@click.option("-o", "--odir", help="output directory for testing ...")
def testdata(odir):
    project_root_path = dirname(dirname(__file__))
    run_cmd(
        f"python3 {project_root_path}/pipelines/main.py run -- workflow --tab {project_root_path}/pipelines/test/test_input.tab --odir {odir} --workers 2 --log-path {odir}/cmd_log.txt",
        dry_run=False)


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)

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
        require_tasks["fasttree"] = fasttree(SE_data=singlereads,
                                             **unify_kwargs)
        require_tasks["seqtk"] = seqtk_summary(SE_data=singlereads,
                                               **unify_kwargs)
        require_tasks["kraken2"] = kraken2_summary(SE_data=singlereads,
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
