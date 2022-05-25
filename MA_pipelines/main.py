
import luigi,click
from os.path import *
import os,sys
__file__ = abspath(realpath(__file__))
sys.path.insert(0, dirname(dirname(__file__)))
import MA_pipelines
from MA_pipelines import run_cmd
from toolkit.parse_file_name import fileparser
from MA_pipelines.HaploCaller import CombineVariants
from MA_pipelines.BootBQSR import CV_4nd
# from MA_pipelines.BootBQSR import 
from collections import defaultdict

project_root_path = dirname(dirname(__file__))
@click.group()
def cli():
    pass


@cli.command()
@click.argument('cmd', nargs=-1)
def run(cmd):
    luigi.run(cmdline_args=cmd)

@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py workflow --tab {project_root_path}/test_set/germline/data_input.tsv --odir {odir}  --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


class workflow(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    bootBQSR = luigi.BoolParameter()
    log_path = luigi.Parameter(default=None)

    def requires(self):
        fp = fileparser(self.tab)
        ############################################################
        
        tasks = defaultdict(list)
        for SN,info in fp.get_full_PE_info().items():
            info_dict = {}
            info_dict['odir'] = self.odir
            info_dict['R1'] = info['R1']
            info_dict['R2'] = info['R2']
            info_dict['SampleID'] = SN
            info_dict['REF'] = info['ref']
            if not self.bootBQSR:
                tasks[f"{SN}_HC"].append(CombineVariants(infodict=info_dict,
                                                            dry_run=self.dry_run))
            else:
                tasks[f"{SN}_HC"].append(CV_4nd(infodict=info_dict,
                                                dry_run=self.dry_run))
            return tasks

            

    def run(self):
        pass

if __name__ == '__main__':
    cli()