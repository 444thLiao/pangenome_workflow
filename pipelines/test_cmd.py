import sys
from os.path import dirname

sys.path.insert(0, dirname(dirname(__file__)))
import click
from pipelines import run_cmd

project_root_path = dirname(dirname(__file__))


@click.group()
def cli():
    pass


@cli.command()
@click.option("-o", "--odir", required=True,help="output directory for testing ...")
def full(odir):
    run_cmd(f"python3 {project_root_path}/pipelines/main.py workflow --tab {project_root_path}/pipelines/test/test_input.tab --odir {odir} --log-path {odir}/cmd_log.txt")


if __name__ == '__main__':
    cli()
