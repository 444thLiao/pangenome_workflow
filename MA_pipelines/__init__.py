import setting as config
from sub_api.cal_Cov_script import bam2info
import os,sys
from subprocess import check_call
from glob import glob
import time
import luigi

class base_luigi_task(luigi.Task):
    infodict = luigi.DictParameter(default=dict())
    dry_run = luigi.BoolParameter(default=False)
    round = luigi.Parameter(default=None)
    
    def get_log_path(self):
        base_log_path = self.infodict.get("log_path",None)
        if base_log_path is not None:
            return base_log_path
        

def record_cmdline(message, default):
    # fixme
    if os.path.isfile(default):
        with open(default, 'a') as f1:
            f1.write(time.ctime() + ' ' * 4 + message + '\n')
    else:
        with open(default, 'w') as f1:
            f1.write('{:#^40}'.format('Starting the somatic pipelines.'))
            f1.write(time.ctime() + ' ' * 4 + message + '\n')


def run_cmd(cmd, dry_run=False, log_file=None, **kwargs):
    outstream = None
    if type(log_file) == str:
        if os.path.isfile(log_file):
            if os.path.getsize(log_file) != 0:
                outstream = open(log_file, 'a')
        if outstream is None:
            valid_path(log_file,check_ofile=True)
            outstream = open(log_file, 'w')
    elif log_file is None:
        outstream = sys.stdout
    else:
        outstream = log_file

    print(cmd, file=outstream)
    outstream.flush()
    if not dry_run:
        check_call(cmd,
                   shell=True,
                   executable="/usr/bin/zsh",
                   stdout=outstream,
                   stderr=outstream,
                   **kwargs)
        outstream.flush()

def valid_path(in_pth,
               check_size=False,
               check_dir=False,
               check_glob=False,
               check_odir=False,
               check_ofile=False):
    if type(in_pth) == str:
        in_pths = [in_pth]
    else:
        in_pths = in_pth[::]
    for in_pth in in_pths:
        in_pth = os.path.abspath(os.path.realpath(in_pth))
        if in_pth is None:
            continue
        if check_glob:
            query_list = glob(in_pth)
            if not query_list:
                raise Exception('Error because of input file pattern %s' % in_pth)
        if check_dir:
            if not os.path.isdir(in_pth):
                raise Exception("Error because %s doesn't exist" % in_pth)
        if check_size:
            if os.path.getsize(in_pth) <= 0:
                raise Exception("Error because %s does not contain content." % in_pth)
        if check_odir:
            if not os.path.isdir(in_pth):
                os.makedirs(in_pth, exist_ok=True)
        if check_ofile:
            odir_file = os.path.dirname(in_pth)
            if not os.path.isdir(odir_file):
                os.makedirs(odir_file, exist_ok=True)
    return True