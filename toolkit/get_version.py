
from .utils import run_cmd
def seqtk_version():
    # todo
    pass


def get_krkn_version(kraken2_p):
    # todo
    '''
    Get the Kraken software version and Kraken Database path.
    '''
    return_text = run_cmd('%s -v' % kraken2_p,
                   get_output=True,
                   dry_run=False)
    version = return_text.rstrip().split('\n')[0].split(' ')[-1]
    return {'softwareKrakenVersion' : version}


def get_mlst_version(mlst_p):
    '''
    Get the MLST software version.
    '''
    return_text = run_cmd('%s -v' % mlst_p,
                   get_output=True,
                   dry_run=False)
    version = return_text.rstrip().split(' ')[1]
    return {'softwareMLSTversion': version}