import os
from os.path import exists, join


def get_project_root():
    from os.path import dirname
    project_root = dirname(dirname(__file__))
    return project_root


def env_exe(name):
    p = os.environ.get('PATH')
    p = p.split(':')
    for _p in p:
        if exists(join(_p, name)):
            return join(_p, name)


project_root = get_project_root()
############################################################
# exe path
############################################################


fastqc_path = env_exe('fastqc') if env_exe('fastqc') else "/tools/anaconda3/envs/pangenome_pipelines/bin/fastqc"

multiqc_path = env_exe('multiqc') if env_exe('multiqc') else "/tools/anaconda3/envs/pangenome_pipelines/bin/multiqc"

trimmomatic_dir = os.environ.get('HOME') + "/anaconda3/envs/wgs/share/trimmomatic-0.39-1"
trimmomatic_setting = "ILLUMINACLIP:%s/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5" % trimmomatic_dir
trimmomatic_jar_path = "%s/trimmomatic.jar" % trimmomatic_dir
fastp_extra_params = ''
fastp_path = env_exe('fastp') if env_exe('fastp') else "/tools/anaconda3/envs/pangenome_pipelines/bin/fastp"


shovill_path = env_exe('shovill') if env_exe('shovill') else "/tools/anaconda3/envs/pangenome_pipelines/bin/shovill"
prokka_path = env_exe('prokka') if env_exe('prokka') else "/tools/prokka/bin/prokka"
roary_path = env_exe('roary') if env_exe('roary') else "/tools/anaconda3/envs/pangenome_pipelines/bin/roary"
quast_path = env_exe('quast.py') if env_exe('quast.py') else "/tools/anaconda3/envs/pangenome_pipelines/bin/quast.py"
pandoo_path = env_exe('pandoo') if env_exe('pandoo') else "/tools/anaconda3/envs/pangenome_pipelines/bin/pandoo"
fasttree_path = "/usr/bin/fasttreeMP"

ISEscan_path = "/tools/ISEScan/isescan.py"
abricate_path = "/tools/abricate/bin/abricate"
abricate_py_path = "%s/toolkit/batch_abricate.py" % project_root
plasmid_detect_path = "%s/toolkit/process_plasmid.py" % project_root

phigaro_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/phigaro"
gubbins_path = "/tools/anaconda3/envs/pangenome_pipelines/bin/run_gubbins.py"

mlst_path = env_exe('mlst') if env_exe('mlst') else '/tools/mlst/bin/mlst'
kraken2_path = env_exe('kraken2') if env_exe('kraken2') else '/home-user/thliao/bin/kraken2'
seqtk_path = env_exe('seqtk') if env_exe('seqtk') else '/usr/bin/seqtk'
roary_plot_path = "/tools/Roary/contrib/roary_plots/roary_plots.py"
mash_path = env_exe('mash') if env_exe('mash') else "/tools/mash-Linux64-v2.1.1/mash"

############################################################
# database path
############################################################
abricate_db = "/tools/abricate/db/"

ariba_db = "/home/db_public/ariba_db"
ariba_db_list = ["CARD", "megares", "plasmidfinder", "resfinder", "srst2_argannot", "vfdb_full", "virulencefinder"]
ariba_str = str({k: "{db}/{name}".format(db=ariba_db, name=k.lower()) for k in ariba_db_list})

mlst_db = '/tools/mlst/db/pubmlst/'
phigaro_config = "/home/db_public/phigaro_db/config.yml"
kraken2_db = "/mnt/home-backup/thliao/kraken2_db"

mash_db = "/home-user/thliao/data/mash_db/RefSeq88n.msh"
# download from https://mash.readthedocs.io/en/latest/data.html
mash_db_summary = '/home-user/thliao/dataassembly_summary_refseq.txt'
# download from
